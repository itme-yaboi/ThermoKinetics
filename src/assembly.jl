# Compilation and the assembler (DESIGN §4.2, §5.1, §10, §11).
#
# `thermal_transfer` (the law side) maps forces to fluxes and heat deliveries.
# `_apply_conductors!` (the bookkeeping side) is the only code that writes state
# derivatives: it moves heat deliveries into entropy at each target's own temperature and
# books the same terms into the ledgers, so the entropy balance closes by construction.

# --- kernel A, thermal row --------------------------------------------------------------

"""
    thermal_transfer(law, Ta, Tb) -> (; J_Q, X_Q, Q_a, Q_b, S_gen)

Kernel A (07 §3.1) for a conductor with only a thermal row. X_Q = 1/T_b − 1/T_a,
J_Q = M_QQ X_Q. Heat deliveries: Q_b = J_Q into b and Q_a = D − Q_b into a, where the released
non-thermal power D is zero for a thermal-only conductor, so Q_a + Q_b == 0 exactly.
S_gen = Q_a/T_a + Q_b/T_b is exactly what is written to the two ends, and it is ≥ 0 in
floating point: J_Q has the sign of T_a − T_b, and so does 1/T_b − 1/T_a.
"""
@inline function thermal_transfer(law::ConductorLaw, Ta, Tb)
    X_Q = (Ta - Tb) / (Ta * Tb)               # = 1/T_b − 1/T_a, without cancellation
    J_Q = conductance_QQ(law, Ta, Tb) * X_Q
    D = zero(J_Q)                             # no non-thermal rows
    Q_b = J_Q
    Q_a = D - Q_b                             # remainder rule: Σ Q_i = D exactly
    S_gen = Q_a / Ta + Q_b / Tb
    return (; J_Q, X_Q, Q_a, Q_b, S_gen)
end

# --- compiled representation ------------------------------------------------------------

struct BodySlot{B<:Body}
    body::B
    index::Int     # node index (bodies are nodes 1..nb)
    offset::Int    # first state slot of this body's charges
    tslot::Int     # state slot of its entropy (0 if it has no thermal form)
end

struct CompiledConductor{L<:ConductorLaw}
    name::Symbol
    law::L
    a::Int         # node index of port :a
    b::Int         # node index of port :b
    ledger::Int    # ledger slot of its Ṡ_gen
end

"""
    CompiledSystem

Flat, type-stable lowering of a `ThermoSystem`. State layout: `u = [u_phys; q_ledger]`.
`u_phys` holds the bodies' charges. `q_ledger` holds, in order, the integrated entropy
production of every conductor, then (entropy, energy) received by every reservoir. Ledgers
never feed back into the dynamics (D-10).
"""
struct CompiledSystem{BG<:Tuple,CG<:Tuple}
    sys::ThermoSystem
    body_groups::BG
    conductor_groups::CG
    # bodies
    body_names::Vector{Symbol}
    bodies::Vector{Body}
    body_offset::Vector{Int}
    body_tslot::Vector{Int}
    # reservoirs
    reservoir_names::Vector{Symbol}
    T_res::Vector{Float64}
    res_ledger::Vector{Int}        # ledger slot of received entropy; energy is the next slot
    # conductors
    conductor_names::Vector{Symbol}
    conductor_ref::Vector{Tuple{Int,Int}}   # (group, position) for accessors
    # lookup
    node_index::Dict{Symbol,Int}
    conductor_index::Dict{Symbol,Int}
    # sizes
    nphys::Int
    nledger::Int
    u0::Vector{Float64}
    Tbuf::Vector{Float64}
end

nbodies(cs::CompiledSystem) = length(cs.body_names)
nnodes(cs::CompiledSystem) = length(cs.body_names) + length(cs.reservoir_names)

body_names(cs::CompiledSystem) = copy(cs.body_names)
reservoir_names(cs::CompiledSystem) = copy(cs.reservoir_names)
conductor_names(cs::CompiledSystem) = copy(cs.conductor_names)

"Group a vector into a tuple of concretely typed vectors, one per concrete element type."
function _group_by_type(items::AbstractVector)
    types = unique(map(typeof, items))
    return Tuple(Vector{T}([x for x in items if typeof(x) === T]) for T in types)
end

"""
    compile(sys) -> CompiledSystem

Validate the model graph and lower it to flat data (DESIGN §11).
"""
function compile(sys::ThermoSystem)
    g = sys.graph
    body_labels = [n for n in sys.nodes if g[n] isa Body]
    res_labels = [n for n in sys.nodes if g[n] isa Reservoir]
    nb = length(body_labels)
    node_index = Dict{Symbol,Int}()
    for (i, n) in enumerate(body_labels)
        node_index[n] = i
    end
    for (j, n) in enumerate(res_labels)
        node_index[n] = nb + j
    end

    # bodies: physical state layout
    slots = BodySlot[]
    body_offset = Int[]
    body_tslot = Int[]
    off = 1
    for (i, n) in enumerate(body_labels)
        b = g[n]
        k = thermal_slot(b.material)
        tslot = k === nothing ? 0 : off + k - 1
        push!(slots, BodySlot(b, i, off, tslot))
        push!(body_offset, off)
        push!(body_tslot, tslot)
        off += ncharges(b.material)
    end
    nphys = off - 1

    # conductors: validation and ledger layout
    cond_labels = [e for e in sys.elements if g[e] isa Conductor]
    nc = length(cond_labels)
    compiled = CompiledConductor[]
    for (k, e) in enumerate(cond_labels)
        el = g[e]
        a, b = port_node(sys, e, :a), port_node(sys, e, :b)
        for d in domains(el), n in (a, b)
            d isa Thermal || throw(ArgumentError("slice 1 supports thermal conductors only (:$e)"))
            has_thermal(g[n]) || throw(ArgumentError("conductor :$e needs a thermal form of motion on :$n"))
        end
        push!(compiled, CompiledConductor(e, el.law, node_index[a], node_index[b], k))
    end
    res_ledger = [nc + 2j - 1 for j in 1:length(res_labels)]
    nledger = nc + 2 * length(res_labels)

    body_groups = _group_by_type([s for s in slots])
    conductor_groups = _group_by_type([c for c in compiled])
    conductor_ref = Vector{Tuple{Int,Int}}(undef, nc)
    for (gi, grp) in enumerate(conductor_groups), (pos, c) in enumerate(grp)
        conductor_ref[c.ledger] = (gi, pos)
    end

    u0 = zeros(nphys + nledger)
    for s in slots
        u0[s.offset:s.offset+ncharges(s.body.material)-1] .= s.body.E0
    end

    return CompiledSystem(sys, body_groups, conductor_groups,
                          body_labels, Body[g[n] for n in body_labels], body_offset, body_tslot,
                          res_labels, Float64[g[n].T for n in res_labels], res_ledger,
                          cond_labels, conductor_ref,
                          node_index, Dict(e => k for (k, e) in enumerate(cond_labels)),
                          nphys, nledger, u0, zeros(nb + length(res_labels)))
end

# --- the right-hand side ----------------------------------------------------------------

@inline function _body_charges(u, s::BodySlot{<:Body{M,N}}) where {M,N}
    return SVector{N}(ntuple(k -> u[s.offset+k-1], Val(N)))
end

function _body_temperatures!(T, group::Vector{<:BodySlot}, uphys)
    @inbounds for s in group
        s.tslot == 0 && continue
        T[s.index] = forces(s.body.material, _body_charges(uphys, s))[s.tslot-s.offset+1]
    end
    return nothing
end

"Deliver heat `Q` [W] to node `i` at its temperature `T`: entropy Q/T, plus reservoir ledgers."
@inline function _deliver!(dphys, dled, cs::CompiledSystem, i::Int, Q, T)
    nb = nbodies(cs)
    if i <= nb
        @inbounds dphys[cs.body_tslot[i]] += Q / T
    else
        j = cs.res_ledger[i-nb]
        @inbounds dled[j] += Q / T        # entropy received
        @inbounds dled[j+1] += Q          # energy received
    end
    return nothing
end

function _apply_conductors!(dphys, dled, T, group::Vector{<:CompiledConductor}, cs)
    @inbounds for c in group
        Ta, Tb = T[c.a], T[c.b]
        tr = thermal_transfer(c.law, Ta, Tb)
        _deliver!(dphys, dled, cs, c.a, tr.Q_a, Ta)
        _deliver!(dphys, dled, cs, c.b, tr.Q_b, Tb)
        dled[c.ledger] += tr.S_gen
    end
    return nothing
end

"""
    assemble!(dphys, dled, uphys, cs, t)

Physical rates and ledger rates in one pass. Bookkeeping is done only here.
"""
function assemble!(dphys, dled, uphys, cs::CompiledSystem, t)
    Tel = promote_type(eltype(uphys), Float64)
    T = Tel === Float64 ? cs.Tbuf : Vector{Tel}(undef, nnodes(cs))
    nb = nbodies(cs)
    for j in eachindex(cs.T_res)
        T[nb+j] = cs.T_res[j]
    end
    foreach(grp -> _body_temperatures!(T, grp, uphys), cs.body_groups)
    fill!(dphys, zero(eltype(dphys)))
    fill!(dled, zero(eltype(dled)))
    foreach(grp -> _apply_conductors!(dphys, dled, T, grp, cs), cs.conductor_groups)
    return nothing
end

"In-place ODE right-hand side; the parameter is the compiled system."
function ode_rhs!(du, u, cs::CompiledSystem, t)
    np = cs.nphys
    assemble!(view(du, 1:np), view(du, np+1:length(du)), view(u, 1:np), cs, t)
    return nothing
end

"""
    jacobian_sparsity(cs) -> SparseMatrixCSC

Structural Jacobian pattern from the graph topology: each conductor couples the charges of
its end bodies; ledger rows depend on the bodies they observe; ledger columns are empty.
"""
function jacobian_sparsity(cs::CompiledSystem)
    n = cs.nphys + cs.nledger
    I, J = collect(1:n), collect(1:n)
    nb = nbodies(cs)
    for grp in cs.conductor_groups, c in grp
        ends = [i for i in (c.a, c.b) if i <= nb]
        cols = [cs.body_tslot[i] for i in ends]
        rows = copy(cols)
        push!(rows, cs.nphys + c.ledger)
        for i in (c.a, c.b)
            if i > nb
                j = cs.res_ledger[i-nb]
                push!(rows, cs.nphys + j, cs.nphys + j + 1)
            end
        end
        for r in rows, col in cols
            push!(I, r); push!(J, col)
        end
    end
    return sparse(I, J, ones(length(I)), n, n)
end
