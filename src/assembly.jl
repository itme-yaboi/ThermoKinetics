# Compilation and the assembler (DESIGN §4.2, §5.1, §10, §11).
#
# Laws map forces to fluxes (`thermal_transfer`, `fluxes`); heat kernels decide where the
# released power goes (`distribute_heat`); `_apply!` (the bookkeeping side) is the only code
# that writes state derivatives. It moves every non-thermal charge out of one port and into
# the other, turns every heat delivery into entropy at the receiving node's own temperature,
# and books the same terms into the ledgers.

# --- heat kernels ----------------------------------------------------------------------------

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

"Heat kernel of a thermal-only conductor (kernel A, thermal row)."
struct ThermalRow end

"""
    WeightsPolicy(w_a, w_b, w_h)

Compiled kernel A for non-thermal rows: heat w_i·D to target i, with Σ w = 1. `Split(w)` is
(1 − w, w, 0) and `Into(h)` is (0, 0, 1).
"""
struct WeightsPolicy
    w_a::Float64
    w_b::Float64
    w_h::Float64
end

"""
    NodePolicy(G_a, G_b, G_h)

Compiled kernel B: the zero-capacity internal node with leg conductances to a, b and h.
"""
struct NodePolicy
    G_a::Float64
    G_b::Float64
    G_h::Float64
end

active(p::WeightsPolicy) = (p.w_a > 0, p.w_b > 0, p.w_h > 0)
active(p::NodePolicy) = (p.G_a > 0, p.G_b > 0, p.G_h > 0)

"""
    distribute_heat(policy, D, Ta, Tb, Th) -> (; Q_a, Q_b, Q_h, S_gen, T_m, parts)

Distribute released power D ≥ 0 as heat over the active targets. Energy closes by the
remainder rule (the last active target takes D minus the others). Temperatures of inactive
targets are never read, so they may be NaN.

`S_gen` is the ledger value, computed in a form that is ≥ 0 in floating point:
- kernel A: Σ w_i D / T_i;
- kernel B: D/T_m (`parts.dissipation`) + Σ G_i (T_m − T_i)²/(T_i T_m) (`parts.conduction`).

It equals Σ Q_i/T_i, the entropy actually written, to round-off.
"""
function distribute_heat(p::WeightsPolicy, D, Ta, Tb, Th)
    z = zero(D)
    Q_a = p.w_a > 0 ? p.w_a * D : z
    Q_h = p.w_h > 0 ? p.w_h * D : z
    Q_b = p.w_b > 0 ? D - Q_a - Q_h : z
    if p.w_b == 0                               # remainder goes to the last active target
        p.w_h > 0 ? (Q_h = D - Q_a) : (Q_a = D - Q_h)
    end
    S_gen = (p.w_a > 0 ? Q_a / Ta : z) + (p.w_b > 0 ? Q_b / Tb : z) + (p.w_h > 0 ? Q_h / Th : z)
    return (; Q_a, Q_b, Q_h, S_gen, T_m = oftype(S_gen, NaN), parts = (; total = S_gen))
end

function distribute_heat(p::NodePolicy, D, Ta, Tb, Th)
    z = zero(D)
    act = active(p)
    Gsum = p.G_a + p.G_b + p.G_h
    T_m = (D + (act[1] ? p.G_a * Ta : z) + (act[2] ? p.G_b * Tb : z) + (act[3] ? p.G_h * Th : z)) / Gsum
    Q_a = act[1] ? p.G_a * (T_m - Ta) : z
    Q_b = act[2] ? p.G_b * (T_m - Tb) : z
    Q_h = act[3] ? p.G_h * (T_m - Th) : z
    if act[3]                                   # remainder rule on the last active leg
        Q_h = D - Q_a - Q_b
    elseif act[2]
        Q_b = D - Q_a
    else
        Q_a = D
    end
    leg(G, on, T) = on ? G * (T_m - T)^2 / (T * T_m) : z
    dissipation = D / T_m
    conduction = leg(p.G_a, act[1], Ta) + leg(p.G_b, act[2], Tb) + leg(p.G_h, act[3], Th)
    return (; Q_a, Q_b, Q_h, S_gen = dissipation + conduction, T_m, parts = (; dissipation, conduction))
end

"Ledger slots per conductor: kernel B reports its production in two parts."
nledger(::ThermalRow) = 1
nledger(::WeightsPolicy) = 1
nledger(::NodePolicy) = 2

# --- compiled representation -------------------------------------------------------------------

struct BodySlot{B<:Body}
    body::B
    index::Int     # node index (bodies are nodes 1..nb)
    offset::Int    # first state slot of this body's charges (= its first force slot)
end

struct CompiledConductor{L<:ConductorLaw,K,N}
    name::Symbol
    law::L
    kernel::K
    a::Int; b::Int; h::Int          # node indices (h = 0: none)
    fa::NTuple{N,Int}               # force-buffer index of each non-thermal row at a
    fb::NTuple{N,Int}
    ca::NTuple{N,Int}               # du index receiving (−) the charge at a
    cb::NTuple{N,Int}
    led::Int                        # du index of the first ledger slot
end

"""
    CompiledSystem

Flat, type-stable lowering of a `ThermoSystem`. State layout: `u = [u_phys; q_ledger]`.
`u_phys` holds the bodies' charges. `q_ledger` holds, in order:
- each conductor's integrated entropy production (two slots for kernel B: dissipation and
  conduction);
- for each reservoir, the charge it received in each of its domains (entropy for
  `Thermal`), then the energy it received.

Ledgers never feed back into the dynamics (D-10). Forces live in a buffer aligned with the
state, `P[1:nphys]` (bodies), followed by the reservoirs' fixed forces.
"""
struct CompiledSystem{BG<:Tuple,CG<:Tuple}
    sys::ThermoSystem
    body_groups::BG
    conductor_groups::CG
    # bodies
    body_names::Vector{Symbol}
    bodies::Vector{Body}
    body_offset::Vector{Int}
    # reservoirs
    reservoir_names::Vector{Symbol}
    reservoirs::Vector{Reservoir}
    res_force::Vector{Float64}          # P[nphys+1:end]
    res_led::Vector{Int}                # du index of the first ledger slot of each reservoir
    # per node (bodies then reservoirs)
    node_T::Vector{Int}                 # force-buffer index of T (0: no thermal DOF)
    node_S::Vector{Int}                 # du index receiving entropy (0: none)
    node_U::Vector{Int}                 # du index of the reservoir energy ledger (0 for bodies)
    # conductors
    conductor_names::Vector{Symbol}
    conductor_ref::Vector{Tuple{Int,Int}}
    conductor_led::Vector{UnitRange{Int}}
    # lookup and sizes
    node_index::Dict{Symbol,Int}
    conductor_index::Dict{Symbol,Int}
    nphys::Int
    nledger::Int
    u0::Vector{Float64}
    Pbuf::Vector{Float64}
end

nbodies(cs::CompiledSystem) = length(cs.body_names)
nnodes(cs::CompiledSystem) = length(cs.body_names) + length(cs.reservoir_names)

body_names(cs::CompiledSystem) = copy(cs.body_names)
reservoir_names(cs::CompiledSystem) = copy(cs.reservoir_names)
conductor_names(cs::CompiledSystem) = copy(cs.conductor_names)

"""
Group a vector into a tuple of concretely typed vectors, one per concrete element type.
Groups are sorted by type name, so models with the same element types share one compiled
`CompiledSystem` type (and one compiled solver) whatever order they were built in.
"""
function _group_by_type(items::AbstractVector)
    types = sort(unique(map(typeof, items)); by = string)
    return Tuple(Vector{T}([x for x in items if typeof(x) === T]) for T in types)
end

function _compile_kernel(el::Conductor, has_h::Bool)
    pol = el.dissipation
    pol === nothing && return ThermalRow()
    pol isa Split && return WeightsPolicy(1 - pol.w, pol.w, 0.0)
    pol isa Into && return WeightsPolicy(0.0, 0.0, 1.0)
    pol isa InternalNode && return NodePolicy(pol.G_a, pol.G_b, pol.G_h)
    error("unknown dissipation policy $(typeof(pol))")
end

"""
    compile(sys) -> CompiledSystem

Validate the model graph and lower it to flat data (DESIGN §11).
"""
function compile(sys::ThermoSystem)
    g = sys.graph
    body_labels = [n for n in sys.nodes if g[n] isa Body]
    res_labels = [n for n in sys.nodes if g[n] isa Reservoir]
    nb, nr = length(body_labels), length(res_labels)
    node_index = Dict{Symbol,Int}()
    for (i, n) in enumerate(body_labels); node_index[n] = i; end
    for (j, n) in enumerate(res_labels); node_index[n] = nb + j; end

    # bodies: state layout (force slots coincide with charge slots)
    slots = BodySlot[]
    body_offset = Int[]
    off = 1
    for (i, n) in enumerate(body_labels)
        b = g[n]
        push!(slots, BodySlot(b, i, off))
        push!(body_offset, off)
        off += ncharges(b.material)
    end
    nphys = off - 1

    # reservoir forces follow the body forces in the force buffer
    reservoirs = Reservoir[g[n] for n in res_labels]
    res_fstart = Int[]
    res_force = Float64[]
    for r in reservoirs
        push!(res_fstart, nphys + length(res_force) + 1)
        append!(res_force, r.forces)
    end

    # conductors: ledger layout first, then reservoir ledgers
    cond_labels = [e for e in sys.elements if g[e] isa Conductor]
    kernels = Any[]
    led = nphys + 1
    conductor_led = UnitRange{Int}[]
    for e in cond_labels
        el = g[e]
        has_h = any(nb_ -> g[e, nb_].role == :h, neighbor_labels(g, e))
        k = _compile_kernel(el, has_h)
        push!(kernels, k)
        push!(conductor_led, led:led+nledger(k)-1)
        led += nledger(k)
    end
    res_led = Int[]
    for r in reservoirs
        push!(res_led, led)
        led += length(r.domains) + 1
    end
    nledger_total = led - 1 - nphys

    # per-node lookup: force index and du target for each domain
    function force_index(i::Int, d::FormOfMotion)
        if i <= nb
            k = findfirst(==(d), collect(domains(slots[i].body)))
            return k === nothing ? 0 : body_offset[i] + k - 1
        else
            k = findfirst(==(d), reservoirs[i-nb].domains)
            return k === nothing ? 0 : res_fstart[i-nb] + k - 1
        end
    end
    function charge_target(i::Int, d::FormOfMotion)
        if i <= nb
            return force_index(i, d)                 # same slot as the charge
        else
            k = findfirst(==(d), reservoirs[i-nb].domains)
            return k === nothing ? 0 : res_led[i-nb] + k - 1
        end
    end
    node_T = [force_index(i, Thermal()) for i in 1:nb+nr]
    node_S = [charge_target(i, Thermal()) for i in 1:nb+nr]
    node_U = [i <= nb ? 0 : res_led[i-nb] + length(reservoirs[i-nb].domains) for i in 1:nb+nr]

    compiled = CompiledConductor[]
    for (k, e) in enumerate(cond_labels)
        el = g[e]
        a, b = port_node(sys, e, :a), port_node(sys, e, :b)
        hl = third_node(el.dissipation)
        ia, ib = node_index[a], node_index[b]
        ih = hl === nothing ? 0 : node_index[hl]
        # internal isolation (n ≤ m): every conducted domain must exist on both ends
        for d in domains(el), (n, i) in ((a, ia), (b, ib))
            force_index(i, d) == 0 &&
                throw(ArgumentError("conductor :$e conducts $d, but :$n has no such form of motion"))
        end
        # heat may only be delivered to thermal nodes
        kern = kernels[k]
        if !(kern isa ThermalRow)
            for (on, n, i) in zip(active(kern), (a, b, hl), (ia, ib, ih))
                on && node_T[i] == 0 &&
                    throw(ArgumentError("conductor :$e sends friction heat to :$n, which has no thermal form of motion"))
            end
        end
        rows = nonthermal_domains(el.law)
        N = length(rows)
        push!(compiled, CompiledConductor(e, el.law, kern, ia, ib, ih,
                                          ntuple(r -> force_index(ia, rows[r]), N), ntuple(r -> force_index(ib, rows[r]), N),
                                          ntuple(r -> charge_target(ia, rows[r]), N), ntuple(r -> charge_target(ib, rows[r]), N),
                                          first(conductor_led[k])))
    end

    body_groups = _group_by_type([s for s in slots])
    conductor_groups = _group_by_type([c for c in compiled])
    conductor_ref = Vector{Tuple{Int,Int}}(undef, length(cond_labels))
    cidx = Dict(e => k for (k, e) in enumerate(cond_labels))
    for (gi, grp) in enumerate(conductor_groups), (pos, c) in enumerate(grp)
        conductor_ref[cidx[c.name]] = (gi, pos)
    end

    u0 = zeros(nphys + nledger_total)
    for s in slots
        u0[s.offset:s.offset+ncharges(s.body.material)-1] .= s.body.E0
    end
    Pbuf = vcat(zeros(nphys), res_force)

    return CompiledSystem(sys, body_groups, conductor_groups,
                          body_labels, Body[g[n] for n in body_labels], body_offset,
                          res_labels, reservoirs, res_force, res_led,
                          node_T, node_S, node_U,
                          cond_labels, conductor_ref, conductor_led,
                          node_index, cidx, nphys, nledger_total, u0, Pbuf)
end

# --- the right-hand side ------------------------------------------------------------------------

@inline function _body_charges(u, s::BodySlot{<:Body{M,N}}) where {M,N}
    return SVector{N}(ntuple(k -> u[s.offset+k-1], Val(N)))
end

function _body_forces!(P, group::Vector{<:BodySlot}, u)
    @inbounds for s in group
        F = forces(s.body.material, _body_charges(u, s))
        for k in eachindex(F)
            P[s.offset+k-1] = F[k]
        end
    end
    return nothing
end

"The force buffer for state `u`: body forces, then reservoir forces."
function force_buffer(cs::CompiledSystem, u::AbstractVector)
    Tel = promote_type(eltype(u), Float64)
    P = Tel === Float64 ? cs.Pbuf : Vector{Tel}(undef, length(cs.Pbuf))
    for j in eachindex(cs.res_force)
        P[cs.nphys+j] = cs.res_force[j]
    end
    foreach(grp -> _body_forces!(P, grp, u), cs.body_groups)
    return P
end

@inline _node_T(P, cs::CompiledSystem, i::Int) = (k = cs.node_T[i]; k == 0 ? oftype(P[1] + 0.0, NaN) : P[k])

"Deliver heat `Q` [W] to node `i` at temperature `T`: entropy Q/T, plus the reservoir energy ledger."
@inline function _heat!(du, cs::CompiledSystem, i::Int, Q, T)
    @inbounds du[cs.node_S[i]] += Q / T
    @inbounds cs.node_U[i] == 0 || (du[cs.node_U[i]] += Q)
    return nothing
end

"""
    transfer(c, P, cs) -> NamedTuple

Evaluate one conductor at force buffer `P`: fluxes, released power, heat deliveries and
production. Pure; used by both the assembler and the `flux` accessor.
"""
function transfer(c::CompiledConductor{<:ConductorLaw,ThermalRow,0}, P, cs::CompiledSystem)
    Ta, Tb = _node_T(P, cs, c.a), _node_T(P, cs, c.b)
    tr = thermal_transfer(c.law, Ta, Tb)
    z = zero(tr.J_Q)
    return (; Q = tr.J_Q, J = SVector{0,typeof(z)}(), D = z, Q_a = tr.Q_a, Q_b = tr.Q_b, Q_h = z,
            T_a = Ta, T_b = Tb, T_h = oftype(Ta, NaN), T_m = oftype(Ta, NaN), S_gen = tr.S_gen,
            parts = (; total = tr.S_gen))
end

function transfer(c::CompiledConductor{<:ConductorLaw,K,N}, P, cs::CompiledSystem) where {K,N}
    ΔP = SVector{N}(ntuple(k -> P[c.fa[k]] - P[c.fb[k]], Val(N)))
    J = fluxes(c.law, ΔP)
    D = sum(ΔP .* J)
    on = active(c.kernel)
    Ta = on[1] ? _node_T(P, cs, c.a) : oftype(D, NaN)
    Tb = on[2] ? _node_T(P, cs, c.b) : oftype(D, NaN)
    Th = on[3] ? _node_T(P, cs, c.h) : oftype(D, NaN)
    h = distribute_heat(c.kernel, D, Ta, Tb, Th)
    return (; Q = zero(D), J, D, h.Q_a, h.Q_b, h.Q_h, T_a = _node_T(P, cs, c.a), T_b = _node_T(P, cs, c.b),
            T_h = c.h == 0 ? oftype(D, NaN) : _node_T(P, cs, c.h), h.T_m, h.S_gen, h.parts)
end

function _apply!(du, P, c::CompiledConductor{<:ConductorLaw,ThermalRow,0}, cs)
    tr = transfer(c, P, cs)
    _heat!(du, cs, c.a, tr.Q_a, tr.T_a)
    _heat!(du, cs, c.b, tr.Q_b, tr.T_b)
    @inbounds du[c.led] += tr.S_gen
    return nothing
end

function _apply!(du, P, c::CompiledConductor{<:ConductorLaw,K,N}, cs) where {K,N}
    tr = transfer(c, P, cs)
    @inbounds for k in 1:N
        du[c.ca[k]] -= tr.J[k]
        du[c.cb[k]] += tr.J[k]
        # a charge moved into or out of a reservoir carries energy P·J
        cs.node_U[c.a] == 0 || (du[cs.node_U[c.a]] -= P[c.fa[k]] * tr.J[k])
        cs.node_U[c.b] == 0 || (du[cs.node_U[c.b]] += P[c.fb[k]] * tr.J[k])
    end
    on = active(c.kernel)
    on[1] && _heat!(du, cs, c.a, tr.Q_a, tr.T_a)
    on[2] && _heat!(du, cs, c.b, tr.Q_b, tr.T_b)
    on[3] && _heat!(du, cs, c.h, tr.Q_h, tr.T_h)
    _ledger!(du, c, tr.parts)
    return nothing
end

@inline _ledger!(du, c, parts::NamedTuple{(:total,)}) = (@inbounds du[c.led] += parts.total; nothing)
@inline function _ledger!(du, c, parts::NamedTuple{(:dissipation, :conduction)})
    @inbounds du[c.led] += parts.dissipation
    @inbounds du[c.led+1] += parts.conduction
    return nothing
end

function _apply_group!(du, P, group::Vector{<:CompiledConductor}, cs)
    for c in group
        _apply!(du, P, c, cs)
    end
    return nothing
end

"""
    assemble!(du, u, cs, t)

Physical and ledger rates in one pass, into the full-length `du`. Only `u[1:nphys]` is read.
All bookkeeping is done here.
"""
function assemble!(du, u, cs::CompiledSystem, t)
    P = force_buffer(cs, u)
    fill!(du, zero(eltype(du)))
    foreach(grp -> _apply_group!(du, P, grp, cs), cs.conductor_groups)
    return nothing
end

"In-place ODE right-hand side; the parameter is the compiled system."
ode_rhs!(du, u, cs::CompiledSystem, t) = assemble!(du, u, cs, t)

"Whether every body's charges in `u` are a valid state of its material."
function all_valid(cs::CompiledSystem, u::AbstractVector)
    for (i, b) in enumerate(cs.bodies)
        N = ncharges(b.material)
        off = cs.body_offset[i]
        valid_state(b.material, view(u, off:off+N-1)) || return false
    end
    return true
end

"""
    jacobian_sparsity(cs) -> SparseMatrixCSC

Structural Jacobian pattern from the graph topology: each conductor couples all charges of
the bodies on its ports; ledger rows depend on the bodies they observe; ledger columns are
empty.
"""
function jacobian_sparsity(cs::CompiledSystem)
    n = cs.nphys + cs.nledger
    I, J = collect(1:n), collect(1:n)
    nb = nbodies(cs)
    body_cols(i) = (off = cs.body_offset[i]; off:off+ncharges(cs.bodies[i].material)-1)
    for grp in cs.conductor_groups, c in grp
        nodes = [i for i in (c.a, c.b, c.h) if i != 0]
        cols = reduce(vcat, [collect(body_cols(i)) for i in nodes if i <= nb]; init = Int[])
        rows = copy(cols)
        append!(rows, cs.conductor_led[cs.conductor_index[c.name]])
        for i in nodes
            i > nb || continue
            j = i - nb
            append!(rows, cs.res_led[j]:cs.res_led[j]+length(cs.reservoirs[j].domains))
        end
        for r in rows, col in cols
            push!(I, r); push!(J, col)
        end
    end
    return sparse(I, J, ones(length(I)), n, n)
end
