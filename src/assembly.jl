# Compilation and the assembler (DESIGN §4.2, §5.1, §10, §11).
#
# Laws map forces to fluxes (`thermal_transfer`, `fluxes`, flow drivers); heat kernels decide
# where the released power goes (`distribute_heat`); `_apply!` (the bookkeeping side) is the
# only code that writes state derivatives. It moves every non-thermal charge out of one port and into
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

"A stream's species rows at both ends: force index of μ, charge target, reservoir y index."
struct CompiledStream{D,N}
    name::Symbol
    driver::D
    a::Int; b::Int
    pa::Int; pb::Int                # force index of −p at each end
    fa::NTuple{N,Int}; fb::NTuple{N,Int}
    ca::NTuple{N,Int}; cb::NTuple{N,Int}
    ya::NTuple{N,Int}; yb::NTuple{N,Int}
    led::Int                        # Ṡ_gen
end

struct CompiledPump{F,M,N}
    name::Symbol
    mdot::F
    mat::M                          # the upstream fluid (for the isentropic delivered state)
    a::Int; b::Int
    pb::Int
    fa::NTuple{N,Int}; fb::NTuple{N,Int}
    ca::NTuple{N,Int}; cb::NTuple{N,Int}
    ya::NTuple{N,Int}
    led::Int                        # Ṡ_gen, then external work
end

struct CompiledBoundary{F,M,N}
    name::Symbol
    mdot::F
    mat::M
    pumped::Bool
    T_in::Float64; p_in::Float64; s_in::Float64
    μ_in::NTuple{N,Float64}; y::NTuple{N,Float64}
    d::Int; pd::Int
    fd::NTuple{N,Int}; cd::NTuple{N,Int}
    led::Int                        # Ṡ_gen, work, received m_1..m_N, S, U
end

struct CompiledOutlet{O,N}
    name::Symbol
    outlet::O
    u::Int; pu::Int
    fu::NTuple{N,Int}; cu::NTuple{N,Int}
    led::Int                        # received m_1..m_N, S, U
end

"""
    CompiledSystem

Flat, type-stable lowering of a `ThermoSystem`. State layout: `u = [u_phys; q_ledger]`.
`u_phys` holds the bodies' charges. `q_ledger` holds, element by element in insertion order:
- conductors: entropy production (two slots for kernel B: dissipation and conduction);
- streams: production;
- pumps: production, then external work supplied;
- flow boundaries: production, work, and the mass of each species, entropy and energy
  *received* (negative: they supply);
- outlets: mass of each species, entropy and energy received.

Then, for each reservoir: the charge received in each of its domains (entropy for
`Thermal`), then energy. Ledgers never feed back into the dynamics (D-10). Forces live in a
buffer aligned with the state: `P[1:nphys]` are the bodies' forces, followed by the
reservoirs' fixed forces.
"""
struct CompiledSystem
    sys::ThermoSystem
    body_groups::Vector{Any}       # concretely typed vectors, one per body type (function barrier)
    element_groups::Vector{Any}    # likewise per element type
    # bodies
    body_names::Vector{Symbol}
    bodies::Vector{Body}
    body_offset::Vector{Int}
    # reservoirs
    reservoir_names::Vector{Symbol}
    reservoirs::Vector{Reservoir}
    res_force::Vector{Float64}
    res_led::Vector{Int}
    # per node (bodies then reservoirs)
    node_T::Vector{Int}                 # force index of T (0: no thermal DOF)
    node_S::Vector{Int}                 # du index receiving entropy (0: none)
    node_U::Vector{Int}                 # du index of a reservoir's energy ledger (0 for bodies)
    # elements
    element_names::Vector{Symbol}
    element_ref::Dict{Symbol,Tuple{Int,Int}}
    conductor_names::Vector{Symbol}
    prod_led::Dict{Symbol,UnitRange{Int}}             # Ṡ_gen slots of every producing element
    work_led::Dict{Symbol,Int}                         # external work of pumps and boundaries
    bound_led::Dict{Symbol,NamedTuple{(:m, :S, :U),Tuple{UnitRange{Int},Int,Int}}}
    ext_charges::Vector{Tuple{FormOfMotion,Int}}       # boundary/outlet received charges
    touch::Vector{Tuple{Vector{Int},UnitRange{Int}}}   # per element: nodes, ledger slots (sparsity)
    moves::Vector{Tuple{Vector{Int},Vector{FormOfMotion},Bool}}  # per element: nodes, charges moved, is boundary
    active::Vector{Int}                                # physical slots some element writes
    # lookup and sizes
    node_index::Dict{Symbol,Int}
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
element_names(cs::CompiledSystem) = copy(cs.element_names)

"""
Group a vector into concretely typed vectors, one per concrete element type, held in a
`Vector{Any}`. Each group's loop is specialized (one dynamic dispatch per group per call, a
function barrier), while `CompiledSystem` keeps a single concrete type. So every model shares
one compiled right-hand side and one compiled solver (D-26).
"""
function _group_by_type(items::AbstractVector)
    types = sort(unique(map(typeof, items)); by = string)
    return Any[Vector{T}([x for x in items if typeof(x) === T]) for T in types]
end

function _compile_kernel(el::Conductor)
    pol = el.dissipation
    pol === nothing && return ThermalRow()
    pol isa Split && return WeightsPolicy(1 - pol.w, pol.w, 0.0)
    pol isa Into && return WeightsPolicy(0.0, 0.0, 1.0)
    pol isa InternalNode && return NodePolicy(pol.G_a, pol.G_b, pol.G_h)
    error("unknown dissipation policy $(typeof(pol))")
end

_nledger(el::Conductor) = nledger(_compile_kernel(el))
_nledger(::Stream) = 1
_nledger(::IdealPump) = 2
_nledger(el::FlowBoundary) = 2 + length(mass_slots(el.material)) + 2
_nledger(el::PassThroughOutlet, nspecies::Int) = nspecies + 2

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
    node_obj(i) = i <= nb ? g[body_labels[i]] : g[res_labels[i-nb]]

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

    reservoirs = Reservoir[g[n] for n in res_labels]
    res_fstart = Int[]
    res_force = Float64[]
    for r in reservoirs
        push!(res_fstart, nphys + length(res_force) + 1)
        append!(res_force, r.forces)
    end

    # per-node lookup
    ndoms(i) = i <= nb ? collect(FormOfMotion, domains(slots[i].body)) : reservoirs[i-nb].domains
    function force_index(i::Int, d::FormOfMotion)
        k = findfirst(==(d), ndoms(i))
        k === nothing && return 0
        return i <= nb ? body_offset[i] + k - 1 : res_fstart[i-nb] + k - 1
    end
    species_of(i) = [d for d in ndoms(i) if d isa Species]
    yindex(i, d) = i <= nb ? 0 : findfirst(==(d), species_of(i))

    # ledger layout, element by element in insertion order, then reservoirs
    elem_labels = copy(sys.elements)
    led = nphys + 1
    elem_led = Dict{Symbol,UnitRange{Int}}()
    for e in elem_labels
        el = g[e]
        n = el isa PassThroughOutlet ? _nledger(el, length(species_of(node_index[port_node(sys, e, :a)]))) : _nledger(el)
        elem_led[e] = led:led+n-1
        led += n
    end
    res_led = Int[]
    for r in reservoirs
        push!(res_led, led)
        led += length(r.domains) + 1
    end
    nledger_total = led - 1 - nphys

    function charge_target(i::Int, d::FormOfMotion)
        i <= nb && return force_index(i, d)                  # a body's charge slot
        k = findfirst(==(d), reservoirs[i-nb].domains)
        return k === nothing ? 0 : res_led[i-nb] + k - 1
    end
    node_T = [force_index(i, Thermal()) for i in 1:nb+nr]
    node_S = [charge_target(i, Thermal()) for i in 1:nb+nr]
    node_U = [i <= nb ? 0 : res_led[i-nb] + length(reservoirs[i-nb].domains) for i in 1:nb+nr]

    compiled = Any[]
    conductor_labels = Symbol[]
    prod_led = Dict{Symbol,UnitRange{Int}}()
    work_led = Dict{Symbol,Int}()
    bound_led = Dict{Symbol,NamedTuple{(:m, :S, :U),Tuple{UnitRange{Int},Int,Int}}}()
    ext_charges = Tuple{FormOfMotion,Int}[]
    touch = Tuple{Vector{Int},UnitRange{Int}}[]
    moves = Tuple{Vector{Int},Vector{FormOfMotion},Bool}[]
    active_set = Set{Int}()
    mark_active!(targets) = foreach(t -> (0 < t <= nphys && push!(active_set, t)), targets)

    function fluid_node(e, i, label)
        node = node_obj(i)
        for d in (Thermal(), Volumetric())
            force_index(i, d) == 0 && throw(ArgumentError(":$e needs a fluid at :$label, which has no $d form of motion"))
        end
        isempty(species_of(i)) && throw(ArgumentError(":$e needs a fluid at :$label, which has no species (mass) charge"))
        node isa Reservoir && node.material === nothing &&
            throw(ArgumentError(":$e: reservoir :$label must carry a material (Reservoir(…; material = …)) to exchange fluid"))
        return nothing
    end
    function same_species(e, i, j, li, lj)
        Set(species_of(i)) == Set(species_of(j)) ||
            throw(ArgumentError(":$e joins :$li and :$lj, whose species differ ($(species_of(i)) vs $(species_of(j)))"))
        return nothing
    end

    for e in elem_labels
        el = g[e]
        r = elem_led[e]
        if el isa Conductor
            a, b = port_node(sys, e, :a), port_node(sys, e, :b)
            hl = third_node(el.dissipation)
            ia, ib = node_index[a], node_index[b]
            ih = hl === nothing ? 0 : node_index[hl]
            for d in domains(el), (n, i) in ((a, ia), (b, ib))
                force_index(i, d) == 0 &&
                    throw(ArgumentError("conductor :$e conducts $d, but :$n has no such form of motion"))
            end
            kern = _compile_kernel(el)
            if !(kern isa ThermalRow)
                for (on, n, i) in zip(active(kern), (a, b, hl), (ia, ib, ih))
                    on && node_T[i] == 0 &&
                        throw(ArgumentError("conductor :$e sends friction heat to :$n, which has no thermal form of motion"))
                end
            end
            rows = nonthermal_domains(el.law)
            N = length(rows)
            c = CompiledConductor(e, el.law, kern, ia, ib, ih,
                                  ntuple(k -> force_index(ia, rows[k]), N), ntuple(k -> force_index(ib, rows[k]), N),
                                  ntuple(k -> charge_target(ia, rows[k]), N), ntuple(k -> charge_target(ib, rows[k]), N),
                                  first(r))
            push!(compiled, c); push!(conductor_labels, e)
            prod_led[e] = r
            nodes = [i for i in (ia, ib, ih) if i != 0]
            push!(touch, (nodes, r))
            push!(moves, ([ia, ib], collect(FormOfMotion, rows), false))
            mark_active!(c.ca); mark_active!(c.cb)
            heated = kern isa ThermalRow ? (ia, ib) : Tuple(i for (on, i) in zip(active(kern), (ia, ib, ih)) if on)
            mark_active!(node_S[i] for i in heated)
        elseif el isa Union{Stream,IdealPump}
            a, b = port_node(sys, e, :a), port_node(sys, e, :b)
            ia, ib = node_index[a], node_index[b]
            fluid_node(e, ia, a); fluid_node(e, ib, b); same_species(e, ia, ib, a, b)
            sp = species_of(ia)
            N = length(sp)
            fa, fb = ntuple(k -> force_index(ia, sp[k]), N), ntuple(k -> force_index(ib, sp[k]), N)
            ca, cb = ntuple(k -> charge_target(ia, sp[k]), N), ntuple(k -> charge_target(ib, sp[k]), N)
            ya = ntuple(k -> yindex(ia, sp[k]), N)
            if el isa Stream
                c = CompiledStream(e, el.driver, ia, ib, force_index(ia, Volumetric()), force_index(ib, Volumetric()),
                                   fa, fb, ca, cb, ya, ntuple(k -> yindex(ib, sp[k]), N), first(r))
            else
                N == 1 || throw(ArgumentError("ideal pump :$e: single-species fluids only (for now)"))
                up = node_obj(ia)
                mat = up isa Body ? up.material : up.material
                c = CompiledPump(e, el.mdot, mat, ia, ib, force_index(ib, Volumetric()), fa, fb, ca, cb, ya, first(r))
                work_led[e] = first(r) + 1
            end
            push!(compiled, c)
            prod_led[e] = first(r):first(r)
            push!(touch, ([ia, ib], r))
            push!(moves, ([ia, ib], collect(FormOfMotion, sp), false))
            mark_active!(ca); mark_active!(cb); mark_active!((node_S[ia], node_S[ib]))
        elseif el isa FlowBoundary
            d = port_node(sys, e, :a)
            id = node_index[d]
            id <= nb || throw(ArgumentError("flow boundary :$e must feed a body, not a reservoir"))
            fluid_node(e, id, d)
            msp = [dd for dd in domains(el.material) if dd isa Species]
            Set(msp) == Set(species_of(id)) ||
                throw(ArgumentError("flow boundary :$e supplies $msp, but :$d holds $(species_of(id))"))
            N = length(msp)
            N == 1 || throw(ArgumentError("flow boundary :$e: single-species fluids only (for now)"))
            E = charges_at(el.material; T = el.T, p = el.p, m = 1.0)
            F = forces(el.material, E)
            ms = mass_slots(el.material)
            mtot = sum(E[i] for i in ms)
            c = CompiledBoundary(e, el.mdot, el.material, el.supply === :pumped, el.T, el.p,
                                 E[thermal_slot(el.material)] / mtot,
                                 ntuple(k -> F[ms[k]], N), ntuple(k -> E[ms[k]] / mtot, N),
                                 id, force_index(id, Volumetric()),
                                 ntuple(k -> force_index(id, msp[k]), N), ntuple(k -> charge_target(id, msp[k]), N),
                                 first(r))
            push!(compiled, c)
            prod_led[e] = first(r):first(r)
            work_led[e] = first(r) + 1
            bound_led[e] = (m = first(r)+2:first(r)+1+N, S = first(r) + 2 + N, U = first(r) + 3 + N)
            for k in 1:N; push!(ext_charges, (msp[k], first(r) + 1 + k)); end
            push!(touch, ([id], r))
            push!(moves, ([id], collect(FormOfMotion, msp), true))
            mark_active!(c.cd); mark_active!((node_S[id],))
        elseif el isa PassThroughOutlet
            uu = port_node(sys, e, :a)
            iu = node_index[uu]
            iu <= nb || throw(ArgumentError("outlet :$e must drain a body, not a reservoir"))
            fluid_node(e, iu, uu)
            sp = species_of(iu)
            N = length(sp)
            c = CompiledOutlet(e, el, iu, force_index(iu, Volumetric()),
                               ntuple(k -> force_index(iu, sp[k]), N), ntuple(k -> charge_target(iu, sp[k]), N), first(r))
            push!(compiled, c)
            bound_led[e] = (m = first(r):first(r)+N-1, S = first(r) + N, U = first(r) + N + 1)
            for k in 1:N; push!(ext_charges, (sp[k], first(r) + k - 1)); end
            push!(touch, ([iu], r))
            push!(moves, ([iu], collect(FormOfMotion, sp), true))
            mark_active!(c.cu); mark_active!((node_S[iu],))
        else
            error("unsupported element $(typeof(el))")
        end
    end

    body_groups = _group_by_type([s for s in slots])
    element_groups = _group_by_type([c for c in compiled])
    element_ref = Dict{Symbol,Tuple{Int,Int}}()
    for (gi, grp) in enumerate(element_groups), (pos, c) in enumerate(grp)
        element_ref[c.name] = (gi, pos)
    end

    u0 = zeros(nphys + nledger_total)
    for s in slots
        u0[s.offset:s.offset+ncharges(s.body.material)-1] .= s.body.E0
    end
    Pbuf = vcat(zeros(nphys), res_force)

    return CompiledSystem(sys, body_groups, element_groups,
                          body_labels, Body[g[n] for n in body_labels], body_offset,
                          res_labels, reservoirs, res_force, res_led,
                          node_T, node_S, node_U,
                          elem_labels, element_ref, conductor_labels, prod_led, work_led, bound_led, ext_charges,
                          touch, moves, sort(collect(active_set)),
                          node_index, nphys, nledger_total, u0, Pbuf)
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
    for grp in cs.body_groups
        _body_forces!(P, grp, u)
    end
    return P
end

@inline _node_T(P, cs::CompiledSystem, i::Int) = (k = cs.node_T[i]; k == 0 ? convert(eltype(P), NaN) : P[k])

"Deliver heat `Q` [W] to node `i` at temperature `T`: entropy Q/T, plus the reservoir energy ledger."
@inline function _heat!(du, cs::CompiledSystem, i::Int, Q, T)
    @inbounds du[cs.node_S[i]] += Q / T
    @inbounds cs.node_U[i] == 0 || (du[cs.node_U[i]] += Q)
    return nothing
end

"Move charge `J` into node `i` (target `t`) at its force `μ`; a reservoir also books the energy μJ."
@inline function _charge!(du, cs::CompiledSystem, i::Int, t::Int, J, μ)
    @inbounds du[t] += J
    @inbounds cs.node_U[i] == 0 || (du[cs.node_U[i]] += μ * J)
    return nothing
end

"Specific entropy and mass fractions of the fluid at node `i` (the upstream state)."
@inline function _upstream(cs::CompiledSystem, u, i::Int, c::NTuple{N,Int}, yi::NTuple{N,Int}) where {N}
    Tel = promote_type(eltype(u), Float64)
    if i <= nbodies(cs)
        m = sum(ntuple(k -> u[c[k]], Val(N)))
        return (convert(Tel, u[cs.node_S[i]] / m), ntuple(k -> convert(Tel, u[c[k]] / m), Val(N)))
    else
        r = cs.reservoirs[i-nbodies(cs)]
        return (convert(Tel, r.s), ntuple(k -> convert(Tel, r.y[yi[k]]), Val(N)))
    end
end

"""
    transfer(c, P, u, cs, t) -> NamedTuple

Evaluate one element at force buffer `P` and state `u`: fluxes, heat deliveries and
production. Pure; used by both the assembler and the accessors.
"""
function transfer(c::CompiledConductor{<:ConductorLaw,ThermalRow,0}, P, u, cs::CompiledSystem, t)
    Ta, Tb = _node_T(P, cs, c.a), _node_T(P, cs, c.b)
    tr = thermal_transfer(c.law, Ta, Tb)
    z = zero(tr.J_Q)
    return (; Q = tr.J_Q, J = SVector{0,typeof(z)}(), D = z, Q_a = tr.Q_a, Q_b = tr.Q_b, Q_h = z,
            T_a = Ta, T_b = Tb, T_h = oftype(Ta, NaN), T_m = oftype(Ta, NaN), S_gen = tr.S_gen,
            parts = (; total = tr.S_gen))
end

function transfer(c::CompiledConductor{<:ConductorLaw,K,N}, P, u, cs::CompiledSystem, t) where {K,N}
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

"""
Stream: upwinded bundle. With J_s = ṁ s_up and D = Σ (μ_a − μ_b) J, the upstream node loses
entropy J_s at its own T, so Q_up = −T_up J_s, and the downstream node takes the rest:
Q_dn = D − Q_up. Production [D + J_s (T_up − T_dn)]/T_dn.
"""
function transfer(c::CompiledStream{D,N}, P, u, cs::CompiledSystem, t) where {D,N}
    mdot = mass_flow(c.driver, -P[c.pa], -P[c.pb])
    fwd = mdot >= 0
    s_up, y = fwd ? _upstream(cs, u, c.a, c.ca, c.ya) : _upstream(cs, u, c.b, c.cb, c.yb)
    J = SVector{N}(ntuple(k -> mdot * y[k], Val(N)))
    Js = mdot * s_up
    Ta, Tb = _node_T(P, cs, c.a), _node_T(P, cs, c.b)
    Dp = sum(ntuple(k -> (P[c.fa[k]] - P[c.fb[k]]) * J[k], Val(N)))
    if fwd
        Q_a = -Ta * Js
        Q_b = Dp - Q_a
        S_gen = (Dp + Js * (Ta - Tb)) / Tb
    else
        Q_b = Tb * Js
        Q_a = Dp - Q_b
        S_gen = (Dp + Js * (Ta - Tb)) / Ta
    end
    return (; mdot, J, Js, D = Dp, Q_a, Q_b, T_a = Ta, T_b = Tb, S_gen, parts = (; total = S_gen))
end

"""
Ideal pump a → b: the fluid leaves a at its own state and is delivered at (s_a, p_b).
External work W = E* − E_a, where E = Σ μ J + T J_s at each state.
"""
function transfer(c::CompiledPump{F,M,N}, P, u, cs::CompiledSystem, t) where {F,M,N}
    mdot = _nonneg_flow(c.mdot(t), "pump")
    s_a, y = _upstream(cs, u, c.a, c.ca, c.ya)
    J = SVector{N}(ntuple(k -> mdot * y[k], Val(N)))
    Js = mdot * s_a
    Ta, Tb = _node_T(P, cs, c.a), _node_T(P, cs, c.b)
    st = isentropic_state(c.mat, s_a, -P[c.pb])
    E_a = sum(ntuple(k -> P[c.fa[k]] * J[k], Val(N))) + Ta * Js
    E_star = st.μ * mdot + st.T * Js
    W = E_star - E_a
    Q_a = -Ta * Js
    Q_b = E_star - sum(ntuple(k -> P[c.fb[k]] * J[k], Val(N)))
    S_gen = (sum(ntuple(k -> (st.μ - P[c.fb[k]]) * J[k], Val(N))) + Js * (st.T - Tb)) / Tb
    return (; mdot, J, Js, W, Q_a, Q_b, T_a = Ta, T_b = Tb, T_star = st.T, S_gen, parts = (; total = S_gen))
end

"""
Flow boundary into node d. Passive: delivered at (T_in, p_in); p_d ≤ p_in is required of every
accepted state (`passive_violations`), not of solver trial points.
Pumped: delivered isentropically at (s_in, p_d), with the work booked.
"""
function transfer(c::CompiledBoundary{F,M,N}, P, u, cs::CompiledSystem, t) where {F,M,N}
    mdot = _nonneg_flow(c.mdot(t), "flow boundary")
    Tel = promote_type(eltype(P), Float64)
    p_d = -P[c.pd]
    if c.pumped
        st = isentropic_state(c.mat, c.s_in, p_d)
        T_del, μ_del = convert(Tel, st.T), ntuple(k -> convert(Tel, st.μ), Val(N))
    else                                   # admissibility (p_d ≤ p_in) is checked on accepted states
        T_del, μ_del = convert(Tel, c.T_in), ntuple(k -> convert(Tel, c.μ_in[k]), Val(N))
    end
    J = SVector{N}(ntuple(k -> mdot * c.y[k], Val(N)))
    Js = mdot * c.s_in
    Td = _node_T(P, cs, c.d)
    E_in = sum(ntuple(k -> c.μ_in[k] * J[k], Val(N))) + c.T_in * Js
    E_del = sum(ntuple(k -> μ_del[k] * J[k], Val(N))) + T_del * Js
    Q_d = E_del - sum(ntuple(k -> P[c.fd[k]] * J[k], Val(N)))
    S_gen = (sum(ntuple(k -> (μ_del[k] - P[c.fd[k]]) * J[k], Val(N))) + Js * (T_del - Td)) / Td
    return (; mdot, J, Js, W = E_del - E_in, E_in, Q_d, T_d = Td, S_gen, parts = (; total = S_gen))
end

"Pass-through outlet: the node's own bundle leaves the model, with no generation."
function transfer(c::CompiledOutlet{O,N}, P, u, cs::CompiledSystem, t) where {O,N}
    mdot = outlet_flow(c.outlet, -P[c.pu], t)
    s_u, y = _upstream(cs, u, c.u, c.cu, ntuple(_ -> 0, Val(N)))
    J = SVector{N}(ntuple(k -> mdot * y[k], Val(N)))
    Js = mdot * s_u
    Tu = _node_T(P, cs, c.u)
    E = sum(ntuple(k -> P[c.fu[k]] * J[k], Val(N))) + Tu * Js
    return (; mdot, J, Js, E, T_u = Tu)
end

function _apply!(du, P, u, c::CompiledConductor{<:ConductorLaw,ThermalRow,0}, cs, t)
    tr = transfer(c, P, u, cs, t)
    _heat!(du, cs, c.a, tr.Q_a, tr.T_a)
    _heat!(du, cs, c.b, tr.Q_b, tr.T_b)
    @inbounds du[c.led] += tr.S_gen
    return nothing
end

function _apply!(du, P, u, c::CompiledConductor{<:ConductorLaw,K,N}, cs, t) where {K,N}
    tr = transfer(c, P, u, cs, t)
    for k in 1:N
        _charge!(du, cs, c.a, c.ca[k], -tr.J[k], P[c.fa[k]])
        _charge!(du, cs, c.b, c.cb[k], tr.J[k], P[c.fb[k]])
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

function _apply!(du, P, u, c::Union{CompiledStream{D,N},CompiledPump{D,M,N}}, cs, t) where {D,M,N}
    tr = transfer(c, P, u, cs, t)
    for k in 1:N
        _charge!(du, cs, c.a, c.ca[k], -tr.J[k], P[c.fa[k]])
        _charge!(du, cs, c.b, c.cb[k], tr.J[k], P[c.fb[k]])
    end
    _heat!(du, cs, c.a, tr.Q_a, tr.T_a)
    _heat!(du, cs, c.b, tr.Q_b, tr.T_b)
    @inbounds du[c.led] += tr.S_gen
    c isa CompiledPump && (@inbounds du[c.led+1] += tr.W)
    return nothing
end

function _apply!(du, P, u, c::CompiledBoundary{F,M,N}, cs, t) where {F,M,N}
    tr = transfer(c, P, u, cs, t)
    @inbounds for k in 1:N
        du[c.cd[k]] += tr.J[k]
        du[c.led+1+k] -= tr.J[k]                 # boundary receives −ṁ_k
    end
    _heat!(du, cs, c.d, tr.Q_d, tr.T_d)
    @inbounds du[c.led] += tr.S_gen
    @inbounds du[c.led+1] += tr.W
    @inbounds du[c.led+2+N] -= tr.Js              # entropy received (negative: supplied)
    @inbounds du[c.led+3+N] -= tr.E_in            # energy received
    return nothing
end

function _apply!(du, P, u, c::CompiledOutlet{O,N}, cs, t) where {O,N}
    tr = transfer(c, P, u, cs, t)
    @inbounds for k in 1:N
        du[c.cu[k]] -= tr.J[k]
        du[c.led+k-1] += tr.J[k]
    end
    @inbounds du[cs.node_S[c.u]] -= tr.Js
    @inbounds du[c.led+N] += tr.Js
    @inbounds du[c.led+N+1] += tr.E
    return nothing
end

function _apply_group!(du, P, u, group::Vector, cs, t)
    for c in group
        _apply!(du, P, u, c, cs, t)
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
    for grp in cs.element_groups
        _apply_group!(du, P, u, grp, cs, t)
    end
    return nothing
end

"In-place ODE right-hand side; the parameter is the compiled system."
ode_rhs!(du, u, cs::CompiledSystem, t) = assemble!(du, u, cs, t)

"""
    passive_violations(cs, u) -> Vector

Passive flow boundaries whose node pressure exceeds their supply pressure at state `u`, as
`(name, p_in, p_node)`. Such a boundary would need work to deliver (label it `:pumped`).
"""
function passive_violations(cs::CompiledSystem, u::AbstractVector)
    out = Tuple{Symbol,Float64,Float64}[]
    P = nothing
    for grp in cs.element_groups, c in grp
        c isa CompiledBoundary && !c.pumped || continue
        P === nothing && (P = copy(force_buffer(cs, u)))
        p_d = -P[c.pd]
        p_d > c.p_in * (1 + 1e-12) && push!(out, (c.name, c.p_in, p_d))
    end
    return out
end

function _throw_passive(v, t)
    name, p_in, p_d = first(v)
    throw(ArgumentError("passive flow boundary :$name delivers at p = $p_in Pa into a node at p = $p_d Pa" *
                        (t === nothing ? "" : " at t = $t") * "; label it supply = :pumped, or raise its pressure"))
end

"Whether `u` is admissible: valid material states and no passive boundary pushing uphill."
admissible(cs::CompiledSystem, u::AbstractVector) = all_valid(cs, u) && isempty(passive_violations(cs, u))

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

Structural Jacobian pattern from the topology: each element couples all charges of the
bodies on its ports, and its ledger rows (and those of the reservoirs it touches) depend on
them; ledger columns are empty.
"""
function jacobian_sparsity(cs::CompiledSystem)
    n = cs.nphys + cs.nledger
    I, J = collect(1:n), collect(1:n)
    nb = nbodies(cs)
    body_cols(i) = (off = cs.body_offset[i]; off:off+ncharges(cs.bodies[i].material)-1)
    for (nodes, led) in cs.touch
        cols = reduce(vcat, [collect(body_cols(i)) for i in nodes if i <= nb]; init = Int[])
        rows = vcat(cols, collect(led))
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
