# Running a compiled system and reading results (DESIGN §10, §12).

"Initial state vector `[u_phys; q_ledger]` (ledgers start at zero)."
initial_state(cs::CompiledSystem) = copy(cs.u0)

"""
    rates(cs, u; t = 0.0) -> du

The full right-hand side (physical and ledger rates) at state `u`.
"""
function rates(cs::CompiledSystem, u::AbstractVector; t::Real = 0.0)
    du = similar(u, promote_type(eltype(u), Float64))
    ode_rhs!(du, u, cs, t)
    return du
end

# --- error norm and tolerances ----------------------------------------------------------------

"""
Integrator error norm over the physical states only (ledgers are excluded, D-10). A
system without bodies has no physical states; its norm falls back to the full vector.
"""
struct PhysNorm <: Function   # a Function, so broadcasting treats it as a scalar callable
    n::Int
end
function (pn::PhysNorm)(u::AbstractArray, t)
    v = pn.n == 0 ? u : view(u, 1:pn.n)
    isempty(v) && return zero(real(eltype(u)))
    return sqrt(sum(abs2, v) / length(v))
end
(pn::PhysNorm)(u::Number, t) = abs(u)

"""
    abstol_vector(cs; atol_scale = 1.0, atol_ledger = 1e-12)

Absolute tolerance per state. Each charge is resolved to the equivalent of its domain's
force tolerance (`default_force_atol`, times `atol_scale`) through its capacity 1/A_kk at
the initial state. Charge magnitudes are not used: entropy has an arbitrary zero (D-20).
"""
function abstol_vector(cs::CompiledSystem; atol_scale::Real = 1.0, atol_ledger::Real = 1e-12)
    tol = fill(Float64(atol_ledger), cs.nphys + cs.nledger)
    for (i, b) in enumerate(cs.bodies)
        off = cs.body_offset[i]
        N = ncharges(b.material)
        E = cs.u0[off:off+N-1]
        A = eos_matrix(b.material, E)
        for (k, d) in enumerate(domains(b.material))
            t = atol_scale * default_force_atol(d) / abs(A[k, k])
            tol[off+k-1] = (isfinite(t) && t > 0) ? t : 1e-12 * max(abs(E[k]), 1.0)
        end
    end
    return tol
end

# --- transient ----------------------------------------------------------------------------------

"""
    ODEProblem(cs, tspan; u0 = initial_state(cs), sparse = cs.nphys > 50)

An `ODEProblem` for the compiled system. With `sparse = true` the Jacobian prototype comes
from the graph topology.
"""
function SciMLBase.ODEProblem(cs::CompiledSystem, tspan; u0 = initial_state(cs),
                              sparse::Bool = cs.nphys > 50, kwargs...)
    f = SciMLBase.ODEFunction{true}(ode_rhs!; jac_prototype = sparse ? jacobian_sparsity(cs) : nothing)
    return SciMLBase.ODEProblem(f, u0, tspan, cs; kwargs...)
end

"""
    simulate(cs, tspan; alg = Rodas5P(), reltol = 1e-8, atol_scale = 1.0, sparse, kwargs...)

Integrate with the ledger states excluded from the error norm. Steps that would leave a
material's valid states (e.g. a gas with V ≤ 0) are rejected. Other keywords (for example
`saveat`) go to `solve`.
"""
function simulate(cs::CompiledSystem, tspan; alg = Rodas5P(), reltol::Real = 1e-8, atol_scale::Real = 1.0,
                  u0 = initial_state(cs), sparse::Bool = cs.nphys > 50, kwargs...)
    prob = SciMLBase.ODEProblem(cs, tspan; u0, sparse)
    sol = SciMLBase.solve(prob, alg; reltol, abstol = abstol_vector(cs; atol_scale),
                          internalnorm = PhysNorm(cs.nphys),
                          isoutofdomain = (u, p, t) -> !all_valid(p, u), kwargs...)
    SciMLBase.successful_retcode(sol) || error("integration failed: retcode $(sol.retcode)")
    return sol
end

# --- steady state -------------------------------------------------------------------------------

"""
    steady_state(cs; u0 = initial_state(cs)) -> u

Solve f_phys(u_phys) = 0 by Newton (ledgers excluded, D-10); the returned ledger block is
zero. Report steady losses with `entropy_production(cs, u)`. Every body must be connected
to a reservoir: isolated components have one singular direction per conserved total, and
the constraint-row treatment (DESIGN §12.2) is not implemented yet.
"""
function steady_state(cs::CompiledSystem; u0 = initial_state(cs))
    _check_anchored(cs)
    cs.nphys == 0 && return copy(u0)
    function resid!(r, x, p)
        ufull = vcat(x, zeros(eltype(x), p.nledger))
        du = similar(ufull)
        assemble!(du, ufull, p, 0.0)
        r .= view(du, 1:p.nphys)
        return nothing
    end
    prob = NonlinearProblem(resid!, u0[1:cs.nphys], cs)
    sol = SciMLBase.solve(prob, NewtonRaphson())
    SciMLBase.successful_retcode(sol) || error("steady-state solve failed: retcode $(sol.retcode)")
    return vcat(sol.u, zeros(cs.nledger))
end

function _check_anchored(cs::CompiledSystem)
    g = cs.sys.graph
    for comp in connected_components(g.graph)
        labels = [label_for(g, v) for v in comp]
        any(l -> g[l] isa Body, labels) || continue
        any(l -> g[l] isa Reservoir, labels) ||
            throw(ArgumentError("steady_state: bodies $(filter(l -> g[l] isa Body, labels)) are not connected " *
                                "to any reservoir (isolated components need constraint rows, not yet implemented)"))
    end
    return nothing
end

# --- accessors ------------------------------------------------------------------------------------

function _body(cs::CompiledSystem, name::Symbol)
    i = get(cs.node_index, name, 0)
    (1 <= i <= nbodies(cs)) || throw(ArgumentError("no body named :$name"))
    return i, cs.bodies[i]
end

function _reservoir(cs::CompiledSystem, name::Symbol)
    i = get(cs.node_index, name, 0)
    i > nbodies(cs) || throw(ArgumentError("no reservoir named :$name"))
    return i - nbodies(cs), cs.reservoirs[i-nbodies(cs)]
end

"""
    charges(cs, v, body) -> SVector

The block of `v` (a state or a rate vector) that belongs to `body`.
"""
function charges(cs::CompiledSystem, v::AbstractVector, name::Symbol)
    i, b = _body(cs, name)
    N = ncharges(b.material)
    off = cs.body_offset[i]
    return SVector{N}(v[off:off+N-1])
end

"""
    force(cs, u, node, domain)

The force of `domain` at a body or reservoir: T, −p, μ, …
"""
function force(cs::CompiledSystem, u::AbstractVector, name::Symbol, d::FormOfMotion)
    i = get(cs.node_index, name, 0)
    i == 0 && throw(ArgumentError("no node named :$name"))
    i > nbodies(cs) && return reservoir_force(cs.reservoirs[i-nbodies(cs)], d)
    b = cs.bodies[i]
    k = findfirst(==(d), collect(domains(b.material)))
    k === nothing && throw(ArgumentError(":$name has no $d form of motion"))
    return forces(b.material, charges(cs, u, name))[k]
end

"Temperature of a body or reservoir [K]."
temperature(cs::CompiledSystem, u::AbstractVector, name::Symbol) = force(cs, u, name, Thermal())

"Pressure of a body or reservoir [Pa] (p = −P_V)."
pressure(cs::CompiledSystem, u::AbstractVector, name::Symbol) = -force(cs, u, name, Volumetric())

"Set a body's temperature in state vector `u` (thermal-only bodies)."
function set_temperature!(cs::CompiledSystem, u::AbstractVector, name::Symbol, T::Real)
    i, b = _body(cs, name)
    ncharges(b.material) == 1 || throw(ArgumentError("set_temperature! needs a thermal-only body"))
    u[cs.body_offset[i]] = charges_at(b.material; T)[1]
    return u
end

function _conductor(cs::CompiledSystem, name::Symbol)
    k = get(cs.conductor_index, name, 0)
    k == 0 && throw(ArgumentError("no conductor named :$name"))
    gi, pos = cs.conductor_ref[k]
    return k, cs.conductor_groups[gi][pos]
end

"""
    flux(cs, u, conductor) -> NamedTuple

One conductor's transfer at state `u`:
- `Q`: thermal-row heat flow from a to b [W] (0 for non-thermal conductors);
- `J`: non-thermal charge fluxes from a to b;
- `D`: released power [W];
- `Q_a`, `Q_b`, `Q_h`: heat delivered to each port [W];
- `T_a`, `T_b`, `T_h`, `T_m`: port temperatures and the internal-node temperature [K]
  (NaN where absent);
- `S_out_of_a = −Q_a/T_a`, `S_into_b = Q_b/T_b` [W/K];
- `S_gen` and its `parts`;
- end labels `a`, `b`.
"""
function flux(cs::CompiledSystem, u::AbstractVector, name::Symbol)
    _, c = _conductor(cs, name)
    P = copy(force_buffer(cs, u))
    tr = transfer(c, P, cs)
    a, b = port_node(cs.sys, name, :a), port_node(cs.sys, name, :b)
    return merge(tr, (; S_out_of_a = -tr.Q_a / tr.T_a, S_into_b = tr.Q_b / tr.T_b, a, b))
end

"Heat flow [W] through a conductor's thermal row, positive from port a to port b."
heat_flow(cs::CompiledSystem, u::AbstractVector, name::Symbol) = flux(cs, u, name).Q

_ledger_total(cs::CompiledSystem, v, k::Int) = sum(v[j] for j in cs.conductor_led[k])

"Entropy production rate [W/K] of every conductor at state `u`."
function entropy_production(cs::CompiledSystem, u::AbstractVector)
    du = rates(cs, u)
    return Dict(n => _ledger_total(cs, du, k) for (k, n) in enumerate(cs.conductor_names))
end

"""
    entropy_production_parts(cs, u)

Per conductor, the production rate split as the ledger records it: `(; total)` for kernel A,
`(; dissipation, conduction)` for an internal node (kernel B).
"""
function entropy_production_parts(cs::CompiledSystem, u::AbstractVector)
    P = copy(force_buffer(cs, u))
    return Dict(n => transfer(_conductor(cs, n)[2], P, cs).parts for n in cs.conductor_names)
end

"Integrated entropy production [J/K] of every conductor, read from the ledger block of `u`."
function generated_entropy(cs::CompiledSystem, u::AbstractVector)
    return Dict(n => _ledger_total(cs, u, k) for (k, n) in enumerate(cs.conductor_names))
end

"""
    reservoir_ledger(cs, v, reservoir) -> (; S, U)

Entropy and energy received by a reservoir (`S = 0` without a thermal DOF), read from the
ledger block of `v`. Integrated values from a state, rates from a rate vector.
"""
function reservoir_ledger(cs::CompiledSystem, v::AbstractVector, name::Symbol)
    j, r = _reservoir(cs, name)
    k = findfirst(d -> d isa Thermal, r.domains)
    S = k === nothing ? zero(eltype(v)) : v[cs.res_led[j]+k-1]
    return (; S, U = v[cs.res_led[j]+length(r.domains)])
end

"Charge of `domain` received by a reservoir (from the ledger block of `v`)."
function reservoir_charge(cs::CompiledSystem, v::AbstractVector, name::Symbol, d::FormOfMotion)
    j, r = _reservoir(cs, name)
    k = findfirst(==(d), r.domains)
    k === nothing && throw(ArgumentError("reservoir :$name has no $d form of motion"))
    return v[cs.res_led[j]+k-1]
end

"""
    total_charge(cs, u, domain)

Σ over bodies of their `domain` charge, plus the charge received by reservoirs. For every
non-thermal domain this is conserved (08 A2).
"""
function total_charge(cs::CompiledSystem, u::AbstractVector, d::FormOfMotion)
    s = 0.0
    for (i, b) in enumerate(cs.bodies)
        k = findfirst(==(d), collect(domains(b.material)))
        k === nothing || (s += u[cs.body_offset[i]+k-1])
    end
    for (j, r) in enumerate(cs.reservoirs)
        k = findfirst(==(d), r.domains)
        k === nothing || (s += u[cs.res_led[j]+k-1])
    end
    return s
end

"Total internal energy of the bodies [J]."
function total_energy(cs::CompiledSystem, u::AbstractVector)
    return sum((internal_energy(b.material, charges(cs, u, b.name)) for b in cs.bodies); init = 0.0)
end

"Total entropy of the bodies [J/K] (thermal charges; arbitrary zero per material)."
function total_entropy(cs::CompiledSystem, u::AbstractVector)
    s = 0.0
    for (i, b) in enumerate(cs.bodies)
        k = thermal_slot(b.material)
        k === nothing || (s += u[cs.body_offset[i]+k-1])
    end
    return s
end

"""
    energy_residual(cs, u) -> (; residual, scale)

Instantaneous energy balance: Σ_bodies P·dE/dt plus the energy rate into reservoirs. It is
zero to round-off at every state (each connection closes exactly). `scale` is the total
energy throughput of the conductors.
"""
function energy_residual(cs::CompiledSystem, u::AbstractVector)
    du = rates(cs, u)
    r = 0.0
    for b in cs.bodies
        r += dot(forces(b.material, charges(cs, u, b.name)), charges(cs, du, b.name))
    end
    for n in cs.reservoir_names
        r += reservoir_ledger(cs, du, n).U
    end
    P = copy(force_buffer(cs, u))
    scale = 0.0
    for n in cs.conductor_names
        _, c = _conductor(cs, n)
        tr = transfer(c, P, cs)
        scale += abs(tr.Q_a) + abs(tr.Q_b) + abs(tr.Q_h)
        for k in eachindex(tr.J)
            scale += abs(tr.J[k]) * (abs(P[c.fa[k]]) + abs(P[c.fb[k]]))
        end
    end
    return (; residual = r, scale)
end

"Non-thermal domains present in the system."
function _nonthermal_domains(cs::CompiledSystem)
    ds = FormOfMotion[]
    for b in cs.bodies, d in domains(b.material)
        d isa Thermal || d in ds || push!(ds, d)
    end
    for r in cs.reservoirs, d in r.domains
        d isa Thermal || d in ds || push!(ds, d)
    end
    return ds
end

"""
    closure(cs, sol) -> NamedTuple

Closure residuals over a transient (DESIGN §10), reported for every run. They are measured
from the problem's initial state `sol.prob.u0`, whatever `saveat` kept, and ledgers count
from their initial values.
- energy: max |ΔU_bodies(t) + ΔU_reservoirs(t)|, which sits at integrator tolerance.
  `energy_rel` divides by the energy moved.
- entropy: max |ΔS_bodies(t) + ΔS_reservoirs(t) − ΔS_gen(t)|, which sits at round-off.
  `entropy_rel` divides by the entropy moved plus produced.
- charge: per non-thermal domain, max |Δ total_charge|, which sits at round-off.
  `charge_rel` is the largest relative to Σ|charge| of the bodies.
"""
function closure(cs::CompiledSystem, sol)
    u0 = sol.prob.u0
    U0 = [internal_energy(b.material, charges(cs, u0, b.name)) for b in cs.bodies]
    S0 = total_entropy(cs, u0)
    led0 = [reservoir_ledger(cs, u0, n) for n in cs.reservoir_names]
    Sgen0 = sum(values(generated_entropy(cs, u0)); init = 0.0)
    doms = _nonthermal_domains(cs)
    Q0 = [total_charge(cs, u0, d) for d in doms]
    Qscale = [max(sum((abs(u0[cs.body_offset[i]+k-1]) for (i, b) in enumerate(cs.bodies)
                       for (k, dd) in enumerate(domains(b.material)) if dd == d); init = 0.0), floatmin()) for d in doms]
    eE = 0.0; eS = 0.0; movedE = 0.0; movedS = 0.0; eQ = 0.0
    for u in sol.u
        Ub = [internal_energy(b.material, charges(cs, u, b.name)) for b in cs.bodies]
        led = [reservoir_ledger(cs, u, n) for n in cs.reservoir_names]
        dU_res = [l.U - l0.U for (l, l0) in zip(led, led0)]
        dS_res = [l.S - l0.S for (l, l0) in zip(led, led0)]
        dSgen = sum(values(generated_entropy(cs, u)); init = 0.0) - Sgen0
        dS = total_entropy(cs, u) - S0
        eE = max(eE, abs(sum(Ub; init = 0.0) - sum(U0; init = 0.0) + sum(dU_res; init = 0.0)))
        eS = max(eS, abs(dS + sum(dS_res; init = 0.0) - dSgen))
        movedE = max(movedE, sum(abs.(Ub .- U0); init = 0.0) + sum(abs, dU_res; init = 0.0))
        movedS = max(movedS, abs(dS) + sum(abs, dS_res; init = 0.0) + abs(dSgen))
        for (q, d) in enumerate(doms)
            eQ = max(eQ, abs(total_charge(cs, u, d) - Q0[q]) / Qscale[q])
        end
    end
    return (; energy = eE, energy_rel = movedE > 0 ? eE / movedE : eE,
            entropy = eS, entropy_rel = movedS > 0 ? eS / movedS : eS,
            charge_rel = eQ)
end
