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

# --- error norm and tolerances ------------------------------------------------------------

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
    abstol_vector(cs; atol_T = 1e-8, atol_ledger = 1e-12)

Absolute tolerance per state. A body's entropy is resolved to the entropy equivalent of
`atol_T` kelvin, using its capacity dS/dT = 1/A_SS at the initial state. Entropy has an
arbitrary zero, so its magnitude is not a meaningful scale.
"""
function abstol_vector(cs::CompiledSystem; atol_T::Real = 1e-8, atol_ledger::Real = 1e-12)
    tol = fill(Float64(atol_ledger), cs.nphys + cs.nledger)
    for (i, b) in enumerate(cs.bodies)
        off = cs.body_offset[i]
        N = ncharges(b.material)
        A = eos_matrix(b.material, cs.u0[off:off+N-1])
        for k in 1:N
            tol[off+k-1] = atol_T / A[k, k]
        end
    end
    return tol
end

# --- transient ----------------------------------------------------------------------------

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
    simulate(cs, tspan; alg = Rodas5P(), reltol = 1e-8, atol_T = 1e-8, sparse, kwargs...)

Integrate with the ledger states excluded from the error norm. Other keywords (for example
`saveat`) go to `solve`.
"""
function simulate(cs::CompiledSystem, tspan; alg = Rodas5P(), reltol::Real = 1e-8, atol_T::Real = 1e-8,
                  u0 = initial_state(cs), sparse::Bool = cs.nphys > 50, kwargs...)
    prob = SciMLBase.ODEProblem(cs, tspan; u0, sparse)
    sol = SciMLBase.solve(prob, alg; reltol, abstol = abstol_vector(cs; atol_T),
                          internalnorm = PhysNorm(cs.nphys), kwargs...)
    SciMLBase.successful_retcode(sol) || error("integration failed: retcode $(sol.retcode)")
    return sol
end

# --- steady state -------------------------------------------------------------------------

"""
    steady_state(cs; u0 = initial_state(cs)) -> u

Solve f_phys(u_phys) = 0 by Newton (ledgers excluded, D-10); the returned ledger block is
zero. Report steady losses with `entropy_production(cs, u)`. Every body must be connected
to a reservoir: isolated components have one singular direction per conserved total, and
the constraint-row treatment (DESIGN §12.2) is not implemented in slice 1.
"""
function steady_state(cs::CompiledSystem; u0 = initial_state(cs))
    _check_anchored(cs)
    cs.nphys == 0 && return copy(u0)
    function resid!(r, x, p)
        dled = similar(x, cs.nledger)
        assemble!(r, dled, x, p, 0.0)
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

# --- accessors ----------------------------------------------------------------------------

function _body(cs::CompiledSystem, name::Symbol)
    i = get(cs.node_index, name, 0)
    (1 <= i <= nbodies(cs)) || throw(ArgumentError("no body named :$name"))
    return i, cs.bodies[i]
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

"Temperature of a body or reservoir [K]."
function temperature(cs::CompiledSystem, u::AbstractVector, name::Symbol)
    i = get(cs.node_index, name, 0)
    i == 0 && throw(ArgumentError("no node named :$name"))
    i > nbodies(cs) && return cs.T_res[i-nbodies(cs)]
    return temperature(cs.bodies[i].material, charges(cs, u, name))
end

"Set a body's temperature in state vector `u` (thermal-only bodies)."
function set_temperature!(cs::CompiledSystem, u::AbstractVector, name::Symbol, T::Real)
    i, b = _body(cs, name)
    ncharges(b.material) == 1 || throw(ArgumentError("set_temperature! needs a thermal-only body"))
    u[cs.body_offset[i]] = charges_at(b.material; T)[1]
    return u
end

"""
    flux(cs, u, conductor) -> NamedTuple

`Q` [W] from a to b; end temperatures `T_a`, `T_b`; entropy leaving a (`S_out_of_a = Q/T_a`)
and arriving at b (`S_into_b = Q/T_b`) [W/K]; production `S_gen`; end labels `a`, `b`.
"""
function flux(cs::CompiledSystem, u::AbstractVector, name::Symbol)
    k = get(cs.conductor_index, name, 0)
    k == 0 && throw(ArgumentError("no conductor named :$name"))
    gi, pos = cs.conductor_ref[k]
    c = cs.conductor_groups[gi][pos]
    a, b = port_node(cs.sys, name, :a), port_node(cs.sys, name, :b)
    Ta, Tb = temperature(cs, u, a), temperature(cs, u, b)
    tr = thermal_transfer(c.law, Ta, Tb)
    return (; Q = tr.Q_b, T_a = Ta, T_b = Tb, S_out_of_a = -tr.Q_a / Ta, S_into_b = tr.Q_b / Tb,
            S_gen = tr.S_gen, a, b)
end

"Heat flow [W] through a conductor, positive from port a to port b."
heat_flow(cs::CompiledSystem, u::AbstractVector, name::Symbol) = flux(cs, u, name).Q

"Entropy production rate [W/K] of every conductor at state `u`."
function entropy_production(cs::CompiledSystem, u::AbstractVector)
    du = rates(cs, u)
    return Dict(n => du[cs.nphys+k] for (k, n) in enumerate(cs.conductor_names))
end

"Integrated entropy production [J/K] of every conductor, read from the ledger block of `u`."
function generated_entropy(cs::CompiledSystem, u::AbstractVector)
    return Dict(n => u[cs.nphys+k] for (k, n) in enumerate(cs.conductor_names))
end

"""
    reservoir_ledger(cs, v, reservoir) -> (; S, U)

Entropy and energy received by a reservoir, read from the ledger block of `v`. Integrated
values from a state, rates from a rate vector.
"""
function reservoir_ledger(cs::CompiledSystem, v::AbstractVector, name::Symbol)
    i = get(cs.node_index, name, 0)
    i > nbodies(cs) || throw(ArgumentError("no reservoir named :$name"))
    j = cs.res_ledger[i-nbodies(cs)]
    return (; S = v[cs.nphys+j], U = v[cs.nphys+j+1])
end

"Total internal energy of the bodies [J]."
function total_energy(cs::CompiledSystem, u::AbstractVector)
    return sum((internal_energy(b.material, charges(cs, u, b.name)) for b in cs.bodies); init = 0.0)
end

"Total entropy of the bodies [J/K] (thermal charges; arbitrary zero per material)."
function total_entropy(cs::CompiledSystem, u::AbstractVector)
    s = 0.0
    for (i, b) in enumerate(cs.bodies)
        cs.body_tslot[i] == 0 || (s += u[cs.body_tslot[i]])
    end
    return s
end

"""
    energy_residual(cs, u) -> (; residual, scale)

Instantaneous energy balance: Σ_bodies P·dE/dt plus the energy rate into reservoirs. It is
zero to round-off for every state (each connection closes exactly). `scale` is the total
heat throughput Σ|Q̇| over conductors.
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
    scale = sum((abs(flux(cs, u, n).Q) for n in cs.conductor_names); init = 0.0)
    return (; residual = r, scale)
end

"""
    closure(cs, sol) -> NamedTuple

Closure residuals over a transient (DESIGN §10), reported for every run. They are measured
from the problem's initial state `sol.prob.u0`, whatever `saveat` kept, and ledgers count
from their initial values.
- energy: max |ΔU_bodies(t) + ΔU_reservoirs(t)|, which sits at integrator tolerance.
  `energy_rel` divides by the energy moved, max Σ|ΔU_i| + Σ|ΔU_res|.
- entropy: max |ΔS_bodies(t) + ΔS_reservoirs(t) − ΔS_gen(t)|, which is a linear invariant and
  so sits at round-off. `entropy_rel` divides by the entropy moved plus produced.
"""
function closure(cs::CompiledSystem, sol)
    u0 = sol.prob.u0
    U0 = [internal_energy(b.material, charges(cs, u0, b.name)) for b in cs.bodies]
    S0 = total_entropy(cs, u0)
    led0 = [reservoir_ledger(cs, u0, n) for n in cs.reservoir_names]
    Sgen0 = sum(values(generated_entropy(cs, u0)); init = 0.0)
    eE = 0.0; eS = 0.0; movedE = 0.0; movedS = 0.0
    for u in sol.u
        Ub = [internal_energy(b.material, charges(cs, u, b.name)) for b in cs.bodies]
        dU_res = [reservoir_ledger(cs, u, n).U - l0.U for (n, l0) in zip(cs.reservoir_names, led0)]
        dS_res = [reservoir_ledger(cs, u, n).S - l0.S for (n, l0) in zip(cs.reservoir_names, led0)]
        dSgen = sum(values(generated_entropy(cs, u)); init = 0.0) - Sgen0
        dS = total_entropy(cs, u) - S0
        eE = max(eE, abs(sum(Ub; init = 0.0) - sum(U0; init = 0.0) + sum(dU_res; init = 0.0)))
        eS = max(eS, abs(dS + sum(dS_res; init = 0.0) - dSgen))
        movedE = max(movedE, sum(abs.(Ub .- U0); init = 0.0) + sum(abs, dU_res; init = 0.0))
        movedS = max(movedS, abs(dS) + sum(abs, dS_res; init = 0.0) + abs(dSgen))
    end
    return (; energy = eE, energy_rel = movedE > 0 ? eE / movedE : eE,
            entropy = eS, entropy_rel = movedS > 0 ? eS / movedS : eS)
end
