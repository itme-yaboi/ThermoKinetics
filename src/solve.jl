# Running a compiled system and reading results (DESIGN §10, §12).

"Initial state vector `[u_phys; q_ledger]` (ledgers start at zero)."
initial_state(cs::CompiledSystem) = copy(cs.u0)

"""
    rates(cs, u; t = 0.0) -> du

The full right-hand side (physical and ledger rates) at state `u`.
"""
function rates(cs::CompiledSystem, u::AbstractVector; t::Real = 0.0)
    v = passive_violations(cs, u)
    isempty(v) || _throw_passive(v, t)
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
material's valid states (e.g. a gas with V ≤ 0), or push a passive boundary uphill, are
rejected. Other keywords (for example
`saveat`) go to `solve`.
"""
function simulate(cs::CompiledSystem, tspan; alg = Rodas5P(), reltol::Real = 1e-8, atol_scale::Real = 1.0,
                  u0 = initial_state(cs), sparse::Bool = cs.nphys > 50, kwargs...)
    prob = SciMLBase.ODEProblem(cs, tspan; u0, sparse)
    sol = SciMLBase.solve(prob, alg; reltol, abstol = abstol_vector(cs; atol_scale),
                          internalnorm = PhysNorm(cs.nphys),
                          isoutofdomain = (u, p, t) -> !admissible(p, u), kwargs...)
    if !SciMLBase.successful_retcode(sol)
        v = passive_violations(cs, sol.u[end])
        isempty(v) || _throw_passive(v, sol.t[end])
        error("integration failed: retcode $(sol.retcode)")
    end
    return sol
end

# --- steady state -------------------------------------------------------------------------------

"""
    steady_state(cs; u0 = initial_state(cs)) -> u

Solve f_phys = 0 by Newton (ledgers excluded, D-10); the returned ledger block is zero.
Report steady losses with `entropy_production(cs, u)`.

Only the *active* charges are unknowns: those some element writes. A rigid vessel's volume,
for instance, is held at its `u0` value. The system must be anchored:
- energy: every connected part has a reservoir, a flow boundary or an outlet;
- each conserved charge: every group of bodies exchanging it reaches a reservoir holding
  it, or a boundary or outlet carrying it.

Isolated groups have one singular direction per conserved total. The constraint-row
treatment (DESIGN §12.2) is not implemented yet; integrate those with `simulate` instead.

Newton (then trust region) needs a reasonable start. For flow networks, warm-start from a
short transient: `steady_state(cs; u0 = simulate(cs, (0, t)).u[end])`.
"""
function steady_state(cs::CompiledSystem; u0 = initial_state(cs))
    _check_anchored(cs)
    act = cs.active
    isempty(act) && return copy(u0)
    base = collect(Float64, u0[1:cs.nphys])
    function resid!(r, x, p)
        ufull = zeros(eltype(x), p.nphys + p.nledger)
        ufull[1:p.nphys] .= base
        ufull[act] .= x
        if !admissible(p, ufull)         # e.g. a trial step to V < 0: reject it, do not throw
            fill!(r, convert(eltype(r), NaN))
            return nothing
        end
        du = similar(ufull)
        assemble!(du, ufull, p, 0.0)
        r .= view(du, act)
        return nothing
    end
    prob = NonlinearProblem(resid!, base[act], cs)
    sol = SciMLBase.solve(prob, NewtonRaphson())
    if !SciMLBase.successful_retcode(sol)            # far from the solution: globalize
        sol = SciMLBase.solve(prob, TrustRegion())
    end
    SciMLBase.successful_retcode(sol) ||
        error("steady-state solve failed: retcode $(sol.retcode). Try a transient warm start: " *
              "steady_state(cs; u0 = simulate(cs, (0, t)).u[end])")
    u = vcat(base, zeros(cs.nledger))
    u[act] .= sol.u
    v = passive_violations(cs, u)
    isempty(v) || _throw_passive(v, nothing)
    return u
end

"Union–find root."
function _root(parent, i)
    while parent[i] != i
        parent[i] = parent[parent[i]]
        i = parent[i]
    end
    return i
end

function _check_anchored(cs::CompiledSystem)
    nb, nn = nbodies(cs), nnodes(cs)
    # energy: any element connects its nodes; reservoirs, boundaries and outlets anchor
    parent = collect(1:nn)
    anchored = falses(nn)
    for j in nb+1:nn; anchored[j] = true; end
    for (nodes, _, isbdry) in cs.moves
        isbdry && foreach(i -> anchored[i] = true, nodes)
    end
    for (nodes, _) in cs.touch, i in nodes[2:end]
        parent[_root(parent, i)] = _root(parent, nodes[1])
    end
    _report_unanchored(cs, parent, anchored, "energy (no reservoir, flow boundary or outlet)")
    # each conserved (non-thermal) charge moved by some element
    doms = unique(d for (_, ds, _) in cs.moves for d in ds)
    for d in doms
        parent = collect(1:nn)
        anchored = falses(nn)
        involved = falses(nn)
        for j in 1:nn-nb
            d in cs.reservoirs[j].domains && (anchored[nb+j] = true)
        end
        for (nodes, ds, isbdry) in cs.moves
            d in ds || continue
            foreach(i -> involved[i] = true, nodes)
            isbdry && foreach(i -> anchored[i] = true, nodes)
            for i in nodes[2:end]
                parent[_root(parent, i)] = _root(parent, nodes[1])
            end
        end
        _report_unanchored(cs, parent, anchored, "$d"; only = involved)
    end
    return nothing
end

function _report_unanchored(cs, parent, anchored, what; only = trues(length(parent)))
    nb = nbodies(cs)
    ok = falses(length(parent))
    for i in eachindex(parent)
        anchored[i] && (ok[_root(parent, i)] = true)
    end
    bad = [cs.body_names[i] for i in 1:nb if only[i] && !ok[_root(parent, i)]]
    isempty(bad) ||
        throw(ArgumentError("steady_state: bodies $bad are not anchored for $what. Isolated groups conserve " *
                            "that total, and the constraint-row treatment is not implemented yet; use simulate"))
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

function _element(cs::CompiledSystem, name::Symbol)
    ref = get(cs.element_ref, name, nothing)
    ref === nothing && throw(ArgumentError("no element named :$name"))
    return cs.element_groups[ref[1]][ref[2]]
end

_eval(cs::CompiledSystem, u::AbstractVector, name::Symbol; t::Real = 0.0) =
    transfer(_element(cs, name), copy(force_buffer(cs, u)), u, cs, t)

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
    c = _element(cs, name)
    c isa CompiledConductor || throw(ArgumentError(":$name is not a conductor; use mass_flow for flow elements"))
    tr = _eval(cs, u, name)
    a, b = port_node(cs.sys, name, :a), port_node(cs.sys, name, :b)
    return merge(tr, (; S_out_of_a = -tr.Q_a / tr.T_a, S_into_b = tr.Q_b / tr.T_b, a, b))
end

"Heat flow [W] through a conductor's thermal row, positive from port a to port b."
heat_flow(cs::CompiledSystem, u::AbstractVector, name::Symbol) = flux(cs, u, name).Q

"""
    mass_flow(cs, u, element; t = 0.0)

Mass flow [kg/s] of a stream (positive from a to b), pump, flow boundary (into its node) or
outlet (out of its node).
"""
function mass_flow(cs::CompiledSystem, u::AbstractVector, name::Symbol; t::Real = 0.0)
    c = _element(cs, name)
    c isa CompiledConductor && throw(ArgumentError(":$name is a conductor; use flux"))
    return _eval(cs, u, name; t).mdot
end

_ledger_total(v, r::UnitRange{Int}) = sum(v[j] for j in r)

"Entropy production rate [W/K] of every producing element at state `u`."
function entropy_production(cs::CompiledSystem, u::AbstractVector; t::Real = 0.0)
    du = rates(cs, u; t)
    return Dict(n => _ledger_total(du, r) for (n, r) in cs.prod_led)
end

"""
    entropy_production_parts(cs, u)

Per producing element, the production rate split as the ledger records it: `(; total)`, or
`(; dissipation, conduction)` for an internal node (kernel B).
"""
function entropy_production_parts(cs::CompiledSystem, u::AbstractVector; t::Real = 0.0)
    return Dict(n => _eval(cs, u, n; t).parts for n in keys(cs.prod_led))
end

"Integrated entropy production [J/K] of every producing element, read from the ledger block of `u`."
generated_entropy(cs::CompiledSystem, u::AbstractVector) =
    Dict(n => _ledger_total(u, r) for (n, r) in cs.prod_led)

"External work [W, or J if integrated] supplied by every pump and flow boundary, read from `v`."
pump_work(cs::CompiledSystem, v::AbstractVector) = Dict(n => v[i] for (n, i) in cs.work_led)

"""
    boundary_ledger(cs, v, element) -> (; m, S, U, W)

What a flow boundary or outlet *received* (mass, entropy, energy; negative for a supplying
boundary) and the external work it took in, read from the ledger block of `v`.
"""
function boundary_ledger(cs::CompiledSystem, v::AbstractVector, name::Symbol)
    r = get(cs.bound_led, name, nothing)
    r === nothing && throw(ArgumentError(":$name is not a flow boundary or outlet"))
    W = haskey(cs.work_led, name) ? v[cs.work_led[name]] : zero(eltype(v))
    return (; m = _ledger_total(v, r.m), S = v[r.S], U = v[r.U], W)
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

Σ over bodies of their `domain` charge, plus what reservoirs, flow boundaries and outlets
received. For every non-thermal domain this is conserved (08 A2).
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
    for (dd, i) in cs.ext_charges
        dd == d && (s += u[i])
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

"Energy terms outside the bodies, from `v`: received by reservoirs, boundaries and outlets, and −(external work)."
function _external_energy_terms(cs::CompiledSystem, v::AbstractVector)
    return vcat([reservoir_ledger(cs, v, n).U for n in cs.reservoir_names],
                [v[r.U] for r in values(cs.bound_led)], [-v[i] for i in values(cs.work_led)])
end

"Entropy received outside the bodies (reservoirs, boundaries, outlets), term by term, from `v`."
function _external_entropy_terms(cs::CompiledSystem, v::AbstractVector)
    return vcat([reservoir_ledger(cs, v, n).S for n in cs.reservoir_names], [v[r.S] for r in values(cs.bound_led)])
end

_external_energy(cs::CompiledSystem, v::AbstractVector) = sum(_external_energy_terms(cs, v); init = 0.0)

"""
    energy_residual(cs, u) -> (; residual, scale)

Instantaneous energy balance: Σ_bodies P·dE/dt, plus the energy rate into reservoirs,
boundaries and outlets, minus external work. It is zero to round-off at every state (each
connection closes exactly). `scale` is the total energy throughput of the elements.
"""
function energy_residual(cs::CompiledSystem, u::AbstractVector)
    du = rates(cs, u)
    r = 0.0
    scale = 0.0
    for b in cs.bodies
        F = forces(b.material, charges(cs, u, b.name))
        dE = charges(cs, du, b.name)
        r += dot(F, dE)
        scale += sum(abs.(F .* dE))
    end
    ext = _external_energy(cs, du)
    r += ext
    scale += abs(ext)
    for n in cs.element_names
        tr = _eval(cs, u, n)
        for f in (:Q_a, :Q_b, :Q_h, :Q_d, :W, :E, :E_in)
            hasproperty(tr, f) && (scale += abs(getproperty(tr, f)))
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
- energy: max |ΔU_bodies + Δ(energy received by reservoirs, boundaries, outlets) −
  Δ(external work)|, which sits at integrator tolerance. `energy_rel` divides by the
  energy moved.
- entropy: max |ΔS_bodies + Δ(entropy received outside) − ΔS_gen|, which sits at round-off.
  `entropy_rel` divides by the entropy moved plus produced.
- charge: per non-thermal domain, max |Δ total_charge|, which sits at round-off.
  `charge_rel` is the largest, relative to the bodies' Σ|charge| (its maximum over the run).
"""
function closure(cs::CompiledSystem, sol)
    u0 = sol.prob.u0
    U0 = [internal_energy(b.material, charges(cs, u0, b.name)) for b in cs.bodies]
    S0 = total_entropy(cs, u0)
    X0, Y0 = _external_energy_terms(cs, u0), _external_entropy_terms(cs, u0)
    Sgen0 = sum(values(generated_entropy(cs, u0)); init = 0.0)
    doms = _nonthermal_domains(cs)
    Q0 = [total_charge(cs, u0, d) for d in doms]
    body_charge(u, d) = sum((abs(u[cs.body_offset[i]+k-1]) for (i, b) in enumerate(cs.bodies)
                             for (k, dd) in enumerate(domains(b.material)) if dd == d); init = 0.0)
    # scale: the largest body total over the run (a tank filled from vacuum starts near zero)
    Qscale = [max(maximum(body_charge(u, d) for u in sol.u), body_charge(u0, d), floatmin()) for d in doms]
    eE = 0.0; eS = 0.0; movedE = 0.0; movedS = 0.0; eQ = 0.0
    for u in sol.u
        Ub = [internal_energy(b.material, charges(cs, u, b.name)) for b in cs.bodies]
        dXt = _external_energy_terms(cs, u) .- X0
        dYt = _external_entropy_terms(cs, u) .- Y0
        dX, dY = sum(dXt; init = 0.0), sum(dYt; init = 0.0)
        dSgen = sum(values(generated_entropy(cs, u)); init = 0.0) - Sgen0
        dS = total_entropy(cs, u) - S0
        eE = max(eE, abs(sum(Ub; init = 0.0) - sum(U0; init = 0.0) + dX))
        eS = max(eS, abs(dS + dY - dSgen))
        movedE = max(movedE, sum(abs.(Ub .- U0); init = 0.0) + sum(abs, dXt; init = 0.0))
        movedS = max(movedS, abs(dS) + sum(abs, dYt; init = 0.0) + abs(dSgen))
        for (q, d) in enumerate(doms)
            eQ = max(eQ, abs(total_charge(cs, u, d) - Q0[q]) / Qscale[q])
        end
    end
    return (; energy = eE, energy_rel = movedE > 0 ? eE / movedE : eE,
            entropy = eS, entropy_rel = movedS > 0 ? eS / movedS : eS,
            charge_rel = eQ)
end
