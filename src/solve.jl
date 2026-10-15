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
    abstol_vector(cs; atol_scale = 1.0, atol_ledger = nothing)

Absolute tolerance per state (DESIGN §12.5, D-31).
- **Charges** are resolved to the equivalent of their domain's force tolerance
  (`default_force_atol`, times `atol_scale`) through their capacity 1/A_kk at the initial state:
  1e-8 K of temperature, 1e-4 Pa of pressure, 1e-8 V, … Charge magnitudes are not used,
  because entropy has an arbitrary zero (D-20).
- **Junction forces:** the domain's force tolerance. **Transformer fluxes:** 1e-10 × `atol_scale`.
- **Ledgers** are scaled like the states they account for. Entropy ledgers (production,
  entropy received) get the median body entropy tolerance; energy ledgers get that times the
  median body temperature; charge ledgers get the median body tolerance of that domain.
  Ledgers start at zero, so a flat tiny value would force small steps for nothing.
  `atol_ledger = x` sets all ledger tolerances to `x` instead.
"""
function abstol_vector(cs::CompiledSystem; atol_scale::Real = 1.0, atol_ledger = nothing)
    tol = fill(1e-12 * atol_scale, nstate(cs) + cs.nledger)
    for (j, jn) in enumerate(cs.junctions), (k, d) in enumerate(jn.domains)   # junction forces
        tol[cs.junction_slots[j][k][2]] = atol_scale * default_force_atol(d)
    end
    for c in cs.transformers                                                   # transformer fluxes
        tol[c.λ] = atol_scale * 1e-10
    end
    by_domain = Dict{FormOfMotion,Vector{Float64}}()
    temperatures = Float64[]
    for (i, b) in enumerate(cs.bodies)
        off = cs.body_offset[i]
        N = ncharges(b.material)
        E = cs.u0[off:off+N-1]
        A = eos_matrix(b.material, E)
        F = forces(b.material, E)
        fl = stiffness_floor(b.material)
        for (k, d) in enumerate(domains(b.material))
            t = atol_scale * default_force_atol(d) / max(abs(A[k, k]), fl[k])
            tol[off+k-1] = (isfinite(t) && t > 0) ? t : 1e-12 * max(abs(E[k]), 1.0)
            push!(get!(by_domain, d, Float64[]), tol[off+k-1])
            d isa Thermal && push!(temperatures, F[k])
        end
    end
    # ledgers: scaled like the states they account for
    med(v) = (w = sort(v); n = length(w); isodd(n) ? w[(n+1)÷2] : (w[n÷2] + w[n÷2+1]) / 2)
    if isempty(temperatures)                       # no thermal bodies: use the reservoirs' temperatures
        for r in cs.reservoirs, (d, f) in zip(r.domains, r.forces)
            d isa Thermal && push!(temperatures, f)
        end
    end
    T_ref = isempty(temperatures) ? 300.0 : med(temperatures)
    ref_S = haskey(by_domain, Thermal()) ? med(by_domain[Thermal()]) : 1e-12 * atol_scale
    for (i, kind) in enumerate(cs.ledger_kind)
        slot = nstate(cs) + i
        tol[slot] = atol_ledger !== nothing ? Float64(atol_ledger) :
                    kind === :S ? ref_S :
                    kind === :U ? ref_S * T_ref :
                    haskey(by_domain, kind) ? med(by_domain[kind]) : 1e-12 * atol_scale
    end
    return tol
end

# --- transient ----------------------------------------------------------------------------------

"""
    ODEProblem(cs, tspan; u0 = initial_state(cs), sparse = cs.nphys > 50)

An `ODEProblem` for the compiled system. With `sparse = true` the Jacobian prototype comes
from the graph topology (ODE systems only). With algebraic unknowns (junctions, transformers),
it is a mass-matrix DAE: 𝕄 = diag(1 for charges, 0 for algebraic rows, 1 for ledgers).
"""
function SciMLBase.ODEProblem(cs::CompiledSystem, tspan; u0 = initial_state(cs),
                              sparse::Bool = cs.nphys > 50, kwargs...)
    if cs.nalg == 0
        f = SciMLBase.ODEFunction{true}(ode_rhs!; jac_prototype = sparse ? jacobian_sparsity(cs) : nothing)
    else
        M = Diagonal(vcat(ones(cs.nphys), zeros(cs.nalg), ones(cs.nledger)))
        f = SciMLBase.ODEFunction{true}(ode_rhs!; mass_matrix = M)
    end
    return SciMLBase.ODEProblem(f, u0, tspan, cs; kwargs...)
end

"""
    consistent_state(cs, u; t = 0.0) -> u

Solve the algebraic unknowns (junction forces, transformer fluxes) for the given charges, so
the DAE starts consistent. Errors if a transformer's constraint P_a = r P_b is violated by the
charges themselves: fix the initial charges instead.
"""
function consistent_state(cs::CompiledSystem, u::AbstractVector; t::Real = 0.0)
    cs.nalg == 0 && return copy(u)
    v = collect(Float64, u)
    P = copy(force_buffer(cs, v, t))
    for c in cs.transformers
        g = P[c.fa] - c.r * P[c.fb]
        abs(g) <= 1e-9 * max(abs(P[c.fa]), abs(c.r * P[c.fb]), 1e-12) ||
            throw(ArgumentError("transformer :$(c.name): initial state violates P_a = r P_b (residual $g); " *
                                "start from consistent charges"))
    end
    alg = cs.nphys+1:cs.nphys+cs.nalg
    function resid!(r, x, p)
        w = similar(x, length(v))
        w .= v
        w[alg] .= x
        du = similar(w)
        assemble!(du, w, p, t)
        r .= view(du, alg)
        return nothing
    end
    sol = SciMLBase.solve(NonlinearProblem(resid!, v[alg], cs), NewtonRaphson(); abstol = 1e-12)
    SciMLBase.successful_retcode(sol) || error("could not initialize the algebraic unknowns: retcode $(sol.retcode)")
    v[alg] .= sol.u
    return v
end

"""
    simulate(cs, tspan; alg = Rodas5P(), reltol = 1e-8, atol_scale = 1.0, ledger_control = true,
             atol_ledger = nothing, sparse, kwargs...)

Integrate the compiled system. Other keywords (`saveat`, `tstops`, `callback`, …) go to `solve`.
- Steps that would leave a material's valid states (e.g. a gas with V ≤ 0), or push a passive
  boundary uphill, are rejected.
- Pass the discontinuities of time-varying reservoir forces and `Switched` gates as `tstops`:
  they are evaluated as left limits up to each stop (D-29).

**Ledgers and tolerances** (DESIGN §12.5, D-31):
- `ledger_control = true` (default): the ledgers (entropy produced per element; charge, energy
  and entropy received by reservoirs and boundaries; external work) are in the integrator's
  error norm, so the losses it reports are as accurate as the states. Ledgers never feed back
  into the dynamics, and steady solves exclude them.
- `ledger_control = false` only when you need just the physical trajectories. Then check that
  `closure(cs, sol).energy_rel` is small: ledgers are uncontrolled quadratures and can be off by
  percent when the states are easy to integrate but the ledger integrands are not (a
  prescribed drive).
- `reltol` applies to everything. Absolute tolerances come from `abstol_vector`: charges from
  force tolerances (`atol_scale` scales them all), ledgers scaled like the states. To compare
  small differences of large ledgers (lost work = Q_h η_C − W), the relative error of the
  difference is about `reltol` × Q/ΔQ, so tighten `reltol` accordingly.
- `atol_ledger = x` overrides the ledger tolerances. Loosen it if steps stay small while ledgers
  are near zero; tighten it if ledgers near zero matter.
"""
function simulate(cs::CompiledSystem, tspan; alg = Rodas5P(), reltol::Real = 1e-8, atol_scale::Real = 1.0,
                  u0 = initial_state(cs), sparse::Bool = cs.nphys > 50 && cs.nalg == 0,
                  ledger_control::Bool = true, atol_ledger = nothing, kwargs...)
    u0 = consistent_state(cs, u0; t = first(tspan))
    prob = SciMLBase.ODEProblem(cs, tspan; u0, sparse)
    extra = cs.nalg == 0 ? (;) : (; initializealg = SciMLBase.NoInit())     # u0 is already consistent
    # discontinuities (`tstops`): evaluate time functions as left limits up to the next stop, and
    # move on to the next interval only once the step landing on the stop is accepted
    stops = sort(collect(Float64, get(kwargs, :tstops, Float64[])))
    nextstop(t) = (k = findfirst(>(t), stops); k === nothing ? Inf : stops[k])
    cs.tbreak[] = nextstop(first(tspan))
    advance = SciMLBase.DiscreteCallback((u, t, integ) -> t in stops, integ -> (cs.tbreak[] = nextstop(integ.t));
                                         save_positions = (false, false))
    cb = haskey(kwargs, :callback) ? SciMLBase.CallbackSet(kwargs[:callback], advance) : advance
    rest = Base.structdiff(values(kwargs), NamedTuple{(:callback,)})
    sol = try
        SciMLBase.solve(prob, alg; reltol, abstol = abstol_vector(cs; atol_scale, atol_ledger),
                        internalnorm = PhysNorm(ledger_control ? length(u0) : nstate(cs)), callback = cb,
                        isoutofdomain = (u, p, t) -> !admissible(p, u), extra..., rest...)
    finally
        cs.tbreak[] = Inf
    end
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
    base = collect(Float64, u0[1:nstate(cs)])
    function resid!(r, x, p)
        ufull = zeros(eltype(x), nstate(p) + p.nledger)
        ufull[1:nstate(p)] .= base
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
    # Scale each residual row by its round-off floor, so the stopping test does not depend on
    # how stiff a connection is: a 1e6 W/K tie has an entropy-rate noise of ~eps·G·T.
    # Row scaling leaves the Newton steps unchanged.
    x0 = base[act]
    J0 = ForwardDiff.jacobian((r, x) -> resid!(r, x, cs), zeros(length(act)), x0)
    rowscale = abs.(J0) * _charge_scales(cs, base)[act]
    for k in eachindex(rowscale)
        (isfinite(rowscale[k]) && rowscale[k] > 0) || (rowscale[k] = 1.0)
    end
    scaled!(r, x, p) = (resid!(r, x, p); r ./= rowscale; nothing)
    prob = NonlinearProblem(scaled!, x0, cs)
    # converged: the scaled residual is at the target, or Newton stalled on the round-off floor
    done(sol) = SciMLBase.successful_retcode(sol) ||
                (all(isfinite, sol.resid) && maximum(abs, sol.resid) <= 1e-10)
    sol = SciMLBase.solve(prob, NewtonRaphson(); abstol = 1e-14)
    if !done(sol)                                    # far from the solution: globalize
        sol = SciMLBase.solve(prob, TrustRegion(); abstol = 1e-14)
    end
    done(sol) ||
        error("steady-state solve failed: retcode $(sol.retcode). Try a transient warm start: " *
              "steady_state(cs; u0 = simulate(cs, (0, t)).u[end])")
    u = vcat(base, zeros(cs.nledger))
    u[act] .= sol.u
    v = passive_violations(cs, u)
    isempty(v) || _throw_passive(v, nothing)
    return u
end

"""
Per physical slot, the charge equivalent of the absolute force: capacity 1/A_kk times
max(|P_k|, 1). It sets each steady residual row's round-off floor (D-27).
"""
function _charge_scales(cs::CompiledSystem, u::AbstractVector)
    s = ones(nstate(cs))
    for (i, b) in enumerate(cs.bodies)
        off = cs.body_offset[i]
        N = ncharges(b.material)
        E = SVector{N}(u[off:off+N-1])
        A = eos_matrix(b.material, E)
        F = forces(b.material, E)
        fl = stiffness_floor(b.material)
        for k in 1:N
            v = max(abs(F[k]), 1.0) / max(abs(A[k, k]), fl[k])
            s[off+k-1] = (isfinite(v) && v > 0) ? v : 1.0
        end
    end
    return s
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
    for j in nb+1:nb+length(cs.reservoirs); anchored[j] = true; end
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
        for j in eachindex(cs.reservoirs)
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
function force(cs::CompiledSystem, u::AbstractVector, name::Symbol, d::FormOfMotion; t::Real = 0.0)
    i = get(cs.node_index, name, 0)
    i == 0 && throw(ArgumentError("no node named :$name"))
    nb, nr = nbodies(cs), length(cs.reservoirs)
    nb < i <= nb + nr && return reservoir_force(cs.reservoirs[i-nb], d, t)
    if i > nb + nr                                        # a junction: its algebraic unknown
        j = i - nb - nr
        k = findfirst(==(d), cs.junctions[j].domains)
        k === nothing && throw(ArgumentError(":$name has no $d form of motion"))
        return u[cs.junction_slots[j][k][2]]
    end
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
    transfer(_element(cs, name), copy(force_buffer(cs, u, t)), u, cs, t)

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
    _unwrap(c) isa CompiledConductor || throw(ArgumentError(":$name is not a conductor; use mass_flow for flow elements"))
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
    _unwrap(c) isa CompiledConductor && throw(ArgumentError(":$name is a conductor; use flux"))
    return _eval(cs, u, name; t).mdot
end

"""
    algebraic_index(cs, junction, domain)

Index in the state (and in `rates`) of a junction's force for `domain`; its `rates` entry is
the junction's net charge inflow (0 on a solution).
"""
function algebraic_index(cs::CompiledSystem, name::Symbol, d::FormOfMotion)
    j = findfirst(==(name), cs.junction_names)
    j === nothing && throw(ArgumentError("no junction named :$name"))
    k = findfirst(==(d), cs.junctions[j].domains)
    k === nothing && throw(ArgumentError("junction :$name has no $d form of motion"))
    return cs.junction_slots[j][k][2]
end

"""
    constraint_residual(cs, u, transformer; t = 0.0)

g = P_a − r P_b for a transformer at state `u` (0 on the constraint manifold). The power
mismatch λ g is the only place a transformer departs from exact power conservation (D-29).
"""
function constraint_residual(cs::CompiledSystem, u::AbstractVector, name::Symbol; t::Real = 0.0)
    _element(cs, name) isa CompiledTransformer || throw(ArgumentError(":$name is not a transformer"))
    return _eval(cs, u, name; t).g
end

_ledger_total(v, r::UnitRange{Int}) = sum((v[j] for j in r); init = zero(eltype(v)))

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

"External work [W, or J if integrated] supplied by every pump, flow boundary and source, read from `v`."
pump_work(cs::CompiledSystem, v::AbstractVector) = Dict(n => v[i] for (n, i) in cs.work_led)

"""
    source_ledger(cs, v, source) -> (; charge, W)

What a `Source` or `GravityField` has received from the model (`charge`; negative when it
supplies charge) and the external work `W` it has supplied, read from the ledger block of `v`.
These are integrated values from a state and rates from a rate vector. For gravity, W = −ΔPE.
"""
function source_ledger(cs::CompiledSystem, v::AbstractVector, name::Symbol)
    c = _unwrap(_element(cs, name))
    c isa CompiledSource || throw(ArgumentError(":$name is not a source"))
    return (; charge = v[c.led], W = v[c.led+1])
end

"""
    boundary_ledger(cs, v, element) -> (; m, S, U, W)

What a flow boundary, outlet, heat source or nuclear source *received* (mass, entropy, energy;
negative for a supplying one) and the external work it took in, read from the ledger block of
`v`. For a heat source, −U is the heat supplied and −S the entropy Q̇/T it brought in.
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

"Element symbols of every species that carries species data."
_elements(cs::CompiledSystem) = sort(unique(a for sd in values(cs.species_info) for a in keys(sd.atoms)))

"""
    element_totals(cs, v) -> Dict{Symbol,Float64}

Amount [mol] of each element (and isotope symbol) over the bodies' species, plus what
reservoirs, boundaries and nuclear sources received: Σ_s a_es m_s/M_s with m_s = `total_charge`.
Reactions and transport conserve it; nuclear sources do too, through their species ledgers.
Linear in `v`, so on a rate vector it gives the rates.
"""
function element_totals(cs::CompiledSystem, v::AbstractVector)
    out = Dict{Symbol,Float64}()
    for e in _elements(cs)
        out[e] = sum((atom_count(sd, e) / sd.M * total_charge(cs, v, Species(s)) for (s, sd) in cs.species_info); init = 0.0)
    end
    return out
end

"""
    reaction_rates(cs, u, reaction; t = 0.0) -> (; ξ̇, A, T, S_gen)

Extent rates ξ̇ [mol/s], affinities A [J/mol] (one per reaction), the body temperature and
the production A·ξ̇/T [W/K] of a `Reaction` element at state `u`.
"""
function reaction_rates(cs::CompiledSystem, u::AbstractVector, name::Symbol; t::Real = 0.0)
    _unwrap(_element(cs, name)) isa CompiledReaction || throw(ArgumentError(":$name is not a reaction"))
    tr = _eval(cs, u, name; t)
    return (; tr.ξ̇, tr.A, tr.T, tr.S_gen)
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
    for slots in cs.junction_slots, (_, ai) in slots           # a junction holds no energy: P_j × net inflow
        r += u[ai] * du[ai]
    end
    for n in cs.element_names
        tr = _eval(cs, u, n)
        for f in (:Q_a, :Q_b, :Q_h, :Q_d, :W, :E, :E_in, :Q)
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
  energy moved: through the connections, and converted inside the bodies (Σ|P·ΔE|, e.g.
  chemical to thermal in a reaction).
- entropy: max |ΔS_bodies + Δ(entropy received outside) − ΔS_gen|, which sits at round-off.
  `entropy_rel` divides by the entropy moved plus produced.
- charge: per non-thermal domain, max |Δ total_charge|, which sits at round-off.
  `charge_rel` is the largest, relative to the bodies' Σ|charge| (its maximum over the run); for
  a species, relative to the total mass of the bodies that hold it.
  Species changed by reactions are not each conserved: their element totals are checked
  instead (`element_totals`), in the same `charge_rel`.
- species: per conserved species (not changed by a reaction), `(; abs, rel)`: the absolute
  drift max |Δ total_charge| [kg] and that drift relative to the species' own largest body
  amount. Unlike `charge_rel`, `rel` is not diluted by the rest of the body, so a real loss of
  a trace species (tritium in water) shows at full size (D-43).
- elements: when reactions act, per element `(; abs, rel)`: drift [mol] and relative to the
  element's largest body amount.
"""
function closure(cs::CompiledSystem, sol)
    u0 = sol.prob.u0
    U0 = [internal_energy(b.material, charges(cs, u0, b.name)) for b in cs.bodies]
    E0b = [charges(cs, u0, b.name) for b in cs.bodies]
    S0 = total_entropy(cs, u0)
    X0, Y0 = _external_energy_terms(cs, u0), _external_entropy_terms(cs, u0)
    Sgen0 = sum(values(generated_entropy(cs, u0)); init = 0.0)
    # converted and reacting totals are not conserved; reacting species conserve their elements
    doms = [d for d in _nonthermal_domains(cs) if !(d in cs.converted) && !(d in cs.reacting)]
    elems = isempty(cs.reacting) ? Symbol[] : _elements(cs)
    E0 = [element_totals(cs, u0)[e] for e in elems]
    Q0 = [total_charge(cs, u0, d) for d in doms]
    body_charge(u, d) = sum((abs(u[cs.body_offset[i]+k-1]) for (i, b) in enumerate(cs.bodies)
                             for (k, dd) in enumerate(domains(b.material)) if dd == d); init = 0.0)
    # scale: the largest body total over the run (a tank filled from vacuum starts near zero). A
    # species is measured against all the mass of the bodies holding it: round-off in the coupled
    # solve scales with that, and a trace species (tritium in water) would otherwise read as drift.
    function holding_mass(u, d)
        s = 0.0
        for (i, b) in enumerate(cs.bodies)
            ds = domains(b.material)
            d in ds || continue
            s += sum((abs(u[cs.body_offset[i]+k-1]) for (k, dd) in enumerate(ds) if dd isa Species); init = 0.0)
        end
        return s
    end
    scale_of(u, d) = d isa Species ? holding_mass(u, d) : body_charge(u, d)
    Qscale = [max(maximum(scale_of(u, d) for u in sol.u), scale_of(u0, d), floatmin()) for d in doms]
    body_elem(u, e) = sum((atom_count(sd, e) / sd.M * body_charge(u, Species(s)) for (s, sd) in cs.species_info); init = 0.0)
    Escale = [max(maximum(body_elem(u, e) for u in sol.u), floatmin()) for e in elems]
    eE = 0.0; eS = 0.0; movedE = 0.0; movedS = 0.0; eQ = 0.0
    dQ = zeros(length(doms)); dEl = zeros(length(elems))       # absolute drifts
    for u in sol.u
        Ub = [internal_energy(b.material, charges(cs, u, b.name)) for b in cs.bodies]
        dXt = _external_energy_terms(cs, u) .- X0
        dYt = _external_entropy_terms(cs, u) .- Y0
        dX, dY = sum(dXt; init = 0.0), sum(dYt; init = 0.0)
        dSgen = sum(values(generated_entropy(cs, u)); init = 0.0) - Sgen0
        dS = total_entropy(cs, u) - S0
        eE = max(eE, abs(sum(Ub; init = 0.0) - sum(U0; init = 0.0) + dX))
        eS = max(eS, abs(dS + dY - dSgen))
        # energy converted between forms inside the bodies (chemical to thermal in a reaction)
        conv = sum((sum(abs.(forces(b.material, charges(cs, u, b.name)) .* (charges(cs, u, b.name) .- E0b[i])))
                    for (i, b) in enumerate(cs.bodies)); init = 0.0)
        movedE = max(movedE, sum(abs.(Ub .- U0); init = 0.0) + sum(abs, dXt; init = 0.0) + conv)
        movedS = max(movedS, abs(dS) + sum(abs, dYt; init = 0.0) + abs(dSgen))
        for (q, d) in enumerate(doms)
            dQ[q] = max(dQ[q], abs(total_charge(cs, u, d) - Q0[q]))
            eQ = max(eQ, dQ[q] / Qscale[q])
        end
        if !isempty(elems)
            el = element_totals(cs, u)
            for (q, e) in enumerate(elems)
                dEl[q] = max(dEl[q], abs(el[e] - E0[q]))
                eQ = max(eQ, dEl[q] / Escale[q])
            end
        end
    end
    # floors: when (almost) nothing moves, measure against the system's own magnitude
    # (round-off in these closures scales with the magnitudes of U and S themselves)
    movedE = max(movedE, 1e-3 * max(sum(abs, U0; init = 0.0), floatmin()))
    movedS = max(movedS, 1e-3 * max(abs(S0), 1.0))
    own(d) = max(maximum(body_charge(u, d) for u in sol.u), body_charge(u0, d), floatmin())
    species = Dict{Symbol,NamedTuple{(:abs, :rel),Tuple{Float64,Float64}}}(
        d.name => (abs = dQ[q], rel = dQ[q] / own(d)) for (q, d) in enumerate(doms) if d isa Species)
    elements = Dict{Symbol,NamedTuple{(:abs, :rel),Tuple{Float64,Float64}}}(
        e => (abs = dEl[q], rel = dEl[q] / Escale[q]) for (q, e) in enumerate(elems))
    return (; energy = eE, energy_rel = eE / movedE,
            entropy = eS, entropy_rel = eS / movedS,
            charge_rel = eQ, species, elements)
end
