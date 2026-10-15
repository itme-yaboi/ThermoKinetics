# Getting started

The workflow is always the same: build a [`ThermoSystem`](@ref) with [`add!`](@ref) and
[`connect!`](@ref), [`compile`](@ref) it, then [`simulate`](@ref) or [`steady_state`](@ref) and
read results with accessors such as [`temperature`](@ref).

## Two solids relaxing to a common temperature

An aluminium-like and a copper-like block are joined by a thermal conductance. They relax to
the capacity-weighted mean temperature, and the entropy produced in the conductor is known in
closed form.

```@example solids
using ThermoKinetics

m = (1.5, 4.0); c = (900.0, 385.0); T = (600.0, 290.0)   # kg, J/(kg K), K
sys = ThermoSystem()
a = add!(sys, Body(:al, ConstantCapacitySolid(m[1], c[1]); T = T[1]))
b = add!(sys, Body(:cu, ConstantCapacitySolid(m[2], c[2]); T = T[2]))
k = connect!(sys, Conductor(ThermalConductance(3.0)), a, b)   # 3 W/K
cs = compile(sys)

C = m .* c                                   # heat capacities
τ = 1 / (3.0 * (1 / C[1] + 1 / C[2]))        # relaxation time
sol = simulate(cs, (0.0, 40τ))
uf = sol.u[end]

Tf = sum(C .* T) / sum(C)
(final = temperature(cs, uf, a), expected = Tf)
```

Entropy production is attributed per element. The total produced over the run,
[`generated_entropy`](@ref), compares with the exact value `Σ Cᵢ ln(T_f / Tᵢ)`; the production
*rate*, [`entropy_production`](@ref), has fallen to almost zero by the end.
[`closure`](@ref) reports the relative energy and entropy balance errors.

```@example solids
Sgen_exact = sum(C .* log.(Tf ./ T))
cl = closure(cs, sol)
(S_generated = generated_entropy(cs, uf)[k], exact = Sgen_exact,
 rate_at_end = entropy_production(cs, uf)[k],
 energy_rel = cl.energy_rel, entropy_rel = cl.entropy_rel)
```

## A reacting gas mixture

An ideal-gas mixture of two isomers in a closed rigid vessel reacts with a mass-action rate
law. Species data come from [`SpeciesData`](@ref) (molar data here, converted
internally); the reaction is declared by its stoichiometry.

```@example reaction
using ThermoKinetics

R̄ = 8.314462618                                   # J/(mol K)
gasdata(name, atoms; cp, s0) = SpeciesData(name; atoms, cv = cp - R̄, s0, basis = :molar)

Δs = 4.0                                           # entropy difference, J/(mol K)
A = gasdata(:nC4, (C = 4, H = 10); cp = 98.0, s0 = 310.0)
B = gasdata(:iC4, (C = 4, H = 10); cp = 98.0, s0 = 310.0 + Δs)
mix = IdealGasMixture(A, B)

kf, T0, V, nA0 = 2.0, 400.0, 1e-3, 0.03
sys = ThermoSystem()
add!(sys, Body(:v, mix; T = T0, V, n = (nC4 = nA0, iC4 = 1e-20 * nA0)))
connect!(sys, Reaction(:nC4 => -1, :iC4 => 1; law = MassAction(kf)), :v; name = :iso)
cs = compile(sys)

sol = simulate(cs, (0.0, 3.0); saveat = [0.05, 0.5, 3.0])
cl = closure(cs, sol)
(T = temperature(cs, sol.u[end], :v), p = pressure(cs, sol.u[end], :v),
 S_produced = generated_entropy(cs, sol.u[end]),
 energy_rel = cl.energy_rel, entropy_rel = cl.entropy_rel)
```

The temperature stays at 400 K because the isomers have equal heat capacities and formation
enthalpies; the reaction is driven purely by the entropy difference, and its entropy
production appears in [`generated_entropy`](@ref).

## A heat source in a rod

[`rod!`](@ref) discretises a conducting rod into thermal cells. Here both ends are cooled by a
coolant reservoir through a film coefficient, a [`HeatSource`](@ref) in every cell heats the
rod uniformly, and [`steady_state`](@ref) solves for the steady profile. The peak temperature
is close to the exact `T∞ + qL/(2h) + qL²/(8λ)` (second-order accurate in the cell size).

```@example rod
using ThermoKinetics

λ, ρ, c, L, area, h = 50.0, 8000.0, 500.0, 0.1, 1e-3, 100.0
n, Tinf, q = 10, 300.0, 1e4                        # cells, K, W/m³
sys = ThermoSystem()
add!(sys, Reservoir(:coolant; T = Tinf))
rod = rod!(sys, :slab; n, L, area, λ, ρ, c, T0 = Tinf, left = (:coolant, h), right = (:coolant, h))
for cell in rod.cells
    connect!(sys, HeatSource(q * area * rod.dx), cell)
end
cs = compile(sys)
u = steady_state(cs)

Tmax = maximum(temperature(cs, u, cell) for cell in rod.cells)
Tc_exact = Tinf + q * L / (2h) + q * L^2 / (8λ)
Sgen = sum(values(entropy_production(cs, u)))      # steady loss rate, W/K
(Tmax = Tmax, exact_centre = Tc_exact, total_Sgen_rate = Sgen)
```
