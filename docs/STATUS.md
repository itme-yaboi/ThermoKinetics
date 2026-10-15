# STATUS

*Rewritten at the end of every slice; this page is not a log. Last updated: 2026-10-14, after slice 3.*

## Where we are

- **Phase:** build. DESIGN.md revision 2 is approved.
- **Done:** slices 1, 2 and 3. The plotting design (DESIGN §17) is written but not built.
- **Next:** slice 4 (coupled conductors, electrical domain).
- **Tests:** 176,167 checks pass in ~4 min (`julia --project=. -e 'using Pkg; Pkg.test()'`). Most of that is one-time compilation of each element type's kernels and of the solver (D-26 §1). The first run after a dependency change also precompiles packages, which takes a few minutes.

## What exists

| Area | In `src/` | Notes |
|---|---|---|
| Domains | `Thermal`, `Volumetric` (force −p), `Species(name)` | traits and default force tolerances (`domains.jl`) |
| Materials | `ConstantCapacitySolid`, `IdealGas(species; R, cv)`, `ElasticVolume`, `Composite` | U(E) closed form; AD forces and Hessian; fixed-mass capacities; `valid_state`; `isentropic_state` (ideal gas). The ideal gas takes T and any two of m, p, V (`materials.jl`) |
| Nodes | `Body`, `Reservoir(name; T, p, material)` or `Reservoir(name, domain => force, …)` | a material reservoir fixes μ(T, p) and s, so it can supply or absorb streams |
| Conductors | `ThermalConductance(G)` (kernel A); `VolumeConductance(L)` with `InternalNode` (kernel B, default), `Split`, `Into` | coupled conductors rejected until slice 4 |
| Flow | `Stream(HydraulicConductance(κ))` (upwinded bundle), `IdealPump(mdot)` (isentropic, work logged), `FlowBoundary(material; T, p, mdot, supply = :passive \| :pumped)`, `PassThroughOutlet(; mdot)` or `(; κ, p_out)` | single species for pumps and boundaries (`elements.jl`, `assembly.jl`) |
| Model | `ThermoSystem`, `add!`, `connect!` (two-port, one-port, third-node), `compile` | bipartite MetaGraphsNext graph lowered to flat arrays; element groups in a `Vector{Any}` (one compiled solver for all models) |
| Running | `simulate` (Rodas5P; ledgers excluded from the error norm; invalid or uphill-passive steps rejected); `steady_state` (active charges only, anchoring checked per conserved quantity, Newton then trust region); `rates` | `solve.jl` |
| Results | `temperature`, `pressure`, `force`, `charges`, `flux`, `heat_flow`, `mass_flow`, `entropy_production(_parts)`, `generated_entropy`, `pump_work`, `reservoir_ledger`, `reservoir_charge`, `boundary_ledger`, `total_charge`, `total_energy`, `total_entropy`, `energy_residual`, `closure` (energy, entropy, charge; counts boundaries and work) | `solve.jl` |
| Spatial | `rod!` | `spatial.jl` |

Validation covered (`test/`):
- **Kernel identities:** thermal row, `Split`/`Into`, internal node, streams in both directions, pumps.
- **08 §A:** A1–A7 on solid, gas and stream networks.
- **08 §B:** B2, B3, B4, B5 (tank filling, including near-evacuated), plus Veinik's §23 and §11 checks.
- **08 §C:** C1, C2, C7 (stream mixing, steady), C8.
- **08 §D:** D1 (thermal and volumetric), D4, D5.
- **Pumped loop:** pump work = heat rejected = T₀ Σ Ṡ_gen (Gouy–Stodola).
- **Application tests (DESIGN §0.1):** composite wall (transient heat transfer); rigid-vessel and loose-piston heating (control mass); tank blowdown and counter-flow heat exchanger vs. ε–NTU (control volume).

## Known shortcuts and open items

- **Not built yet:**
  - coupled conductors and mass conductors — slice 4;
  - steady state for isolated groups (constraint rows);
  - turbulent (√Δp) flow driver;
  - multi-species pumps and boundaries;
  - inertance-driven flows;
  - real-fluid backends.
- **Steady state** of strongly coupled flow networks needs a transient warm start: `steady_state(cs; u0 = simulate(cs, (0, t)).u[end])` (D-26 §3).
- **Tolerances:**
  - Abstol is per domain, scaled by capacity (D-20 §1, D-22 §5). The relative tolerance still multiplies |S|, which depends on the reference.
  - Linear invariants hold to round-off amplified by implicit-solve conditioning (D-22 §6).
  - Stream production is ≥ 0 exactly only in exact arithmetic; floating-point tolerance is about −1e-12 (D-26 §4).
- **Compile time:** first use of each element type compiles its kernels; the solver compiles once (D-26 §1). Further cuts (precompile workloads) are possible later.
- **Network plotting** is designed (DESIGN §17, D-23), not built. It must use the Freivalds–Glagoļevs layout in `_reference/graphing/`, with Lengauer (local only) for routing. Separate track, scheduled when the user asks.
- **The old library is not consulted** (D-21).

## Next: slice 4

- Coupled conductors with a thermal row (kernel A, Onsager matrix data): a Knudsen barrier in the (mass, entropy) basis; a thermoelectric couple.
- `Electrical` domain, `LinearCapacitor`, resistor (kernel B); mass conductors on μ.
- Tests: C3, C4, C5, C6 (with basis transformations), D1 (electrical), D2; application test for electrical circuits (DESIGN §0.1).
