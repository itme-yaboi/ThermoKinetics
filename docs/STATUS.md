# STATUS

*Rewritten at the end of every slice; this page is not a log. Last updated: 2026-10-14, after slice 4.*

## Where we are

- **Phase:** build. DESIGN.md revision 2 is approved.
- **Done:** slices 1–4. The plotting design (DESIGN §17) is written but not built.
- **Next:** slice 5, the mass-matrix DAE path (transformer converters, zero-capacity non-thermal nodes, kinetic DOFs, gyrator converters). It must come before balance-of-plant cycle work.
- **Tests:** 184,436 checks pass in ~3 min (`julia --project=. -e 'using Pkg; Pkg.test()'`). Most of that is one-time compilation of each element type's kernels and of the solver (D-26 §1). The first run after a dependency change also precompiles packages, which takes a few minutes.

## What exists

| Area | In `src/` | Notes |
|---|---|---|
| Domains | `Thermal`, `Volumetric` (force −p), `Species(name)`, `Electrical` (force φ) | traits and default force tolerances |
| Materials | `ConstantCapacitySolid`, `IdealGas`, `ElasticVolume`, `LinearCapacitor`, `Composite` | U(E) closed form; AD forces and Hessian; fixed-mass capacities; `valid_state`; `isentropic_state` (ideal gas) |
| Nodes | `Body`, `Reservoir(name; T, p, material)` or `Reservoir(name, domain => force, …)` | material reservoirs supply and absorb streams |
| Conductors | `ThermalConductance(G)`; `LinearConductance` via `VolumeConductance`, `Resistor`, `MassConductance` with `InternalNode` / `Split` / `Into`; coupled `ThermoelectricConductance`, `KnudsenBarrier`, `OnsagerConductance(rows, M; W, magnetic)` (kernel A, `Split`, ‖Bᵀx‖² ledger) | basis maps: `thermoelectric_from_veinik`, `barrier_mass_basis`, `electrokinetic_mass_basis` |
| Flow | `Stream(HydraulicConductance(κ))`, `IdealPump`, `FlowBoundary(…; supply)`, `PassThroughOutlet` | single species for pumps and boundaries |
| Model | `ThermoSystem`, `add!`, `connect!`, `compile` | bipartite MetaGraphsNext graph lowered to flat arrays; element groups in a `Vector{Any}` |
| Running | `simulate` (Rodas5P; ledgers excluded from the error norm; inadmissible steps rejected); `steady_state` (active charges, anchoring per conserved quantity, row-scaled Newton then trust region); `rates` | |
| Results | `temperature`, `pressure`, `force`, `charges`, `flux`, `heat_flow`, `mass_flow`, `entropy_production(_parts)`, `generated_entropy`, `pump_work`, `reservoir_ledger`, `reservoir_charge`, `boundary_ledger`, `total_charge`, `total_energy`, `total_entropy`, `energy_residual`, `closure` | |
| Spatial | `rod!` | |

Validation covered (`test/`):
- **Kernel identities:** thermal row, `Split`/`Into`, internal node, streams, pumps, coupled kernel at ΔT/T ~ 1.
- **08 §A:** A1–A7 on solid, gas and stream networks. A5 also covers Onsager validation and Onsager–Casimir (W does no dissipation).
- **08 §B:** B2–B5.
- **08 §C:**
  - C1, C2, C3 (Joule, also through an internal node);
  - C4 (Knudsen: the basis transformation, and simulated p′/p″ = √(T′/T″));
  - C5 (Seebeck (278) exactly, Kelvin 2, Kelvin 1 for constant α);
  - C6 (Saxén in both bases);
  - C7, C8.
- **08 §D:** D1 (thermal, volumetric, electrical), D2, D4, D5. The μ contact rule for mass conductors.
- **Application tests (DESIGN §0.1):**
  - composite wall (transient heat transfer);
  - rigid-vessel and loose-piston heating (control mass);
  - tank blowdown and counter-flow heat exchanger (control volume);
  - resistor self-heating with its loss map, and thermoelectric generator efficiency including η_max (electrical circuits).

## Known shortcuts and open items

- **Not built yet:**
  - DAE path, converters, kinetic DOFs, inductors — slice 5;
  - steady state for isolated groups;
  - turbulent flow driver;
  - multi-species pumps and boundaries;
  - inertance-driven flows;
  - real-fluid backends;
  - temperature-dependent Seebeck α(T);
  - electrochemical potential at material boundaries;
  - the membrane-permeation application test.
- **Steady state** of strongly coupled flow networks needs a transient warm start (D-26 §3). Residual rows are scaled by their round-off floor (D-27 §7).
- **Tolerances:**
  - Abstol is per domain, scaled by capacity.
  - Linear invariants hold to round-off amplified by implicit-solve conditioning.
  - Stream production is ≥ 0 exactly only in exact arithmetic (D-26 §4).
- **Compile time:** first use of each element type compiles its kernels; the solver compiles once (D-26 §1).
- **Network plotting** is designed (DESIGN §17, D-23), not built. It must use the Freivalds–Glagoļevs layout in `_reference/graphing/`. Every phase's reference is now there (DESIGN §17.1, D-28).
- **The old library is not consulted** (D-21).

## Next: slice 5

- Mass-matrix DAE form (𝕄 du/dt = f with singular 𝕄) on Rodas5P/FBDF; algebraic rows with λ flux variables.
- `Momentum`/`AngularMomentum`/`FluxLinkage` domains and inertial materials; gyrator converters (piston, motor); transformer converters (DAE); zero-capacity non-thermal nodes.
- Tests: E1 (Carnot with reversible converters), E2 (Carnot with bridges), DAE consistency (constraint residuals at round-off, Pᵀ F ≡ 0); control-mass piston–cylinder compression application test.
