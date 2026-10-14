# STATUS

*Rewritten at the end of every slice; this page is not a log. Last updated: 2026-10-14, after slice 2.*

## Where we are

- **Phase:** build. DESIGN.md revision 2 is approved.
- **Done:** slices 1 and 2.
- **Next:** slice 3 (streams).
- **Tests:** 174,545 checks pass in ~75 s (`julia --project=. -e 'using Pkg; Pkg.test()'`). Most of that time is solver compilation for each new combination of element types (D-22 §9). The first run after a dependency change also precompiles packages, which takes a few minutes.

## What exists

| Area | In `src/` | Notes |
|---|---|---|
| Domains | `Thermal`, `Volumetric` (force −p), `Species(name)` | traits and default force tolerances (`domains.jl`) |
| Materials | `ConstantCapacitySolid`, `IdealGas(species; R, cv)`, `ElasticVolume(C; V0, p0)`, `Composite(parts...)` | U(E) closed form; forces and Hessian by AD; `capacity_matrix` and `is_stable` on the fixed-mass block; `valid_state` (`materials.jl`) |
| Nodes | `Body`, `Reservoir(name; T, p)` or `Reservoir(name, domain => force, …)` | reservoirs ledger received charge per domain and energy |
| Conductors | `ThermalConductance(G)` (kernel A); `VolumeConductance(L)` with `InternalNode(G_a, G_b[, G_h, h])` (kernel B, default), `Split(w)` or `Into(node)` (kernel A) | coupled conductors rejected until slice 4 (`elements.jl`, `assembly.jl`) |
| Model | `ThermoSystem`, `add!`, `connect!` (third-node port for `Into`/`InternalNode(h)`), `compile` | bipartite MetaGraphsNext graph lowered to flat arrays; groups sorted by type |
| Running | `simulate` (Rodas5P; ledgers excluded from the error norm; invalid steps rejected), `steady_state` (anchored networks only), `rates` | `solve.jl` |
| Results | `temperature`, `pressure`, `force`, `charges`, `flux`, `heat_flow`, `entropy_production(_parts)`, `generated_entropy`, `reservoir_ledger`, `reservoir_charge`, `total_charge`, `total_energy`, `total_entropy`, `energy_residual`, `closure` (energy, entropy, charge) | `solve.jl` |
| Spatial | `rod!` (insulated, fixed-T or convective ends) | `spatial.jl` |

Validation covered (`test/`):
- **Kernel identities:** thermal row, `Split`/`Into`, internal node.
- **08 §A:** A1, A3, A6 on solid and gas networks; A2 (volume and species mass); A4 (raw Hessians); A5; A7.
- **08 §B:** B2, B3, B4, plus Veinik's §23 equal-pressure case and §11 (an uncoupled body is a pure accumulator; an ideal gas compressed with heat removed elsewhere is exactly isentropic).
- **08 §C:** C1, C2, C8.
- **08 §D:** D1 (thermal analogue and volumetric proper), D4, D5.
- **Application tests (DESIGN §0.1):** composite wall (transient heat transfer); rigid-vessel and loose-piston heating (control mass).

## Known shortcuts and open items

- **Not built yet:**
  - coupled conductors (thermal + non-thermal rows) and mass conductors — slice 4;
  - steady state for isolated networks (constraint rows).
- **Tolerances:**
  - Abstol is per domain, scaled by capacity (D-20 §1, D-22 §5). The relative tolerance still multiplies |S|, which depends on the reference; revisit with real-fluid reference states.
  - Linear invariants (charge, entropy closure) hold to round-off amplified by implicit-solve conditioning, ≲1e-11 relative (D-22 §6).
- **Compile latency:** several seconds for each new combination of element types (D-22 §9).
- **Sparse Jacobian:** used only above 50 states.
- **The old library is not consulted** (D-21). Specific correlations or validated cases are requested from the author when needed.

## Next: slice 3

- Streams (an upwinded charge bundle), hydraulic drivers (ṁ = κΔp), `IdealPump`, labeled `FlowBoundary` (`:passive` or `:pumped`), `PassThroughOutlet`.
- Reservoirs that carry a material (so s, h and μ are defined for inflow).
- Tests: B5 (tank filling), C7 (stream mixing), a pumped-loop Gouy–Stodola closure, and application tests for control-volume systems (tank blowdown; counter-flow heat exchanger vs. ε–NTU, if it fits the slice).
