# STATUS

*Rewritten at the end of every slice; this page is not a log. Last updated: 2026-10-14, after slice 1.*

## Where we are

- **Phase:** build. DESIGN.md revision 2 is approved.
- **Done:** slice 1 (thermal domain).
- **Next:** slice 2 (volumetric and species domains, ideal gas).
- **Tests:** 51,219 checks pass in ~35 s (`julia --project=. -e 'using Pkg; Pkg.test()'`). The first run also precompiles dependencies, which takes a few minutes.

## What exists (slice 1)

| Area | In `src/` | Notes |
|---|---|---|
| Domains | `Thermal` | traits only (`domains.jl`) |
| Materials | `ConstantCapacitySolid(m, c; T_ref)` | U(S) closed form; forces and Hessian by ForwardDiff (`materials.jl`) |
| Nodes | `Body`, `Reservoir(name; T)` | reservoir is thermal-only (D-20 §3) |
| Elements | `Conductor(ThermalConductance(G))` | kernel A, thermal row only (D-20 §2) |
| Model | `ThermoSystem`, `add!`, `connect!`, `compile` | bipartite MetaGraphsNext graph lowered to flat arrays (`system.jl`, `assembly.jl`) |
| Running | `simulate` (Rodas5P, ledgers excluded from the error norm), `steady_state` (Newton, anchored networks only), `rates` | `solve.jl` |
| Results | `temperature`, `flux`, `heat_flow`, `entropy_production`, `generated_entropy`, `reservoir_ledger`, `total_energy`, `total_entropy`, `energy_residual`, `closure` | `solve.jl` |
| Spatial | `rod!` (insulated, fixed-T or convective ends) | `spatial.jl` |

Validation covered (`test/`):
- kernel identities;
- A1, A3, A6 (A2 vacuous; A4 and A5 trivial or constructor-level);
- B2, C1, C2, C8;
- D1 (thermal analogue), D4 (covers D3 for thermal), D5 (slab).

## Known shortcuts and open items

All are logged in DECISIONS D-20:
- **Not built yet:**
  - kernel B (internal node) and non-thermal rows — slice 4;
  - `valid_state` / `isoutofdomain` — needed by the ideal gas in slice 2;
  - steady state for isolated networks (constraint rows).
- **Abstol is scaled by capacity** (DESIGN §11 updated). The relative tolerance still multiplies |S|, which depends on the reference; revisit with real-fluid reference states.
- **Sparse Jacobian** is used only above 50 states.
- **The old-library review phase** (`docs/EXISTING_LIBRARY.md`) has not started. The user decides when.

## Next: slice 2

- `Volumetric` and `Species` domains, `IdealGas` (closed form), `Composite` material, volume conductor.
- Fixed-mass capacity functions; `valid_state`.
- A general `Reservoir(name, forces)`.
- Tests: A4, A7, B3, B4, plus A2 (now meaningful).
