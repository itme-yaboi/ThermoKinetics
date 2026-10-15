# STATUS

*Rewritten at the end of every slice; this page is not a log. Last updated: 2026-10-14, after slice 5.*

## Where we are

- **Phase:** build. DESIGN.md revision 2 is approved.
- **Done:** slices 1–5. The plotting design (DESIGN §17) is written but not built.
- **Next:** not yet chosen; see DESIGN §15 "Later". Candidates:
  - the real-fluid backend (Clapeyron) for Rankine and Brayton cycles;
  - E3/E4;
  - the plotting track;
  - reactions and nuclear sources.
- **Open question for the user:** should `ledger_control` default to `true`? (D-29 §11; it touches D-10.)
- **Tests:** 184,615 checks pass in ~1.5 min (`julia --project=. -e 'using Pkg; Pkg.test()'`). The first run after a dependency change also precompiles packages, which takes a few minutes.

## What exists

| Area | In `src/` | Notes |
|---|---|---|
| Domains | `Thermal`, `Volumetric` (−p), `Species`, `Electrical` (φ), `Momentum` (one axis, w), `AngularMomentum` (ω), `FluxLinkage` (i) | traits and default force tolerances |
| Materials | `ConstantCapacitySolid`, `IdealGas`, `ElasticVolume`, `LinearCapacitor`, `InertialMass`, `RotationalInertia`, `LinearInductor`, `Composite` | U(E) closed form; AD forces and Hessian |
| Nodes | `Body`; `Reservoir` (constant or time-varying forces `domain => t -> P(t)`, or a material); `Junction` (zero capacity, algebraic forces) | |
| Conductors | `ThermalConductance`; `LinearConductance` (`VolumeConductance`, `Resistor`, `MassConductance`) with `InternalNode` / `Split` / `Into`; coupled `ThermoelectricConductance`, `KnudsenBarrier`, `OnsagerConductance` | basis maps from Veinik's data |
| Converters | `Converter(W, domains)` (antisymmetric W), `Gyrator`, `InductorBranch`; `Transformer` (index-reduced DAE, Baumgarte τ) | power-exact; converted totals are not conserved (by design) |
| Flow | `Stream`, `IdealPump`, `FlowBoundary`, `PassThroughOutlet` | single species for pumps and boundaries |
| Switching | `Switched(element; on = t -> Bool)` | discontinuities as `tstops`; evaluated as left limits |
| Running | `simulate` (Rodas5P; mass-matrix DAE when junctions or transformers exist; `consistent_state`; inadmissible steps rejected, NaN on invalid trial stages; `ledger_control`); `steady_state` | |
| Results | as before, plus `algebraic_index`, `constraint_residual` | |
| Spatial | `rod!` | |

Validation covered (`test/`):
- **Kernel identities**, all element kinds.
- **08 §A:** A1–A7.
- **08 §B:** B2–B5.
- **08 §C:** C1–C8.
- **08 §D:** D1 (thermal, volumetric, electrical), D2, D4, D5.
- **08 §E:**
  - E1 (Carnot: the deficit and Ṡ_gen fall tenfold for a tenfold G);
  - E2 (η = 1 − T₂/T₁ at the internal temperatures, exactly; Gouy–Stodola).
- **Converters and DAE:**
  - an LC oscillator (exact period and energy);
  - RLC through a junction;
  - a frictionless gas spring (isentropic to the ulp);
  - a DC motor with a junction (steady speed V/(Rb/k + k));
  - a transformer (RC charging sees C₁ + C₂/r²; constraint to 1e-8);
  - switching with a time-varying reservoir (exact piecewise exponentials).
- **Application tests (DESIGN §0.1):**
  - composite wall;
  - rigid-vessel and loose-piston heating;
  - piston–cylinder compression by the atmosphere with seal friction (gas exactly on its isentrope; friction heat and entropy);
  - tank blowdown;
  - counter-flow heat exchanger;
  - resistor self-heating;
  - thermoelectric generator.

## Known shortcuts and open items

- **Ledger error control** is opt-in (`ledger_control = true`). Without it, ledgers are uncontrolled quadratures when the physical states are easy to integrate, as with a prescribed drive; `closure` reveals it (D-29 §11). Awaiting the user's decision on the default.
- **Transformers** conserve power on the constraint manifold. The off-manifold mismatch stays at integrator tolerance, and the initial charges must satisfy the constraint (D-29 §1). DAE systems use a dense Jacobian.
- **Not built yet:**
  - vector momentum;
  - steady state for isolated groups;
  - turbulent flow driver;
  - multi-species pumps and boundaries;
  - inertance-driven flows;
  - real-fluid backends;
  - α(T) thermoelectrics;
  - contact potentials;
  - the membrane application test;
  - E3, E4;
  - reactions and nuclear sources;
  - two-phase flow.
- **Steady state** of strongly coupled flow networks needs a transient warm start (D-26 §3).
- **Tolerances:** linear invariants hold to round-off amplified by implicit-solve conditioning; "exact" zero-rate charges move by ~1 ulp (D-29 §8); stream production is ≥ 0 exactly only in exact arithmetic (D-26 §4).
- **Network plotting** is designed (DESIGN §17, D-23, D-28), not built. All its reference PDFs are local only (D-30).
- **The old library is not consulted** (D-21).
