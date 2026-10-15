# STATUS

*Rewritten at the end of every slice; this page is not a log. Last updated: 2026-10-15, after slice 6.*

## Where we are

- **Phase:** build. DESIGN.md revision 2 is approved. Of the revision 3 additions, mechanics (DESIGN §18) is accepted and built; rockets (§19) are still proposed.
- **Done:** slices 1–6. Slice 6 was mechanical and rotational systems (D-32, D-34, D-35). The plotting design (DESIGN §17) is written but not built.
- **Next:** slice 7, reactions and nuclear sources. Not started; scope and tests are to be detailed at the start of the slice.
  - It includes a mass-action rate law. Only if that presents an issue does it move to future work, with the issue recorded (D-34).
  - Slice 8 needs an ideal-gas mixture with formation enthalpies from it.
- **After that:** slice 8, rockets (DESIGN §19, D-33; proposed, not reviewed). The framework can model them; the limits are in §19.6.
- **Later:** the real-fluid backend (Clapeyron) for Rankine and Brayton; E3/E4; the plotting track; see DESIGN §15.
- **Tests:** 186,336 checks pass in ~1.6 min (`julia --project=. -e 'using Pkg; Pkg.test()'`). The first run after a dependency change also precompiles packages, which takes a few minutes.

## What exists

| Area | In `src/` | Notes |
|---|---|---|
| Domains | `Thermal`, `Volumetric` (−p), `Species`, `Electrical` (φ), `Momentum` (one axis, w), `AngularMomentum` (ω), `FluxLinkage` (i), `Displacement` (x, F), `Angle` (θ, M) | traits and default force tolerances |
| Materials | `ConstantCapacitySolid`, `IdealGas`, `ElasticVolume`, `LinearCapacitor`, `InertialMass`, `RotationalInertia`, `LinearInductor`, `LinearSpring`, `TorsionSpring`, `PendulumPotential`, `Composite` | U(E) closed form; AD forces and Hessian; `stiffness_floor` for non-convex U; `inertial_mass` |
| Nodes | `Body`; `Reservoir` (constant or time-varying forces `domain => t -> P(t)`, or a material); `Junction` (zero capacity, algebraic forces) | |
| Conductors | `ThermalConductance`; `LinearConductance` (`VolumeConductance`, `Resistor`, `MassConductance`), `QuadraticDrag`, `CoulombFriction`, with `InternalNode` / `Split` / `Into`; coupled `ThermoelectricConductance`, `KnudsenBarrier`, `OnsagerConductance` | basis maps from Veinik's data; monotone laws are D ≥ 0 by type |
| Converters | `Converter(W, domains)` (antisymmetric W), `Gyrator`, `InductorBranch`, `SpringBranch`, `TorsionBranch`; `Transformer` (index-reduced DAE, Baumgarte τ: gears, rack and pinion) | power-exact; converted totals are not conserved (by design) |
| Sources | `Source(domain => rate)` (force, torque, current; constant or t -> rate); `GravityField(g)` (−Mg; special-purpose, flagged, D-34) | charge received and external work ledgered (`source_ledger`) |
| Flow | `Stream`, `IdealPump`, `FlowBoundary`, `PassThroughOutlet` | single species for pumps and boundaries |
| Switching | `Switched(element; on = t -> Bool)` (also clutches and brakes) | discontinuities as `tstops`; evaluated as left limits |
| Running | `simulate` (Rodas5P; mass-matrix DAE when junctions or transformers exist; `consistent_state`; inadmissible steps rejected, NaN on invalid trial stages; ledgers error-controlled by default with scaled tolerances, `atol_ledger` override); `steady_state` | tolerances: DESIGN §12.5 |
| Results | as before, plus `algebraic_index`, `constraint_residual`, `source_ledger` | |
| Spatial | `rod!` | |

Validation covered (`test/`):
- **Kernel identities**, all element kinds (spring branches Pᵀ W P = 0; friction laws D ≥ 0 and odd).
- **08 §A:** A1–A7, plus A1/A3 on random mechanical chains.
- **08 §B:** B2–B5.
- **08 §C:** C1–C8.
- **08 §D:** D1 (thermal, volumetric, electrical), D2, D4, D5.
- **08 §E:**
  - E1 (Carnot: the deficit and Ṡ_gen fall tenfold for a tenfold G);
  - E2 (η = 1 − T₂/T₁ at the internal temperatures, exactly; Gouy–Stodola).
- **08 §G:** G1–G9.
  - **Translational:**
    - mass–spring–damper (free and forced);
    - two-mass normal modes;
    - spring under gravity, including the half-energy law and a direct steady state;
    - pendulum period vs. the elliptic integral, at θ₀ = 0.01, 1, 2.5 and π/2, plus damped heat;
    - quadratic drag and Coulomb friction.
  - **Rotational:**
    - gear-train spin-up with the reflected inertia;
    - torsional two-inertia oscillation;
    - clutch loss independent of the conductance.
- **Converters and DAE:**
  - an LC oscillator (exact period and energy);
  - RLC through a junction;
  - a frictionless gas spring (isentropic to the ulp);
  - a DC motor with a junction (steady speed V/(Rb/k + k)), and the same motor driving a load through a gearbox;
  - rack and pinion;
  - a transformer (RC charging sees C₁ + C₂/r²; constraint to 1e-8);
  - switching with a time-varying reservoir (exact piecewise exponentials).
- **Application tests (DESIGN §0.1):**
  - composite wall;
  - rigid-vessel and loose-piston heating;
  - piston–cylinder compression by the atmosphere with seal friction (gas exactly on its isentrope; friction heat and entropy);
  - tank blowdown;
  - counter-flow heat exchanger;
  - resistor self-heating;
  - thermoelectric generator;
  - translational mechanics (08 G1, G2, G5);
  - rotational mechanics (08 G6–G8).

## Known shortcuts and open items

- **`GravityField`** is a special-purpose momentum source, flagged to revisit (D-34).
  - Revisit when hydrostatic or elevation heads (μ + gz in streams), altitude-dependent gravity, or height as a state are needed.
  - It ledgers work only; the ∫ w dt ledger for variable mass is due in slice 8.
- **Mechanics limits** (DESIGN §18.7): one axis; no 3-D rigid-body rotation (Euler's equations); Coulomb friction regularized, with no true sticking; no lossy (mesh-efficiency) transformer; no linkages or contact.
- **Ledger error control** is on by default, with ledger tolerances scaled like the states (D-31). Overhead is 0–10 % except where it buys accuracy. When it is required vs. optional, and how to set tolerances: DESIGN §12.5 and the `simulate` docstring.
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
  - heat sources;
  - the membrane application test;
  - E3, E4;
  - reactions and nuclear sources;
  - two-phase flow;
  - rockets (§19).
- **Steady state** of strongly coupled flow networks needs a transient warm start (D-26 §3).
- **Tolerances:** linear invariants hold to round-off amplified by implicit-solve conditioning; "exact" zero-rate charges move by ~1 ulp (D-29 §8); stream production is ≥ 0 exactly only in exact arithmetic (D-26 §4).
- **Network plotting** is designed (DESIGN §17, D-23, D-28), not built. All its reference PDFs are local only (D-30).
- **The old library is not consulted** (D-21).
