# STATUS

*Rewritten at the end of every slice; this page is not a log. Last updated: 2026-10-15, after slice 7.*

## Where we are

- **Repository:** the package is `ThermoKinetics` (D-38), on the private GitHub repository `itme-yaboi/ThermoKinetics` (`origin`). `docs/veinik/` and the Veinik section notes are local only and must never be committed (D-36, D-37).
- **Phase:** build. DESIGN.md revision 2 is approved. Of the revision 3 additions, mechanics (DESIGN §18) is accepted and built; rockets (§19) are still proposed.
- **Done:** slices 1–7. Slice 7 was reactions, nuclear and heat sources (D-39 to D-43, all accepted). The plotting design (DESIGN §17) is written but not built.
- **Next:** slice 8, rockets (DESIGN §19, D-33; proposed, not reviewed). Slice 7 supplies the gas mixture with formation enthalpies and the reaction element; the limits are in §19.6.
- **Later:** the real-fluid backend (Clapeyron) for Rankine and Brayton; E3/E4; the plotting track; see DESIGN §15.
- **Tests:** 192,952 checks pass in ~2.3 min (`julia --project=. -e 'using Pkg; Pkg.test()'`). The first run after a dependency change also precompiles packages, which takes a few minutes.

## What exists

| Area | In `src/` | Notes |
|---|---|---|
| Domains | `Thermal`, `Volumetric` (−p), `Species`, `Electrical` (φ), `Momentum` (one axis, w), `AngularMomentum` (ω), `FluxLinkage` (i), `Displacement` (x, F), `Angle` (θ, M) | traits and default force tolerances |
| Materials | `ConstantCapacitySolid`, `IdealGas`, `ElasticVolume`, `LinearCapacitor`, `InertialMass`, `RotationalInertia`, `LinearInductor`, `LinearSpring`, `TorsionSpring`, `PendulumPotential`, `Composite`; `IdealGasMixture`, `IdealSolution` from `SpeciesData` | U(E) closed form; AD forces and Hessian; `stiffness_floor` for non-convex U; `inertial_mass`; mixtures: formation data, molar masses from atoms, constant c_v, trace floor (D-41) |
| Nodes | `Body`; `Reservoir` (constant or time-varying forces `domain => t -> P(t)`, or a material); `Junction` (zero capacity, algebraic forces) | |
| Conductors | `ThermalConductance`; `LinearConductance` (`VolumeConductance`, `Resistor`, `MassConductance`), `QuadraticDrag`, `CoulombFriction`, with `InternalNode` / `Split` / `Into`; coupled `ThermoelectricConductance`, `KnudsenBarrier`, `OnsagerConductance` | basis maps from Veinik's data; monotone laws are D ≥ 0 by type |
| Converters | `Converter(W, domains)` (antisymmetric W), `Gyrator`, `InductorBranch`, `SpringBranch`, `TorsionBranch`; `Transformer` (index-reduced DAE, Baumgarte τ: gears, rack and pinion) | power-exact; converted totals are not conserved (by design) |
| Sources | `Source(domain => rate)` (force, torque, current; constant or t -> rate); `GravityField(g)` (−Mg; special-purpose, flagged, D-34); `HeatSource(Q̇)`; `NuclearSource(ν; rate = FirstOrder(λ or half_life) or t -> ξ̇, E_dep, external)` | charge received and external work ledgered (`source_ledger`); heat and nuclear sources: energy, entropy and species received (`boundary_ledger`); nucleons and charge checked at compile; dU = E_dep ξ̇ |
| Reactions | `Reaction(ν...; law)` or `Reaction([ν₁, ν₂]; law)` with `AffinityLaw(L)` (scalar or PSD matrix) or `MassAction(k0; b, Ea, orders)` | inside one body, dU = 0, A·ξ̇ ≥ 0 exactly; element balance checked at compile; `reaction_rates`, `element_totals` |
| Flow | `Stream`, `IdealPump`, `FlowBoundary`, `PassThroughOutlet` | single species for pumps and boundaries |
| Switching | `Switched(element; on = t -> Bool)` (also clutches and brakes) | discontinuities as `tstops`; evaluated as left limits |
| Running | `simulate` (Rodas5P; mass-matrix DAE when junctions or transformers exist; `consistent_state`; inadmissible steps rejected, NaN on invalid trial stages; ledgers error-controlled by default with scaled tolerances, `atol_ledger` override); `steady_state` | tolerances: DESIGN §12.5 |
| Results | as before, plus `algebraic_index`, `constraint_residual`, `source_ledger`, `reaction_rates`, `element_totals` | closure checks element totals when reactions act; per-species and per-element absolute drift (`closure(...).species`, `.elements`, D-43) |
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
- **08 §H:** H1–H10 (isomerization, equilibrium by Gibbs minimization and entropy maximization, adiabatic combustion, calorimeter, Onsager-coupled reactions, tritium decay, Bateman chain, heat source, irradiated water, approach to zero); mixture and reaction kernel identities; A1/A3 and A → 0 on random reacting networks.
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
  - rotational mechanics (08 G6–G8);
  - slab with internal heating (parabolic profile, second-order convergence);
  - water–gas shift equilibrium at 1000 K vs. Gibbs minimization, heat of reaction to the furnace;
  - tritium storage bed: inventory, decay heat to the wall, quasi-steady bed temperature;
  - tritium generation in irradiated water (D(n,γ)T, then decay) vs. Bateman, to saturation.

## Known shortcuts and open items

- **`GravityField`** is a special-purpose momentum source, flagged to revisit (D-34).
  - Revisit when hydrostatic or elevation heads (μ + gz in streams), altitude-dependent gravity, or height as a state are needed.
  - It ledgers work only; the ∫ w dt ledger for variable mass is due in slice 8.
- **Mechanics limits** (DESIGN §18.7): one axis; no 3-D rigid-body rotation (Euler's equations); Coulomb friction regularized, with no true sticking; no lossy (mesh-efficiency) transformer; no linkages or contact.
- **Mixtures and reactions** (DESIGN §8.1, D-41):
  - constant c_v per species, so high-temperature accuracy is limited (NASA polynomials come with the property backends);
  - the trace floor: below 1e-13 of a body's mass the mixing term is smoothed, so a depleted species can read about −1e-12 of the body's mass instead of failing (`trace` keyword). When it matters and how to handle it: D-42. **Check D-44 for any reaction or nuclear process whose low-mass-fraction species you care about** (discharge-level tritium, ~10 kBq/L, sits at the floor);
  - mass action needs a body with a volume (gas); a solution uses `AffinityLaw`;
  - no mixture reservoirs or multi-species flow boundaries yet; reacting streams are out of scope;
  - `steady_state` refuses closed reacting groups (their element totals are conserved: the unimplemented constraint-row case); use `simulate`.
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
  - the membrane application test;
  - E3, E4;
  - temperature-dependent c_p (NASA polynomials);
  - radiation conductance;
  - two-phase flow;
  - rockets (§19).
- **Steady state** of strongly coupled flow networks needs a transient warm start (D-26 §3).
- **Tolerances:** linear invariants hold to round-off amplified by implicit-solve conditioning; "exact" zero-rate charges move by ~1 ulp (D-29 §8); stream production is ≥ 0 exactly only in exact arithmetic (D-26 §4).
- **Network plotting** is designed (DESIGN §17, D-23, D-28), not built. All its reference PDFs are local only (D-30).
- **The old library is not consulted** (D-21).
