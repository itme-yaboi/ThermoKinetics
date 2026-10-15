# FAILURE MODES

*Known edge cases: where the framework as built may give a wrong result, or will fail, and we proceed anyway. Last updated: 2026-10-15, after slice 10b.*

## What belongs here

Only scenarios in which the current code may (or will) mishandle a model: a silently wrong number, a misattributed loss, a solve that fails, a refusal a user might not expect. Each was accepted knowingly, usually as the price of a decision recorded in `DECISIONS.md`.

Not here: features that are simply not built (DESIGN §15 and `STATUS.md` list those), and the rationale for decisions (`DECISIONS.md`).

**Rules for this file**
- One entry per failure mode, with a stable ID (`FM-n`). IDs are never reused.
- Entries are kept current, not appended to: when the behaviour changes, rewrite the entry. When a fix removes the failure, delete the entry and record the fix in `DECISIONS.md`, citing the ID.
- When a decision is accepted knowing it can fail in some cases, add the entry in the same commit.
- **Silent** means a result comes back with no error; the dangerous kind. **Loud** means an error, a refusal or a failed solve.

Each entry gives: what happens, when, why it was accepted (decision), how to detect it, what to do, and when to revisit.

---

## 1. Real fluids and phases (Clapeyron backend)

### FM-1 Phase test trusts the single-phase inversion first — silent, believed rare
- **What:** for a pure fluid the uniform-fluid inversion S(T; V, n) = S runs first. Its state is accepted as single-phase when V/n lies outside (v_l, v_v) at its own temperature T₁; only otherwise does the lever-rule flash run (D-53 item 1).
- **Why it is normally safe:** a stable uniform state is the equilibrium, which is unique at given (S, V, m). Any spurious root of the inversion inside the unstable region lies below the isochore's dome temperature, so V/n is inside the dome there and the state goes to the flash anyway.
- **When it can fail:** the argument assumes a well-behaved equation of state. A model with extra loops or artifacts (spurious high-density phases of multiparameter equations, see FM-2; cubic or SAFT models far from their fitted range; below the triple point) could give an inversion root T₁ at which V/n is outside the computed (v_l, v_v) although the stable state is two-phase. The body is then reported as a metastable uniform fluid: wrong T, p and μ, `quality` returns `nothing`, and no error is raised.
- **Detect:** for a doubtful state, compare with Clapeyron's own flash at the same (u, v), or check that p is not below p_sat(T) for a liquid-like or above it for a vapour-like state.
- **Do:** stay within the models validated in `test_twophase.jl` (IAPWS-95 water, Span–Wagner CO₂) and their normal range. Test a new model in the dome before trusting it.
- **Revisit when:** a new pure-fluid model is used for two-phase work, or a case runs deep below T_c near the triple point.

### FM-2 Spurious saturation pairs in multiparameter equations — silent if the guards miss
- **What:** the Span–Wagner CO₂ equation has equal-p, equal-μ volume pairs at high density that are not the real liquid–vapour equilibrium (D-53 item 3). A saturation Newton solve started too far away converges to one, and the phase test then misreads dome states.
- **Guards in place:** a warm start only within 0.1 % in T of the last solve, otherwise Clapeyron's `saturation_pressure` (ancillary-equation start); the polished pair must stay near its start, be distinct, and have ∂p/∂v < 0 in both phases.
- **When it can still fail:** if Clapeyron's own `saturation_pressure` returns a spurious pair, the guards accept it (they only check closeness to that start).
- **Detect:** v_l, v_v or p_sat far from tables or from the model's ancillary equations.
- **Revisit when:** a new model shows misread dome states, or Clapeyron's saturation solver changes.

### FM-3 Near-critical states inside the dome — silent
- **What:** very close to T_c the two saturated volumes merge (the polish needs ln v_v − ln v_l > 1e-6), and the saturation solve can fail. The flash treats a temperature where saturation fails as "above the dome", so it finds no root in that thin cap, and the phase test keeps the uniform-fluid state: wrong T and p for a state that is really two-phase.
- **Detect:** T within a fraction of a percent of T_c and `quality` returns `nothing`.
- **Do:** avoid operating points inside the dome right at the critical point; sCO₂ cycles above T_c are unaffected.
- **Revisit when:** a case operates in the near-critical two-phase cap.

### FM-4 Below the triple point: no solid phase — silent
- **What:** the saturation range starts at T_lo = 0.4 T_c, not at the triple point (D-53 item 2). Between T_lo and the triple point, dome states are liquid–vapour mixtures of a metastable (supercooled) liquid; solids are never formed.
- **Where it bites:** water between 258.8 K and 273.16 K (supercooled water, no ice); CO₂ between 121.7 K and 216.6 K (liquid CO₂ where dry ice would form).
- **Below T_lo:** a dome state is treated as single-phase (wrong), or the inversion refuses it (loud).
- **Do:** keep pure fluids above their triple point.
- **Revisit when:** freezing, sublimation or cryogenic CO₂ matters.

### FM-5 Mixtures are always single-phase — silent, will fail
- **What:** a real-fluid mixture (e.g. GERG-2008 air) is evaluated as one uniform phase everywhere, with no phase check at all (D-51 item 5, D-52). Inside its two-phase envelope it is a metastable or unstable uniform fluid: wrong T, p, μ, with no error.
- **Where it bites:** cryogenic air or nitrogen–oxygen separation, condensing flue gas or wet gas, any compressed or subcooled mixture near its bubble point.
- **Mixture isentropic states** may converge on a metastable branch when the end state is two-phase (or fail loudly with a `DomainError`).
- **Do:** check doubtful mixture states against Clapeyron's flash outside the model. A mixture phase check is a planned option (D-51 review); a mixture flash is future work (DESIGN §15).
- **Revisit when:** any process or cryogenic model uses a mixture near or inside its envelope.

### FM-6 Phase-selective draws keep pumping when the selected phase runs out — silent
- **What:** `IdealPump(mdot; draw = :liquid)` takes saturated liquid while the body is two-phase. If the body dries out (x → 1) it becomes superheated vapour, and from a single-phase body every draw takes the body's own state, so the prescribed ṁ continues as vapour. Likewise `draw = :vapor` from a body that fills with liquid. A real trap or a dry pump would stop or cavitate; here nothing flags it, and what is drawn jumps at the boundary.
- **Detect:** watch `quality` on the drawing body; a value at 0 or 1, or `nothing`, during a run means the draw is no longer phase-selective.
- **Do:** size inventories so the body stays two-phase, or switch the pump off (`Switched`) on a quality condition.
- **Revisit when:** a model needs a trap that closes, level control, or pump cavitation.

### FM-7 Steady state of a two-phase body with a free volume is indeterminate — loud (untested)
- **What:** inside the dome the Hessian is singular along the tie line: at fixed T and p, any quality is an equilibrium. If a two-phase body's volume is also free (on a piston or a pressure reservoir) and nothing else fixes its inventory split, `steady_state` faces a singular Jacobian. This follows from the phase rule; it has not been run.
- **Do:** fix the volume (a rigid vessel, as in the Rankine test) or the mass, or use `simulate` to a steady state.
- **Revisit when:** a steady-state solve of such a body is needed.

### FM-8 Entropy datum and the relative tolerance — may degrade accuracy
- **What:** the integrator's relative tolerance multiplies |S|, whose zero is arbitrary (D-20 §1). With a real-fluid datum (IAPWS-95: s = 0 for saturated liquid at the triple point) a body with large |S| gets a looser absolute error in S than one near the datum. D-20 flagged this for review "when real-fluid reference states arrive"; they now have.
- **Detect:** closure and results that change when `reltol` is tightened.
- **Do:** tighten `reltol` for long real-fluid runs where entropy accuracy matters.
- **Revisit:** now due; consider scaling the relative tolerance by capacity as the absolute one is.

---

## 2. Mixtures, reactions and trace species

### FM-9 Trace floor: μ of very dilute species is not ideal — silent
- **What:** below 1e-13 of a body's mass (the `trace` keyword), the mixing term is smoothed (D-41, D-42, D-44). Conservation, energy, entropy and first-order nuclear rates stay exact. But μ of the dilute species is too high (about 4.6 RT at 100× below the floor), so any flux or rate driven by it is wrong: permeation and diffusion (about 25× too small for a 1e-15 vs 1e-14 pair), partition equilibria, mass-action rates. Its amount can read about −1e-12 of the body's mass.
- **Where it bites:** discharge-level or environmental tritium (~10 kBq/L ≈ 2e-13 of the mass of water), any trace species whose transport matters.
- **Detect:** a species' mass fraction near or below 1e-13 of its body. `closure` dilutes a trace species' loss in the body's scale; check `closure(cs, sol).species[:X].abs` instead (D-43).
- **Do:** lower `trace` on that body to ~1e-15 (not lower: round-off), or put the dilute species in its own smaller body; read slightly negative amounts as zero.
- **Revisit when:** a fusion case needs μ below 1e-15 of a body (that needs log-amount states, a design change).

### FM-10 Constant c_v per species — silent at high temperature
- **What:** ideal-gas mixtures use constant heat capacities (D-39 item 1). Flame temperatures, rocket chambers and other high-temperature states carry the error of that approximation; dissociation is absent unless species are added.
- **Do:** use the real-fluid backend where it covers the species, or treat high-temperature mixture results as approximate.
- **Revisit when:** NASA polynomials are added to the mixtures.

### FM-11 Equilibrium chambers book the combustion loss on the injectors — silent misattribution
- **What:** in a well-stirred lumped reactor at A ≈ 0, the reaction element produces almost nothing and the combustion irreversibility appears on the feeds (their generation includes the jump to the equilibrium μ) (D-48 item 1). The total is exact; the per-element loss map is not what a designer expects.
- **Do:** read feed generation into a reacting body as combustion loss.
- **Revisit when:** a finite-rate or two-zone chamber is built (D-48 future work).

---

## 3. Flow, nozzles and rockets

### FM-12 Nozzles are always choked — silent
- **What:** `Nozzle` has no back-pressure check (D-48 item 4). An unchoked or strongly over-expanded nozzle (shocks inside) is computed as choked with the isentropic exit state.
- **Do:** check p_e against the ambient pressure yourself; treat sea-level results of high-expansion nozzles with care.
- **Revisit when:** a nozzle can unchoke in a model (start-up, shutdown, throttling).

### FM-13 No elevation head in streams — silent, will fail for vertical loops
- **What:** gravity is a momentum source on moving bodies (`GravityField`, D-34, flagged). Streams carry μ, not μ + gz, so hydrostatic and elevation heads in flow networks are absent.
- **Where it bites:** natural circulation, tall balance-of-plant loops, any pump sized against a static head.
- **Do:** add the head as a pressure offset by hand, or model the column explicitly.
- **Revisit when:** elevation matters in a flow network (D-34's revisit list).

### FM-14 Strongly coupled flow networks may not solve to steady state from a cold start — loud
- **What:** `steady_state` can fail on a cold start of a strongly coupled network (e.g. a 2 × N-cell heat exchanger) (D-26 §3).
- **Do:** warm start with a transient `simulate`, then call `steady_state` from that state.

### FM-15 Stream entropy production is not sign-exact in floating point — reporting noise
- **What:** a stream's Ṡ_gen is ≥ 0 in exact arithmetic but can read about −1e-12 (relative) per element, because the advected terms carry |s| with an arbitrary reference (D-26 §4).
- **Do:** treat tiny negative stream production as zero; it is not a second-law violation.

---

## 4. Mechanics and electrical

### FM-16 Coulomb friction has no true sticking — silent creep
- **What:** `CoulombFriction` is tanh-regularized (D-35 item 5). A body that should be held by static friction creeps at w = w_ε atanh(F_applied/F_c), which is small but never zero.
- **Do:** make w_ε small relative to the speeds that matter; do not use it to hold a load indefinitely.
- **Revisit when:** a model needs stick–slip or true holding.

### FM-17 No contact potentials between electrical materials — silent
- **What:** the electrical force is φ (D-27 item 9). A circuit joining different conductor materials misses their contact potentials (the jump in electrochemical potential at a material boundary).
- **Do:** keep electrical connections within one material, or add the contact potential as a source by hand.

### FM-18 Mass conductors are isothermal — silent for non-isothermal transport
- **What:** `MassConductance` moves mass at μ without entropy (D-27 item 6). Between bodies at different temperatures it misses the transported heat and thermal diffusion.
- **Do:** use a coupled law (`KnudsenBarrier`, or `OnsagerConductance` with a thermal row) when the two ends differ in temperature.

---

## 5. Solver and numerics

### FM-19 Transformer constraint drift — small, loud only through closure
- **What:** transformers conserve power on the constraint manifold; off it the mismatch is held at integrator tolerance (Baumgarte τ), and the initial charges must satisfy the constraint (D-29 item 1). Inconsistent initial charges put a power error into the run.
- **Detect:** `constraint_residual`, and energy closure.
- **Do:** start from consistent charges (`consistent_state`).

### FM-20 Rodas5P overshoots fast decay to a bound — may push a charge past zero
- **What:** Rodas5P's stability function is negative for hλ ≳ 6 on linear decay (D-41 finding). A charge relaxing quickly to a bound (zero mass, zero amount) can overshoot it on a large step. Harmless with the trace floor; elsewhere it can trigger step rejections or a slightly negative value.
- **Do:** if a run stalls on rejections near a bound, tighten tolerances or limit `dtmax`.

### FM-21 Turning off ledger error control — silent if done
- **What:** with `ledger_control = false`, loss figures can be wrong by percent when the physical states are easy to integrate but a ledger is not (D-29 §11, D-31). Only `closure` shows it.
- **Do:** keep the default (`true`).

---

## 6. Concurrency

### FM-22 Real-fluid materials are not thread-safe — silent
- **What:** a `RealFluid` keeps its warm starts and saturation cache in `Ref`s shared by every body made of it (D-51 item 3; slice 10b caches). Evaluating one material from two threads at once can mix a cached saturation pair from one temperature with a lookup at another, giving wrong saturation volumes without an error.
- **Do:** create one `RealFluid` per thread or per concurrently running model.
- **Revisit when:** models are run in parallel threads.
