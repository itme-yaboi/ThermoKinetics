# DESIGN — VeinikThermo.jl

Status: **revision 2, reviewed** (2026-10-14). Revision 1's open questions have been answered and folded in (§16 records the answers). Slices 1–5 are built (see `STATUS.md`). Decisions with alternatives are logged in `DECISIONS.md` (referenced as D-n).

**Revision 3 additions (2026-10-15):**
- **Mechanical and rotational systems** (§18, D-32): reviewed and accepted with the user's decisions in D-34, and **built as slice 6** (D-35).
- **Rockets** (§19, D-33): still proposed and not built (slice 8).

Between them they add §0.1 rows, two domains (§1), materials (§2.3), monotone friction laws (§5.2), spring branches (§6.1), momentum-carrying streams and a nozzle element (§7), force, torque and gravity sources (§8), and slices 6–8 (§15).

## 0. Design in one paragraph

The state vector is the concatenation of every body's **charges**. Each body's material provides **one potential U(E)**; forces are P = ∂U/∂E and nothing else. Physics elements (conductor laws, converters, stream drivers, reactions, sources) are **pure functions from forces to fluxes**. They never write to the state derivative. A single, small, audited **assembler** turns those fluxes into charge and entropy rates using a fixed bookkeeping rule per element *kind*. The rule conserves every non-thermal charge and energy identically, and it makes entropy production non-negative through the types of the laws. The model is held as a **graph** (Graphs.jl + MetaGraphsNext.jl) for building, validation and queries. It is then compiled to flat, type-stable arrays for the right-hand side. **No ModelingToolkit** anywhere (D-2). Adding a form of motion means declaring a new domain and writing materials and laws that use it. The assembler and solver do not change.

```
Domain (declaration) ─┐
                      ├─> Material: U(E), groups ─> Body (charges = state)
Property backend ─────┘                              │
                                                     ▼
Model graph (bodies, reservoirs, elements; ports as edges) ──compile──> flat layout
                                                     │
Laws: forces ──> fluxes ──> Transfer ──> Assembler (bookkeeping) ──> du/dt (physical) + ledger rates
```

## 0.1 Intended functionality (the end goal)

The library exists so that users can build **practically useful thermodynamic simulations** from a small set of element kinds, with every result carrying its energy and entropy ledger. The target uses are listed below. Each slice should move at least one row forward, and a row counts as delivered only when its **application test** passes. An application test is a model a user would actually build, checked against an independent engineering result.

| Use | Typical models | Built from | Slice | Application test |
|---|---|---|---|---|
| **Transient heat transfer** | multilayer walls, components heating up or cooling down, fins, heat sinks, convection and radiation boundaries, internal heating | solid bodies, thermal conductors (G, radiation), reservoirs, heat sources, `rod!`/mesh helpers | 1 ✅ (radiation and heat sources still to add) | composite wall vs. series resistances, steady and warm-up (`test_applications.jl`) ✅ |
| **Control-mass systems** | gas in a closed vessel, piston–cylinder compression and expansion, calorimetry, closed heating and cooling | coupled bodies U(S, V, m), volume conductor, converters (piston + momentum) | 2 ✅ (isochoric, isobaric), 5 ✅ (converter piston) | rigid-vessel and loose-piston heating vs. closed-system first law (`test_applications.jl`) ✅; piston–cylinder compression by the atmosphere: gas exactly on its isentrope, friction heat = W_atm − ΔU, entropy = Q_f/T ✅ |
| **Control-volume systems** | tank filling and blowdown, mixing chambers, pipes, heat exchangers, open loops | streams, hydraulic drivers, `IdealPump`, labeled flow boundaries, pass-through outlets | 3 ✅ | tank blowdown (remaining gas exactly isentropic) and counter-flow heat exchanger vs. ε–NTU (first-order convergence, Richardson-extrapolated) in `test_applications.jl` ✅; 08 B5 tank filling and C7 mixing in `test_streams.jl` |
| **Power and refrigeration cycles** | Rankine, Brayton, Carnot, Otto; fusion balance-of-plant loops | real-fluid backend (Clapeyron) with two-phase flash, pumps and turbines (converters, including transformer form via DAE), heat exchangers, switched connections | 3, 5 (Carnot ✅: E1, E2 in `test_cycles.jl`), later (real fluids: Rankine, Brayton) | Rankine and Brayton efficiency and net work vs. property-table hand calculations; second-law efficiency and per-component exergy destruction summing to the lost work |
| **Electrical circuits** | RC and RLC networks, Joule heating of conductors, motors and generators, thermoelectric devices | `Electrical`, `FluxLinkage`, capacitors, resistors (kernel B), gyrator converters, coupled conductors | 4 ✅ (R, C, thermoelectric), 5 ✅ (inductors, LC/RLC, DC motor, transformer) | RC charging and the two-capacitor paradox (D1, D2, `test_electrical.jl`) ✅; resistor self-heating with its heat path and loss map ✅; thermoelectric generator efficiency vs. the classical TEG balance and η_max ✅ (`test_applications.jl`) |
| **Mass transfer** | diffusion and permeation through walls and membranes (e.g. tritium permeation), drying, osmosis, thermal diffusion | `Species` domains, mass conductors on μ, coupled conductors | 2, 4 (mass conductors on μ ✅, Knudsen barrier ✅), later (membrane application test) | steady permeation flux through a membrane vs. Fick/Sieverts-type closed form, with the partition jump handled by μ |
| **Process systems** | mixers, splitters and recycle loops, flash drums and separators, compressor trains with intercooling, gas processing and purification, cryogenic plants (e.g. a fusion fuel-cycle tritium plant or cryoplant) | multicomponent streams, real-fluid backend with phase-equilibrium flash, `IdealPump` and compressor/turbine converters, heat exchangers, reactions, switched connections | 3, later | mixer–splitter recycle loop vs. hand mass and energy balance; adiabatic flash drum vs. Rachford–Rice with the backend's K-values; intercooled compressor train work vs. hand calculation, with exergy destruction per stage |
| **Chemical and nuclear** | reacting mixtures near equilibrium, decay heat, tritium decay inventories | internal transformations (reactions), nuclear sources | later | equilibrium composition vs. Gibbs minimization; decay inventory and heat vs. Bateman |
| **Coupled and multi-physics** | thermoelectric, electrokinetic, magnetocaloric, MHD pressure drop | coupled conductors (Onsager–Casimir when B ≠ 0), coupled materials | 4 ✅, later (MHD) | C4, C5, C6 ✅ (`test_coupled.jl`); Hartmann-flow MHD pressure drop |
| **Translational mechanics** (§18) | mass–spring–damper, multi-mass chains, force or base excitation, vibration isolation, pendulums, viscous, drag and Coulomb friction | `Momentum` bodies (`InertialMass`), `Displacement` springs (`LinearSpring`), `PendulumPotential`, spring branches (converters), dampers (momentum conductors, kernel A or B), force sources, gravity, monotone friction laws | 5 (masses, gyrators, prescribed drives), 6 ✅ | free mass–spring–damper vs. the analytic response for ζ < 1, = 1, > 1, with the heat to the environment equal to the initial mechanical energy and entropy Q/T; forced response amplitude and phase vs. the transfer function, and mean power in = mean dissipation; pendulum period vs. 4√(ℓ/g)·K(sin θ₀/2) from small (→ 2π√(ℓ/g)) to large amplitude, and the damped pendulum's heat (08 G1, G2, G5; `test_applications.jl`) ✅ |
| **Rotational mechanics** (§18) | flywheels, shafts, gear trains, rack and pinion, clutches and brakes, motor-driven loads, torsional vibration | `AngularMomentum` bodies (`RotationalInertia`), `Angle` torsion springs, `Transformer` gears, bearing friction, torque sources, switched clutches | 5 (inertias, transformers, DC motor), 6 ✅ | spin-up through a two-stage gear train vs. the reflected inertia J₁ + J₂/r₁² + J₃/(r₁r₂)²; two-inertia torsional frequency; clutch engagement: angular momentum conserved, loss ½J₁J₂Δω²/(J₁+J₂) for any clutch conductance (08 G6–G8; `test_applications.jl`) ✅ |
| **Rockets and propulsion** (proposed, §19) | the rocket equation, cold-gas thrusters, a simplified bipropellant engine (feed, combustion chamber, nozzle) on a test stand and in flight | moving variable-mass bodies, momentum-carrying streams, reactions (slice 7) in an ideal-gas mixture with formation enthalpies, an ideal nozzle element, gravity, drag | 8 (proposed; needs 6 ✅ and 7) | cold-gas thruster Δv vs. ∫ c_eff dm/M along the isentropic blowdown (closed form); static fire: T_c, p_c, c*, C_F and I_sp vs. a hand calculation with the same property model; flight: Δv vs. Tsiolkovsky c_eff ln(m₀/m_f) − g t_b, with momentum and energy closure |

Cross-cutting, for every use:
- **network plotting**: a schematic diagram of any model (bodies, reservoirs and elements laid out with the compact orthogonal algorithm of §17), optionally coloured by results such as temperature, pressure or entropy production;
- transient and steady solutions;
- **loss attribution** (Ṡ_gen and exergy destruction per element);
- energy and entropy closure on every run;
- SI throughout.

A model mixing several rows (a fusion blanket loop: nuclear heating, conduction, coolant streams, heat exchanger, power cycle) must need no new solver code. That is the test of the architecture.

## 1. Forms of motion (domains)

A domain is a **declaration with no physics** (`07` §2):

```julia
abstract type FormOfMotion end
struct Thermal     <: FormOfMotion end              # S [J/K], T [K]
struct Volumetric  <: FormOfMotion end              # V [m³], −p [Pa]
struct Electrical  <: FormOfMotion end              # q [C],  φ [V]
struct Species     <: FormOfMotion; name::Symbol end  # m_k [kg], μ_k [J/kg]
struct Momentum     <: FormOfMotion end              # p [kg·m/s], w [m/s] (one axis; vector form later, D-29)
struct AngularMomentum <: FormOfMotion end
struct FluxLinkage <: FormOfMotion end              # Λ [Wb], i [A] (inductor)
struct Displacement <: FormOfMotion end             # x [m], F [N] (spring; §18.1)
struct Angle       <: FormOfMotion end              # θ [rad], M [N·m] (torsion spring; §18.1)
# … polarization, magnetic, surface, gravitational height

# traits (the whole "domain API")
charge_symbol(d), charge_unit(d), force_symbol(d), force_unit(d)
ncomponents(d)              # 1, or N for vector charges
producible(d) = false       # only producible(::Thermal) = true
```

- **No sign flags exist.** The volumetric force is −p because ∂U/∂V = −p for every material (non-negotiable 3).
- `producible` is used in one place: the assembler may write a production term only to a `Thermal` slot.
- Domains are values, so `Species(:T2)` and `Species(:He3)` are distinct domains with no new types needed.
- A body has **at most one** `Thermal` domain (`07` §2).
- Kinetic and inertial forms are never both present. There is only `Momentum`.

## 2. Materials: the state function

### 2.1 Interface

```julia
abstract type Material end
domains(mat)::Tuple                 # ordered charges of a body made of this material
internal_energy(mat, E)             # U(E), E an AbstractVector (generic in eltype for AD)
forces(mat, E)       = ∇U           # default: ForwardDiff.gradient
eos_matrix(mat, E)   = ∇²U          # default: ForwardDiff.hessian, then symmetrized
parts(mat)                          # coupling groups (§2.2)
mass_slots(mat)                     # indices of Species charges (for the fixed-mass block)
valid_state(mat, E)::Bool           # e.g. m > m_min, T within backend range
```

Derived quantities, written once for all materials:
- `capacity_matrix(mat, E; hold=:forces | :charges)` works on the **fixed-mass block** only. The full Hessian is never inverted (D-8).
- `is_stable(mat, E)` checks that the fixed-mass block is positive definite.

Parameters that are not charges stay as parameters. The mass of an incompressible solid is the main example: such a body has no `Species` domain.

### 2.2 Coupled, uncoupled, coupling groups

Coupling is a property of the potential (`07` §1). It is expressed by **composition** (D-12):

```julia
gas_and_wall_charge = Composite(IdealGas(He; ...),      # part 1: (Thermal, Volumetric, Species(:He))
                                LinearCapacitor(C))     # part 2: (Electrical,)
```

- A `Composite` is a list of parts. Each part is a material over a **disjoint** subset of domains. U = Σ U_part.
- Each part is one **coupling group**. ∇²U is block-diagonal with one block per part, and the Jacobian sparsity pattern follows it.
- `forces(::Composite)` evaluates each part only on its own slice of E. Uncoupled separability (A7) therefore holds by construction.
- At most one part may contain `Thermal`.
- A body with no coupled part is a pure accumulator (Veinik §11).
- Optional DOFs (e.g. a gas that gains `Electrical` as a plasma) are declared by the material and enabled when the model is built. The DOF set is fixed for the life of a compiled system.

### 2.3 Analytic materials (closed-form U(E), no inversion)

| Material | Charges | U | Notes |
|---|---|---|---|
| `ConstantCapacitySolid(m, c; T_ref)` | S | m c T_ref·exp(S/(m c)) | T = T_ref e^{S/mc}; S₀ = m c ln(T₀/T_ref) |
| `IdealGas(R, c_v; refs)` | S, V, m | m c_v T(s, v), T = T_r (v/v_r)^{−R/c_v} exp((s − s_r)/c_v) | homogeneous of degree 1; μ = g = h − Ts with consistent references |
| `ElasticVolume(C; V0, p0)` | V | −p0(V − V0) + (V − V0)²/(2C) | uncoupled volumetric store (Veinik's linear one-DOF body); not a gas (D-22) |
| `LinearCapacitor(C)` | q | q²/(2C) | |
| `InertialMass(m)` | p⃗ | \|p⃗\|²/(2m) | force = velocity |
| `LinearSpring(k; x0)`, `TorsionSpring(k; θ0)` (§18.2) | x or θ | ½k(x − x₀)² | force = tension or torque |
| `PendulumPotential(m, g, ℓ)` (§18.2) | θ | mgℓ(1 − cos θ) | not convex for \|θ\| ≥ π/2 (an inverted pendulum is unstable); `stiffness_floor` = mgℓ |
| `Moving(material; m_dry)` (proposed, §19.2) | the material's charges and p | U_int(E) + p²/(2M), M = m_dry + Σ m_k | variable-mass kinetic body; μ_k = μ_k,int − w²/2 |

### 2.4 Property backends and the Newton inversion

Real-fluid models are explicit in (T, ρ). **Clapeyron.jl is the first backend** (AD-friendly Helmholtz); CoolProp comes later through a wrapper (D-7).

1. **Reduce to one unknown using homogeneity.** From the charges (S, V, m₁…m_n) form v = V/m, s = S/m and the composition. Only **T** is unknown.
2. **Solve s(T; v, x) = s_target** with safeguarded Newton.
   - ∂s/∂T|_v = c_v/T > 0, so the function is monotone in a stable single phase.
   - The bracket comes from the backend's validity range.
   - Warm start from the last converged T (a per-body cache).
   - Fall back to bisection if a Newton step leaves the bracket.
   - Converge to |ΔT| ≤ 4 eps·T.
3. **Forces from the Helmholtz function** a(T, v, x): T; −p = ∂a/∂v; μ_k = ∂A/∂m_k; U = m(a + Ts).
4. **Derivatives via the implicit-function theorem.** `forces` on Dual inputs solves on the values, then propagates the partials through T = T(s, v, x). The Jacobian is exact, and AD never runs through the iterations.
5. **∇²U** comes from the second derivatives of a, through the change of variables. A4 checks its symmetry *before* symmetrizing, which tests the backend's consistency.
6. **Guards.** `valid_state` rejects m < m_min (B5 starts from p₀ = 1 Pa) and out-of-range states. The integrator's `isoutofdomain` hook uses it to reject steps.
7. **Two-phase** (water/steam, later): backends must provide an (s, v, x)-flash behind the same `forces` interface.
8. **Isentropic state at a given pressure** (needed by the ideal pump, §7.2): solve p(s, v) = p_target for v at fixed s. This is a 1-D monotone root find (∂p/∂v|_s < 0); it is closed-form for the ideal gas.

## 3. Bodies and reservoirs

```julia
Body(name, material, E₀)                    # E₀ in charges; helpers accept (T₀, p₀, …) and invert
Reservoir(name, forces; material = nothing) # fixed forces (Dirichlet), infinite capacity
```

- A **body** owns a contiguous slice of the physical state: its charges, in `domains(material)` order.
- A **reservoir** has no dynamic state. If a stream draws from it, it needs a material so that s, h and μ at its (T, p) are consistent.
- Every reservoir has **ledger** entries (§10) for the charge, energy and entropy it receives. "Bodies + reservoir ledgers" then form a closed account.

## 4. Model graph, element contract, and Transfer kinds

### 4.1 The model graph (D-17)

The model is a `MetaGraphsNext.MetaGraph` over an undirected `Graphs.SimpleGraph`, with **two vertex kinds**:
- **node vertices:** bodies and reservoirs;
- **element vertices:** conductors, converters, streams, pumps, sources, reactions, switches.

Each **edge is a port**: it links an element to a node and carries metadata (role `:a`, `:b`, `:h`, `:in`, `:out`, …; domains; direction).

The bipartite element–port form is used because:
- several elements can join the same pair of nodes (a conductor *and* a stream), which a simple graph cannot store as edges;
- some elements have three or more ports (a conductor with an internal node to a wall, a multi-port converter).

The graph serves building (`add!`, `connect!`), validation, topology queries and export (plotting, reports). Examples of queries:
- connected components (isolated sub-networks for steady-state constraints and for invariant tests);
- "which elements touch this body";
- the Jacobian sparsity pattern.

The same graph also carries **display metadata** for network plotting (labels, component kind, box size, optional pinned positions; §17). Display data never reaches `compile` or the physics.

**The graph is not used in the hot loop.** `compile` lowers it to flat arrays (§11). The compiled form is plain data with no generated code, so the core never depends on a symbolic layer and can always be inspected (the lesson from the earlier library's ModelingToolkit dependency).

### 4.2 Element contract

```julia
evaluate(element, forces_view, t) -> Transfer     # pure: forces in, fluxes out
apply!(du, ledger_rates, transfer, layout)        # assembler-owned; one method per Transfer kind
```

| Transfer kind | Produced by | Each non-thermal charge | Energy | Entropy written |
|---|---|---|---|---|
| `ConductorTransfer` | conductors (thermal, charge, volume, mass, coupled), surface exchange | conserved | exact: Σ heat deliveries = released power | Σ Q_i/T_i ≥ 0 (§5) |
| `ConverterTransfer` | ideal converters | **not conserved** (the one exception) | exact: Pᵀ W P ≡ 0 | none |
| `StreamTransfer` | streams, ideal pumps, flow boundaries, pass-through outlets; nozzles (proposed, §19.4) | conserved (momentum too, when both ends carry it) | exact; pump work booked as external | advected + generation ≥ 0 |
| `InternalTransfer` | reactions (inside one body) | per stoichiometry (elements) | exact (dU = 0) | A·J/T ≥ 0 |
| `SourceTransfer` | external sources (heat, charge, nuclear; force, torque and gravity, §18.6) | external, ledgered | external, ledgered | Q̇/T booked as **boundary** flow |

The five `apply!` methods are the *only* code that writes `du` (D-3). A wrong law gives a wrong flux, but it cannot create charge, break energy conservation, or produce negative entropy.

## 5. Conductors

### 5.1 One bookkeeping rule, two heat-distribution kernels

A conductor joins a and b across a declared tuple of domains, all of which must exist on both ends (n ≤ m).
- Non-thermal rows carry fluxes J_k from a to b.
- They release power D = Σ_k (P_k^a − P_k^b) J_k.
- Heat is delivered to a set of thermal targets i ∈ {a, b, h} as Q_i (positive into i).

The **bookkeeping** (`apply!` for `ConductorTransfer`) is the same for every conductor:

- E_k: a −= J_k, b += J_k
- S_i += Q_i / T_i for every thermal target
- ledger: Ṡ_gen = Σ_i Q_i / T_i
- the last Q_i is computed as D − Σ(others), so energy closes to round-off by construction

**Energy:** Σ_i Q_i − Σ_k P_k^a J_k + Σ_k P_k^b J_k = D − D = 0.

What differs between conductors is *how the Q_i are distributed*. Two kernels exist, and each is proven to give Σ Q_i/T_i ≥ 0.

**Kernel A — direct exact form** (`07` §3.1, generalized; D-4). Weights w_i ≥ 0 with Σw_i = 1 over the targets.
- Forces: X_Q = 1/T_b − 1/T_a and Y_k = ΔP_k · Σ_i w_i/T_i. Law: J = M x with x = (X_Q, Y).
- Heat: Q_a = −J_Q + w_a D, Q_b = J_Q + w_b D, Q_h = w_h D.
- Proof: Σ Q_i/T_i = J_Q X_Q + Σ_k J_k Y_k = xᵀ M_sym x ≥ 0 at any ΔT.
- `Split(w)` is (w_a, w_b) = (1 − w, w). `Into(h)` is w_h = 1 (deposition straight into a wall or reservoir).
- **Used for:** pure thermal conductors (the energy form Q̇ = GΔT), and coupled conductors with a thermal row (thermoelectric, thermodiffusion, Knudsen barrier).

**Kernel B — internal thermal node, eliminated analytically** (`07` §3.2). This is the **default for conductors with non-thermal dissipation** (Joule, viscous, diffusion).
- A zero-capacity node m has thermal legs of conductance G_i ≥ 0 (at least one > 0) to its targets: a, b and/or a wall h.
- The non-thermal law is J = L ΔP, where L has a PSD symmetric part. Engineering laws map exactly (Ohm J = Δφ/R, Hagen–Poiseuille mass conductance on μ). L does not depend on T_m, so D = ΔPᵀ L_sym ΔP ≥ 0 is explicit.
- Node balance: Σ_i G_i (T_m − T_i) = D, which gives **T_m = (Σ G_i T_i + D)/Σ G_i** (explicit, no DAE). Then Q_i = G_i (T_m − T_i).
- Proof: Σ Q_i/T_i = D/T_m + Σ_i Q_i (1/T_i − 1/T_m). Both terms are ≥ 0, because Q_i has the sign of T_m − T_i.
- The ledger reports the two parts separately: "dissipation at the node" and "conduction from the node". Losses are thus attributed where they physically occur.
- The legs also conduct heat between a and b through the element (a resistor conducts heat). That is real physics, and it is accounted for exactly.
- Nonlinear legs (radiation) would make the node balance a monotone scalar equation solved by safeguarded 1-D Newton. That is allowed later; slice 4 uses linear legs.

**Choice rules:**
- Pure thermal → kernel A, energy form.
- Non-thermal dissipative → kernel B by default. Its leg conductances are required arguments; omitting them is an error that names the alternatives. `Split(w)` or `Into(h)` (kernel A) can be chosen explicitly when the node's conductances are unknown or irrelevant. `Into(h)` is the G_h → ∞ limit of kernel B.
- Coupled with a thermal row → kernel A with weights on a and b only, since D can be negative (thermoelectric generation).
- A target without `Thermal` is a build error.

### 5.2 Conductance types (second law by type)

```julia
Conductance(B, W)     # M = B Bᵀ + W;  W = −Wᵀ exactly; W ≠ 0 only via a MagneticField context
ScalarConductance(g)  # g ≥ 0 asserted at construction (errors, never clamps)
```

- **Exact single-row mappings** (kernel A):
  - M_QQ = G·T_a T_b gives Q̇ = G ΔT.
  - Radiation: M_QQ = εσF·T_aT_b(T_a + T_b)(T_a² + T_b²) gives Q̇ = εσF(T_a⁴ − T_b⁴).
- **Coupled data:** a law supplies a factor B, so that M = BBᵀ + W, at the current state. Ṡ_gen = ‖Bᵀx‖² ≥ 0 holds exactly whatever state dependence is used.
  - **Calibrated laws** (`ThermoelectricConductance`, `KnudsenBarrier`) choose that dependence so their engineering law holds exactly at any finite ΔT: J_e = (Δφ + αΔT)/R with Π = αT̄; J_m = 0 exactly at p/√T equal.
  - **Constant matrices** (`OnsagerConductance`) are validated symmetric PSD. Coefficient data from other bases are mapped by **congruence** M = A L Aᵀ (not only diagonal), which preserves symmetry, PSD and entropy production: `barrier_mass_basis` ((J_v, J_s) → (heat, mass)), `electrokinetic_mass_basis` ((J_e, J_v) → (charge, mass)), `thermoelectric_from_veinik` ((277) → R, K, α). See D-27.
- **Gauge invariance:** only differences ΔP_k enter.
- **Monotone single-row laws (§18.4):** J = f(ΔP) with f odd, nondecreasing and f(0) = 0, so D = ΔP·f(ΔP) ≥ 0 by type. Kernel B's proof needs only D ≥ 0 and D independent of T_m, so these use kernel B unchanged. They come only as closed forms, so the sign stays a property of the type: `QuadraticDrag(c)` (J = c\|Δ\|Δ) and `CoulombFriction(F_c; w_ε)` (J = F_c tanh(Δ/w_ε)).

### 5.3 Special conductors (same kernels)

- **Thermal conductor** (`G`, radiation): kernel A, thermal row only.
- **Surface exchange** (Robin): a conductor to a reservoir.
- **Volume conductor** (damped piston or membrane): `Volumetric` row on −p, kernel B (viscous loss at the seal or membrane).
- **Mass conductor** (barrier diffusion): `Species` rows on μ. Fick-type data are converted at construction.
- **Resistor**: `Electrical` row, kernel B. C3, D1 and D2 use `Into(reservoir)` or a node whose only leg goes to the reservoir.

## 6. Converters

### 6.1 Ideal gyrator-form converters (now)

An ideal converter is a constant antisymmetric coupling: F = W P with W = −Wᵀ, so Pᵀ W P ≡ 0. That means power is conserved and nothing is produced (D-6).
- Two-port form: J₁ = r·P₂ and J₂ = −r·P₁.
- **Piston:** `Volumetric` ↔ `Momentum`, with r = A_p.
- **Motor:** `FluxLinkage` ↔ `AngularMomentum`.
- **Spring branch (§18.3):** `SpringBranch()` (`Momentum`, `Momentum`, `Displacement`) and `TorsionBranch()` (`AngularMomentum`, `AngularMomentum`, `Angle`). These are the mechanical twins of `InductorBranch`.

Converters are explicit (an ODE) because one side carries a kinetic DOF. Real converters are an ideal converter plus loss conductors.

### 6.2 Algebraic unknowns: junctions and transformers (mass-matrix DAE; built in slice 5)

The state becomes `[charges; algebraic unknowns; ledgers]` with mass matrix
𝕄 = diag(1, 0, 1), solved by Rodas5P (or FBDF) with a dense Jacobian. The algebraic rows and
unknowns are physical states for the solver (in its error norm, in steady solves), not ledgers.
`consistent_state` solves the algebraic unknowns before integration starts.

- **Zero-capacity junctions** (`Junction(name, domains…)`): a massless mechanical node, or an
  electrical node with no capacitance.
  - Its forces are algebraic unknowns. Its charge-target rows accumulate the net inflow, so its
    row is Kirchhoff's current law, Σ J = 0.
  - This is index 1 (as originally designed) whenever the connected laws depend on the junction
    force.
  - No thermal or species junctions: a zero-capacity node cannot hold entropy, and mixing needs a
    thermal DOF.
- **Transformers** (`Transformer(r; a, b, τ)`): F_a = λ and F_b = −rλ under P_a = r P_b.
  - **Correction (D-29):** between two storing bodies this constraint is **index 2**, not 1. It
    involves only the bodies' states, like a rigid link between two masses, so λ appears only
    after one differentiation.
  - It is solved in **index-reduced form**: λ's row is d/dt(P_a − r P_b) + (P_a − r P_b)/τ = 0.
    The derivative uses the Hessian rows of the end bodies' U (by AD) dotted with their charge
    rates, formed in a second assembly pass, with Baumgarte stabilization at time constant τ.
  - Power is λ (P_a − r P_b): zero on the constraint. Off it, the mismatch is bounded by the
    stabilized drift (integrator tolerance). `constraint_residual` reports it.
  - The initial charges must satisfy the constraint; `simulate` refuses them otherwise.
- **Converters with a kinetic port** stay explicit (§6.1). An ideal gear between two inertias is
  index 2 for the same reason; model one combined inertia, or use a `Transformer`.

## 7. Streams, pumps and flow boundaries

### 7.1 Streams between bodies

```julia
Stream(u, d, driver)   # driver: ṁ = κ(state)·(p_u − p_d), κ ≥ 0 by type (laminar, turbulent)
```

`StreamTransfer` `apply!` (`07` §4):
- Upwind by sign(ṁ).
- The bundle is ṁ_k = ṁ y_k^u, the entropy ṁ s_u, and momentum if kept. V is never advected into a fixed control volume.
- u loses ṁ_k and ṁ s_u; its energy loss is ṁ h_u (Euler's relation).
- d gains ṁ_k and entropy (ṁ h_u − Σ μ_k^d ṁ_k)/T_d.
- Ṡ_gen = ṁ[s_u(T_u − T_d) + Σ y_k(μ_k^u − μ_k^d)]/T_d ≥ ṁ v_u (p_u − p_d)/T_d ≥ 0.

Hydraulic drivers **only set ṁ** (no double counting).

### 7.2 Ideal pump (work-supplied flow, D-13)

```julia
IdealPump(u, d, ṁ_spec)   # ṁ_spec: ṁ(t) or a law; flow from u to d regardless of Δp
```

- The upstream fluid is taken isentropically from p_u to p_d. Its delivered state is (s_u, p_d), found by inversion §2.4-8.
- The bundle delivered to d carries h* = h(s_u, p_d).
- External shaft work Ẇ = ṁ(h* − h_u) is booked in the ledger as **external work**. For a liquid this is ṁ v Δp; for a gas it is the exact isentropic enthalpy rise.
- Inlet-face generation is computed against the delivered state at p_d, so it is ≥ 0 by concavity.
- Ẇ < 0 (flow down the gradient) is an ideal turbine returning work. It is equally valid, and equally ledgered.
- This makes pumped loops possible now, without the shaft/momentum converter. The physical pump (a converter from a motor plus loss conductors) remains the long-term model; `IdealPump` is its ledgered idealization.

### 7.3 Flow boundaries (labeled; unlabeled is an error)

```julia
FlowBoundary(d, material, T_in, p_in, ṁ(t); supply = :passive | :pumped)   # no default
PassThroughOutlet(u, ṁ_spec)    # external outflow at the upstream state, no generation
```

- `:passive`: the fluid arrives at (T_in, p_in).
  - Every evaluation checks p_in ≥ p_d. A violation is an **error**, not a warning, and it names the boundary and the time.
  - Inlet-face generation is internal production (the C7 mixing entropy).
- `:pumped`: an ideal pump between the boundary state and d, exactly as in §7.2. The external work is ledgered.
- Advected entropy at boundaries is booked as **boundary** flow.

### 7.4 Notes

- The ṁ = 0 crossing: the right-hand side is continuous, but its Jacobian is not.
- Inertance-driven ṁ is deferred until the momentum domain exists (accepted).
- With kinetic or gravitational DOFs kept, μ = ∂U/∂m at fixed momentum and height. This is automatic.
- **Momentum-carrying streams (proposed, §19.3).** When both ends hold `Momentum`, the bundle adds ṁ w_u. Σ P·J then gives ṁ(h + w²/2), and the generation gains ṁ(w_u − w_d)²/(2T_d) ≥ 0, which is zero between rigidly linked bodies.
- **Ideal nozzle (proposed, §19.4):** a stream-kind outlet that turns chamber enthalpy into jet velocity. Its bundle carries momentum, so thrust follows from momentum conservation.

## 8. Internal transformations and sources

- **Reactions** (`InternalTransfer`):
  - J_r = M_r·(A/T) with M_r PSD; dm = νJ; dS = A·J/T ≥ 0; dU = 0.
  - ν is validated against the atom matrix.
- **Nuclear processes** (`SourceTransfer`):
  - A prescribed rate with ν validated for nucleon number and charge.
  - Q̇ is deposited as Q̇/T at the local T and booked as external.
- **External sources** (heating, cooling, charge injection): `SourceTransfer`, booked as boundary flow. Their entropy can have either sign, and they never touch Ṡ_gen.
- **Mechanical sources (§18.6, built in slice 6):**
  - `Source(domain => rate(t))` prescribes a non-thermal charge inflow into a body: a force into `Momentum`, a torque into `AngularMomentum`, a current into `Electrical`. Its energy P·rate is booked as external work, and the charge as received (negatively) by the source. Heat sources are not this element and are not built yet.
  - `GravityField(g)` is a momentum source, −M g on a body, where M is the body's own inertial mass. Its work is ledgered. It is a special-purpose element, flagged to revisit (D-34).

## 9. Switched connections (built in slice 5: Boolean gates, D-29)

`Switched(element, schedule)` scales the element's `Transfer` by a gate g ∈ {0, 1}, held as a parameter.
- g·M stays PSD and g·W stays antisymmetric, so every guarantee survives switching.
- Switch times use `PresetTimeCallback` or `ContinuousCallback`, with reinitialization at each switch.

## 10. Ledgers (separable) and diagnostics

Ledgers are pure quadratures that never feed back into the dynamics (D-10):
- Ṡ_gen per element, with the parts split for kernel B and streams;
- per reservoir: charge, energy and entropy received;
- per source, pump and boundary: energy, work and entropy delivered.

**Separability:**
- The ODE state is u = [u_phys; q_ledger]. The right-hand side computes physical rates and ledger rates in the same pass, into separate blocks.
- **Error norm:** by default the integrator's `internalnorm` also covers the ledgers (`ledger_control = true`, D-31, superseding the exclusion in D-10). Their absolute tolerances are scaled like the states they account for (§12.5).
- **Steady state:** solves only f_phys(u_phys) = 0. Ledger *rates* are reported at the solution; Ṡ_gen [W/K] per element is the steady-state loss map.
- **Jacobian:** ledger columns are structurally zero (no feedback), and the sparsity prototype says so.

**Diagnostics reported on every transient run:**
- **Energy-balance closure residual**: r_E(t) = U(t) − U(0) − ∫(boundary energy + external work) dt, plus its value relative to the energy throughput. It sits at integrator tolerance (§12.3). It is reported, not assumed.
- **Entropy closure**: S(t) − S(0) − ∫ boundary entropy dt − Σ S_gen.
- **Charge closure** per domain (round-off).
- The loss ranking T₀ S_gen per element, and contact-rule checks.

## 11. Compilation and the right-hand side

`compile(model_graph)`:
1. **Validation** at build time:
   - port domains exist on both ends;
   - at most one `Thermal` per body;
   - conductor kernel choice and leg conductances given;
   - weights sum to 1 on thermal targets;
   - stoichiometry is conserved;
   - W is antisymmetric;
   - flow boundaries are labeled;
   - stream ends are fluid bodies.
2. **Layout:** physical slots per body and charge, then DAE algebraic slots (slice 5+), then ledger slots.
3. **Grouping:** one concretely typed vector per element type, held in a `Vector{Any}`. Each group's loop is specialized (a function barrier, one dynamic dispatch per group per call), while `CompiledSystem` keeps a single concrete type, so all models share one compiled solver (D-26).
4. **Sparsity prototype** from the graph's incidence and the materials' coupling groups.
5. **Scaling:** a per-state abstol vector: each charge is resolved to the equivalent of a force tolerance, through its capacity at the initial state (for entropy, `atol_T`/A_SS). Charge magnitudes are not used, because entropy has an arbitrary zero (D-20).
6. **Isolated components**, from `connected_components` on the graph after removing reservoirs and boundaries. These identify the conserved totals for steady-state constraint rows.

```
rhs!(du, u, p, t):
  du .= 0
  for each body: P_i = forces(mat_i, E_i)    # once per call
  for each element group, each element:
      tr = evaluate(element, P, t)
      apply!(du_phys, du_ledger, tr)
```

Dual-capable buffers use `PreallocationTools.DiffCache`. The Newton warm-start cache stores values only.

## 12. Solver strategy

### 12.1 Transient

- **Default:** `Rodas5P`, a Rosenbrock method with an exact AD Jacobian and a sparse prototype.
- **Large models:** `FBDF` or `QNDF` with a sparse LU (KLU).
- **Both handle mass-matrix DAEs** (§6.2).
- **Linear invariants are preserved to round-off:** 1ᵀf ≡ 0 implies 1ᵀJ = 0, and that carries through every Rosenbrock stage and every BDF Newton update. Charge conservation is therefore structural through the integrator. In practice "round-off" means eps times the conditioning of each implicit solve, accumulated over the steps: measured up to ~1e-11 relative, so tests use 1e-10, still far below integrator tolerance (D-22).

### 12.2 Steady state

- NonlinearSolve on f_phys = 0, using an AD Jacobian and a transient warm start. The fallback is pseudo-transient continuation.
- In each isolated component, one row per conserved total is replaced by its constraint.

### 12.3 What "energy conserved" means numerically (accepted)

- Exact per connection and per evaluation: the instantaneous residual is at round-off.
- The integrated drift is at integrator tolerance, because U is nonlinear in the state.
- It is reported on every transient run as the closure residual (§10).
- An energy-exact discrete-gradient integrator remains a possible later addition. It would not change the core.

### 12.4 No ModelingToolkit (firm)

The core is built directly on Graphs.jl/MetaGraphsNext.jl, OrdinaryDiffEq and NonlinearSolve. There is no symbolic layer and no generated code. No MTK backend or export is planned (D-2).

### 12.5 Tolerances and ledger control (guidance, D-31)

**What is controlled.** The integrator's error norm covers the bodies' charges, the algebraic
unknowns, and, by default, the ledgers (`ledger_control = true`). `reltol` (default 1e-8)
applies to all of them. The absolute tolerances (`abstol_vector`) are:

| State | Absolute tolerance | Why |
|---|---|---|
| A body's charge (S, V, m, q, p, …) | its domain's force tolerance (`default_force_atol`: 1e-8 K, 1e-4 Pa, 1e-6 J/kg, 1e-8 V, 1e-9 m/s, …) divided by A_kk, i.e. times its capacity, at the initial state | resolves what matters (the force) and ignores charge magnitudes, since entropy has an arbitrary zero (D-20) |
| Junction force | the domain's force tolerance | it *is* a force |
| Transformer flux λ | 1e-10 | |
| Entropy ledger (Ṡ_gen per element, entropy received) | the median body entropy tolerance | same resolution as the states it accounts for |
| Energy ledger (energy received, external work) | that × the median body temperature | the energy equivalent |
| Charge ledger (mass, volume, charge received) | the median body tolerance of that domain | |

Ledgers start at zero, so a flat tiny absolute tolerance would force small steps for nothing.
In the benchmark that cost a 300-cell slab 3× the steps; scaled, the overhead is gone.
`atol_scale` multiplies every tolerance; `atol_ledger = x` replaces all ledger tolerances
with x.

**When ledger control is required.** Keep the default (on) whenever ledger values are used:
losses per element, exergy destruction, heat or work delivered, efficiencies, closure. It is
indispensable when the physical states are easy to integrate but a ledger integrand is not:
- prescribed drives and time-varying reservoirs: V linear in t, while the work integrand p·w
  varies;
- isentropic or adiabatic legs, where S is constant but p changes;
- switched cycles, and long runs that accumulate ledgers;
- differences of large ledgers, such as the lost work Q_h η_C − W.

Measured on the Carnot cycle, turning it off left the work 0.24 % wrong (closure 2.4e-3); on, it
is exact to 4e-9.

**When it may be turned off.** Only when you need just the physical trajectories
(temperatures, pressures, flows) and not the ledgers, or for a very large model where the
ledger integrands track the states (conduction, diffusion) and every millisecond matters. Then
always check `closure(cs, sol)`. `energy_rel` within about 10× `reltol` means the ledgers are
still sound.

**Measured cost** (on vs off, after warm-up):
- conduction, diffusion, gas networks, heat exchangers, RLC: 0–10 %;
- prescribed-drive Carnot: 2.3×, which is exactly the cost of getting its work ledger right.

**Handling tolerances in practice:**
- *Loss of a small difference of large quantities* (lost work, small heat leaks): its relative
  error is about `reltol` × (large quantity)/(difference). Tighten `reltol` so that is small
  enough.
- *Many steps while ledgers are still near zero:* loosen `atol_ledger`, or check that the bodies'
  force tolerances suit the problem (`atol_scale`).
- *A ledger near zero that matters* (a tiny leak): tighten `atol_ledger`.
- *Closure worse than expected:* check that `ledger_control` is on and that discontinuities are
  passed as `tstops`. Then tighten `reltol`.

## 13. Structural guarantees for the non-negotiables

| # | Non-negotiable | How it is guaranteed by construction |
|---|---|---|
| 1 | Non-thermal charges conserved by connections; converters the only exception; internal transformations explicit; production ≥ 0; external sources separate; entropy never conserved | States are charges. Only `apply!` writes `du`. Conductor and stream kernels move the *same* J out of one port and into the other. Converters are a separate kind whose only freedom is an antisymmetric W (or, in the DAE form, a constraint with zero power). Reactions and nuclear processes are separate kinds with stoichiometry validated at build time. Entropy writes come only from: heat deliveries Q_i/T_i, proven ≥ 0 per kernel (§5.1); stream generation (≥ 0 by concavity, with no prescribed uphill flow except through a ledgered pump); reaction A·J/T with PSD M_r; and external sources booked as boundary flow. Integrators preserve linear invariants. |
| 2 | Energy across every connection; exact discrete conductor | Conductors: the last heat delivery is D − Σ(others), and Σ Q_i = D exactly. Converters: Pᵀ W P ≡ 0, or λ(P₁ − rP₂) = 0. Streams: Euler's relation. Pumps: work booked explicitly. Kernel A is `07` §3.1, valid at any ΔT; kernel B is exact by node elimination. Integrated drift is reported (§12.3). |
| 3 | P = ∂U/∂E, volumetric −p, no sign exceptions | Forces come only from `forces(mat, E)` = ∇U (AD, or IFT for backends). Domains have no sign field. |
| 4 | Symmetric static Hessian; L with PSD symmetric part; symmetric unless B ≠ 0 | ∇²U is the Hessian of one scalar. `Conductance(B, W)` forms BBᵀ + W, and W is reachable only through a `MagneticField` context. Kernel B's L uses the same type. |
| 5 | True forces across material boundaries | Laws receive forces from the assembler, never concentrations. The `Species` force is μ = ∂U/∂m. Fick data are converted at construction. |
| 6 | SI only | Float64 SI internally. Unit constants do not exist in the code. Unitful.jl only at constructor boundaries. |

## 14. Testing strategy

- **Kernel identity tests** (before any bodies exist): random forces, temperatures, weights, PSD M, antisymmetric W and leg conductances.
  - Every `apply!` closes energy to round-off.
  - Kernel A: Ṡ_gen = xᵀMx ≥ 0.
  - Kernel B: Ṡ_gen = D/T_m + Σ Q_i(1/T_i − 1/T_m) ≥ 0.
  - Cases include ΔT/T ~ 1 and coupled rows.
- **08 §A** as seeded property tests on random networks:
  - A1: the instantaneous residual at round-off; the closure residual within a tolerance-scaled bound.
  - A2 at round-off.
  - A3: per-element Ṡ_gen ≥ 0 at every evaluation.
  - A4–A7 as the relevant domains arrive.
- **Analytic cases B → C → D → E** in slice order. Each test states its tolerance and why.
- **Closure checks** (§10) at the end of every integration test.
- **Basis-transformation tests** for C4 and C6: transform Veinik's volume-basis coefficients to the (mass, entropy) and (charge, mass) bases. Verify that the transformed matrix stays symmetric PSD and that it reproduces the volume-basis predictions (300) and (323).

## 15. Vertical slices

1. **Slice 1 — ✅ done 2026-10-14** (deviations: D-20; status: `STATUS.md`):
   - `Thermal` domain; `ConstantCapacitySolid`; `Body`, `Reservoir`.
   - Thermal conductor (kernel A, energy form; the general kernel A/B code with identity tests).
   - Model graph + `compile`, `rhs!`, separable ledgers, abstol, sparsity.
   - `Rodas5P` transient; NonlinearSolve steady state.
   - `rod(...)` helper; closure diagnostics.
   - Tests: kernel identity; A1, A3, A6; **B2, C1, C2**; C8; D4.
2. **Slice 2 — ✅ done 2026-10-14** (deviations: D-22): `Volumetric`, `Species`, `IdealGas`, `Composite`, volume conductor, fixed-mass capacities. Tests: A4, A7, B3, B4.
3. **Slice 3 — ✅ done 2026-10-14** (deviations: D-26): streams, hydraulic drivers, `IdealPump`, labeled flow boundaries, pass-through outlet. Tests: B5, C7, and a pumped-loop closure test: a loop with an `IdealPump` and a cooler to a reservoir at T₀. At steady state, pump work = heat rejected, and pump work = T₀ · Σ Ṡ_gen over every element (Gouy–Stodola).
4. **Slice 4 — ✅ done 2026-10-14** (deviations: D-27): coupled conductors (Knudsen barrier, thermoelectric couple); `Electrical` and `LinearCapacitor`; resistor with kernel B. Tests: C3, C4, C5, C6 (with the basis transformations), D1, D2.
5. **Slice 5 — ✅ done 2026-10-14** (deviations: D-29): the mass-matrix DAE path — transformer converters, zero-capacity non-thermal nodes, kinetic DOFs and gyrator converters. Tests: E1, E2, plus DAE consistency (constraint residuals at round-off, Pᵀ F ≡ 0).
- **Plotting track P (independent of the physics slices; not started).** Network plotting per §17: display metadata on the graph, the Freivalds–Glagoļevs layout, edge routing, a Makie-based renderer. It can be scheduled between any two slices when the user asks for it.
6. **Slice 6 — ✅ done 2026-10-15** (deviations: D-35): mechanical and rotational systems (§18, D-32, D-34).
   - Added `Displacement` and `Angle`, `LinearSpring`, `TorsionSpring` and `PendulumPotential` (with the `stiffness_floor` tolerance hook), `SpringBranch` and `TorsionBranch`, `QuadraticDrag` and `CoulombFriction`, `Source` and `GravityField`.
   - Gears, rack and pinion and clutches use the existing `Transformer` and `Switched`.
   - Tests: 08 §G G1–G9, a random-chain A1/A3 property test, rack and pinion, and a motor with a gearbox.
7. **Slice 7 — next: reactions and nuclear sources** (§8, D-11). Scope and tests are to be detailed when the slice starts.
   - It includes a sign-preserving nonlinear (mass-action) rate law alongside the linear affinity law. Only if it presents an issue does it move to future work, with the issue recorded (D-34).
   - Slice 8 needs two things from it: an ideal-gas mixture material with per-species formation enthalpy and standard entropy (so the affinity carries the heat of reaction), and the reaction element.
8. **Slice 8 — proposed (revision 3, not reviewed): rockets** (§19, D-33). It adds `Moving`, momentum-carrying streams and pumps, `Nozzle`, and a displacement (∫ w dt) ledger on `GravityField` for variable-mass bodies. Tests are in §19.5. It needs slices 6 and 7.
9. **Later** (order to be agreed with the user): E3 (rectangular cycle; switching and time-varying reservoirs already exist), E4 (Otto), the Clapeyron backend with Newton inversion, two-phase flash, the momentum vector domain, inertance-driven streams, 3-D rigid-body rotation (§18.7).

## 16. Review record (revision 1 questions)

| Q | Topic | Answer | Where |
|---|---|---|---|
| 1 | Energy-conservation wording | Accepted. Report the closure residual on transient runs. 08 A1 and 07 §7 reworded. | §10, §12.3 |
| 2 | Third-node deposition | Accepted. Combined with the reviewer's suggestion: internal node (kernel B) is the default for non-thermal dissipation; `Into(h)`/`Split(w)` remain explicit options. | §5.1 |
| 3 | Gyrator-form converters; transformer deferred | Accepted, but DAE/transformer converters are required before balance-of-plant cycle work, so scheduled as slice 5. | §6.2, §15 |
| 4 | Flow-boundary inlet against pressure | Unlabeled boundaries are an error; passive violations are errors; an explicit work-supplied `IdealPump` / `:pumped` boundary logs the work. | §7.2–7.3 |
| 5 | Inertance-driven streams | Deferred (accepted). | §15 |
| 6 | C6 basis | Convert to (charge, mass) and test the transformation. 08 C6 updated. | §14 |
| 7 | Ledgers in the state | Accepted, separable: excluded from the error norm and from steady-state solves. | §10 |
| 8 | First backend | Clapeyron.jl. | §2.4 |
| — | ModelingToolkit | **Not used, anywhere.** Graph packages for topology instead. | §4.1, §12.4 |

## 17. Network plotting (planned; not implemented)

A basic network-plotting ability is part of the intended functionality (§0.1). It turns a model into a schematic diagram that can be checked by eye and put in a report. Results can optionally be overlaid: node colour by T or p, element colour or width by Q̇ or Ṡ_gen, which gives a picture of the loss map. It is a separate track (§15 P) and touches no physics.

### 17.1 Reference algorithm (required)

The layout **must** implement the node-placement algorithm in:

> K. Freivalds and J. Glagoļevs, *Graph Compact Orthogonal Layout Algorithm*, Institute of Mathematics and Computer Science, University of Latvia.
> Local copy: `_reference/graphing/Freivalds, Graph Compact Orthogonal Layout Algorith.pdf` (12 pages).

What the paper specifies, and what the implementation follows:
- **Drawing model.** Nodes are rectangles of prescribed minimum size (w_i, h_i) with a minimum gap δ. Edges are orthogonal polylines. Only point-wise edge crossings are allowed.
- **Three phases** (Fig. 1): node placement, then edge routing, then normalization (overlap removal with minimal node movement). The paper's contribution, and the part specified in detail, is node placement.
- **Grid cell size c** (eq. 1), from L_min and L_max of the boxes plus δ. Each node occupies ⌈(w_i + δ)/c⌉ × ⌈(h_i + δ)/c⌉ cells.
- **Edge-length function** (eq. 2): d(v_i, v_j) = d_e(border to border) + (1/20)·min(|Δx_c|/(w′_i + w′_j), |Δy_c|/(h′_i + h′_j)). The second term favours aligned centres.
- **Compaction** (§3.1, eq. 3). It works in one direction at a time:
  - minimize Σ_(i,j)∈E (z_i + w′_i/2 − z_j − w′_j/2)² subject to z_j − z_i ≥ γ·w′_i for every edge of the visibility graph S;
  - round with x′_i = ⌊z_i⌋;
  - γ ≥ 1 leaves slack; γ = 1 gives the tightest packing.
- **Main loop** (Algorithm 1, p. 7):
  - grid 5√|V| × 5√|V|; iterationCount = 90√|V|;
  - temperature T₀ = 2√|V|, cooled by k = (0.2/T₀)^(1/iterationCount);
  - each node moves to the median of its neighbours plus a random displacement ∝ T, placed at the nearest free cell (searching Manhattan radius d + 1), and is swapped with a neighbouring node when it did not move;
  - compaction every 9th iteration, alternating directions, with γ = 3;
  - stage 1 uses unit-size nodes; compaction with the real sizes (`expand = true`) switches to stage 2, where γ decreases to 1 over the last compactions.
- **Starting layouts** (§3.3): random (default in Algorithm 1), force-directed (T = 3, 100 iterations, compaction every 3rd), or BFS (iterationCount = 10√|V|, T = 0.2√|V|).
  - The paper finds force-directed best. Use it once a force-directed placement is available; keep random as the fallback.
  - Use a seeded RNG so layouts are reproducible.

What the paper **does not** specify (it cites other work), and the plan for it.

The other phases, and the compaction solver, come from the paper's references, all now in `_reference/graphing/`:

| Ref. | Paper | File | Used for |
|---|---|---|---|
| [10] | K. Freivalds, P. Kikusts, *Optimum Layout Adjustment Supporting Ordering Constraints in Graph-Like Diagram Drawing*, Proc. Latvian Academy of Sciences B, 2001 (12 pp.) | `layopt.pdf` | **normalization** (the third phase): segment-based QP with obstacle-graph ordering constraints; the solver the authors used for eq. (3) |
| — | T. Dwyer, Y. Koren, K. Marriott, *IPSep-CoLa: An Incremental Procedure for Separation Constraint Layout of Graphs*, IEEE TVCG 12(5), 2006 (8 pp.) | `An Incremental Procedure for.pdf` | **the QP solver**: `solve_QPSC`, gradient projection (Fig. 8) with block merge/split projection and Lagrange multipliers (Fig. 9) |
| [7] | T. Dwyer, K. Marriott, P. J. Stuckey, *Fast Node Overlap Removal*, GD 2005, LNCS 3843, pp. 153–164 (12 pp.) | `Fast node overlap removal.pdf` | **non-overlap constraints**: O(n log n) scan-line generation (`generate_Cxno`, `generate_Cyno`, Fig. 1); `satisfy_VPSC` / `solve_VPSC` |
| [3] | T. C. Biedl, B. P. Madden, I. G. Tollis, *The Three-Phase Method: A Unified Approach to Orthogonal Graph Drawing*, IJCGA 10(6), 2000, pp. 553–580 (journal version of the GD 1997 paper; 28 pp.) | `biedl2000.pdf` | **port assignment** (§3.3): ordering the edges on each side of a box (right side clockwise, left side counter-clockwise); the non-overlapping-trees conditions |
| [14] | T. Lengauer, *Combinatorial Algorithms for Integrated Circuit Layout*, 1990/1992 (book) | see below (local only) | **edge routing**: maze running, track graph, bend minimization |

**Lengauer** (the paper's [14]). Local copy, **not in git** (copyrighted book, D-25): `_reference/graphing/(Applicable Theory in Computer Science) Thomas Lengauer (auth.) - Combinatorial Algorithms for Integrated Circuit Layout-Vieweg+Teubner Verlag (1992).pdf` (714 PDF pages; **PDF page = printed page + 18**).

| Need | Lengauer section | Printed pages (PDF pages) |
|---|---|---|
| Shortest paths, Dijkstra, goal-directed (A*-type) single-pair search | §3.8, esp. §3.8.4–3.8.5 | 61–78 (79–96) |
| Maze running (Lee), coding techniques, goal-directed search, suboptimal paths | §8.4.1 | 395–400 (413–418) |
| Line searching (Hightower, Mikami–Tabuchi) and the **track graph**: few-bend routes without a full grid | §8.4.2 | 400–409 (418–427) |
| **Minimizing bends** via the edge-cost transformation of Fig. 8.12 (Dijkstra on an expanded grid) | §8.4.4 | 411–413 (429–431) |
| Sequential routing of many nets (order, rip-up ideas) | §8.3–8.4 | 392–413 (410–431) |
| Graph-based 1-D compaction (constraint graph, longest paths), **wire balancing** | §10.2.2, §10.2.2.3 | 596–619 (614–637) |
| Quadratic optimization, gradient projection (background to the QP solver) | §4.8.1, §4.8.3 | 204–210 (222–228) |
| Simulated annealing and cooling schedules (background to Algorithm 1) | §4.3 | 144–160 (162–178) |

The plan, by phase:
1. **Node placement:** Freivalds–Glagoļevs Algorithm 1 (above). This is the required core.
2. **Compaction inside node placement** (eq. (3)): build the visibility graph from the grid, then solve with **IPSep-CoLa's `solve_QPSC`**:
   - repeat: a steepest-descent step with exact line search (step s = gᵀg / gᵀAg); `project` onto the separation constraints by merging blocks along violated constraints and expanding a block at its smallest-Lagrange-multiplier active constraint; then a second step β = argmin over [0, 1];
   - then `split_blocks` wherever a Lagrange multiplier is negative;
   - until the step is small and nothing split.

   A is the Laplacian of the drawn edges (eq. (3)) plus a small diagonal regularization, so the problem is strictly convex. This is the same algorithm as Freivalds–Kikusts' rigid-component shift/merge with a spanning-tree cut, so it is the solver the layout paper used. It replaces option (b), a generic QP (D-28).
3. **Edge routing** (Lengauer §8.4):
   - a shortest orthogonal route per port edge on the grid left free by the boxes, with the bend-cost transformation of §8.4.4 (Fig. 8.12), so that length plus a bend penalty is minimized exactly by Dijkstra;
   - goal-directed search (§3.8.5.1) for speed;
   - the track graph (§8.4.2.3) if grids get large.
4. **Port assignment** (Biedl et al. §3.3): spread the edges that share a box side over distinct ports, ordering them by where their routes go (clockwise on the right, counter-clockwise on the left, and likewise for top and bottom) to avoid crossings at the box. Add grid lines where a side has more edges than ports.
5. **Normalization** (Freivalds–Kikusts; Fig. 1(c) of the layout paper):
   - one coordinate per *segment*, meaning box sides and every route segment, one direction at a time;
   - ordering constraints from the segment obstacle (visibility) graph, plus minimum box size and the gap δ;
   - objective: node drift Σ ((x_l + x_r)/2 − x_c)² plus w Σ (x_r − x_l)², with w = 10 as the paper suggests;
   - solved with the same `solve_QPSC`; vertical segments first, then horizontal.

   This also separates overlapping route segments from different edges. Box-only overlap constraints for the node stage come from the Dwyer–Marriott–Stuckey scan line (`generate_Cxno`, `generate_Cyno`), which gives O(n) constraints in O(n log n).

**Solver tests (when implemented):**
- **Separation-constraint QP:** every constraint is satisfied, the KKT conditions hold (active multipliers ≥ 0, the gradient is balanced by the active constraints), and the result matches small problems solved by hand.
- **`project` alone:** it is exactly the VPSC least-squares projection of Dwyer–Marriott–Stuckey.

### 17.2 What gets drawn

The model graph is bipartite (§4.1): node vertices (bodies, reservoirs) and element vertices (conductors, later streams, converters, pumps…), with ports as edges.
- **Default: lay out the full bipartite graph.** Bodies and reservoirs are larger labelled boxes. Elements are small boxes (the paper handles mixed node sizes directly). Ports are edges.
  - This draws multi-port elements naturally: an `InternalNode` leg or `Into(h)` is just a third edge.
  - It also keeps parallel conductors visually distinct.
- **Option: collapse two-port elements into labelled edges** for a cleaner schematic. Multi-port elements stay boxes.

### 17.3 Display metadata on the graph

Display data lives with the model graph and is never read by `compile`:
- a display label;
- the component kind (body / reservoir / conductor kind / …), which picks the default box size and style;
- the box size (w, h);
- an optional pinned position (excluded from moves);
- an optional group (e.g. the cells of one `rod!` drawn as a row).

Implementation options:
- (a) a display side-table `Dict{Symbol,DisplayInfo}` in `ThermoSystem`, filled with defaults by `add!`/`connect!`;
- (b) wrapping vertex data as `(component, display)`.

(a) is preferred: it leaves the physics data untouched and is trivially optional.

Layout output: a `NetworkLayout` with box rectangles per vertex and polyline routes per port edge. It is a plain-data result, so it can be tested without any renderer and reused across renders.

### 17.4 Rendering and dependencies

- The layout (pure Julia: grid, moves, compaction, routing) lives in the core package, because it is deterministic and testable.
- Rendering uses Makie (CairoMakie for files, GLMakie for interactive) through a **package extension** (weak dependency), so users of the physics never load a plotting stack. `plot_network(sys_or_cs, u = nothing; overlay = …)` draws the layout, with optional result overlays.
- No ModelingToolkit and no symbolic layer (D-2). GraphMakie's own layouts are not used; it may be used for drawing primitives only if that turns out simpler.

### 17.5 Tests (when implemented)

- **Layout invariants on random and rod-shaped models:**
  - no two boxes overlap (gap ≥ δ);
  - every route is axis-aligned and connects the right ports;
  - routes do not pass through boxes;
  - layouts are reproducible for a fixed seed.
- **Compaction:**
  - the QP solution satisfies every separation constraint;
  - total edge length does not increase across a compaction step;
  - γ = 1 gives no slack along the critical chain.
- **Quality smoke tests:** grid-like models (a rod, a 2-D slab) come out planar with straight edges, like the partial-grid results in the paper's Fig. 4(a).
- **Renderer:** a CairoMakie smoke test that writes a PNG for a small model with an overlay (in the extension's tests).

## 18. Mechanical and rotational systems (revision 3; built in slice 6)

*Proposed 2026-10-15, accepted by the user with the decisions in D-34, built as slice 6 (D-35). Decisions: D-32, D-34, D-35.*

Most of what mechanics needs already exists:
- the `Momentum` and `AngularMomentum` forms of motion, `InertialMass` and `RotationalInertia`;
- gyrators, transformers and junctions;
- momentum conductors: the piston's seal friction and the DC motor's bearing;
- time-varying reservoirs (prescribed velocity drives) and switching.

What is missing:
- the elastic forms of motion (springs);
- prescribed forces and torques;
- gravity;
- nonlinear friction.

None of these needs a new solver or a new `Transfer` kind.

### 18.1 Two new forms of motion (Veinik (5), (6), (12), (13))

```julia
struct Displacement <: FormOfMotion end   # x [m],   force F [N]     (translation, tension and compression: X dx)
struct Angle        <: FormOfMotion end   # θ [rad], torque M [N·m]  (rotation, torsion: M dφ)
```

- The force is ∂U/∂x, the spring's tension. As everywhere, there are no sign flags.
- Default force tolerances are 1e-8 N and 1e-8 N·m, to be checked on the test cases. `atol_scale` applies as usual.
- Veinik's inertial forms (9) and (10) are not added. Kinetic and inertial forms are one phenomenon (`04`), and the library keeps the kinetic one.

### 18.2 Materials (pure accumulators)

| Material | Charge | U | Force |
|---|---|---|---|
| `LinearSpring(k; x0 = 0)` | x | ½k(x − x₀)² | k(x − x₀) |
| `TorsionSpring(k; θ0 = 0)` | θ | ½k(θ − θ₀)² | k(θ − θ₀) |
| `PendulumPotential(m, g, ℓ)` | θ | mgℓ(1 − cos θ) | mgℓ sin θ |

- A nonlinear spring, for example Duffing's ½kx² + ¼k₃x⁴, is just another U.
- Each is a one-DOF accumulator. Kinetic ↔ elastic conversion happens in a converter *between* two bodies, as Veinik §11 and §2.2 require: a body with no coupled part never converts.
- `PendulumPotential` holds gravity reduced to the pendulum's single coordinate, so the pendulum's energy stays inside the model.
- It is not convex: A = mgℓ cos θ ≤ 0 for |θ| ≥ π/2. That is physics (an inverted pendulum is unstable), so `is_stable` is correctly false there.
- It breaks the 1/A_kk tolerance scaling (D-20 §1). A zero A currently falls back to 1e-12·max(|E|, 1), and a small A gives a huge tolerance.
  - **Built:** the material hook `stiffness_floor(mat)` (per charge, default 0). Tolerances and steady-state scaling use max(|A_kk|, floor). For `PendulumPotential` the floor is mgℓ, the value at θ = 0.

### 18.3 Spring branches (converters)

A spring is a body, joined to its two moving nodes by a 3-port converter. It is the mechanical twin of `InductorBranch` (the force–current analogy).

- **`SpringBranch()`**
  - Ports (`Momentum` a, `Momentum` b, `Displacement` s), with forces (w_a, w_b, F).
  - dp_a/dt = +F, dp_b/dt = −F, dx/dt = w_b − w_a.
  - W = [0 0 1; 0 0 −1; −1 1 0], which is antisymmetric, so the branch is power-exact.
- **`TorsionBranch()`:** the same, with `AngularMomentum` and `Angle` (a compliant shaft or gear mesh).
- **A spring to a fixed frame:** end b is a `Momentum` reservoir at w = 0, as for the piston.
- **A pendulum** is `RotationalInertia(mℓ²)` plus `PendulumPotential` plus `TorsionBranch` to a frame reservoir. Equivalently, it is the two bodies joined by `Gyrator(-1; a = AngularMomentum(), b = Angle())`, which gives dL/dt = −M and dθ/dt = ω.
- Both sides store energy, so the system stays an ODE.

### 18.4 Dampers and friction

- **Viscous friction** already exists: `LinearConductance(Momentum(), c)` or `(AngularMomentum(), b)`.
  - J = cΔw and D = cΔw² ≥ 0.
  - It uses kernel B by default: heat at the damper's node, with legs to a wall or the environment. `Into(h)` is the alternative.
- **Monotone nonlinear laws** (§5.2) are new conductance types: J = f(ΔP) with f odd, nondecreasing and f(0) = 0, so D ≥ 0 by type.
  - `QuadraticDrag(c)`: J = c|Δw|Δw. This is aerodynamic drag, with c = ½ρC_dA.
  - `CoulombFriction(F_c; w_ε)`: J = F_c tanh(Δw/w_ε). Regularized; true stick–slip needs a switched or complementarity model (later).

### 18.5 Gears and other kinematic ratios (already built)

- **Gear pair:** `Transformer(r; a = AngularMomentum(), b = AngularMomentum())`.
  - ω_a = r ω_b; the torques are τ_a = λ and τ_b = −rλ.
  - A gear train is transformers in series. Reflected inertia (J_a + J_b/r²) follows without being coded.
- **Rack and pinion, capstan, lead screw:** a `Transformer` between `AngularMomentum` and `Momentum` (ω = v/R). A lever is `Momentum`–`Momentum`.
- **Index 2 (D-29):** the initial speeds must satisfy the ratio. A compliant mesh (`TorsionBranch` with the mesh stiffness) is the ODE alternative, and it also models mesh vibration.
- **Losses:** viscous bearing conductors on each shaft for now. A load-dependent mesh efficiency (loss ∝ transmitted torque) is not a function of ΔP between two nodes; it needs a lossy transformer (later).
- **Clutch:** `Switched(Conductor(LinearConductance(AngularMomentum(), G); …))` between two shafts. A brake is the same switched conductor to the frame.

### 18.6 Forces, torques and gravity (external sources, §8)

- **`Source(domain => rate(t))` on a body** prescribes a non-thermal charge inflow:
  - a force F(t) into `Momentum`, or a torque into `AngularMomentum`;
  - its energy P·rate is booked as external work.
- **Already available alternatives:**
  - base excitation: a prescribed velocity drive (a time-varying `Momentum` reservoir);
  - a dead weight: a constant-force `Displacement` reservoir through a `Gyrator`.
- **`GravityField(g)` on bodies with `Momentum`** (chosen, D-34; **flagged: may need to change**):
  - a momentum source −M g, where M is the body's own inertial mass (`inertial_mass(mat)`: `InertialMass` m, or a `Composite`'s parts). For a `Moving` body (§19.2) M will be its total mass, a function of its state;
  - the momentum axis points up;
  - its work (−ΔPE) is ledgered. For a constant mass the height change is −W/(Mg). The displacement ledger ∫ w dt, which a variable mass needs, is deferred to slice 8;
  - it is one-axis, consistent with D-29 §4.
- **The general rule behind this choice** (D-34, CLAUDE.md): an explicit, ledgered, special-purpose element is preferred to a change of the framework's rules.
- **Not chosen: height as a charge with U ∋ Mgh** (Veinik (11)).
  - It is the same form as dQ = F dx. For a constant mass it already works with no rule change: a separate weight body joined by a converter. `PendulumPotential` is that pattern.
  - For a variable mass, h must live in the body that holds the mass. dh/dt = w then needs a converter between two domains of the *same* body, which the converter rule forbids (ports on distinct nodes, D-29 §3). It would also convert inside a body whose p–h part is uncoupled.
  - **Revisit when** hydrostatic or elevation heads (μ + gz in streams), altitude-dependent gravity, or height as a state are needed.

### 18.7 Out of scope (would need a design change)

- **3-D rigid-body rotation** (Euler's equations, gyroscopic precession, nutation).
  - It needs vector angular momentum (D-29 §4) and a state-dependent antisymmetric operator *inside one body*: dL/dt = L × ω with ω = ∂U/∂L.
  - That operator is power-exact (ωᵀ(L × ω) = 0). It is a new element kind, an internal reversible coupling (the L(x) of GENERIC). Not proposed now.
- **Planar or spatial linkages:** position constraints are index 3.
- **Contact, impact and true stick–slip.**

### 18.8 Tests

**Invariants and kernels:**
- `SpringBranch` and `TorsionBranch` give Pᵀ W P ≡ 0.
- The new conductance types give D ≥ 0 at random ΔP, including ΔP near 0.
- A4 holds for the new materials.
- Source ledgers close energy to round-off.

**Analytic and application cases** (each with closure):

1. **Free mass–spring–damper:** x(t) vs. the analytic solution for ζ = 0.1, 1 and 2.
   - Undamped, the energy is conserved to integrator tolerance.
   - Damped, the heat into the reservoir equals ½kx₀² + p₀²/2m and the entropy equals Q/T_res, for any c (as D1).
2. **Forced oscillator:** a sinusoidal force source.
   - After the transient, the amplitude is F₀/√((k − mΩ²)² + (cΩ)²) and the phase is atan(cΩ/(k − mΩ²)).
   - Mean input power equals mean dissipation, ½cΩ²X².
3. **Two masses with springs:** the normal-mode frequencies.
4. **Spring under gravity:** static deflection Mg/k, and oscillation about it at the unchanged frequency.
5. **Pendulum:**
   - period 2π√(ℓ/g) at θ₀ = 0.01;
   - 4√(ℓ/g)·K(sin(θ₀/2)) at θ₀ = 2.5 rad, with K computed by the arithmetic–geometric mean (no new dependency);
   - with pivot friction, it dissipates mgℓ(1 − cos θ₀) in total.
6. **Gear train:**
   - a constant torque τ on shaft 1 through two stages r₁, r₂ with inertias J₁, J₂, J₃ gives ω₁ = τt/(J₁ + J₂/r₁² + J₃/(r₁r₂)²);
   - the constraint residual stays at integrator tolerance.
7. **Two inertias on a shaft:** ω_n = √(k(J₁ + J₂)/(J₁J₂)).
8. **Clutch engagement** (the rotational two-capacitor problem, D2):
   - the final ω is (J₁ω₁ + J₂ω₂)/(J₁ + J₂);
   - the loss is ½J₁J₂(ω₁ − ω₂)²/(J₁ + J₂) for any clutch conductance, and the entropy is loss/T.
9. **DC motor + gearbox + load:** steady speed and current vs. a hand calculation.

## 19. Rockets: the rocket equation and a simplified engine (proposed, revision 3; slice 8)

*Proposed 2026-10-15, not reviewed. Feasibility study and design. Decisions: D-33.*

**Verdict: possible within the framework.** It needs reactions (slice 7) and the mechanics additions (§18, built in slice 6), plus three new pieces:
- a moving variable-mass material;
- momentum carried by streams;
- an ideal nozzle element.

There is no new solver and no new `Transfer` kind: the nozzle is a stream-kind element and gravity is a source. The limits are listed in §19.6.

### 19.1 The physics in Veinik's terms

A rocket chains three conversions:
1. **Chemical → thermal:** combustion in the chamber, an internal transformation (§8) at dU = 0.
2. **Thermal → kinetic:** the nozzle, an advected charge bundle whose enthalpy becomes jet velocity.
3. **Momentum to the vehicle:** the momentum the jet carries backward is what the vehicle gains forward. Momentum is a charge and is conserved, so no special "thrust" element exists.

The momentum balance of the vehicle (mass M, velocity w, exhaust velocity v_e relative to it, exit area A_e) is:

  d(Mw)/dt = −ṁ(w − v_e) + (p_e − p_a)A_e − Mg,  with dM/dt = −ṁ,

  so M dw/dt = ṁ v_e + (p_e − p_a)A_e − Mg ≡ ṁ c_eff − Mg.

With constant c_eff and no gravity this is Tsiolkovsky: Δv = c_eff ln(m₀/m_f).

### 19.2 Moving bodies of variable mass

`Moving(material; m_dry = 0)` adds a `Momentum` charge to any material that has mass:

  U(E, p) = U_int(E) + p²/(2M),  M = m_dry + Σ_k m_k.

- **Forces:** the force on p is w = p/M. Mass and momentum are coupled (one coupling group), so μ_k = μ_k,int − w²/2 at fixed momentum (`07` §4; §7.4 "automatic"). T and p do not depend on w (Galilean invariance).
- **Rigid links:** every body that travels with the vehicle is `Moving` and carries its own momentum: tanks, the chamber, and the dry structure (lumped into a tank's m_dry).
  - They are held at one velocity by `Transformer(1; a = Momentum(), b = Momentum())`, which already exists (index 2, D-29).
  - Its λ is the structural force between the bodies; for example, the thrust that the chamber transmits to the tank.
- **Why not one vehicle momentum whose inertia sums the other bodies' masses:** U must be a function of the body's own charges (§2.1).
- **Why the tanks must move too:** a stream from a moving tank to a non-moving chamber would book ½ṁw² as dissipation, which is wrong. Every body that exchanges mass with a moving body moves with it.
- **Homogeneity:** U is first-order homogeneous only for m_dry = 0. Capacities use the fixed-mass block, as always (D-8).

### 19.3 Streams and pumps that carry momentum

When both ends hold `Momentum`, the bundle adds momentum ṁ w_u (§7.1 "momentum if kept"; not yet built).
- **Energy:** the upstream body loses ṁ(g_u − w_u²/2) + T_u ṁ s_u + w_u ṁ w_u = ṁ(h_u + w_u²/2). The downstream body receives the same, which leaves:

    T_d Ṡ_d = ṁ(h_u − g_d) + ½ṁ(w_u − w_d)².

- **Generation** is today's stream generation plus ṁ(w_u − w_d)²/(2T_d) ≥ 0. This is the inelastic-mixing loss, and it is zero between rigidly linked bodies.
- **`IdealPump`** computes its work from the internal material, h(s_u, p_d) − h_u, using `U_int`. On linked bodies the w-terms cancel.

### 19.4 The ideal nozzle (a new stream-kind outlet)

`Nozzle(A_t, ε; η = 1)` runs from a chamber body c to outside the model, like a `PassThroughOutlet`. Its parameters are the throat area A_t, the area ratio ε = A_e/A_t, and an optional kinetic-energy efficiency η.

**Mass flow.** The throat is choked:

  ṁ = A_t p_c √(k/(R T_c)) (2/(k+1))^((k+1)/(2(k−1))),

with the frozen chamber composition setting k and R (Veinik §59 (448), (453)). The slice may support choked flow only, and raise an error when p_a/p_c rises above the critical ratio β_cr, as a passive-boundary violation does.

**Exit state.**
- Isentropic, frozen expansion to the exit Mach number set by ε. The supersonic root of the area–Mach relation (Veinik (454)) is found by a safeguarded 1-D Newton.
- From the material: T_e, p_e, h_e, and v_e = √(2(h_c − h_e)).
- With η < 1: v_e² = 2η(h_c − h_e,s), h_e = h_c − v_e²/2, and s_e = s(h_e, p_e) > s_c.

**Bundle taken from c** (c is `Moving`; on a test stand, see below):
- species ṁ y_k at the chamber composition;
- entropy ṁ s_c;
- momentum J_p = ṁ(w − v_e) − p_e A_e, the stream thrust (the exit pressure pushes the vehicle forward).

**Exact energy closure.**
- Σ P·J removes ṁ(h_c + w²/2 − w v_e) − w p_e A_e from c.
- That equals the energy the jet carries across the moving exit plane, ṁ(h_e + (w − v_e)²/2) − w p_e A_e. The last term is the flow work at an exit plane moving at w.
- The two are equal exactly when h_c = h_e + v_e²/2, which is how v_e is defined.
- The jet's mass, energy, entropy ṁ s_e and momentum are booked as boundary flows leaving the model.
- The generation, ṁ(s_e − s_c) ≥ 0, is booked on the nozzle. It is zero for η = 1.
- The gas held up in the nozzle is neglected (quasi-steady nozzle).

**Ambient pressure p_a.** The atmosphere acts on the vehicle's projected exit area. This is the existing piston–atmosphere converter, `Gyrator(-A_e; a = Volumetric(), b = Momentum())`, between an atmosphere reservoir and the vehicle. Thrust is then ṁ v_e + (p_e − p_a)A_e. In vacuum, leave it out.

**Test stand (static fire).**
- The chamber is not `Moving`. Its momentum port is a `Momentum` reservoir at w = 0, which receives −J_p = ṁ v_e + p_e A_e.
- That reservoir's momentum-ledger rate *is* the thrust (a load cell), and it receives no work (w = 0).

**Terminology.** The diverging part of a rocket nozzle accelerates supersonic flow; it is not a diffuser. A subsonic diffuser (deceleration and pressure recovery) is the same relations on the subsonic branch, and can be added later if needed.

### 19.5 The simplified engine and the tests

**Engine model.**
- **Feed:**
  - Gaseous propellants (H₂, O₂) are held in `Moving` `IdealGas` tanks, one species each.
  - Pump-fed: an `IdealPump` per propellant (single species, as built), at prescribed ṁ, which sets O/F.
  - Pressure-fed: a `Stream` with `HydraulicConductance` (blowdown).
- **Chamber:**
  - A `Moving` body with an ideal-gas mixture (H₂, O₂, H₂O; optionally dissociation species) and one reaction, 2H₂ + O₂ ⇌ 2H₂O, from slice 7.
  - The reaction is run fast (its time constant ≪ the residence time), so the chamber sits at chemical equilibrium, as standard performance calculations assume.
  - The nozzle is frozen: no reaction inside it.
- **Flight:**
  - `GravityField(g)` on every moving body.
  - Optionally `QuadraticDrag` on `Momentum` between the vehicle and an air reservoir at w = 0, with `Into` the air's thermal reservoir.

**Tests.**
1. **Kernel-level checks:**
   - nozzle Σ P·J = jet energy flux at random states and random w (round-off);
   - generation ṁ(s_e − s_c) ≥ 0;
   - momentum closure: vehicle + jet ledger + gravity ledger = 0.
2. **Galilean invariance:** the same engine at w = 0 and w = 3 km/s gives identical chamber T and p, thrust and generation.
3. **Cold-gas thruster:** a `Moving` gas tank with dry mass, through an ideal nozzle to vacuum.
   - The tank stays on its isentrope (as in the tank-blowdown test), so c_eff(m) = c₀(m/m₀)^((k−1)/2) in closed form.
   - Δv vs. ∫ c_eff dm/(m_dry + m) by quadrature.
   - Thrust vs. ṁ c_eff at every saved time.
4. **Static fire:** the bipropellant engine on a stand, at steady state, vs. a hand calculation with the same property model, solved independently in the test:
   - T_c: the adiabatic flame temperature at equilibrium composition;
   - p_c = ṁ c*/A_t, with c* = √(R T_c/k)/Γ;
   - the exit Mach number from ε;
   - C_F and I_sp = F/(ṁ g₀).
   - **Loss map:** generation at injection (pump and inlet mixing), in combustion (the reaction, the dominant term), and in the nozzle (zero for η = 1, ṁ(s_e − s_c) otherwise), plus the exergy the jet carries away.
5. **Flight:** the engine on a vehicle with gravity, starting from the steady chamber state.
   - Δv vs. c̄_eff ln(m₀/m_f) − g t_b, within the bound (max − min c_eff)·ln(m₀/m_f). The tank states drift during the burn, so c_eff is not exactly constant.
   - Momentum and energy closure.
   - Burnout altitude from the gravity ledger's ∫ w dt.

### 19.6 Limits of the simplified model

- **Lumped chamber; quasi-1-D isentropic nozzle** in closed form.
  - The nozzle flow is not resolved: there are no internal shocks (over-expanded separation), and boundary-layer and divergence losses enter only through η.
  - Resolving it with a chain of bodies would need inertance-driven compressible streams, which is deferred and not worth it here.
- **Constant c_p per species** (calorically perfect), which overpredicts T_c.
  - Real performance needs c_p(T) (NASA polynomials) and dissociation species.
  - c_p(T) needs the 1-D Newton inversion of §2.4, planned with the Clapeyron backend. The element interfaces do not change.
- **Frozen nozzle flow;** no shifting equilibrium.
- **Not real combustion kinetics.** The linear affinity law J = M_r A/T is right for an equilibrium chamber, but it does not model ignition or kinetics-limited combustion. Those need a sign-preserving nonlinear rate law (mass action: J ∝ 1 − e^(−A/RT), with A·J ≥ 0). That law is planned for slice 7 (D-34).
- **Gaseous propellants only,** until the real-fluid backend exists (liquid and cryogenic tanks).
- **One axis:** vertical ascent, or straight-line flight in free space. Gravity turns and attitude dynamics need vector momentum and rigid-body rotation (§18.7).
- **No feed-line or injector dynamics, and no combustion instability.** These need inertance (deferred).
