# Architecture and structure

This page describes how the code is organised and how a model travels from construction to
results. The *why* behind each choice is recorded in the project's design records (DESIGN,
DECISIONS); entry numbers such as D-26 refer to DECISIONS.

## In one paragraph

A model is a graph of **nodes** (bodies, reservoirs, junctions) and **elements** (conductors,
streams, converters, sources, reactions, …). The state is the bodies' **charges** (entropy,
volume, species masses, momentum, …). Every force is the gradient of one material potential,
``P = \partial U / \partial E``. Elements map forces to fluxes, and a single assembler does
*all* the bookkeeping: it moves charges, delivers heat, and records entropy production and
external exchanges in **ledgers**. [`compile`](@ref) lowers the graph to flat, typed data
once; the solvers then call one right-hand side that never touches the graph. Nothing is
generated symbolically (no ModelingToolkit, D-2), so every step can be inspected.

## Source layout

Files are included in this order by `src/ThermoKinetics.jl`:

| File | Role | Main contents |
|---|---|---|
| `domains.jl` | Forms of motion | [`FormOfMotion`](@ref) and its subtypes ([`Thermal`](@ref), [`Volumetric`](@ref), [`Species`](@ref), [`Electrical`](@ref), [`Momentum`](@ref), …): names, units, default force tolerances. No physics. |
| `materials.jl` | State functions | [`Material`](@ref), [`forces`](@ref) and [`eos_matrix`](@ref) by automatic differentiation, the analytic materials ([`ConstantCapacitySolid`](@ref), [`IdealGas`](@ref), springs, inertias, …) and [`Composite`](@ref). |
| `species.jl` | Multi-species materials | [`SpeciesData`](@ref) with the atomic-mass table, [`IdealGasMixture`](@ref), [`IdealSolution`](@ref), the trace floor (D-41). |
| `elements.jl` | Nodes and laws | [`Body`](@ref), [`Reservoir`](@ref), [`Junction`](@ref); conductor laws and dissipation policies; streams, pumps, boundaries, outlets; converters, transformers; [`Source`](@ref), [`GravityField`](@ref), [`Switched`](@ref). Elements hold parameters only. |
| `reactions.jl` | Internal transformations and sources | [`Reaction`](@ref) with [`AffinityLaw`](@ref) / [`MassAction`](@ref), [`NuclearSource`](@ref) with [`FirstOrder`](@ref), [`HeatSource`](@ref). |
| `system.jl` | The model graph | [`ThermoSystem`](@ref) (a `MetaGraphsNext` bipartite graph), [`add!`](@ref), [`connect!`](@ref), port roles. |
| `assembly.jl` | Compilation and the right-hand side | [`compile`](@ref) → [`CompiledSystem`](@ref); the compiled element structs; `transfer` and `_apply!` per element kind; the conductor kernels; `assemble!`; the Jacobian sparsity pattern. |
| `solve.jl` | Running and results | [`simulate`](@ref), [`steady_state`](@ref), tolerances (`abstol_vector`), accessors ([`flux`](@ref), [`mass_flow`](@ref), [`reaction_rates`](@ref), …), ledgers, [`closure`](@ref). |
| `spatial.jl` | Mesh helpers | [`rod!`](@ref): a 1-D finite-volume chain of bodies and conductors. |

## From model to results

```text
 build                         compile (once)                         run
 ─────                         ──────────────                         ───
 ThermoSystem()                validate ports, domains, kernels,      simulate(cs, tspan)
 add!(sys, Body / Reservoir)   stoichiometry, nucleon balance  ──►    │  SciML ODEProblem (Rodas5P),
 connect!(sys, element, …)     lay out the state vector               │  mass-matrix DAE if junctions
        │                      allocate ledger slots                  │  or transformers exist
        │  graph of nodes      group elements by concrete type        ▼
        └─ and elements ─────► CompiledSystem (flat, typed) ───────► solution u(t)
                                                                      │
                                steady_state(cs) ◄──────────── Newton │
                                                                      ▼
                                       results: force, temperature, flux, entropy_production,
                                       reservoir/boundary/source ledgers, closure, element_totals
```

### 1. Building

[`ThermoSystem`](@ref) is a bipartite graph. Nodes are added with [`add!`](@ref); elements
with [`connect!`](@ref), which records each port's role (`:a`, `:b`, `:h` for a conductor's
third thermal node, `:p1 … :pN` for a converter). The graph is used for building, validation
and topology queries only; it is never touched in the right-hand side.

### 2. Compiling

[`compile`](@ref) (in `assembly.jl`) validates the model and lowers it to a
[`CompiledSystem`](@ref):

- **Validation:** every port's domain exists on its node; flow elements see fluids; reactions
  balance every element; nuclear sources balance nucleons and charge; converter matrices are
  antisymmetric; boundaries are labelled.
- **State layout:** `u = [u_phys; u_alg; q_ledger]`.
  - `u_phys` is the bodies' charges, body by body, in each material's `domains` order.
  - `u_alg` holds the algebraic unknowns of the DAE: junction forces, then transformer fluxes.
  - `q_ledger` holds the ledgers, element by element in insertion order, then per reservoir.
- **Force buffer:** `P = [body forces; reservoir forces; junction forces]`. For a body, a
  charge's slot and its force's slot have the same index.
- **Grouping:** compiled elements are stored as one concretely typed vector per element type,
  held in a `Vector{Any}`. Each group's loop is specialised (a function barrier), while
  `CompiledSystem` itself has one concrete type, so every model shares one compiled solver
  (D-26).

### 3. The right-hand side

`assemble!(du, u, cs, t)` is the whole model:

1. If any body's charges are not a valid state of its material, fill `du` with NaN (the
   integrator rejects the step).
2. Fill the force buffer: `forces(material, E)` for every body (automatic differentiation of
   U), reservoir forces (possibly time-varying), junction forces from `u_alg`.
3. For every element group and element: `tr = transfer(c, P, u, cs, t)` computes fluxes,
   heat deliveries and production from forces (pure); `_apply!` books them into `du`.
4. For transformers, a second pass forms the index-reduced constraint rows.

**Only the `_apply!` methods write `du`.** A wrong law gives a wrong flux, but it cannot create
charge, break energy conservation or write negative production. That is how the
non-negotiables are enforced structurally rather than by care in each law.

### The element kinds and their bookkeeping

| Compiled kind | Built from | Each non-thermal charge | Energy | Entropy written |
|---|---|---|---|---|
| `CompiledConductor` (kernels `ThermalRow`, `WeightsPolicy`, `NodePolicy`, `CoupledPolicy`) | [`Conductor`](@ref) with a law | the same J out of one port and into the other | heat deliveries sum to the released power exactly | Σ Q_i/T_i ≥ 0 (exact discrete kernels) |
| `CompiledStream`, `CompiledPump`, `CompiledBoundary`, `CompiledOutlet` | [`Stream`](@ref), [`IdealPump`](@ref), [`FlowBoundary`](@ref), [`PassThroughOutlet`](@ref) | advected bundle, conserved | exact; pump and boundary work ledgered | advected entropy plus generation ≥ 0 |
| `CompiledConverter`, `CompiledTransformer` | [`Converter`](@ref), [`Gyrator`](@ref), branches, [`Transformer`](@ref) | **not** conserved (the one exception) | Pᵀ W P ≡ 0, or λ(P_a − rP_b) on the constraint | none |
| `CompiledSource` | [`Source`](@ref), [`GravityField`](@ref) | external, ledgered | external work ledgered | none |
| `CompiledReaction` | [`Reaction`](@ref) | per stoichiometry; elements conserved | dU = 0 inside the body | A·ξ̇/T ≥ 0, exact in floating point |
| `CompiledNuclear` | [`NuclearSource`](@ref) | prescribed, ledgered per species | dU = E_dep ξ̇, external | boundary entropy, either sign |
| `CompiledHeat` | [`HeatSource`](@ref) | none | Q̇, external | Q̇/T as boundary entropy |
| `CompiledSwitched` | [`Switched`](@ref) | as the gated element, times a gate in {0, 1} | | |

## Materials

Each material is **one potential** U(E). [`forces`](@ref) is its gradient and
[`eos_matrix`](@ref) its Hessian, both by ForwardDiff, so ``P = \partial U/\partial E`` and
Maxwell symmetry hold by construction, for every material. Coupling between forms of motion is
expressed by what U depends on: a single material couples everything it contains;
[`Composite`](@ref) joins independent parts (block-diagonal Hessian).

Materials are closed-form (no inversion). The multi-species materials carry
[`SpeciesData`](@ref) (atoms, formation enthalpy, standard entropy), derive molar masses from
the atoms so reactions conserve mass by construction, and smooth the ideal-mixing term below a
trace floor so a depleted species cannot stall the integrator (D-41, D-42, D-44).

## Ledgers and diagnostics

Ledgers are pure quadratures in `q_ledger`; they never feed back into the dynamics. They record:

- entropy production per producing element ([`entropy_production`](@ref),
  [`generated_entropy`](@ref), [`entropy_production_parts`](@ref));
- per reservoir: the charge, entropy and energy received ([`reservoir_ledger`](@ref),
  [`reservoir_charge`](@ref));
- per pump, boundary and source: work, energy, entropy and mass exchanged with the outside
  ([`pump_work`](@ref), [`boundary_ledger`](@ref), [`source_ledger`](@ref)).

[`closure`](@ref) checks a finished run: energy (to integrator tolerance), entropy (to
round-off: the change in the bodies' entropy equals what crossed the boundary plus what was
produced), and every conserved charge, with element totals ([`element_totals`](@ref)) replacing
species that reactions change. It also reports each species' and element's absolute drift.
[`energy_residual`](@ref) checks the instantaneous balance at a single state.

## Solvers

- **Transient:** [`simulate`](@ref) builds a SciML `ODEProblem` and integrates with `Rodas5P`
  (Rosenbrock, exact AD Jacobian, sparse pattern for large ODE models). With junctions or
  transformers it is a mass-matrix DAE, started consistent by [`consistent_state`](@ref).
- **Tolerances** (`abstol_vector`): each charge is resolved to the equivalent of a force
  tolerance through its capacity (1e-8 K, 1e-4 Pa, …), never by its magnitude, since entropy
  has an arbitrary zero. Ledgers are in the error norm by default, scaled like the states.
- **Admissibility:** steps that leave a material's valid states, or push a passive boundary
  uphill, are rejected. Discontinuities of time-varying inputs are passed as `tstops` and
  evaluated as left limits.
- **Steady state:** [`steady_state`](@ref) solves the physical rates = 0 for the charges some
  element writes, by Newton and then trust region, with residual rows scaled by their
  round-off floor. Every group of bodies must be anchored to a reservoir or boundary for energy
  and for each conserved charge; isolated groups are refused (use `simulate`).

## Where the non-negotiables live

| Non-negotiable | Enforced in |
|---|---|
| Charges conserved by connections; converters the only exception; production ≥ 0 | `_apply!` methods (the same J leaves one port and enters the other); converter kernels (antisymmetric W only); reaction and nuclear stoichiometry checked in `compile`; production computed in sign-exact forms |
| Energy conserved across every connection | conductor kernels (the last heat delivery is the released power minus the others); Pᵀ W P ≡ 0; stream bundles by Euler's relation; work and external energy ledgered |
| P = ∂U/∂E, volumetric force −p | `forces` is the only way a force is computed; domains carry no sign flags |
| Symmetric static Hessian; PSD kinetic matrices | U is one scalar per material; `OnsagerConductance` and `AffinityLaw` check symmetry and PSD and factor L = BBᵀ |
| True forces (μ) across boundaries | laws receive forces from the assembler, never concentrations |
| SI only | Float64 SI throughout; no unit constants in the code |

## Extending the library

- **A new form of motion:** a subtype of [`FormOfMotion`](@ref) with its symbols, units and a
  default force tolerance (`domains.jl`). No solver change.
- **A new material:** a subtype of [`Material`](@ref) with `domains` and a closed-form
  `internal_energy`; optionally `valid_state`, `charges_at`, `stiffness_floor`. Forces and
  the Hessian come for free.
- **A new conductor law:** a `ConductorLaw` with `nonthermal_domains` and `fluxes`; the
  existing kernels do the bookkeeping. The law must make the released power non-negative by
  type.
- **A new element kind** (rare; a design change, D-34): its parameter struct (`elements.jl`
  or `reactions.jl`), a `connect!` method (`system.jl`), its ledger count, ledger kinds and
  compile branch, a compiled struct with `transfer` and `_apply!` (`assembly.jl`), and its
  accessors (`solve.jl`). Prefer an explicit special-purpose element within the existing kinds
  over bending a rule.

## Tests

| File | Covers |
|---|---|
| `test_build.jl` | construction and validation; every export has a docstring |
| `test_materials.jl`, `test_kernel.jl` | state functions; conductor bookkeeping identities |
| `test_invariants.jl` | invariants (energy, second law, equilibrium) on random networks |
| `test_thermostatics.jl`, `test_thermokinetics.jl`, `test_unsteady.jl` | analytic cases: equilibrium, transport, transients |
| `test_streams.jl`, `test_electrical.jl`, `test_coupled.jl` | streams and boundaries; circuits; coupled conductors |
| `test_kinetic.jl`, `test_cycles.jl`, `test_mechanics.jl` | converters, DAEs and switching; Carnot cycles; mechanics |
| `test_reactions.jl` | mixtures, reactions, nuclear and heat sources |
| `test_applications.jl` | one application test per delivered use (walls, vessels, tanks, heat exchangers, circuits, mechanics, internal heating, equilibrium, tritium) |

Every integration test also checks [`closure`](@ref).

## Design records

- **DESIGN**: the target architecture and the vertical slices; what is built and what is
  planned.
- **DECISIONS**: append-only, one entry per decision with alternatives (D-1, D-2, …).
- **STATUS**: where the project stands, rewritten after every slice.

They are included under *Project records* in this site.
