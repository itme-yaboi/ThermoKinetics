# ThermoKinetics.jl

ThermoKinetics is a Julia library for multi-domain thermodynamic modelling, built on
A. I. Veinik's generalized thermodynamics. It is intended for thermal and balance-of-plant
analysis, where correctness, traceability and loss attribution (where entropy is produced)
matter more than raw speed. All quantities are SI.

## Core ideas

- **Bodies and charges.** A system is assembled from *bodies*. Each body holds a set of
  *forms of motion* (domains), each a conjugate pair of a generalized charge `E` and a
  generalized force `P`: entropy and temperature (thermal), volume and pressure
  (volumetric), species amounts and chemical potentials, electric charge and potential,
  momentum and velocity, and so on.
- **One state function.** A body's material defines its energy `U(E)`. Forces are always
  `P = ∂U/∂E`, so the volumetric force is `-p`. Coupling between domains and static reciprocity
  come from `U` being a single potential.
- **Element kinds.** Bodies are connected by *conductors* (dissipative transport, possibly
  with coupled fluxes), *streams* (advection), *converters* (reversible transformation of one
  charge into another at a fixed ratio), *reactions* and *nuclear processes* (internal
  transformations), and *sources* and *reservoirs* (boundaries). Adding a form of motion means
  declaring a new `(E, P)` pair, not writing a new solver.
- **Ledgers and closure.** One assembler does all conservation and entropy bookkeeping.
  Energy is conserved across every connection, other charges are conserved except through
  converters and explicit internal transformations, and entropy production is never negative.
  [`closure`](@ref) checks the energy and entropy balance of a run, and
  [`entropy_production`](@ref) attributes entropy production to each element.

## Pages

- [Getting started](getting_started.md): short runnable examples.
- [Architecture and structure](architecture.md): how the library is organised.
- API reference: [Domains](api/domains.md), [Materials](api/materials.md),
  [Nodes](api/nodes.md), [Conductors and dissipation](api/conductors.md),
  [Converters and transformers](api/converters.md), [Flow](api/flow.md),
  [Sources, reactions and nuclear processes](api/sources.md), [Switching](api/switching.md),
  [System building and running](api/system.md), [Results and diagnostics](api/results.md),
  [Spatial helpers](api/spatial.md).
- Project records: [Design](project/DESIGN.md), [Status](project/STATUS.md) and
  [Decisions](project/DECISIONS.md) (copies of the repository documents).

## Module

```@docs
ThermoKinetics
```
