"""
    ThermoKinetics

Multi-domain thermodynamic modeling on A. I. Veinik's generalized thermodynamics.
The state is the bodies' charges; forces are P = ∂U/∂E; elements map forces to
fluxes; one assembler does all conservation and entropy bookkeeping. See
`docs/DESIGN.md`. SI units throughout.

Implemented so far (see docs/STATUS.md): thermal, volumetric and species domains; solids,
ideal gases, elastic volumes and composites; reservoirs; thermal and volume conductors
with their dissipation policies; streams, ideal pumps, flow boundaries and outlets; the
electrical domain, resistors, mass conductors and coupled (thermoelectric, Knudsen, Onsager)
conductors; kinetic forms of motion, ideal converters, junctions and transformers (mass-matrix
DAE), switched connections and time-varying reservoirs; elastic forms of motion (springs, the
pendulum potential), spring branches, monotone friction laws, force and torque sources and gravity;
ideal-gas mixtures and ideal solutions with formation data, reactions (affinity and mass-action
laws), nuclear sources and heat sources; moving variable-mass bodies, momentum-carrying streams
and pumps, and an ideal nozzle (rockets).
"""
module ThermoKinetics

using ForwardDiff
using StaticArrays
using Graphs
using MetaGraphsNext
using SparseArrays
using LinearAlgebra
using SciMLBase
using OrdinaryDiffEqRosenbrock: Rodas5P
using NonlinearSolve: NonlinearProblem, NewtonRaphson, TrustRegion

include("domains.jl")
include("materials.jl")
include("species.jl")
include("moving.jl")
include("elements.jl")
include("reactions.jl")
include("system.jl")
include("assembly.jl")
include("solve.jl")
include("spatial.jl")

# domains
export FormOfMotion, Thermal, Volumetric, Species, Electrical, Momentum, AngularMomentum, FluxLinkage,
       Displacement, Angle, charge_symbol, charge_unit, force_symbol, force_unit, producible
# materials
export Material, ConstantCapacitySolid, IdealGas, ElasticVolume, LinearCapacitor, InertialMass, RotationalInertia,
       LinearInductor, LinearSpring, TorsionSpring, PendulumPotential, Composite,
       domains, internal_energy, forces, eos_matrix, temperature, pressure, charges_at,
       capacity_matrix, is_stable, valid_state, mass_slots,
       SpeciesData, IdealGasMixture, IdealSolution, molar_mass, nucleons, species_names, species_data, Moving
# elements
export Body, Reservoir, ThermalConductance, VolumeConductance, LinearConductance, Resistor, MassConductance,
       ThermoelectricConductance, KnudsenBarrier, OnsagerConductance, Conductor, Split, Into, InternalNode,
       HydraulicConductance, Stream, IdealPump, FlowBoundary, PassThroughOutlet,
       Junction, Converter, Gyrator, InductorBranch, Transformer, Switched,
       SpringBranch, TorsionBranch, QuadraticDrag, CoulombFriction, Source, GravityField,
       Reaction, AffinityLaw, MassAction, NuclearSource, FirstOrder, HeatSource, Nozzle
# basis transformations of Veinik's coefficient data
export thermoelectric_from_veinik, barrier_mass_basis, electrokinetic_mass_basis
# system
export ThermoSystem, add!, connect!, compile, CompiledSystem, body_names, reservoir_names, conductor_names, element_names
# running
export initial_state, consistent_state, rates, simulate, steady_state
# results and diagnostics
export charges, force, set_temperature!, flux, heat_flow, entropy_production, entropy_production_parts,
       generated_entropy, reservoir_ledger, reservoir_charge, boundary_ledger, pump_work, source_ledger, mass_flow,
       total_charge, total_energy, total_entropy, reaction_rates, element_totals, nozzle_state,
       energy_residual, closure, algebraic_index, constraint_residual
# spatial helpers
export rod!

end # module ThermoKinetics
