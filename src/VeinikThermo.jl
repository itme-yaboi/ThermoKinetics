"""
    VeinikThermo

Multi-domain thermodynamic modeling on A. I. Veinik's generalized thermodynamics.
The state is the bodies' charges; forces are P = ∂U/∂E; elements map forces to
fluxes; one assembler does all conservation and entropy bookkeeping. See
`docs/DESIGN.md`. SI units throughout.

Implemented so far (see docs/STATUS.md): thermal, volumetric and species domains; solids,
ideal gases, elastic volumes and composites; reservoirs; thermal and volume conductors
with their dissipation policies; streams, ideal pumps, flow boundaries and outlets.
"""
module VeinikThermo

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
include("elements.jl")
include("system.jl")
include("assembly.jl")
include("solve.jl")
include("spatial.jl")

# domains
export FormOfMotion, Thermal, Volumetric, Species, charge_symbol, charge_unit, force_symbol, force_unit, producible
# materials
export Material, ConstantCapacitySolid, IdealGas, ElasticVolume, Composite,
       domains, internal_energy, forces, eos_matrix, temperature, pressure, charges_at,
       capacity_matrix, is_stable, valid_state, mass_slots
# elements
export Body, Reservoir, ThermalConductance, VolumeConductance, Conductor, Split, Into, InternalNode,
       HydraulicConductance, Stream, IdealPump, FlowBoundary, PassThroughOutlet
# system
export ThermoSystem, add!, connect!, compile, CompiledSystem, body_names, reservoir_names, conductor_names, element_names
# running
export initial_state, rates, simulate, steady_state
# results and diagnostics
export charges, force, set_temperature!, flux, heat_flow, entropy_production, entropy_production_parts,
       generated_entropy, reservoir_ledger, reservoir_charge, boundary_ledger, pump_work, mass_flow,
       total_charge, total_energy, total_entropy,
       energy_residual, closure
# spatial helpers
export rod!

end # module VeinikThermo
