"""
    VeinikThermo

Multi-domain thermodynamic modeling on A. I. Veinik's generalized thermodynamics.
The state is the bodies' charges; forces are P = ∂U/∂E; elements map forces to
fluxes; one assembler does all conservation and entropy bookkeeping. See
`docs/DESIGN.md`. SI units throughout.

Slice 1: the thermal domain, constant-capacity solids, reservoirs and thermal conductors.
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
using NonlinearSolve: NonlinearProblem, NewtonRaphson

include("domains.jl")
include("materials.jl")
include("elements.jl")
include("system.jl")
include("assembly.jl")
include("solve.jl")
include("spatial.jl")

# domains
export FormOfMotion, Thermal, charge_symbol, charge_unit, force_symbol, force_unit, producible
# materials
export Material, ConstantCapacitySolid, domains, internal_energy, forces, eos_matrix, temperature, charges_at
# elements
export Body, Reservoir, ThermalConductance, Conductor
# system
export ThermoSystem, add!, connect!, compile, CompiledSystem, body_names, reservoir_names, conductor_names
# running
export initial_state, rates, simulate, steady_state
# results and diagnostics
export charges, set_temperature!, flux, heat_flow, entropy_production, generated_entropy, reservoir_ledger,
       total_energy, total_entropy, energy_residual, closure
# spatial helpers
export rod!

end # module VeinikThermo
