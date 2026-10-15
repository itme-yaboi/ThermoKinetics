using ThermoKinetics
using Test
using Random
using LinearAlgebra
using StaticArrays
using ForwardDiff

include("helpers.jl")

@testset "ThermoKinetics" begin
    include("test_build.jl")          # construction, validation, units hygiene
    include("test_materials.jl")      # state functions: A4, A7, B3
    include("test_kernel.jl")         # conductor bookkeeping identities (kernel A, thermal row)
    include("test_invariants.jl")     # 08 §A: A1 energy, A3 second law, A6 equilibrium (A2, A4, A5 where applicable)
    include("test_thermostatics.jl")  # 08 §B: B2
    include("test_thermokinetics.jl") # 08 §C: C1, C2, C8
    include("test_unsteady.jl")       # 08 §D: D1 (thermal analogue), D4 (= D3 thermal), D5 (slab)
    include("test_streams.jl")        # slice 3: B5, C7, streams, pumps, boundaries
    include("test_electrical.jl")     # slice 4: C3, D1, D2, resistors, mass conductors
    include("test_coupled.jl")        # slice 4: C4, C5, C6, coupled kernel, Onsager–Casimir
    include("test_kinetic.jl")        # slice 5: kinetic DOFs, converters, junctions, transformers, switching
    include("test_cycles.jl")         # slice 5, 9: 08 §E cycles (E1, E2 Carnot; E3 rectangular)
    include("test_mechanics.jl")      # slice 6: 08 §G G3, G4, G9; springs, friction laws, sources, gears
    include("test_reactions.jl")      # slice 7: 08 §H H1–H8, H10; mixtures, reactions, nuclear and heat sources
    include("test_rockets.jl")        # slice 8: moving bodies, momentum streams, nozzle, gravity on variable mass
    include("test_applications.jl")   # DESIGN §0.1 application tests
end
