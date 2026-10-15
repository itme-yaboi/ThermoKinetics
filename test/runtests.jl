using VeinikThermo
using Test
using Random
using LinearAlgebra
using StaticArrays

include("helpers.jl")

@testset "VeinikThermo" begin
    include("test_build.jl")          # construction, validation, units hygiene
    include("test_materials.jl")      # state functions: A4, A7, B3
    include("test_kernel.jl")         # conductor bookkeeping identities (kernel A, thermal row)
    include("test_invariants.jl")     # 08 §A: A1 energy, A3 second law, A6 equilibrium (A2, A4, A5 where applicable)
    include("test_thermostatics.jl")  # 08 §B: B2
    include("test_thermokinetics.jl") # 08 §C: C1, C2, C8
    include("test_unsteady.jl")       # 08 §D: D1 (thermal analogue), D4 (= D3 thermal), D5 (slab)
    include("test_streams.jl")        # slice 3: B5, C7, streams, pumps, boundaries
    include("test_applications.jl")   # DESIGN §0.1 application tests
end
