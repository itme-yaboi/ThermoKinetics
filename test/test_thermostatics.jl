@testset "08 §B thermostatics" begin
    @testset "B2 calorimetry" begin
        m = (1.5, 4.0); c = (900.0, 385.0); T = (600.0, 290.0)   # aluminium-ish and copper-ish blocks
        sys = ThermoSystem()
        a = add!(sys, Body(:al, ConstantCapacitySolid(m[1], c[1]); T = T[1]))
        b = add!(sys, Body(:cu, ConstantCapacitySolid(m[2], c[2]); T = T[2]))
        k = connect!(sys, Conductor(ThermalConductance(3.0)), a, b)
        cs = compile(sys)
        C = m .* c
        τ = 1 / (3.0 * (1 / C[1] + 1 / C[2]))
        sol = simulate(cs, (0.0, 40τ))
        uf = sol.u[end]

        Tf = sum(C .* T) / sum(C)                         # energy-weighted
        Tf_entropy_weighted = sum((C ./ T) .* T) / sum(C ./ T)
        @test !(Tf ≈ Tf_entropy_weighted)                 # the two candidates really differ here
        @test temperature(cs, uf, a) ≈ Tf atol = 1e-6
        @test temperature(cs, uf, b) ≈ Tf atol = 1e-6

        Sgen_exact = sum(C .* log.(Tf ./ T))
        @test Sgen_exact > 0
        @test generated_entropy(cs, uf)[k] ≈ Sgen_exact rtol = 1e-7
        @test total_entropy(cs, uf) - total_entropy(cs, initial_state(cs)) ≈ Sgen_exact rtol = 1e-7
        @test total_energy(cs, uf) ≈ total_energy(cs, initial_state(cs)) rtol = 1e-9
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10
    end
end
