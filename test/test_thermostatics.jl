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

@testset "08 §B thermostatics: gases" begin
    He, N2 = gas_he(), gas_n2()

    function two_gases(g1, T1, p1, V1, g2, T2, p2, V2)
        sys = ThermoSystem()
        a = add!(sys, Body(:a, g1; T = T1, V = V1, m = p1 * V1 / (g1.R * T1)))
        b = add!(sys, Body(:b, g2; T = T2, V = V2, m = p2 * V2 / (g2.R * T2)))
        piston = connect!(sys, Conductor(VolumeConductance(5e-9); dissipation = InternalNode(G_a = 0.05, G_b = 0.05)), a, b)
        wall = connect!(sys, Conductor(ThermalConductance(0.05)), a, b)
        return compile(sys), a, b, piston, wall
    end

    @testset "B4 two-gas joining (piston + diathermal wall)" begin
        cs, a, b, piston, wall = two_gases(He, 600.0, 3e5, 1e-3, N2, 300.0, 1e5, 2e-3)
        u0 = initial_state(cs)
        sol = simulate(cs, (0.0, 3000.0))
        uf = sol.u[end]
        m = (charges(cs, u0, a)[3], charges(cs, u0, b)[3])
        Tf = (m[1] * He.cv * 600.0 + m[2] * N2.cv * 300.0) / (m[1] * He.cv + m[2] * N2.cv)
        pf = (m[1] * He.R + m[2] * N2.R) * Tf / 3e-3
        for x in (a, b)
            @test temperature(cs, uf, x) ≈ Tf rtol = 1e-7       # T equalizes
            @test pressure(cs, uf, x) ≈ pf rtol = 1e-7          # p equalizes
        end
        # both charges of each body moved
        @test abs(charges(cs, uf, a)[1] - charges(cs, u0, a)[1]) > 1e-3 * abs(m[1] * He.cv)
        @test abs(charges(cs, uf, a)[2] - charges(cs, u0, a)[2]) > 1e-2 * 1e-3
        Sgen = sum(values(generated_entropy(cs, uf)))
        @test Sgen > 0
        @test total_entropy(cs, uf) - total_entropy(cs, u0) ≈ Sgen rtol = 1e-9
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10 && cl.charge_rel <= 1e-10
    end

    @testset "B4/§23: equal pressures at the start, yet volume still moves" begin
        p = 1e5
        cs, a, b, _, _ = two_gases(N2, 500.0, p, 1e-3, N2, 300.0, p, 1e-3)
        u0 = initial_state(cs)
        uf = simulate(cs, (0.0, 3000.0)).u[end]
        m = (charges(cs, u0, a)[3], charges(cs, u0, b)[3])
        Tf = (m[1] * 500.0 + m[2] * 300.0) / (m[1] + m[2])
        @test pressure(cs, uf, a) ≈ p rtol = 1e-7                # same gas: p ends where it began …
        Va_f = m[1] * N2.R * Tf / p
        @test charges(cs, uf, a)[2] ≈ Va_f rtol = 1e-7          # … but the hot gas gave up volume
        @test Va_f < 1e-3 * (1 - 0.05)
    end

    @testset "an uncoupled body is a pure accumulator (§11): no conversion" begin
        # solid + elastic volume, compressed by a high-pressure reservoir; friction heat goes
        # to a separate sink, so nothing may reach the body's thermal DOF
        body = Composite(ConstantCapacitySolid(1.0, 500.0), ElasticVolume(1e-9; V0 = 1e-3, p0 = 1e5))
        sys = ThermoSystem()
        x = add!(sys, Body(:x, body, charges_at(body, (T = 300.0,), (p = 1e5,))))
        hp = add!(sys, Reservoir(:hp; p = 3e5))
        sink = add!(sys, Reservoir(:sink; T = 300.0))
        connect!(sys, Conductor(VolumeConductance(1e-9); dissipation = Into(sink)), hp, x)
        cs = compile(sys)
        u0 = initial_state(cs)
        uf = simulate(cs, (0.0, 20.0)).u[end]
        @test charges(cs, uf, x)[1] == charges(cs, u0, x)[1]     # entropy untouched, exactly
        @test pressure(cs, uf, x) ≈ 3e5 rtol = 1e-7
        # contrast: an ideal gas compressed the same way heats up (coupling converts)
        sys2 = ThermoSystem()
        g = add!(sys2, Body(:g, N2; T = 300.0, p = 1e5, m = 1e-3))
        hp2 = add!(sys2, Reservoir(:hp; p = 3e5))
        sink2 = add!(sys2, Reservoir(:sink; T = 300.0))
        connect!(sys2, Conductor(VolumeConductance(1e-9); dissipation = Into(sink2)), hp2, g)
        cs2 = compile(sys2)
        # its entropy is untouched too, so the compression is exactly isentropic
        Tg = temperature(cs2, simulate(cs2, (0.0, 20.0)).u[end], g)
        @test Tg ≈ 300.0 * 3.0^(N2.R / (N2.cv + N2.R)) rtol = 1e-7
    end
end
