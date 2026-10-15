# Slice 3: streams (advected charge bundles), hydraulic drivers, ideal pumps, labeled flow
# boundaries and pass-through outlets (DESIGN §7). 08 B5, C7; A1–A3, A6 with streams;
# a pumped loop checked against Gouy–Stodola.

@testset "streams, pumps and flow boundaries" begin
    N2, He = gas_n2(), gas_he()
    cp = N2.cv + N2.R
    k = cp / N2.cv

    "Specific entropy and chemical potential of an ideal gas at (T, p), via the material."
    function sμ(g, T, p)
        E = charges_at(g; T, p, m = 1.0)
        return E[1], forces(g, E)[3]
    end

    @testset "construction and validation" begin
        @test_throws ArgumentError HydraulicConductance(-1e-8)
        @test_throws ArgumentError FlowBoundary(N2; T = 300.0, p = 2e5, mdot = t -> 1e-3)          # no supply label
        @test_throws ArgumentError FlowBoundary(N2; T = 300.0, p = 2e5, mdot = t -> 1e-3, supply = :weird)
        @test_throws ArgumentError PassThroughOutlet()                                           # needs mdot or (κ, p_out)

        air = Reservoir(:air; T = 300.0, p = 1e5, material = N2)
        @test domains(air) == (Thermal(), Volumetric(), Species(:N2))
        @test ThermoKinetics.reservoir_force(air, Species(:N2)) ≈ sμ(N2, 300.0, 1e5)[2] rtol = 1e-12

        sys = ThermoSystem()
        a = add!(sys, Body(:a, N2; T = 300.0, p = 2e5, m = 1e-3))
        b = add!(sys, Body(:b, N2; T = 300.0, p = 1e5, m = 1e-3))
        h = add!(sys, Body(:h, He; T = 300.0, p = 1e5, m = 1e-4))
        s = add!(sys, Body(:s, ConstantCapacitySolid(1.0, 500.0); T = 300.0))
        bare = add!(sys, Reservoir(:bare; T = 300.0, p = 1e5))              # no material: cannot supply fluid
        connect!(sys, Stream(HydraulicConductance(1e-8)), a, b)
        @test compile(sys) isa CompiledSystem
        for (x, y) in ((a, h), (a, s), (a, bare))
            sys2 = deepcopy(sys)
            connect!(sys2, Stream(HydraulicConductance(1e-8)), x, y)
            @test_throws ArgumentError compile(sys2)
        end
    end

    @testset "stream bookkeeping identities (both flow directions)" begin
        rng = MersenneTwister(5)
        for _ in 1:200
            sys = ThermoSystem()
            Ta, Tb = 200 + 600rand(rng), 200 + 600rand(rng)
            pa, pb = 0.5e5 + 3e5rand(rng), 0.5e5 + 3e5rand(rng)
            a = add!(sys, Body(:a, N2; T = Ta, p = pa, m = 1e-3))
            b = add!(sys, Body(:b, N2; T = Tb, p = pb, m = 2e-3))
            st = connect!(sys, Stream(HydraulicConductance(1e-8)), a, b)
            cs = compile(sys)
            u = initial_state(cs)
            du = rates(cs, u)
            mdot = 1e-8 * (pa - pb)
            @test charges(cs, du, a)[3] ≈ -mdot rtol = 1e-12          # mass leaves the upstream …
            @test charges(cs, du, b)[3] ≈ mdot rtol = 1e-12           # … and arrives downstream
            @test charges(cs, du, a)[2] == 0 && charges(cs, du, b)[2] == 0   # V is not advected
            r = energy_residual(cs, u)
            @test abs(r.residual) <= 1e-11 * r.scale
            # generation: ṁ[s_u(T_u − T_d) + μ_u − μ_d]/T_d, upwinded, ≥ 0 (flow runs down Δp)
            (u_, d_, T_u, T_d, p_u, p_d) = mdot >= 0 ? (a, b, Ta, Tb, pa, pb) : (b, a, Tb, Ta, pb, pa)
            s_u, μ_u = sμ(N2, T_u, p_u)
            _, μ_d = sμ(N2, T_d, p_d)
            Sgen = entropy_production(cs, u)[st]
            @test Sgen ≈ abs(mdot) * (s_u * (T_u - T_d) + μ_u - μ_d) / T_d rtol = 1e-8
            @test Sgen >= -1e-14 * abs(mdot * s_u)
            @test Sgen >= abs(mdot) * (N2.R * T_u / p_u) * (p_u - p_d) / T_d * (1 - 1e-9)   # concavity bound
        end
    end

    @testset "ideal pump: isentropic delivery, logged work" begin
        sys = ThermoSystem()
        a = add!(sys, Body(:a, N2; T = 300.0, p = 1e5, m = 1e-2))
        b = add!(sys, Body(:b, N2; T = 350.0, p = 3e5, m = 1e-2))
        mdot = 2e-4
        pmp = connect!(sys, IdealPump(t -> mdot), a, b)                  # uphill: work is supplied
        cs = compile(sys)
        u = initial_state(cs)
        du = rates(cs, u)
        W = pump_work(cs, du)[pmp]
        @test W ≈ mdot * cp * 300.0 * (3.0^(N2.R / cp) - 1) rtol = 1e-10   # ṁ (h(s_a, p_b) − h_a)
        r = energy_residual(cs, u)
        @test abs(r.residual) <= 1e-11 * r.scale                        # residual counts the work
        T_star = 300.0 * 3.0^(N2.R / cp)
        Sgen = entropy_production(cs, u)[pmp]
        @test Sgen ≈ mdot * cp * (T_star - 350.0 - 350.0 * log(T_star / 350.0)) / 350.0 rtol = 1e-8
        @test Sgen >= 0
    end

    @testset "passive boundary below the body pressure is an error" begin
        sys = ThermoSystem()
        cv = add!(sys, Body(:cv, N2; T = 300.0, p = 2e5, V = 1e-3))
        connect!(sys, FlowBoundary(N2; T = 300.0, p = 1e5, mdot = t -> 1e-4, supply = :passive), cv)
        cs = compile(sys)
        @test_throws ArgumentError rates(cs, initial_state(cs))
        # the same inflow, labeled pumped, is allowed and logs its (here negative) work
        sys2 = ThermoSystem()
        cv2 = add!(sys2, Body(:cv, N2; T = 300.0, p = 2e5, V = 1e-3))
        fb = connect!(sys2, FlowBoundary(N2; T = 300.0, p = 1e5, mdot = t -> 1e-4, supply = :pumped), cv2)
        cs2 = compile(sys2)
        @test pump_work(cs2, rates(cs2, initial_state(cs2)))[fb] > 0     # pumped up from 1 to 2 bar
    end

    @testset "B5 adiabatic tank filling" begin
        V, T_in, p_in = 1e-2, 300.0, 5e5
        for p0 in (1.5e5, 1.0)                              # partly filled, and near-evacuated (p₀ = 1 Pa)
            T0 = 300.0
            sys = ThermoSystem()
            line = add!(sys, Reservoir(:line; T = T_in, p = p_in, material = N2))
            tank = add!(sys, Body(:tank, N2; T = T0, p = p0, V))
            connect!(sys, Stream(HydraulicConductance(1e-7)), line, tank)
            cs = compile(sys)
            u0 = initial_state(cs)
            m1 = charges(cs, u0, tank)[3]
            sol = simulate(cs, (0.0, 200.0))
            uf = sol.u[end]
            # first law for the control volume: m₂u₂ − m₁u₁ = (m₂ − m₁) h_in, with p₂ = p_in
            m2 = m1 + (p_in * V / N2.R - m1 * T0) * N2.cv / (cp * T_in)
            T2 = p_in * V / (N2.R * m2)
            @test pressure(cs, uf, tank) ≈ p_in rtol = 1e-7
            @test charges(cs, uf, tank)[3] ≈ m2 rtol = 1e-6
            @test temperature(cs, uf, tank) ≈ T2 rtol = 1e-6
            p0 == 1.0 && @test temperature(cs, uf, tank) / T_in ≈ k rtol = 1e-4   # evacuated limit T₂ = k T_in
            Sgen = sum(values(generated_entropy(cs, uf)))
            @test Sgen > 0
            cl = closure(cs, sol)
            @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9 && cl.charge_rel <= 1e-10
        end
    end

    @testset "C7 stream mixing (steady)" begin
        T1, T2, m1, m2, p_in, p_out = 400.0, 300.0, 2e-3, 3e-3, 2e5, 1e5
        sys = ThermoSystem()
        cv = add!(sys, Body(:cv, N2; T = 350.0, p = 1.5e5, V = 1e-3))
        f1 = connect!(sys, FlowBoundary(N2; T = T1, p = p_in, mdot = t -> m1, supply = :passive), cv)
        f2 = connect!(sys, FlowBoundary(N2; T = T2, p = p_in, mdot = t -> m2, supply = :passive), cv)
        out = connect!(sys, PassThroughOutlet(κ = 1e-7, p_out = p_out), cv)
        cs = compile(sys)
        u = steady_state(cs)
        T = (m1 * T1 + m2 * T2) / (m1 + m2)                  # energy balance (same c_p)
        p = p_out + (m1 + m2) / 1e-7                         # outlet driver at steady state
        @test temperature(cs, u, cv) ≈ T rtol = 1e-9
        @test pressure(cs, u, cv) ≈ p rtol = 1e-9
        # entropy generated at the inlet faces = closed-form mixing (and throttling) entropy rate
        Sgen = entropy_production(cs, u)
        expected = m1 * (cp * log(T / T1) - N2.R * log(p / p_in)) + m2 * (cp * log(T / T2) - N2.R * log(p / p_in))
        @test Sgen[f1] + Sgen[f2] ≈ expected rtol = 1e-8
        @test Sgen[f1] > 0 && Sgen[f2] > 0
        @test !haskey(Sgen, out) || Sgen[out] == 0           # a pass-through outlet produces nothing
        # the steady ledgers balance: entropy out = entropy in + generation
        du = rates(cs, u)
        bin = boundary_ledger(cs, du, f1).S + boundary_ledger(cs, du, f2).S
        bout = boundary_ledger(cs, du, out).S
        @test bin + bout ≈ Sgen[f1] + Sgen[f2] rtol = 1e-8
    end

    @testset "pumped loop: pump work = T₀ Σ Ṡ_gen (Gouy–Stodola)" begin
        T0 = 300.0
        sys = ThermoSystem()
        a = add!(sys, Body(:a, N2; T = T0, p = 1e5, V = 5e-3))
        b = add!(sys, Body(:b, N2; T = T0, p = 1e5, V = 5e-3))
        env = add!(sys, Reservoir(:env; T = T0))
        pmp = connect!(sys, IdealPump(t -> 5e-4), a, b)
        connect!(sys, Stream(HydraulicConductance(5e-9)), b, a)
        connect!(sys, Conductor(ThermalConductance(0.5)), b, env)
        connect!(sys, Conductor(ThermalConductance(0.5)), a, env)
        cs = compile(sys)
        sol = simulate(cs, (0.0, 2000.0))
        uf = sol.u[end]
        du = rates(cs, uf)
        W = pump_work(cs, du)[pmp]
        Qout = reservoir_ledger(cs, du, env).U
        @test W ≈ Qout rtol = 1e-7                           # first law: work in = heat out
        @test W ≈ T0 * sum(values(entropy_production(cs, uf))) rtol = 1e-7
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9 && cl.charge_rel <= 1e-10
    end

    @testset "08 §A with streams: random gas networks" begin
        rng = MersenneTwister(91)
        for trial in 1:10
            sys = ThermoSystem()
            n = rand(rng, 2:5)
            bodies = [add!(sys, Body(Symbol(:g, i), N2; T = 250 + 500rand(rng), p = 0.5e5 + 2e5rand(rng),
                                     V = 1e-3 * (0.5 + rand(rng)))) for i in 1:n]
            pairs = vcat([(rand(rng, 1:i-1), i) for i in 2:n], [(1, n)])
            links = Symbol[]
            for (i, j) in pairs
                i == j && continue
                push!(links, connect!(sys, Stream(HydraulicConductance(1e-9 * (1 + 9rand(rng)))), bodies[i], bodies[j]))
                connect!(sys, Conductor(ThermalConductance(0.01 + 0.09rand(rng))), bodies[i], bodies[j])
            end
            cs = compile(sys)
            u0 = initial_state(cs)
            sol = simulate(cs, (0.0, 3000.0); saveat = 30.0)
            cl = closure(cs, sol)
            @test cl.energy_rel <= 1e-6                                        # A1
            @test cl.charge_rel <= 1e-10                                       # A2: mass and volume
            @test cl.entropy_rel <= 1e-9
            S = [total_entropy(cs, u) for u in sol.u]
            @test all(diff(S) .>= -1e-9 * maximum(abs, S))                     # A3
            @test all(minimum(values(entropy_production(cs, u))) >= -1e-12 for u in sol.u[1:10:end])
            uf = sol.u[end]
            for x in bodies                                                    # A6: T and p equal
                @test temperature(cs, uf, x) ≈ temperature(cs, uf, bodies[1]) rtol = 1e-6
                @test pressure(cs, uf, x) ≈ pressure(cs, uf, bodies[1]) rtol = 1e-6
            end
        end
    end
end
