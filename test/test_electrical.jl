# Slice 4: the electrical domain, capacitors, resistors (kernel B or Into), mass conductors on μ.
# 08 C3, D1 (electrical), D2; the contact rule for mass (μ, not p or concentration).

@testset "electrical domain, resistors, mass conductors" begin
    @testset "construction and validation" begin
        @test charge_symbol(Electrical()) == :q && charge_unit(Electrical()) == "C"
        @test force_symbol(Electrical()) == :φ && force_unit(Electrical()) == "V"
        cap = LinearCapacitor(2e-6)
        E = charges_at(cap; φ = 5.0)
        @test E[1] ≈ 1e-5 rtol = 1e-14
        @test forces(cap, E)[1] ≈ 5.0 rtol = 1e-14                   # φ = ∂U/∂q = q/C
        @test internal_energy(cap, E) ≈ 0.5 * 2e-6 * 25 rtol = 1e-14
        @test_throws ArgumentError LinearCapacitor(0.0)
        @test_throws ArgumentError Resistor(-1.0)
        @test_throws ArgumentError Resistor(0.0)                       # an ideal short is not a conductor
        @test_throws ArgumentError MassConductance(:N2, -1.0)
        @test_throws ArgumentError Conductor(Resistor(10.0))           # dissipates: needs a policy
    end

    @testset "C3 Joule heating (Into a sink, and through an internal node)" begin
        V, R, T = 12.0, 4.0, 300.0
        for (label, policy) in ((:into, env -> Into(env)), (:node, env -> InternalNode(G_h = 0.5, h = env)))
            sys = ThermoSystem()
            hi = add!(sys, Reservoir(:hi, Electrical() => V))
            lo = add!(sys, Reservoir(:lo, Electrical() => 0.0))
            env = add!(sys, Reservoir(:env; T))
            k = connect!(sys, Conductor(Resistor(R); dissipation = policy(env)), hi, lo)
            cs = compile(sys)
            u = initial_state(cs)
            du = rates(cs, u)
            @test reservoir_charge(cs, du, lo, Electrical()) ≈ V / R rtol = 1e-14      # current
            @test reservoir_ledger(cs, du, env).U ≈ V^2 / R rtol = 1e-12             # Q̇ = ΔV²/R
            @test entropy_production(cs, u)[k] ≈ V^2 / (R * T) rtol = 1e-12          # Ṡ_gen = ΔV²/(RT)
            r = energy_residual(cs, u)
            @test abs(r.residual) <= 1e-12 * r.scale
            if label === :node                     # the resistor's core runs hot: T_m = T + P/G
                @test flux(cs, u, k).T_m ≈ T + (V^2 / R) / 0.5 rtol = 1e-12
                parts = entropy_production_parts(cs, u)[k]
                @test parts.dissipation ≈ (V^2 / R) / (T + (V^2 / R) / 0.5) rtol = 1e-12
                @test parts.conduction > 0
            end
        end
    end

    @testset "D1 RC charging (electrical): half the energy is dissipated, for any R" begin
        C, V0, T = 1e-3, 10.0, 300.0
        for R in (1.0, 100.0)
            sys = ThermoSystem()
            cap = add!(sys, Body(:cap, LinearCapacitor(C); φ = 0.0))
            src = add!(sys, Reservoir(:src, Electrical() => V0))
            env = add!(sys, Reservoir(:env; T))
            k = connect!(sys, Conductor(Resistor(R); dissipation = Into(env)), src, cap)
            cs = compile(sys)
            τ = R * C
            sol = simulate(cs, (0.0, 40τ); saveat = τ / 2)
            for (t, u) in zip(sol.t, sol.u)
                @test force(cs, u, cap, Electrical()) ≈ V0 * (1 - exp(-t / τ)) atol = 1e-7 * V0
            end
            uf = sol.u[end]
            @test reservoir_ledger(cs, uf, env).U ≈ 0.5 * C * V0^2 rtol = 1e-7       # ½ C V₀²
            @test generated_entropy(cs, uf)[k] ≈ 0.5 * C * V0^2 / T rtol = 1e-7
            cl = closure(cs, sol)
            @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10 && cl.charge_rel <= 1e-10
        end
    end

    @testset "D2 two-capacitor paradox" begin
        C, V0, T = 1e-3, 10.0, 300.0
        for R in (0.5, 50.0)
            sys = ThermoSystem()
            a = add!(sys, Body(:a, LinearCapacitor(C); φ = V0))
            b = add!(sys, Body(:b, LinearCapacitor(C); φ = 0.0))
            env = add!(sys, Reservoir(:env; T))
            k = connect!(sys, Conductor(Resistor(R); dissipation = Into(env)), a, b)
            cs = compile(sys)
            sol = simulate(cs, (0.0, 40 * R * C / 2))
            uf = sol.u[end]
            @test force(cs, uf, a, Electrical()) ≈ V0 / 2 rtol = 1e-7
            @test force(cs, uf, b, Electrical()) ≈ V0 / 2 rtol = 1e-7
            @test reservoir_ledger(cs, uf, env).U ≈ 0.25 * C * V0^2 rtol = 1e-7      # half of ½CV₀²
            @test generated_entropy(cs, uf)[k] ≈ 0.25 * C * V0^2 / T rtol = 1e-7
            @test total_charge(cs, uf, Electrical()) ≈ C * V0 rtol = 1e-12           # charge conserved
            cl = closure(cs, sol)
            @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10 && cl.charge_rel <= 1e-10
        end
    end

    @testset "mass conductor on μ: equilibrium equalizes μ, not p (contact rule)" begin
        N2 = gas_n2()
        sys = ThermoSystem()
        a = add!(sys, Body(:a, N2; T = 400.0, p = 2e5, V = 1e-3))
        b = add!(sys, Body(:b, N2; T = 300.0, p = 1e5, V = 1e-3))
        env = add!(sys, Reservoir(:env; T = 300.0))
        k = connect!(sys, Conductor(MassConductance(:N2, 1e-9); dissipation = Into(env)), a, b)
        cs = compile(sys)
        u0 = initial_state(cs)
        sol = simulate(cs, (0.0, 2e4))
        uf = sol.u[end]
        μa, μb = force(cs, uf, a, Species(:N2)), force(cs, uf, b, Species(:N2))
        @test μa ≈ μb rtol = 1e-7
        @test !isapprox(pressure(cs, uf, a), pressure(cs, uf, b); rtol = 1e-2)      # p differs at equilibrium
        @test total_charge(cs, uf, Species(:N2)) ≈ total_charge(cs, u0, Species(:N2)) rtol = 1e-13
        @test all(≥(0), values(entropy_production(cs, uf)))
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9 && cl.charge_rel <= 1e-10
    end
end
