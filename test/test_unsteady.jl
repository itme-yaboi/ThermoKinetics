using SpecialFunctions: erfc

@testset "08 §D unsteady" begin
    @testset "D1 lumped relaxation (thermal analogue)" begin
        # 03 §3: for the thermal DOF nothing is "dissipated"; the analogous check is the
        # entropy produced, which (like ½KX₀²) is independent of the conductance.
        m, c, T0, Tr = 2.0, 400.0, 500.0, 300.0
        C = m * c
        Sgen_exact = C * log(Tr / T0) + C * (T0 - Tr) / Tr
        for G in (5.0, 50.0)
            sys = ThermoSystem()
            b = add!(sys, Body(:b, ConstantCapacitySolid(m, c); T = T0))
            r = add!(sys, Reservoir(:r; T = Tr))
            k = connect!(sys, Conductor(ThermalConductance(G)), b, r)
            cs = compile(sys)
            τ = C / G
            sol = simulate(cs, (0.0, 40τ); saveat = τ / 4)
            for (t, u) in zip(sol.t, sol.u)
                @test temperature(cs, u, b) - Tr ≈ (T0 - Tr) * exp(-t / τ) atol = 1e-6
            end
            @test generated_entropy(cs, sol.u[end])[k] ≈ Sgen_exact rtol = 1e-7
            cl = closure(cs, sol)
            @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10
        end
    end

    @testset "D4 entropy-form diffusivity (D3, thermal)" begin
        # semi-infinite slab, step change of surface temperature; the state is entropy,
        # and the solution must follow the erfc profile with D = λ/(ρc) (= χ/(ρκ), (371))
        λ, ρ, c = 50.0, 8000.0, 500.0
        D = λ / (ρ * c)
        Ts, Ti, t_end = 400.0, 300.0, 100.0
        L, n = 0.3, 300                                  # 2√(Dt) ≈ 0.07 m ≪ L
        sys = ThermoSystem()
        r = add!(sys, Reservoir(:surface; T = Ts))
        rd = rod!(sys, :slab; n, L, area = 1.0, λ, ρ, c, T0 = Ti, left = r, right = nothing)
        cs = compile(sys)
        sol = simulate(cs, (0.0, t_end); saveat = [25.0, t_end])
        for (t, u) in zip(sol.t, sol.u)
            err = maximum(abs(temperature(cs, u, cell) - (Ti + (Ts - Ti) * erfc(x / (2 * sqrt(D * t)))))
                          for (cell, x) in zip(rd.cells, rd.x))
            @test err <= 0.005 * (Ts - Ti)               # spatial discretization, Δx = 1 mm
        end
        cl = closure(cs, sol)                            # saveat omits t = 0: closure must use prob.u0
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10
        # the entropy formulation's diffusivity is the energy formulation's: χ/(ρκ) = λ/(ρc)
        T = 350.0
        χ, κ = λ / T, c / T
        @test χ / (ρ * κ) ≈ D rtol = 1e-14
    end

    @testset "D5 Biot criterion (slab)" begin
        # slab of thickness L, insulated on one face, convective (Robin) on the other;
        # Bi = hL/λ. Lumped model: θ = exp(-Bi·Fo).
        λ, ρ, c, L, n = 50.0, 8000.0, 500.0, 0.05, 50
        T0, Tf = 500.0, 300.0
        α = λ / (ρ * c)
        maxdiff = Dict{Float64,Float64}()
        for Bi in (0.05, 0.1, 10.0)
            h = Bi * λ / L
            sys = ThermoSystem()
            fluid = add!(sys, Reservoir(:fluid; T = Tf))
            rd = rod!(sys, :slab; n, L, area = 1.0, λ, ρ, c, T0, left = (fluid, h), right = nothing)
            cs = compile(sys)
            Fo_end = 3 / min(Bi, 1.0)
            t_end = Fo_end * L^2 / α
            sol = simulate(cs, (0.0, t_end); saveat = t_end / 100)
            d = 0.0
            for (t, u) in zip(sol.t, sol.u)
                θmean = (sum(temperature(cs, u, cell) for cell in rd.cells) / n - Tf) / (T0 - Tf)
                d = max(d, abs(θmean - exp(-Bi * α * t / L^2)))
            end
            maxdiff[Bi] = d
            cl = closure(cs, sol)
            @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10
        end
        @test maxdiff[0.05] < 0.05
        @test maxdiff[0.1] < 0.05
        @test maxdiff[10.0] > 0.05      # the lumped model fails where the criterion says it should
    end
end

@testset "08 §D unsteady: non-thermal" begin
    @testset "D1 lumped relaxation (volumetric): half-energy loss for any conductance" begin
        # an elastic volume (capacity K = C) charged through a damped port from a pressure
        # reservoir; friction heat goes to a thermal sink at T
        C, V0, p0, pr, Tsink = 2e-9, 1e-3, 1e5, 4e5, 300.0
        X0 = pr - p0
        for L in (1e-10, 1e-8)
            body = ElasticVolume(C; V0, p0)
            sys = ThermoSystem()
            x = add!(sys, Body(:x, body, charges_at(body; p = p0)))
            r = add!(sys, Reservoir(:r; p = pr))
            sink = add!(sys, Reservoir(:sink; T = Tsink))
            k = connect!(sys, Conductor(VolumeConductance(L); dissipation = Into(sink)), r, x)
            cs = compile(sys)
            τ = C / L
            sol = simulate(cs, (0.0, 40τ); saveat = τ / 4)
            for (t, u) in zip(sol.t, sol.u)
                @test pr - pressure(cs, u, x) ≈ X0 * exp(-t / τ) atol = 1e-6 * X0
            end
            led = reservoir_ledger(cs, sol.u[end], sink)
            @test led.U ≈ 0.5 * C * X0^2 rtol = 1e-7          # dissipated energy ½ K X₀²
            @test generated_entropy(cs, sol.u[end])[k] ≈ 0.5 * C * X0^2 / Tsink rtol = 1e-7
            cl = closure(cs, sol)
            @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10 && cl.charge_rel <= 1e-10
        end
    end
end
