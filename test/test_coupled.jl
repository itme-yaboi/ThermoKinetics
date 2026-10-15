# Slice 4: coupled conductors with a thermal row (kernel A, 07 §3.1), Onsager–Casimir, and the
# Veinik → library basis transformations. 08 C4 (Knudsen), C5 (thermoelectric), C6 (Saxén), A5.

@testset "coupled conductors" begin
    N2 = gas_n2()
    cp = N2.cv + N2.R

    @testset "A5: Onsager matrices are validated at construction" begin
        rows = (Thermal(), Electrical())
        @test OnsagerConductance(rows, [2.0 1.0; 1.0 3.0]) isa OnsagerConductance
        @test_throws ArgumentError OnsagerConductance(rows, [2.0 1.0; 0.5 3.0])          # asymmetric, no field
        @test_throws ArgumentError OnsagerConductance(rows, [1.0 2.0; 2.0 1.0])          # indefinite
        W = [0.0 0.4; -0.4 0.0]
        @test_throws ArgumentError OnsagerConductance(rows, [2.0 1.0; 1.0 3.0]; W)       # W needs a field
        @test OnsagerConductance(rows, [2.0 1.0; 1.0 3.0]; W, magnetic = true) isa OnsagerConductance
        @test_throws ArgumentError OnsagerConductance(rows, [2.0 1.0; 1.0 3.0]; W = [0.0 0.4; 0.4 0.0], magnetic = true)
        @test_throws ArgumentError OnsagerConductance((Electrical(), Thermal()), [2.0 1.0; 1.0 3.0])  # thermal row first
        @test_throws ArgumentError ThermoelectricConductance(R = 1.0, K = -0.1, α = 1e-3)
        @test_throws ArgumentError KnudsenBarrier(N2; g = -1.0)
    end

    @testset "coupled kernel identities at large ΔT (thermoelectric between reservoirs)" begin
        rng = MersenneTwister(41)
        for _ in 1:2000
            Ta, Tb = 10.0^(1 + 3rand(rng)), 10.0^(1 + 3rand(rng))
            φa, φb = 10randn(rng), 10randn(rng)
            law = ThermoelectricConductance(R = 10.0^(2rand(rng) - 1), K = 10.0^(2rand(rng) - 2), α = 1e-3 * randn(rng))
            sys = ThermoSystem()
            a = add!(sys, Reservoir(:a, Thermal() => Ta, Electrical() => φa))
            b = add!(sys, Reservoir(:b, Thermal() => Tb, Electrical() => φb))
            k = connect!(sys, Conductor(law; dissipation = Split(rand(rng))), a, b)
            cs = compile(sys)
            u = initial_state(cs)
            r = energy_residual(cs, u)
            @test abs(r.residual) <= 1e-11 * r.scale
            f = flux(cs, u, k)
            @test f.S_gen >= 0                                              # ‖Bᵀx‖², sign-exact
            @test f.S_gen ≈ f.Q_a / Ta + f.Q_b / Tb rtol = 1e-8 atol = 1e-14 * (abs(f.Q_a) / Ta + abs(f.Q_b) / Tb)
            # the engineering law is reproduced exactly at any ΔT
            @test f.J[1] ≈ ((φa - φb) + law.α * (Ta - Tb)) / law.R rtol = 1e-10 atol = 1e-12 * abs(law.α * (Ta - Tb) / law.R)
        end
    end

    @testset "Onsager–Casimir: the antisymmetric part does no dissipation" begin
        M = [2.0 0.5; 0.5 1.0]
        W = [0.0 0.8; -0.8 0.0]
        for (rows, forces_a, forces_b) in (((Electrical(), Volumetric()), (Electrical() => 3.0, Volumetric() => -2e5),
                                                                     (Electrical() => 0.0, Volumetric() => -1e5)),)
            sys = ThermoSystem()
            a = add!(sys, Reservoir(:a, Thermal() => 300.0, forces_a...))
            b = add!(sys, Reservoir(:b, Thermal() => 300.0, forces_b...))
            env = add!(sys, Reservoir(:env; T = 300.0))
            sym = connect!(sys, Conductor(OnsagerConductance(rows, M); dissipation = Into(env)), a, b)
            hall = connect!(sys, Conductor(OnsagerConductance(rows, M; W, magnetic = true); dissipation = Into(env)), a, b)
            cs = compile(sys)
            u = initial_state(cs)
            fs, fh = flux(cs, u, sym), flux(cs, u, hall)
            ΔP = SA[3.0, -1e5]
            @test fh.J ≈ (M + W) * ΔP rtol = 1e-12                       # a transverse (Hall-like) flux …
            @test fh.D ≈ fs.D rtol = 1e-10                               # … that releases no extra power
            @test fh.S_gen ≈ fs.S_gen rtol = 1e-10
        end
    end

    @testset "C5 thermoelectric: Seebeck, Peltier (Kelvin 2), no Thomson heat for constant α (Kelvin 1)" begin
        T0 = 400.0
        # per-leg data in Veinik's (277) variables: J_E = L11 δφ + L12 ΔT/T, J_Q = L21 δφ + L22 ΔT/T
        legs = ((L11 = 2.0, L12 = 2.0 * 2e-4 * T0, L22 = T0 * (0.8 + 2.0 * (2e-4)^2 * T0)),     # p-type, α = +2e-4
                (L11 = 1.5, L12 = -1.5 * 1.5e-4 * T0, L22 = T0 * (0.6 + 1.5 * (1.5e-4)^2 * T0)))   # n-type, α = −1.5e-4
        laws = [thermoelectric_from_veinik(l.L11, l.L12, l.L22; T = T0) for l in legs]
        @test laws[1].α ≈ 2e-4 rtol = 1e-12
        @test laws[1].R ≈ 0.5 rtol = 1e-12
        @test laws[1].K ≈ 0.8 rtol = 1e-12
        @test laws[2].α ≈ -1.5e-4 rtol = 1e-12
        @test_throws ArgumentError thermoelectric_from_veinik(1.0, 10.0, 1.0; T = T0)         # L11 L22 < L12²

        # (278) open circuit, a couple: hot junction H, two cold electrodes held at Tc
        Th, Tc = 500.0, 300.0
        sys = ThermoSystem()
        H = add!(sys, Reservoir(:H, Thermal() => Th, Electrical() => 0.0))
        cold = add!(sys, Reservoir(:cold; T = Tc))
        electrode(name) = add!(sys, Body(name, Composite(ConstantCapacitySolid(1.0, 400.0), LinearCapacitor(1e-3)),
                                         charges_at(Composite(ConstantCapacitySolid(1.0, 400.0), LinearCapacitor(1e-3)),
                                                    (T = Tc,), (φ = 0.0,))))
        C1, C2 = electrode(:C1), electrode(:C2)
        connect!(sys, Conductor(laws[1]), H, C1)                           # coupled: Split(½) by default
        connect!(sys, Conductor(laws[2]), H, C2)
        connect!(sys, Conductor(ThermalConductance(1e6)), C1, cold)
        connect!(sys, Conductor(ThermalConductance(1e6)), C2, cold)
        cs = compile(sys)
        u = steady_state(cs)
        for (C, l, law) in ((C1, legs[1], laws[1]), (C2, legs[2], laws[2]))
            δφ = 0.0 - force(cs, u, C, Electrical())                      # φ_H − φ_C
            ΔT = Th - temperature(cs, u, C)
            @test δφ / ΔT ≈ -l.L12 / (l.L11 * T0) rtol = 1e-9             # (278), exactly at finite ΔT
        end
        couple_V = force(cs, u, C1, Electrical()) - force(cs, u, C2, Electrical())
        @test couple_V ≈ (laws[1].α - laws[2].α) * (Th - Tc) rtol = 1e-5  # electrodes sit within ~0.2 mK of Tc

        # Kelvin 2: isothermal leg carrying current: Peltier heat Π = J_Q/J_e = αT = −(δφ/ΔT)·T
        for law in laws
            sys2 = ThermoSystem()
            a = add!(sys2, Reservoir(:a, Thermal() => T0, Electrical() => 0.1))
            b = add!(sys2, Reservoir(:b, Thermal() => T0, Electrical() => 0.0))
            k = connect!(sys2, Conductor(law), a, b)
            cs2 = compile(sys2)
            f = flux(cs2, initial_state(cs2), k)
            Π = f.Q / f.J[1]
            @test Π ≈ law.α * T0 rtol = 1e-12
            @test Π ≈ -(-law.α) * T0 rtol = 1e-12                         # −(δφ/ΔT)_open · T
        end

        # Kelvin 1 (ε₂ − ε₁ = dΠ/dT − Π/T = 0 for constant α): a discretized leg carrying the same
        # current has exactly the temperature profile of a non-thermoelectric leg (Joule only)
        function leg_profile(α; N = 12, I = 2.0, Rtot = 0.5, Ktot = 0.8)
            sys3 = ThermoSystem()
            body = Composite(ConstantCapacitySolid(0.1, 400.0), LinearCapacitor(1e-4))
            hot = add!(sys3, Reservoir(:hot, Thermal() => Th, Electrical() => I * Rtot - α * (Th - Tc)))
            cld = add!(sys3, Reservoir(:cld, Thermal() => Tc, Electrical() => 0.0))
            cells = [add!(sys3, Body(Symbol(:c, i), body, charges_at(body, (T = 400.0,), (φ = 0.0,)))) for i in 1:N]
            seg = ThermoelectricConductance(R = Rtot / (N + 1), K = Ktot * (N + 1), α = α)
            nodes = vcat([hot], cells, [cld])
            for i in 1:N+1
                connect!(sys3, Conductor(seg), nodes[i], nodes[i+1])
            end
            cs3 = compile(sys3)
            u3 = steady_state(cs3)
            return [temperature(cs3, u3, c) for c in cells], mass_or_current(cs3, u3)
        end
        mass_or_current(cs3, u3) = flux(cs3, u3, first(ThermoKinetics.conductor_names(cs3))).J[1]
        T_te, I_te = leg_profile(3e-4)
        T_0, I_0 = leg_profile(0.0)
        @test I_te ≈ 2.0 rtol = 1e-9
        @test I_0 ≈ 2.0 rtol = 1e-9
        @test T_te ≈ T_0 rtol = 1e-10
    end

    @testset "C4 Knudsen transpiration" begin
        T, p = 350.0, 1e5
        v = N2.R * T / p
        s = charges_at(N2; T, p, m = 1.0)[1]
        # Veinik's barrier, (J_v, J_s) basis, (291)–(301): L₁₂ = −ξp/(2T)
        ξ, χ = 1e-9, 1e-4
        LV = [ξ -ξ*p/(2T); -ξ*p/(2T) χ]
        M = barrier_mass_basis(LV; T, v, s)                               # → (heat, mass) rows of kernel A
        @test M ≈ M' rtol = 1e-14
        @test isposdef(Symmetric(M))
        # the transformed matrix makes Veinik's predictions:
        # zero mass flux ⇒ Δp/ΔT = p/(2T)   (299)
        # transport heat at uniform T: J_Q/J_m − T s = Q* = −RT/2   (298)
        M_QQ, M_Qm, M_mm = M[1, 1], M[1, 2], M[2, 2]
        Δμ_over_ΔT = -(M_Qm / M_mm) / T                                    # from J_m = M_mm Δμ/T + M_mQ ΔT/T² = 0
        @test (Δμ_over_ΔT + s) / v ≈ p / (2T) rtol = 1e-10                 # Δμ = vΔp − sΔT
        @test M_Qm / M_mm - T * s ≈ -N2.R * T / 2 rtol = 1e-10
        # the exact KnudsenBarrier law coincides with it near equilibrium (ξ = g v²/T)
        g = ξ * T / v^2
        kb = KnudsenBarrier(N2; g, K = 0.0)
        μ = forces(N2, charges_at(N2; T, p, m = 1.0))[3]
        Mk = ThermoKinetics.coupled_matrix(kb, T, T, SA[μ], SA[μ], 0.5)
        @test Mk[2, 2] ≈ M_mm rtol = 1e-12
        @test Mk[1, 2] ≈ M_Qm rtol = 1e-10

        # simulation: two vessels held at T′ = 400 K and T″ = 300 K, joined by the barrier
        sys = ThermoSystem()
        a = add!(sys, Body(:a, N2; T = 400.0, p = 1e5, V = 1e-3))
        b = add!(sys, Body(:b, N2; T = 300.0, p = 1e5, V = 1e-3))
        ra = add!(sys, Reservoir(:ra; T = 400.0))
        rb = add!(sys, Reservoir(:rb; T = 300.0))
        connect!(sys, Conductor(ThermalConductance(100.0)), ra, a)     # the barrier conducts K ΔT = 0.1 W
        connect!(sys, Conductor(ThermalConductance(100.0)), rb, b)
        k = connect!(sys, Conductor(KnudsenBarrier(N2; g = 1e-8, K = 1e-3)), a, b)
        cs = compile(sys)
        sol = simulate(cs, (0.0, 2e4))
        uf = sol.u[end]
        Ta, Tb = temperature(cs, uf, a), temperature(cs, uf, b)        # within 1 mK of 400 and 300 K
        @test Ta ≈ 400.0 rtol = 1e-5
        @test Tb ≈ 300.0 rtol = 1e-5
        @test pressure(cs, uf, a) / pressure(cs, uf, b) ≈ sqrt(Ta / Tb) rtol = 1e-8    # (300), exactly
        @test abs(flux(cs, uf, k).J[1]) < 1e-14
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9 && cl.charge_rel <= 1e-10
    end

    @testset "C6 electrokinetic barrier: Saxén in the (charge, mass) basis" begin
        v = 1e-3                                                          # liquid, m³/kg
        L = [2e-3 3e-9; 3e-9 5e-12]                                        # Veinik (J_e, J_v) vs (Δφ, Δp)
        M = electrokinetic_mass_basis(L; v)
        @test M ≈ [L[1, 1] L[1, 2]/v; L[1, 2]/v L[2, 2]/v^2] rtol = 1e-14
        @test isposdef(Symmetric(M))
        # streaming potential (Δφ/Δp) at J_e = 0, and electro-osmosis (J_v/J_e) at Δp = 0, in both bases
        stream_V = -L[1, 2] / L[1, 1]                                     # (319)
        osmo_V = L[2, 1] / L[1, 1]                                        # (320)
        stream_M = (-M[1, 2] / M[1, 1]) * v                               # Δφ/Δμ · dμ/dp, dμ = v dp at fixed T
        osmo_M = (M[2, 1] / M[1, 1]) * v                                  # J_m/J_e · v
        @test stream_M ≈ stream_V rtol = 1e-12
        @test osmo_M ≈ osmo_V rtol = 1e-12
        @test stream_M ≈ -osmo_M rtol = 1e-12                             # Saxén (323)
        # and the assembled conductor produces those fluxes
        sys = ThermoSystem()
        a = add!(sys, Reservoir(:a, Thermal() => 300.0, Electrical() => 1.0, Species(:w) => 50.0))
        b = add!(sys, Reservoir(:b, Thermal() => 300.0, Electrical() => 0.0, Species(:w) => 50.0))
        k = connect!(sys, Conductor(OnsagerConductance((Electrical(), Species(:w)), M); dissipation = Split(0.5)), a, b)
        cs = compile(sys)
        f = flux(cs, initial_state(cs), k)
        @test (f.J[2] / f.J[1]) * v ≈ osmo_V rtol = 1e-12                 # electro-osmotic mass flow at Δμ = 0
    end
end
