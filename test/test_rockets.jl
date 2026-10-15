# Slice 8: rockets (DESIGN §19). Moving variable-mass bodies, momentum-carrying streams and pumps,
# species-subset feeds, the ideal choked nozzle (chamber on a vehicle or on a test stand), and
# gravity on a moving body. Kernel and property tests at single states (no long transients).

@testset "rockets: moving bodies, momentum streams, nozzle" begin
    H2, O2, H2O, N2s = H2g(), O2g(), H2Og(), N2g()
    mix4() = IdealGasMixture(H2, O2, H2O, N2s)
    g0 = 9.81

    "Per-kg k, R, c_p of a mixture from species data, independent of the library."
    function frozen_props(sds, y)
        cv = sum(y[i] * sds[i].cv for i in eachindex(sds))
        R = sum(y[i] * R̄ / molar_mass(sds[i]) for i in eachindex(sds))
        return (; cv, R, cp = cv + R, k = (cv + R) / cv)
    end
    "Area–Mach function A/A* (Veinik (454))."
    areamach(M, k) = (1 / M) * ((2 / (k + 1)) * (1 + (k - 1) / 2 * M^2))^((k + 1) / (2 * (k - 1)))

    "A random Moving mixture chamber state: (body kwargs, species masses, inner material, m_dry, w)."
    function random_chamber(rng)
        mix = mix4()
        sds = (H2, O2, H2O, N2s)
        T = 300 + 3200rand(rng); p = 10.0^(5 + 2rand(rng))
        m = Tuple(10.0^(-3 + 2rand(rng)) for _ in 1:4)
        V = sum(m[i] * R̄ / molar_mass(sds[i]) for i in 1:4) * T / p
        w = 6000 * (rand(rng) - 0.5)
        m_dry = 5rand(rng)
        return (; mix, sds, T, p, m, V, w, m_dry, At = 1e-3 * (0.5 + rand(rng)),
                  ε = 1.5 + 48.5rand(rng), η = rand(rng) < 0.5 ? 1.0 : 0.9)
    end

    function chamber_system(c; gravity = false, nozzle = true)
        sys = ThermoSystem()
        add!(sys, Body(:chamber, Moving(c.mix; m_dry = c.m_dry); T = c.T, V = c.V,
                       m = (H2 = c.m[1], O2 = c.m[2], H2O = c.m[3], N2 = c.m[4]), w = c.w))
        nozzle && connect!(sys, Nozzle(c.At, c.ε; η = c.η), :chamber; name = :nozzle)
        gravity && connect!(sys, GravityField(g0), :chamber; name = :gravity)
        return compile(sys)
    end

    @testset "Moving material: state function" begin
        rng = MersenneTwister(81)
        for _ in 1:200, inner in (gas_n2(), mix4())
            m_dry = 5rand(rng); w = 6000 * (rand(rng) - 0.5)
            T = 300 + 3000rand(rng); V = 10.0^(-3 + 2rand(rng))
            st = inner isa IdealGas ? (; T, V, m = 10.0^(-3 + 2rand(rng))) :
                 (; T, V, m = (H2 = 10.0^(-3 + 2rand(rng)), O2 = 10.0^(-3 + 2rand(rng)),
                               H2O = 10.0^(-3 + 2rand(rng)), N2 = 10.0^(-3 + 2rand(rng))))
            Ei = charges_at(inner; st...)
            mv = Moving(inner; m_dry)
            M = m_dry + sum(Ei[3:end])
            E = vcat(Ei, M * w)                                   # charges (inner..., p)
            @test collect(charges_at(mv; w, st...)) ≈ collect(E) rtol = 1e-12
            F, Fi = forces(mv, E), forces(inner, Ei)
            @test internal_energy(mv, E) ≈ internal_energy(inner, Ei) + (M * w)^2 / (2M) rtol = 1e-12
            @test F[end] ≈ w rtol = 1e-12                          # w = p/M
            @test F[1] ≈ Fi[1] rtol = 1e-12                        # T and −p unchanged
            @test F[2] ≈ Fi[2] rtol = 1e-12
            for k in 3:length(Fi)                                  # μ_k = μ_k,int − w²/2
                @test F[k] ≈ Fi[k] - w^2 / 2 rtol = 1e-10 atol = 1e-10 * (abs(Fi[k]) + w^2)
            end
            A = eos_matrix(mv, E)
            @test maximum(abs, A - A') <= 1e-9 * maximum(abs, A)
            @test is_stable(mv, E)
            mv0 = Moving(inner)                                    # homogeneous for m_dry = 0
            E0 = vcat(Ei, sum(Ei[3:end]) * w)
            λ = 0.1 + 10rand(rng)
            @test internal_energy(mv0, λ * E0) ≈ λ * internal_energy(mv0, E0) rtol = 1e-10
        end
    end

    @testset "momentum-carrying stream" begin
        N2 = gas_n2()
        rng = MersenneTwister(82)
        κ = 1e-8
        for _ in 1:100
            Ta, Tb = 200 + 600rand(rng), 200 + 600rand(rng)
            pa, pb = 0.5e5 + 3e5rand(rng), 0.5e5 + 3e5rand(rng)
            ma, mb = 1e-3 * (0.5 + rand(rng)), 2e-3 * (0.5 + rand(rng))
            da, db = 3rand(rng), 3rand(rng)
            wa, wb = 6000 * (rand(rng) - 0.5), 6000 * (rand(rng) - 0.5)
            function mk(w1, w2, moving)
                sys = ThermoSystem()
                if moving
                    add!(sys, Body(:a, Moving(N2; m_dry = da); T = Ta, p = pa, m = ma, w = w1))
                    add!(sys, Body(:b, Moving(N2; m_dry = db); T = Tb, p = pb, m = mb, w = w2))
                else
                    add!(sys, Body(:a, N2; T = Ta, p = pa, m = ma))
                    add!(sys, Body(:b, N2; T = Tb, p = pb, m = mb))
                end
                connect!(sys, Stream(HydraulicConductance(κ)), :a, :b; name = :st)
                cs = compile(sys)
                return cs, initial_state(cs)
            end
            cs0, u0 = mk(0.0, 0.0, false)
            S0 = entropy_production(cs0, u0)[:st]
            # (a) equal velocities: same production as the static pair
            cs1, u1 = mk(wa, wa, true)
            @test entropy_production(cs1, u1)[:st] ≈ S0 rtol = 1e-9
            r = energy_residual(cs1, u1)
            @test abs(r.residual) <= 1e-11 * r.scale
            # (b) different velocities: extra ṁ(w_u − w_d)²/(2 T_d); momentum ṁ w_u moves downstream
            cs2, u2 = mk(wa, wb, true)
            mdot = κ * (pa - pb)
            Td = mdot >= 0 ? Tb : Ta
            wu = mdot >= 0 ? wa : wb
            S2 = entropy_production(cs2, u2)[:st]
            @test S2 - S0 ≈ abs(mdot) * (wa - wb)^2 / (2Td) rtol = 1e-8 atol = 1e-12 * abs(S0)
            du = rates(cs2, u2)
            if mdot >= 0
                @test charges(cs2, du, :a)[end] ≈ -mdot * wu rtol = 1e-10
                @test charges(cs2, du, :b)[end] ≈ mdot * wu rtol = 1e-10
            else
                @test charges(cs2, du, :b)[end] ≈ mdot * wu rtol = 1e-10     # upstream is b: loses |ṁ| w_b
                @test charges(cs2, du, :a)[end] ≈ -mdot * wu rtol = 1e-10
            end
            r = energy_residual(cs2, u2)
            @test abs(r.residual) <= 1e-11 * r.scale
            @test S2 >= 0
        end
    end

    @testset "validation" begin
        N2 = gas_n2()
        function sys_with(f)
            sys = ThermoSystem()
            add!(sys, Body(:a, Moving(N2); T = 300.0, p = 2e5, m = 1e-3, w = 5.0))
            add!(sys, Body(:b, N2; T = 300.0, p = 1e5, m = 1e-3))
            f(sys)
            return sys
        end
        @test_throws ArgumentError compile(sys_with(s -> connect!(s, Stream(HydraulicConductance(1e-8)), :a, :b)))
        @test_throws ArgumentError compile(sys_with(s -> connect!(s, IdealPump(t -> 1e-5), :b, :a)))
        @test_throws ArgumentError compile(sys_with(s -> connect!(s, PassThroughOutlet(κ = 1e-7, p_out = 1e5), :a)))
        @test_throws ArgumentError compile(sys_with(s -> connect!(s, FlowBoundary(N2; T = 300.0, p = 3e5, mdot = t -> 1e-5, supply = :pumped), :a)))
        @test_throws ArgumentError compile(sys_with(s -> connect!(s, Nozzle(1e-4, 5.0), :b)))   # not moving, no stand
        @test compile(sys_with(s -> connect!(s, Nozzle(1e-4, 5.0), :a))) isa CompiledSystem
        @test_throws ArgumentError Nozzle(0.0, 5.0)
        @test_throws ArgumentError Nozzle(-1e-4, 5.0)
        @test_throws ArgumentError Nozzle(1e-4, 1.0)
        @test_throws ArgumentError Nozzle(1e-4, 0.5)
        @test_throws ArgumentError Nozzle(1e-4, 5.0; η = 0.0)
        @test_throws ArgumentError Nozzle(1e-4, 5.0; η = 1.1)
        @test Nozzle(1e-4, 5.0; η = 1.0) isa Nozzle
        # one species, two datums (IdealGas vs species data): refused
        sysd = ThermoSystem()
        add!(sysd, Body(:tank, IdealGas(:H2; R = R̄ / molar_mass(H2), cv = H2.cv); T = 300.0, p = 1e5, m = 1e-3))
        add!(sysd, Body(:chamber, IdealGasMixture(H2, O2, H2O); T = 300.0, V = 1e-3, n = (H2 = 0.01, O2 = 0.01, H2O = 0.01)))
        connect!(sysd, IdealPump(t -> 1e-5), :tank, :chamber)
        @test_throws ArgumentError compile(sysd)
    end

    @testset "pumps between moving bodies and into richer mixtures" begin
        N2 = gas_n2()
        mdot = 2e-4
        function pumped(moving, wa, wb)
            sys = ThermoSystem()
            if moving
                add!(sys, Body(:a, Moving(N2; m_dry = 1.0); T = 300.0, p = 1e5, m = 1e-2, w = wa))
                add!(sys, Body(:b, Moving(N2; m_dry = 2.0); T = 350.0, p = 3e5, m = 1e-2, w = wb))
            else
                add!(sys, Body(:a, N2; T = 300.0, p = 1e5, m = 1e-2))
                add!(sys, Body(:b, N2; T = 350.0, p = 3e5, m = 1e-2))
            end
            connect!(sys, IdealPump(t -> mdot), :a, :b; name = :pump)
            cs = compile(sys)
            return cs, initial_state(cs)
        end
        cs0, u0 = pumped(false, 0.0, 0.0)
        W0 = pump_work(cs0, rates(cs0, u0))[:pump]
        for w in (-2500.0, 0.0, 1234.5, 2999.0)
            cs, u = pumped(true, w, w)
            du = rates(cs, u)
            @test pump_work(cs, du)[:pump] ≈ W0 rtol = 1e-10
            @test charges(cs, du, :a)[end] ≈ -mdot * w atol = 1e-12 rtol = 1e-10
            @test charges(cs, du, :b)[end] ≈ mdot * w atol = 1e-12 rtol = 1e-10
            r = energy_residual(cs, u)
            @test abs(r.residual) <= 1e-11 * r.scale
        end
        # one-species tank into a richer mixture (upstream species a subset of downstream)
        tank = IdealGasMixture(H2)
        chamber = IdealGasMixture(H2, O2, H2O)
        rng = MersenneTwister(83)
        for _ in 1:20
            sys = ThermoSystem()
            add!(sys, Body(:tank, tank; T = 100 + 300rand(rng), p = 1e5 + 1e5rand(rng), m = (H2 = 0.05,)))
            Tc = 500 + 2000rand(rng); pc = 2e5 + 2e6rand(rng)
            mc = (H2 = 0.01 * (1 + rand(rng)), O2 = 0.08 * (1 + rand(rng)), H2O = 0.05 * (1 + rand(rng)))
            Vc = (mc.H2 * R̄ / molar_mass(H2) + mc.O2 * R̄ / molar_mass(O2) + mc.H2O * R̄ / molar_mass(H2O)) * Tc / pc
            add!(sys, Body(:chamber, chamber; T = Tc, V = Vc, m = mc))
            connect!(sys, IdealPump(t -> 1e-4), :tank, :chamber; name = :pump)
            cs = compile(sys)
            u = initial_state(cs)
            @test entropy_production(cs, u)[:pump] >= -1e-12
            r = energy_residual(cs, u)
            @test abs(r.residual) <= 1e-11 * r.scale
            du = rates(cs, u)
            @test charges(cs, du, :tank)[3] ≈ -1e-4 rtol = 1e-12
        end
        # the same pump with Moving bodies (rocket feed)
        sys = ThermoSystem()
        add!(sys, Body(:tank, Moving(tank; m_dry = 1.0); T = 150.0, p = 3e5, m = (H2 = 0.05,), w = 40.0))
        add!(sys, Body(:chamber, Moving(chamber; m_dry = 3.0); T = 1000.0, V = 1e-3,
                       m = (H2 = 1e-4, O2 = 8e-4, H2O = 1e-4), w = 40.0))
        connect!(sys, IdealPump(t -> 1e-4), :tank, :chamber; name = :pump)
        cs = compile(sys); u = initial_state(cs)
        @test entropy_production(cs, u)[:pump] >= -1e-12
        r = energy_residual(cs, u)
        @test abs(r.residual) <= 1e-11 * r.scale
    end

    @testset "nozzle kernel identities (random states)" begin
        rng = MersenneTwister(84)
        for _ in 1:250
            c = random_chamber(rng)
            cs = chamber_system(c)
            u = initial_state(cs)
            ns = nozzle_state(cs, u, :nozzle)
            mtot = sum(c.m)
            y = collect(c.m) ./ mtot
            fp = frozen_props(c.sds, y)
            k, R, cp = fp.k, fp.R, fp.cp
            Tc, pc, At, ε, η = c.T, c.p, c.At, c.ε, c.η
            @test ns.k ≈ k rtol = 1e-10
            @test ns.R ≈ R rtol = 1e-10
            @test ns.cp ≈ cp rtol = 1e-10
            # area–Mach relation, supersonic root
            @test ns.M_e > 1
            @test areamach(ns.M_e, k) ≈ ε rtol = 1e-10
            # choked mass flow, recomputed from the species data
            mdot = At * pc * sqrt(k / (R * Tc)) * (2 / (k + 1))^((k + 1) / (2 * (k - 1)))
            @test ns.mdot ≈ mdot rtol = 1e-12
            @test mass_flow(cs, u, :nozzle) ≈ mdot rtol = 1e-12
            # exit state and velocity
            Tes = Tc / (1 + (k - 1) / 2 * ns.M_e^2)
            @test ns.T_es ≈ Tes rtol = 1e-10
            @test ns.p_e ≈ pc * (Tes / Tc)^(k / (k - 1)) rtol = 1e-10
            @test ns.v_e^2 ≈ 2 * cp * (Tc - ns.T_e) rtol = 1e-12
            @test ns.v_e^2 ≈ 2 * η * cp * (Tc - ns.T_es) rtol = 1e-12
            @test ns.A_e ≈ ε * At rtol = 1e-14
            @test ns.thrust ≈ ns.mdot * ns.v_e + ns.p_e * ns.A_e rtol = 1e-12
            @test ns.J_p ≈ ns.mdot * (c.w - ns.v_e) - ns.p_e * ns.A_e rtol = 1e-12
            # energy: the chamber's Σ P dE/dt = −E_jet = −[ṁ(h_e + (w − v_e)²/2) − w p_e A_e],
            # with h_e = h_c − v_e²/2 and h_c from the inner material at the inner charges
            du = rates(cs, u)
            mat = Moving(c.mix; m_dry = c.m_dry)
            E = charges(cs, u, :chamber)
            dE = charges(cs, du, :chamber)
            Edot = sum(forces(mat, E) .* dE)
            Ei = charges_at(c.mix; T = Tc, V = c.V, m = (H2 = c.m[1], O2 = c.m[2], H2O = c.m[3], N2 = c.m[4]))
            hc = (internal_energy(c.mix, Ei) + c.p * c.V) / mtot
            Ejet = ns.mdot * (hc - ns.v_e^2 / 2 + (c.w - ns.v_e)^2 / 2) - c.w * ns.p_e * ns.A_e
            escale = ns.mdot * (abs(hc) + c.w^2 + ns.v_e^2)
            @test Edot ≈ -ns.E_jet rtol = 1e-10 atol = 1e-12 * escale
            @test ns.E_jet ≈ Ejet rtol = 1e-10 atol = 1e-10 * escale
            # generation
            S = entropy_production(cs, u)[:nozzle]
            @test S >= -1e-12 * ns.mdot * cp
            if η == 1.0
                @test abs(S) <= 1e-12 * ns.mdot * cp
            else
                @test S ≈ ns.mdot * (cp * log(ns.T_e / Tc) - R * log(ns.p_e / pc)) rtol = 1e-9
                @test S > 0
            end
            @test ns.S_gen ≈ S rtol = 1e-10 atol = 1e-14 * ns.mdot * cp
            # momentum: vehicle + jet = 0, then with gravity on the chamber too
            pscale = ns.mdot * (abs(c.w) + ns.v_e) + ns.p_e * ns.A_e
            @test charges(cs, du, :chamber)[end] + boundary_ledger(cs, du, :nozzle).p ≈ 0 atol = 1e-11 * pscale
            @test boundary_ledger(cs, du, :nozzle).p ≈ ns.J_p rtol = 1e-10 atol = 1e-12 * pscale
            @test boundary_ledger(cs, du, :nozzle).m ≈ ns.mdot rtol = 1e-12
            csg = chamber_system(c; gravity = true)
            ug = initial_state(csg)
            dug = rates(csg, ug)
            Mtot = c.m_dry + mtot
            @test source_ledger(csg, dug, :gravity).charge ≈ g0 * Mtot rtol = 1e-12
            resid = charges(csg, dug, :chamber)[end] + boundary_ledger(csg, dug, :nozzle).p +
                    source_ledger(csg, dug, :gravity).charge
            @test abs(resid) <= 1e-11 * (pscale + g0 * Mtot)
            r = energy_residual(csg, ug)
            @test abs(r.residual) <= 1e-11 * r.scale
        end
        # single-species IdealGas chamber
        N2 = gas_n2()
        for _ in 1:50
            Tc = 300 + 3200rand(rng); pc = 10.0^(5 + 2rand(rng)); w = 6000 * (rand(rng) - 0.5)
            ε = 1.5 + 48.5rand(rng); η = rand(rng) < 0.5 ? 1.0 : 0.9; At = 1e-3
            sys = ThermoSystem()
            add!(sys, Body(:chamber, Moving(N2; m_dry = 2.0); T = Tc, p = pc, m = 0.5, w))
            connect!(sys, Nozzle(At, ε; η), :chamber; name = :nozzle)
            cs = compile(sys); u = initial_state(cs)
            ns = nozzle_state(cs, u, :nozzle)
            k = (N2.cv + N2.R) / N2.cv
            @test areamach(ns.M_e, k) ≈ ε rtol = 1e-10
            @test ns.M_e > 1
            @test ns.mdot ≈ At * pc * sqrt(k / (N2.R * Tc)) * (2 / (k + 1))^((k + 1) / (2 * (k - 1))) rtol = 1e-12
            @test ns.v_e^2 ≈ 2 * (N2.cv + N2.R) * (Tc - ns.T_e) rtol = 1e-12
            du = rates(cs, u)
            Edot = sum(forces(Moving(N2; m_dry = 2.0), charges(cs, u, :chamber)) .* charges(cs, du, :chamber))
            @test Edot ≈ -ns.E_jet rtol = 1e-10
            S = entropy_production(cs, u)[:nozzle]
            if η == 1.0
                @test abs(S) <= 1e-12 * ns.mdot * ns.cp
            else
                @test S > 0
            end
        end
    end

    @testset "Galilean invariance" begin
        rng = MersenneTwister(85)
        for _ in 1:30
            c = random_chamber(rng)
            cs0 = chamber_system((; c..., w = 0.0)); cs1 = chamber_system((; c..., w = 3000.0))
            u0, u1 = initial_state(cs0), initial_state(cs1)
            a, b = nozzle_state(cs0, u0, :nozzle), nozzle_state(cs1, u1, :nozzle)
            for f in (:mdot, :T_e, :p_e, :v_e, :thrust)
                @test getfield(a, f) ≈ getfield(b, f) rtol = 1e-12
            end
            @test entropy_production(cs0, u0)[:nozzle] ≈ entropy_production(cs1, u1)[:nozzle] rtol = 1e-12 atol = 1e-15
            @test temperature(cs0, u0, :chamber) ≈ temperature(cs1, u1, :chamber) rtol = 1e-12
            @test pressure(cs0, u0, :chamber) ≈ pressure(cs1, u1, :chamber) rtol = 1e-12
            @test temperature(cs1, u1, :chamber) ≈ c.T rtol = 1e-10
            @test pressure(cs1, u1, :chamber) ≈ c.p rtol = 1e-10
        end
    end

    @testset "test stand: thrust through a momentum reservoir" begin
        rng = MersenneTwister(86)
        for _ in 1:50
            c = random_chamber(rng)
            sys = ThermoSystem()
            add!(sys, Body(:chamber, c.mix; T = c.T, V = c.V,
                           m = (H2 = c.m[1], O2 = c.m[2], H2O = c.m[3], N2 = c.m[4])))
            add!(sys, Reservoir(:stand, Momentum() => 0.0))
            connect!(sys, Nozzle(c.At, c.ε; η = c.η), :chamber, :stand; name = :nozzle)
            cs = compile(sys); u = initial_state(cs)
            ns = nozzle_state(cs, u, :nozzle)
            du = rates(cs, u)
            F = ns.mdot * ns.v_e + ns.p_e * ns.A_e                    # vacuum thrust
            @test ns.thrust ≈ F rtol = 1e-12
            @test reservoir_charge(cs, du, :stand, Momentum()) ≈ F rtol = 1e-10
            @test abs(reservoir_ledger(cs, du, :stand).U) <= 1e-12 * ns.mdot * ns.cp * c.T   # w = 0: no work
            r = energy_residual(cs, u)
            @test abs(r.residual) <= 1e-11 * r.scale
            @test entropy_production(cs, u)[:nozzle] >= -1e-12 * ns.mdot * ns.cp
        end
    end

    @testset "Çengel 9e Example 17-7" begin
        # Çengel 9e Example 17-7: air (k = 1.4, R = 287), p0 = 1 MPa, T0 = 800 K, A* = 20 cm²,
        # exit Mach 2 (A_e/A* = 1.6875). Table values: p_e = 0.1278 MPa, T_e = 444.5 K,
        # v_e = 845.1 m/s, A_e = 33.75 cm², ṁ = 2.86 kg/s.
        air = IdealGas(:air; R = 287.0, cv = 287.0 / 0.4)
        sys = ThermoSystem()
        add!(sys, Body(:chamber, air; T = 800.0, p = 1.0e6, m = 1.0))
        add!(sys, Reservoir(:stand, Momentum() => 0.0))
        connect!(sys, Nozzle(20e-4, 1.6875), :chamber, :stand; name = :nozzle)
        cs = compile(sys); u = initial_state(cs)
        ns = nozzle_state(cs, u, :nozzle)
        @test ns.M_e ≈ 2.0 rtol = 2e-3
        @test ns.p_e ≈ 0.1278e6 rtol = 2e-3
        @test ns.T_e ≈ 444.5 rtol = 2e-3
        @test ns.v_e ≈ 845.1 rtol = 2e-3
        @test ns.A_e ≈ 33.75e-4 rtol = 2e-3
        @test ns.mdot ≈ 2.86 rtol = 2e-3
    end

    @testset "gravity on a moving body" begin
        # the field acts on the CURRENT mass; ledger: charge +gM, W = −gMw, x' = w
        N2 = gas_n2()
        m_dry = 2.0
        for (mf, w) in ((0.5, 10.0), (0.1, -7.0))
            sys = ThermoSystem()
            add!(sys, Body(:tank, Moving(N2; m_dry); T = 300.0, p = 2e5, m = mf, w))
            connect!(sys, GravityField(g0), :tank; name = :gravity)
            cs = compile(sys); u = initial_state(cs)
            du = rates(cs, u)
            sl = source_ledger(cs, du, :gravity)
            @test sl.charge ≈ g0 * (m_dry + mf) rtol = 1e-12
            @test sl.W ≈ -g0 * (m_dry + mf) * w rtol = 1e-12
            @test sl.x ≈ w rtol = 1e-12
            @test charges(cs, du, :tank)[end] ≈ -g0 * (m_dry + mf) rtol = 1e-12
            r = energy_residual(cs, u)
            @test abs(r.residual) <= 1e-11 * r.scale
        end
        # an InertialMass keeps working, with the x ledger
        sys = ThermoSystem()
        add!(sys, Body(:mass, InertialMass(3.0); w = 4.0))
        connect!(sys, GravityField(g0), :mass; name = :gravity)
        cs = compile(sys); u = initial_state(cs)
        sl = source_ledger(cs, rates(cs, u), :gravity)
        @test sl.charge ≈ 3g0 rtol = 1e-14
        @test sl.W ≈ -3g0 * 4.0 rtol = 1e-14
        @test sl.x ≈ 4.0 rtol = 1e-14
        # after the body loses mass through a nozzle, the force follows the current mass
        sys = ThermoSystem()
        add!(sys, Body(:tank, Moving(N2; m_dry); T = 300.0, p = 5e5, m = 0.5, w = 0.0))
        connect!(sys, Nozzle(1e-5, 4.0), :tank; name = :nozzle)
        connect!(sys, GravityField(g0), :tank; name = :gravity)
        cs = compile(sys); u0 = initial_state(cs)
        sol = simulate(cs, (0.0, 0.05))                                # one short run
        u1 = sol.u[end]
        m1 = charges(cs, u1, :tank)[3]
        @test m1 < 0.5
        for (u, mf) in ((u0, 0.5), (u1, m1))
            sl = source_ledger(cs, rates(cs, u), :gravity)
            @test sl.charge ≈ g0 * (m_dry + mf) rtol = 1e-10
        end
        r = energy_residual(cs, u1)
        @test abs(r.residual) <= 1e-8 * r.scale
    end
end
