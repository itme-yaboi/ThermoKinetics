# Slice 10b: two-phase pure fluids through the Clapeyron backend (DESIGN §2.6, D-52).
# Dome states are built from Clapeyron's own saturation properties (an independent route) and compared
# with the library's flash; draws are checked against the saturated phases' specific entropies.
# All Clapeyron functions are qualified (names such as `internal_energy` clash with ThermoKinetics).

@testset "two-phase pure fluids (Clapeyron backend)" begin
    Rbar = 8.314462618                                   # J/(mol K), only used for tolerances

    function make_fluid(name::Symbol)
        if name == :co2
            mdl = Clapeyron.SingleFluid("carbon dioxide"); return (RealFluid(mdl; species = :CO2), mdl)
        else
            mdl = Clapeyron.IAPWS95(); return (RealFluid(mdl; species = :H2O), mdl)
        end
    end

    molar_mass(mdl) = collect(Clapeyron.mw(mdl))[1] / 1000          # kg/mol

    "Saturated properties at T (molar): everything from Clapeyron's saturation and (V, T) functions."
    function sat_props(mdl, T)
        ps, vl, vv = Clapeyron.saturation_pressure(mdl, T)
        return (; ps, vl, vv,
                sl = Clapeyron.VT_entropy(mdl, vl, T), sv = Clapeyron.VT_entropy(mdl, vv, T),
                ul = Clapeyron.VT_internal_energy(mdl, vl, T), uv = Clapeyron.VT_internal_energy(mdl, vv, T),
                gl = Clapeyron.VT_gibbs_free_energy(mdl, vl, T))
    end

    "Charges (S, V, m) of a wet mixture of quality x at T with mass m."
    function dome_charges(mdl, T, x, m)
        M = molar_mass(mdl); n = m / M; sp = sat_props(mdl, T)
        return [n * ((1 - x) * sp.sl + x * sp.sv), n * ((1 - x) * sp.vl + x * sp.vv), m]
    end

    "A random (T, x, m) in the dome of the named fluid."
    function random_dome(rng, name)
        T = name == :co2 ? 220 + 83rand(rng) : 280 + 360rand(rng)
        return T, 0.02 + 0.96rand(rng), 10.0^(-3 + 3rand(rng))
    end

    rng = MersenneTwister(1011)

    @testset "dome states vs Clapeyron's saturation" begin
        for name in (:water, :co2)
            f, mdl = make_fluid(name)
            M = molar_mass(mdl)
            for _ in 1:40
                T, x, m = random_dome(rng, name)
                n = m / M; sp = sat_props(mdl, T)
                E = dome_charges(mdl, T, x, m)
                P = forces(f, E)
                @test valid_state(f, E)
                @test P[1] ≈ T rtol = 1e-8
                @test -P[2] ≈ sp.ps rtol = 1e-8
                @test P[3] * M ≈ sp.gl rtol = 1e-8 atol = 1e-6 * Rbar * T
                q = quality(f, E)
                @test q ≈ x rtol = 1e-7 atol = 1e-9
                U = ThermoKinetics.internal_energy(f, E)
                @test U ≈ n * ((1 - x) * sp.ul + x * sp.uv) rtol = 1e-9 atol = 1e-9 * n * Rbar * T
            end
        end
    end

    @testset "continuity at the dome boundary" begin
        for (name, T) in ((:water, 450.0), (:water, 600.0), (:co2, 260.0), (:co2, 295.0))
            f, mdl = make_fluid(name)
            M = molar_mass(mdl); m = 0.7; n = m / M
            sp = sat_props(mdl, T)
            for (x, v) in ((1e-7, sp.vl), (1 - 1e-7, sp.vv))
                Ew = dome_charges(mdl, T, x, m)
                Es = charges_at(f; T, V = n * v, m)             # the single-phase saturated state
                Pw = forces(f, Ew); Ps = forces(f, Es)
                @test Pw[1] ≈ Ps[1] rtol = 1e-5
                @test Pw[2] ≈ Ps[2] rtol = 1e-5
                @test Pw[3] ≈ Ps[3] rtol = 1e-5 atol = 1e-5 * Rbar * T / M
            end
            # just outside the dome at the same T: compressed liquid below v_l, superheated vapour above v_v
            # (not at the saturated entropy: expanding saturated vapour isentropically enters the dome
            # for a wet fluid such as water)
            El = charges_at(f; T, V = n * sp.vl * (1 - 1e-4), m)
            Ev = charges_at(f; T, V = n * sp.vv * (1 + 1e-4), m)
            @test valid_state(f, El)
            @test valid_state(f, Ev)
            @test quality(f, El) === nothing
            @test quality(f, Ev) === nothing
            # a single-phase state built from (T, p) away from the dome has no quality
            @test quality(f, charges_at(f; T = T + 30, p = sp.ps / 3, m)) === nothing
        end
    end

    @testset "exact derivatives inside the dome" begin
        for name in (:water, :co2)
            f, mdl = make_fluid(name)
            for _ in 1:4
                T, x, m = random_dome(rng, name)
                x = 0.2 + 0.6 * x                                # keep away from the boundary
                E = dome_charges(mdl, T, x, m)
                J = ForwardDiff.jacobian(e -> forces(f, e), E)
                Jfd = zeros(3, 3)
                for k in 1:3
                    h = 1e-6 * abs(E[k])
                    Ep = copy(E); Em = copy(E)
                    Ep[k] += h; Em[k] -= h
                    Jfd[:, k] = (Vector(forces(f, Ep)) - Vector(forces(f, Em))) ./ (2h)
                end
                for i in 1:3, j in 1:3
                    atol = 1e-5 * sum(abs(J[i, k] * E[k]) for k in 1:3) / abs(E[j])
                    @test J[i, j] ≈ Jfd[i, j] rtol = 1e-4 atol = atol
                end
                A = Matrix(eos_matrix(f, E))
                @test maximum(abs, A - A') ≤ 1e-6 * maximum(abs, A)
                λ = abs.(eigvals(Symmetric((A + A') / 2)))
                @test minimum(λ) ≤ 1e-6 * maximum(λ)             # the tie line
            end
        end
    end

    @testset "isentropic state into the dome" begin
        # Çengel 9e Example 10-1: steam at 3 MPa, 350 °C expanded to 75 kPa
        f, mdl = make_fluid(:water)
        M = molar_mass(mdl)
        s1 = charges_at(f; T = 623.15, p = 3e6, m = 1.0)[1]
        p2 = 75e3
        lo, hi = 330.0, 420.0                                     # bisection for T_sat(p2)
        for _ in 1:80
            mid = (lo + hi) / 2
            Clapeyron.saturation_pressure(mdl, mid)[1] < p2 ? (lo = mid) : (hi = mid)
        end
        Tsat = (lo + hi) / 2
        st = ThermoKinetics.isentropic_state(f, s1, p2)
        @test st.T ≈ Tsat rtol = 1e-7
        sp = sat_props(mdl, Tsat)
        sl, sv = sp.sl / M, sp.sv / M
        @test sl < s1 < sv
        x = (s1 - sl) / (sv - sl)
        @test x ≈ 0.8861 rtol = 2e-3                              # Çengel 9e Example 10-1
        @test st.μ ≈ sp.gl / M rtol = 1e-7 atol = 1e-6 * Rbar * Tsat / M
    end

    @testset "phase-selective draws" begin
        f, mdl = make_fluid(:water)
        M = molar_mass(mdl)
        T0, x0, m0 = 360.0, 0.5, 1.0
        sp = sat_props(mdl, T0)
        sl, sv = sp.sl / M, sp.sv / M
        Eup = dome_charges(mdl, T0, x0, m0)
        Edn = charges_at(f; T = 400.0, p = 1e4, m = 1.0)         # superheated vapour, lower pressure
        mdot = 1e-3

        function build(el)
            sys = ThermoSystem()
            add!(sys, Body(:up, f, Eup))
            add!(sys, Body(:dn, f, Edn))
            connect!(sys, el, :up, :dn; name = :link)
            cs = compile(sys)
            return cs, initial_state(cs)
        end

        function check(el, expect_s; mdot_fixed = true)
            cs, u = build(el)
            du = rates(cs, u)
            dU = charges(cs, du, :up)
            md = -dU[3]
            @test md > 0
            mdot_fixed && @test md ≈ mdot rtol = 1e-12
            @test dU[1] ≈ -md * expect_s rtol = 1e-8
            r = energy_residual(cs, u)
            @test abs(r.residual) <= 1e-12 * r.scale
            sg = entropy_production(cs, u)
            @test sg[:link] >= -1e-12 * (abs(dU[1]) + 1e-300)
            Enew = Eup .+ 1.0 .* Vector(dU)                       # one second of drawing
            return quality(f, Enew)
        end

        qliq = check(IdealPump(t -> mdot; draw = :liquid), sl)
        qvap = check(IdealPump(t -> mdot; draw = :vapor), sv)
        qmix = check(IdealPump(t -> mdot; draw = :mixture), Eup[1] / m0)
        @test qliq > x0                                           # drawing liquid raises the quality
        @test qvap < x0                                           # drawing vapour lowers it
        @test qmix !== nothing
        @test check(IdealPump(t -> mdot), Eup[1] / m0) ≈ qmix rtol = 1e-12   # :mixture is the default

        # a stream out of the same body behaves the same way
        qs = check(Stream(HydraulicConductance(1e-9); draw = :liquid), sl; mdot_fixed = false)
        @test qs > x0
        qs = check(Stream(HydraulicConductance(1e-9); draw = :vapor), sv; mdot_fixed = false)
        @test qs < x0

        @test_throws ArgumentError IdealPump(t -> 1.0; draw = :steam)
        @test_throws ArgumentError Stream(HydraulicConductance(1e-9); draw = :steam)

        # from a single-phase (superheated) body every draw takes the body's own state
        Esh = charges_at(f; T = 500.0, p = 1e5, m = 1.0)
        function rates_sh(draw)
            sys = ThermoSystem()
            add!(sys, Body(:up, f, Esh))
            add!(sys, Body(:dn, f, Edn))
            connect!(sys, IdealPump(t -> mdot; draw), :up, :dn; name = :link)
            cs = compile(sys)
            return Vector(rates(cs, initial_state(cs)))
        end
        dm = rates_sh(:mixture)
        for draw in (:liquid, :vapor)
            d = rates_sh(draw)
            @test maximum(abs, d .- dm) <= 1e-12 * maximum(abs, dm)
        end
    end
end
