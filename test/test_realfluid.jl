# Slice 10a: the single-phase real-fluid backend (RealFluid, through Clapeyron.jl; DESIGN §2.4, §2.5, D-50).
# Every check compares the material (V-based Helmholtz function + temperature inversion) with
# Clapeyron's own pressure-based property functions, an independent code path.
# All Clapeyron functions are qualified: ThermoKinetics and Clapeyron both export names such as
# `internal_energy`, `pressure` and `temperature`.

@testset "real fluids (Clapeyron backend)" begin
    Rbar = 8.314462618                                   # J/(mol K), only used for tolerances

    "Returns (material, Clapeyron model, species names) for the named fluid."
    function make_fluid(name::Symbol)
        if name == :helium
            mdl = Clapeyron.GERG2008(["helium"]); return (RealFluid(mdl; species = :He), mdl, (:He,))
        elseif name == :co2
            mdl = Clapeyron.SingleFluid("carbon dioxide"); return (RealFluid(mdl; species = :CO2), mdl, (:CO2,))
        elseif name == :water
            mdl = Clapeyron.IAPWS95(); return (RealFluid(mdl; species = :H2O), mdl, (:H2O,))
        else
            mdl = Clapeyron.GERG2008(["nitrogen", "oxygen", "argon"])
            return (RealFluid(mdl; species = (:N2, :O2, :Ar)), mdl, (:N2, :O2, :Ar))
        end
    end

    "Molar masses [kg/mol] as a vector."
    molar_masses(mdl) = collect(Clapeyron.mw(mdl)) ./ 1000

    "Mass argument for `charges_at`: a Real for one species, a NamedTuple for several."
    mass_arg(sp, mvec) = length(sp) == 1 ? mvec[1] : NamedTuple{sp}(Tuple(mvec))

    "Clapeyron's view of a state: per-mole composition zn, total moles, mean molar mass."
    function cl_view(mdl, mvec)
        n = mvec ./ molar_masses(mdl)
        ntot = sum(n); zn = n ./ ntot
        return (; n, ntot, zn, Mbar = sum(zn .* molar_masses(mdl)))
    end

    "A random single-phase (T, p, masses) for the fluid."
    function random_state(rng, name, mdl)
        mtot = 10.0^(-3 + 3rand(rng))
        if name == :helium
            return 20 + 1480rand(rng), 10.0^(4 + 3.3rand(rng)), [mtot]
        elseif name == :co2
            Tc = Clapeyron.crit_pure(mdl)[1]
            T = 220 + 880rand(rng)
            p = T < Tc ? Clapeyron.saturation_pressure(mdl, T)[1] * 0.9rand(rng) + 1e3 :
                         10.0^(4 + 3.4rand(rng))
            return T, p, [mtot]
        elseif name == :water
            if rand(rng) < 0.5                                       # compressed liquid
                T = 280 + 320rand(rng)
                return T, Clapeyron.saturation_pressure(mdl, T)[1] + 1e6 + 2e7rand(rng), [mtot]
            else                                                     # superheated vapor
                T = 400 + 700rand(rng)
                p = T < 640 ? 0.9rand(rng) * Clapeyron.saturation_pressure(mdl, T)[1] + 1e3 :
                              10.0^(4 + 3.3rand(rng))
                return T, p, [mtot]
            end
        else
            w = [0.755, 0.232, 0.013] .* (1 .+ 0.1 .* (rand(rng, 3) .- 0.5))
            return 150 + 1350rand(rng), 10.0^(4 + 2.7rand(rng)), mtot .* w ./ sum(w)
        end
    end

    names = (:helium, :co2, :water, :air)
    rng = MersenneTwister(1010)

    @testset "forces and state properties vs Clapeyron's (p, T) route" begin
        for name in names
            f, mdl, sp = make_fluid(name)
            M = molar_masses(mdl)
            @test domains(f) == (Thermal(), Volumetric(), (Species(s) for s in sp)...)
            for _ in 1:60
                T, p, mvec = random_state(rng, name, mdl)
                c = cl_view(mdl, mvec)
                m = sum(mvec)
                E = charges_at(f; T, p, m = mass_arg(sp, mvec))
                @test length(E) == 2 + length(sp)
                P = forces(f, E)
                @test valid_state(f, E)
                @test P[1] ≈ T rtol = 1e-9
                @test -P[2] ≈ p rtol = 1e-8
                μcl = Clapeyron.chemical_potential(mdl, p, T, c.zn)
                for k in eachindex(sp)
                    @test P[2 + k] * M[k] ≈ μcl[k] rtol = 1e-8 atol = 1e-6 * Rbar * T
                end
                vcl = Clapeyron.volume(mdl, p, T, c.zn) / c.Mbar            # m³/kg
                @test E[2] / m ≈ vcl rtol = 1e-9
                cvkg = Clapeyron.isochoric_heat_capacity(mdl, p, T, c.zn) / c.Mbar
                scl = Clapeyron.entropy(mdl, p, T, c.zn) / c.Mbar            # J/(kg K)
                @test E[1] / m ≈ scl rtol = 1e-9 atol = 1e-9 * cvkg
                ucl = m * Clapeyron.internal_energy(mdl, p, T, c.zn) / c.Mbar
                @test ThermoKinetics.internal_energy(f, E) ≈ ucl rtol = 1e-9 atol = 1e-9 * m * cvkg * T
            end
        end
    end

    @testset "homogeneity" begin
        for name in names
            f, mdl, sp = make_fluid(name)
            for _ in 1:10
                T, p, mvec = random_state(rng, name, mdl)
                E = charges_at(f; T, p, m = mass_arg(sp, mvec))
                U = ThermoKinetics.internal_energy(f, E)
                for λ in (0.37, 2.5, 11.0)
                    @test ThermoKinetics.internal_energy(f, λ .* E) ≈ λ * U rtol = 1e-10
                end
            end
        end
    end

    @testset "Hessian: symmetry, stability, T/C_v" begin
        for name in names
            f, mdl, sp = make_fluid(name)
            for _ in 1:15
                T, p, mvec = random_state(rng, name, mdl)
                c = cl_view(mdl, mvec)
                E = charges_at(f; T, p, m = mass_arg(sp, mvec))
                A = Matrix(eos_matrix(f, E))
                @test size(A) == (length(E), length(E))
                @test maximum(abs, A - A') ≤ 1e-7 * maximum(abs, A)            # before any symmetrizing
                @test is_stable(f, E)
                Cv = c.ntot * Clapeyron.isochoric_heat_capacity(mdl, p, T, c.zn)   # J/K
                @test A[1, 1] ≈ T / Cv rtol = 1e-7
            end
        end
    end

    @testset "exact derivatives (implicit-function Jacobian)" begin
        for name in names
            f, mdl, sp = make_fluid(name)
            for _ in 1:8
                T, p, mvec = random_state(rng, name, mdl)
                E = Vector(charges_at(f; T, p, m = mass_arg(sp, mvec)))
                J = ForwardDiff.jacobian(e -> forces(f, e), E)
                Jfd = zeros(length(E), length(E))
                for k in eachindex(E)
                    h = 1e-6 * max(abs(E[k]), 1e-12)
                    Ep = copy(E); Em = copy(E)
                    Ep[k] += h; Em[k] -= h
                    Jfd[:, k] = (Vector(forces(f, Ep)) - Vector(forces(f, Em))) ./ (2h)
                end
                for i in axes(J, 1), j in axes(J, 2)
                    # scale of the (i, j) entry: the size of the terms that make up row i
                    atol = 1e-5 * sum(abs(J[i, k] * E[k]) for k in eachindex(E)) / abs(E[j])
                    @test J[i, j] ≈ Jfd[i, j] rtol = 1e-5 atol = atol
                end
            end
        end
    end

    @testset "robustness: near-critical CO₂ and cold starts" begin
        f, mdl, sp = make_fluid(:co2)
        M = molar_masses(mdl)
        for _ in 1:20
            T = 304.5 + 5.5rand(rng); p = 7.4e6 + 1.1e6rand(rng)
            c = cl_view(mdl, [0.7])
            E = charges_at(f; T, p, m = 0.7)
            P = forces(f, E)
            @test P[1] ≈ T rtol = 1e-7
            @test -P[2] ≈ p rtol = 1e-7
            @test P[3] * M[1] ≈ Clapeyron.chemical_potential(mdl, p, T, c.zn)[1] rtol = 1e-7 atol = 1e-6 * Rbar * T
            @test E[2] / 0.7 ≈ Clapeyron.volume(mdl, p, T, c.zn) / c.Mbar rtol = 1e-7
        end
        for name in (:helium, :water, :air)
            f, mdl, sp = make_fluid(name)
            mv = name == :air ? [0.75, 0.23, 0.02] : [1.0]
            pp = name == :water ? 5e6 : 1e6
            Ecold = charges_at(f; T = 300.0, p = pp, m = mass_arg(sp, mv))
            Ehot = charges_at(f; T = 1400.0, p = pp, m = mass_arg(sp, mv))
            @test forces(f, Ecold)[1] ≈ 300.0 rtol = 1e-9
            @test forces(f, Ehot)[1] ≈ 1400.0 rtol = 1e-9        # far from the cached 300 K
            @test forces(f, Ecold)[1] ≈ 300.0 rtol = 1e-9         # and back
            @test forces(f, Ehot)[1] ≈ 1400.0 rtol = 1e-9
        end
    end

    @testset "guards" begin
        f, mdl, sp = make_fluid(:water)
        ps, vl, vv = Clapeyron.saturation_pressure(mdl, 400.0)
        Mw = molar_masses(mdl)[1]
        vmid = sqrt(vl * vv) / Mw                                  # m³/kg, between saturated liquid and vapor
        Eliq = charges_at(f; T = 400.0, p = ps * 1.5, m = 0.5)     # compressed liquid, for a plausible S
        Edome = [Eliq[1], vmid * 0.5, 0.5]                         # 400 K is below T_c; v inside the dome
        # since slice 10b a pure fluid in the dome is a valid saturated mixture (test_twophase.jl)
        @test valid_state(f, Edome)
        @test 0 < quality(f, Edome) < 1
        @test valid_state(f, Eliq)
        @test valid_state(f, charges_at(f; T = 400.0, p = ps * 0.5, m = 0.5))        # superheated vapor
        Ev = Vector(Eliq)
        Ev[2] = 0.0;   @test !valid_state(f, Ev)
        Ev[2] = -1e-3; @test !valid_state(f, Ev)
        Em = Vector(Eliq); Em[3] = -1e-3
        @test !valid_state(f, Em)
    end

    @testset "isentropic state at a given pressure" begin
        # (T₁, p₁, p₂) per fluid: compression and expansion, all single-phase
        cases = Dict(
            :helium => ((300.0, 1e6, 4e6), (700.0, 5e6, 1e6), (350.0, 1e5, 3e6), (1000.0, 1e6, 2e5)),
            :co2    => ((400.0, 1e6, 6e6), (700.0, 2e7, 3e6), (350.0, 1e5, 2e6), (1000.0, 1e6, 2e5)),
            # (700 K, 3 MPa) → 0.1 MPa would end inside the dome (s₁ < s_g at 0.1 MPa): stop at 1 MPa
            :water  => ((700.0, 3e6, 1e6), (600.0, 1e5, 1e6), (800.0, 5e6, 5e5)),
            :air    => ((300.0, 1e5, 2e6), (900.0, 2e6, 1e5), (420.0, 2e5, 1.2e6)))
        for name in names
            f, mdl, sp = make_fluid(name)
            M = molar_masses(mdl)
            mvec = name == :air ? [0.75, 0.23, 0.02] : [1.0]
            c = cl_view(mdl, mvec)
            for (T1, p1, p2) in cases[name]
                E = charges_at(f; T = T1, p = p1, m = mass_arg(sp, mvec))
                s1 = E[1] / sum(mvec)
                # a mixture's isentropic state needs its composition (mass fractions)
                st = name == :air ? ThermoKinetics.isentropic_state(f, s1, p2, Tuple(mvec ./ sum(mvec))) :
                                    ThermoKinetics.isentropic_state(f, s1, p2)
                cvkg = Clapeyron.isochoric_heat_capacity(mdl, p2, st.T, c.zn) / c.Mbar
                s2 = Clapeyron.entropy(mdl, p2, st.T, c.zn) / c.Mbar
                @test s2 ≈ s1 rtol = 1e-9 atol = 1e-9 * cvkg
                μcl = Clapeyron.chemical_potential(mdl, p2, st.T, c.zn)
                μ = st.μ isa Number ? (st.μ,) : st.μ
                for k in eachindex(sp)
                    @test μ[k] * M[k] ≈ μcl[k] rtol = 1e-8 atol = 1e-6 * Rbar * st.T
                end
                @test (p2 > p1) == (st.T > T1)            # compression heats, expansion cools
            end
        end
    end

    @testset "datum rule" begin
        # the same species held by two materials with different datums is refused
        water = RealFluid(Clapeyron.IAPWS95(); species = :H2O)
        steam = IdealGas(:H2O; R = 461.5, cv = 1410.0)
        sys = ThermoSystem()
        a = add!(sys, Body(:a, water; T = 600.0, p = 2e6, m = 0.1))
        b = add!(sys, Body(:b, steam; T = 600.0, p = 1e6, m = 0.1))
        connect!(sys, Stream(HydraulicConductance(1e-8)), a, b)
        @test_throws ArgumentError compile(sys)

        # two RealFluids with different backends for the same species (IAPWS-95 vs a Peng–Robinson cubic;
        # Clapeyron's IAPWS95() and SingleFluid("water") are the same equation, so those two agree)
        w2 = RealFluid(Clapeyron.PR(["water"]); species = :H2O)
        sys2 = ThermoSystem()
        a2 = add!(sys2, Body(:a, water; T = 600.0, p = 2e6, m = 0.1))
        b2 = add!(sys2, Body(:b, w2; T = 600.0, p = 1e6, m = 0.1))
        connect!(sys2, Stream(HydraulicConductance(1e-8)), a2, b2)
        @test_throws ArgumentError compile(sys2)

        # the same material on both ends is accepted
        sys3 = ThermoSystem()
        a3 = add!(sys3, Body(:a, water; T = 600.0, p = 2e6, m = 0.1))
        b3 = add!(sys3, Body(:b, water; T = 600.0, p = 1e6, m = 0.1))
        connect!(sys3, Stream(HydraulicConductance(1e-8)), a3, b3)
        @test compile(sys3) isa CompiledSystem
    end

    @testset "reference values" begin
        # Çengel 9e Example 7-5 / Table A-6 (superheated steam), 5 MPa, 450 °C:
        # h = 3317.2 kJ/kg, s = 6.8210 kJ/(kg K); the IAPWS-95 datum is the steam-table datum
        f, mdl, sp = make_fluid(:water)
        E = charges_at(f; T = 450 + 273.15, p = 5e6, m = 2.0)
        h = (ThermoKinetics.internal_energy(f, E) - forces(f, E)[2] * E[2]) / 2.0     # (U + pV)/m; the force is −p
        @test h / 1e3 ≈ 3317.2 rtol = 1e-3
        @test E[1] / 2.0 / 1e3 ≈ 6.8210 rtol = 1e-3

        # helium at 300 K, 1 MPa: density 1.5970 kg/m³ (NIST)
        fh, _, _ = make_fluid(:helium)
        Eh = charges_at(fh; T = 300.0, p = 1e6, m = 1.0)
        @test 1.0 / Eh[2] ≈ 1.5970 rtol = 1e-3

        # CO₂ critical point from the model: 304.13 K, 7.3773 MPa
        _, mdc, _ = make_fluid(:co2)
        Tc, pc, _ = Clapeyron.crit_pure(mdc)
        @test Tc ≈ 304.13 rtol = 1e-3
        @test pc ≈ 7.3773e6 rtol = 1e-3
    end
end
