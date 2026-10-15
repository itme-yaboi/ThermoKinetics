# Slice 7: ideal-gas mixtures and ideal solutions with formation data, reactions (affinity and
# mass-action laws), nuclear sources and heat sources (DESIGN §8.1–8.2, 08 §H).

@testset "reactions, nuclear and heat sources" begin
    H2, O2, H2O, N2, I2, HI = H2g(), O2g(), H2Og(), N2g(), I2g(), HIg()

    @testset "species data and atomic masses" begin
        @test molar_mass(H2) ≈ 2 * 1.008e-3 rtol = 1e-12
        @test molar_mass(H2O) ≈ (2 * 1.008 + 15.999) * 1e-3 rtol = 1e-12
        T_atom = SpeciesData(:T; atoms = (T = 1,), cv = 1.0)
        @test molar_mass(T_atom) ≈ 3.01604928e-3 rtol = 1e-9
        @test nucleons(T_atom) == 3 && nucleons(H2O) == 18
        @test H2O.cv ≈ (33.588 - R̄) / molar_mass(H2O) rtol = 1e-14           # molar data stored per kg
        @test_throws ArgumentError SpeciesData(:X; atoms = (Unobtainium = 1,), cv = 1.0)
        @test_throws ArgumentError SpeciesData(:X; atoms = (H = 0,), cv = 1.0)
        @test_throws ArgumentError SpeciesData(:X; atoms = (H = 1,), cv = -1.0)
        @test_throws ArgumentError SpeciesData(:X; atoms = (H = 1,), cv = 1.0, basis = :weird)
        # nuclides outside the table: atomic data given with the species
        Xa = SpeciesData(:Xa; atoms = (X100 = 1,), atom_data = (X100 = (M = 0.0999, A = 100),), cv = 1.0)
        @test molar_mass(Xa) == 0.0999 && nucleons(Xa) == 100
    end

    @testset "ideal-gas mixture: state function" begin
        mix = IdealGasMixture(H2, O2, H2O, N2)
        @test domains(mix) == (Thermal(), Volumetric(), Species(:H2), Species(:O2), Species(:H2O), Species(:N2))
        @test species_names(mix) == (:H2, :O2, :H2O, :N2)
        rng = MersenneTwister(71)
        gs = MolarGas.((H2, O2, H2O, N2))
        for _ in 1:200
            T = 200 + 2800rand(rng); V = 10.0^(-4 + 3rand(rng))
            n = Tuple(10.0^(-6 + 6rand(rng)) for _ in 1:4)
            E = charges_at(mix; T, V, n = (H2 = n[1], O2 = n[2], H2O = n[3], N2 = n[4]))
            F = forces(mix, E)
            @test F[1] ≈ T rtol = 1e-12
            p = sum(n) * R̄ * T / V
            @test -F[2] ≈ p rtol = 1e-12                                          # Dalton
            for k in 1:4                                                          # μ_i = h_i − T s_i (per kg)
                pk = n[k] * R̄ * T / V
                μ_molar = mh_(gs[k], T) - T * ms_(gs[k], T, pk)
                @test F[2+k] * molar_mass(mix, species_names(mix)[k]) ≈ μ_molar rtol = 1e-10 atol = 1e-7 * R̄ * T
            end
            @test internal_energy(mix, E) ≈ sum(n[k] * mu_(gs[k], T) for k in 1:4) rtol = 1e-11 atol = 1e-9 * sum(n) * R̄ * T
            λ = 0.1 + 10rand(rng)
            @test internal_energy(mix, λ * E) ≈ λ * internal_energy(mix, E) rtol = 1e-10 atol = 1e-9 * sum(n) * R̄ * T * λ
            A = eos_matrix(mix, E)
            @test maximum(abs, A - A') <= 1e-10 * maximum(abs, A)
            @test is_stable(mix, E)
        end
        # one species reproduces IdealGas (h_f = c_p T_ref, s0 = 0 at (T_ref, p_ref))
        Rn2 = R̄ / molar_mass(N2)
        g = IdealGas(:N2; R = Rn2, cv = 2.5Rn2)
        one = IdealGasMixture(SpeciesData(:N2; atoms = (N = 2,), cv = g.cv, h_f = (g.cv + g.R) * 298.15, s0 = 0.0))
        for _ in 1:50
            E = random_charges(rng, g)
            @test forces(one, E) ≈ forces(g, E) rtol = 1e-11
            @test internal_energy(one, E) ≈ internal_energy(g, E) rtol = 1e-12
        end
        @test_throws ArgumentError IdealGasMixture(H2, H2)                         # duplicate species
        # trace floor (D-41): an absent species is a valid state with finite μ; above the floor, exactly ideal
        E0 = charges_at(mix; T = 500.0, V = 1e-3, n = (H2 = 0.0, O2 = 0.01, H2O = 0.01, N2 = 0.01))
        @test valid_state(mix, E0) && all(isfinite, forces(mix, E0))
        Eneg = E0 - SVector(0.0, 0.0, 1e-3 * 1e-13 * sum(E0[3:end]), 0.0, 0.0, 0.0)   # slightly negative, round-off size
        @test valid_state(mix, Eneg) && all(isfinite, forces(mix, Eneg)) && is_stable(mix, Eneg)
        @test internal_energy(mix, 3Eneg) ≈ 3internal_energy(mix, Eneg) rtol = 1e-12          # still homogeneous
        @test !valid_state(mix, E0 - SVector(0.0, 0.0, 1e-5 * sum(E0[3:end]), 0.0, 0.0, 0.0))  # not a real negative amount
        @test_throws ArgumentError IdealGasMixture(H2, O2; trace = 0.0)
    end

    @testset "ideal solution: state function" begin
        w = liqdata(:H2O, (H = 2, O = 1); c = 75.385, h_f = -285830.0, s0 = 69.95)
        hdo = liqdata(:HDO, (H = 1, D = 1, O = 1); c = 75.385, h_f = -285830.0, s0 = 69.95)
        sol = IdealSolution(w, hdo)
        @test domains(sol) == (Thermal(), Species(:H2O), Species(:HDO))
        rng = MersenneTwister(72)
        for _ in 1:100
            T = 280 + 100rand(rng); n = (55.5 * rand(rng), 0.02 * rand(rng) + 1e-12)
            E = charges_at(sol; T, n = (H2O = n[1], HDO = n[2]))
            F = forces(sol, E)
            @test F[1] ≈ T rtol = 1e-12
            x = n ./ sum(n)
            for (k, sd) in enumerate((w, hdo))
                M = molar_mass(sd)
                h = sd.h_f * M + sd.cv * M * (T - 298.15)
                s = sd.s0 * M + sd.cv * M * log(T / 298.15) - R̄ * log(x[k])
                @test F[1+k] * M ≈ h - T * s rtol = 1e-10
            end
            @test internal_energy(sol, 2E) ≈ 2internal_energy(sol, E) rtol = 1e-12
            @test is_stable(sol, E)
        end
    end

    @testset "validation" begin
        mix = IdealGasMixture(H2, O2, H2O, N2)
        function build(el; mat = mix)
            sys = ThermoSystem()
            n = NamedTuple{species_names(mat)}(Tuple(1e-2 for _ in species_names(mat)))
            body = Volumetric() in domains(mat) ? Body(:gas, mat; T = 300.0, V = 1e-3, n) : Body(:gas, mat; T = 300.0, n)
            add!(sys, body)
            connect!(sys, el, :gas; name = :el)
            return sys
        end
        @test compile(build(Reaction(:H2 => -2, :O2 => -1, :H2O => 2; law = AffinityLaw(1e-6)))) isa CompiledSystem
        @test_throws ArgumentError compile(build(Reaction(:H2 => -1, :O2 => -1, :H2O => 1; law = AffinityLaw(1e-6))))  # atoms
        @test_throws ArgumentError compile(build(Reaction(:H2 => -1, :I2 => -1, :HI => 2; law = AffinityLaw(1e-6))))   # not in the body
        @test_throws ArgumentError AffinityLaw(-1.0)
        @test_throws ArgumentError AffinityLaw([1.0 2.0; 2.0 1.0])                         # not PSD
        @test_throws ArgumentError AffinityLaw([1.0 0.1; 0.2 1.0])                         # not symmetric
        @test_throws ArgumentError MassAction(-1.0)
        @test_throws ArgumentError Reaction(:H2 => -2, :H2 => 2; law = AffinityLaw(1.0))  # repeated species
        # mass action needs concentrations: a volume
        w = liqdata(:A, (C = 4, H = 10); c = 100.0); b = liqdata(:B, (C = 4, H = 10); c = 100.0)
        sol = IdealSolution(w, b)
        @test_throws ArgumentError compile(build(Reaction(:A => -1, :B => 1; law = MassAction(1.0)); mat = sol))
        @test compile(build(Reaction(:A => -1, :B => 1; law = AffinityLaw(1.0)); mat = sol)) isa CompiledSystem
        # nuclear: nucleons and charge
        Tg = gasdata(:T, (T = 1,); cp = 20.786); He3g = gasdata(:He3, (He3 = 1,); cp = 20.786)
        He4g = gasdata(:He4, (He4 = 1,); cp = 20.786)
        nm = IdealGasMixture(Tg, He3g, He4g)
        @test compile(build(NuclearSource(:T => -1, :He3 => 1; rate = FirstOrder(; half_life = 12.32YEAR)); mat = nm)) isa CompiledSystem
        @test_throws ArgumentError compile(build(NuclearSource(:T => -1, :He4 => 1; rate = FirstOrder(1e-9)); mat = nm))
        @test compile(build(NuclearSource(:T => -1, :He4 => 1; rate = FirstOrder(1e-9), external = (nucleons = 1, charge = 0)); mat = nm)) isa CompiledSystem
        @test_throws ArgumentError NuclearSource(:T => -2, :He3 => 2; rate = FirstOrder(1e-9))   # no single parent
        @test_throws ArgumentError FirstOrder(-1.0)
        @test_throws ArgumentError HeatSource(Inf)
        err = try; Source(Thermal() => 1.0); catch e; e; end
        @test err isa ArgumentError && occursin("HeatSource", err.msg)
    end

    @testset "reaction kernel: dU = 0, A·ξ̇ ≥ 0, elements and mass conserved" begin
        mix = IdealGasMixture(H2, O2, H2O, N2)
        rng = MersenneTwister(73)
        laws = (AffinityLaw(1e-3), MassAction(10.0; orders = (H2 = 1, O2 = 1)),
                MassAction(1e3; b = 0.5, Ea = 4e4), MassAction(1.0; orders = (H2 = 2, O2 = 1)))
        for law in laws, _ in 1:100
            sys = ThermoSystem()
            n = (H2 = 10.0^(-12 + 11rand(rng)), O2 = 10.0^(-12 + 11rand(rng)), H2O = 10.0^(-12 + 11rand(rng)), N2 = rand(rng))
            add!(sys, Body(:gas, mix; T = 300 + 2700rand(rng), V = 1e-3 * (0.1 + rand(rng)), n))
            rx = connect!(sys, Reaction(:H2 => -2, :O2 => -1, :H2O => 2; law), :gas)
            cs = compile(sys)
            u = initial_state(cs)
            du = rates(cs, u)
            r = reaction_rates(cs, u, rx)
            @test r.S_gen >= 0                                                     # exactly, in floating point
            @test sign(r.ξ̇[1]) == sign(r.A[1]) || r.ξ̇[1] == 0
            @test r.S_gen ≈ r.A[1] * r.ξ̇[1] / r.T rtol = 1e-12 atol = 1e-300
            @test entropy_production(cs, u)[rx] == r.S_gen
            F = forces(mix, charges(cs, u, :gas)); dE = charges(cs, du, :gas)
            @test abs(dot(F, dE)) <= 1e-12 * sum(abs.(F .* dE)) + 1e-300           # dU = 0
            el = element_totals(cs, du)
            @test abs(el[:H]) <= 1e-13 * abs(r.ξ̇[1]) * 4 + 1e-300 && abs(el[:O]) <= 1e-13 * abs(r.ξ̇[1]) * 2 + 1e-300
            dm = dE[3:end]
            @test abs(sum(dm)) <= 1e-14 * sum(abs, dm) + 1e-300                     # mass, by construction
            @test abs(energy_residual(cs, u).residual) <= 1e-12 * energy_residual(cs, u).scale + 1e-300
        end
        # no affinity, no rate (law level)
        @test ThermoKinetics.extent_rate(MassAction(1.0), 0.0, 300.0, 1e-3, (1.0, 1.0)) == 0.0
        @test ThermoKinetics.extent_rate(AffinityLaw(1.0), 0.0, 300.0, 1e-3, (1.0, 1.0)) == 0.0
    end

    @testset "H1 reversible isomerization: exact exponential at constant T and p" begin
        Δs = 4.0                                                                    # J/(mol K)
        A_ = gasdata(:nC4, (C = 4, H = 10); cp = 98.0, s0 = 310.0)
        B_ = gasdata(:iC4, (C = 4, H = 10); cp = 98.0, s0 = 310.0 + Δs)
        mix = IdealGasMixture(A_, B_)
        K = exp(Δs / R̄)
        kf, T0, V = 2.0, 400.0, 1e-3
        nA0, nB0 = 0.03, 1e-20 * 0.03
        sys = ThermoSystem()
        add!(sys, Body(:v, mix; T = T0, V, n = (nC4 = nA0, iC4 = nB0)))
        connect!(sys, Reaction(:nC4 => -1, :iC4 => 1; law = MassAction(kf)), :v; name = :iso)
        cs = compile(sys)
        ts = [0.05, 0.2, 0.5, 1.0, 3.0]
        sol = simulate(cs, (0.0, 3.0); saveat = ts)
        N = nA0 + nB0
        neq = N / (1 + K)
        for (t, u) in zip(sol.t, sol.u)
            @test moles(cs, u, :v).nC4 ≈ neq + (nA0 - neq) * exp(-kf * (1 + 1 / K) * t) rtol = 1e-7
            @test temperature(cs, u, :v) ≈ T0 rtol = 1e-9
            @test pressure(cs, u, :v) ≈ N * R̄ * T0 / V rtol = 1e-9
        end
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-7 && cl.entropy_rel <= 1e-9 && cl.charge_rel <= 1e-12
    end

    @testset "H2 equilibrium composition: Gibbs minimization and entropy maximization" begin
        mix = IdealGasMixture(H2, I2, HI)
        gH2, gI2, gHI = MolarGas(H2), MolarGas(I2), MolarGas(HI)
        V, a, b, ε = 1e-3, 0.010, 0.008, 1e-14
        n_of(ξ) = (a - ξ, b - ξ, ε + 2ξ)
        @testset "at fixed T (vessel on a reservoir)" begin
            T = 700.0
            sys = ThermoSystem()
            add!(sys, Body(:v, mix; T, V, n = (H2 = a, I2 = b, HI = ε)))
            add!(sys, Reservoir(:bath; T))
            connect!(sys, Conductor(ThermalConductance(1.0)), :v, :bath)
            connect!(sys, Reaction(:H2 => -1, :I2 => -1, :HI => 2; law = MassAction(1.0)), :v; name = :rx)
            cs = compile(sys)
            sol = simulate(cs, (0.0, 200.0))
            uf = sol.u[end]
            # independent: minimize G(ξ) = Σ n μ(T, p_i) (Δn = 0, so the same minimizer as F at fixed V)
            G(ξ) = (n = n_of(ξ); sum(n[k] * (mh_(g, T) - T * ms_(g, T, n[k] * R̄ * T / V)) for (k, g) in enumerate((gH2, gI2, gHI))))
            ξeq = golden_max(ξ -> -G(ξ), -ε / 2 + 1e-18, min(a, b) - 1e-18)
            nf = moles(cs, uf, :v)
            @test nf.HI ≈ n_of(ξeq)[3] rtol = 1e-7
            @test nf.H2 ≈ n_of(ξeq)[1] rtol = 1e-7
            @test abs(reaction_rates(cs, uf, :rx).A[1]) <= 1e-6 * R̄ * T                # A → 0 (contact rule)
            cl = closure(cs, sol)
            @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9 && cl.charge_rel <= 1e-12
        end
        @testset "adiabatic at fixed V: entropy maximum at fixed U, V and elements" begin
            T0 = 700.0
            sys = ThermoSystem()
            add!(sys, Body(:v, mix; T = T0, V, n = (H2 = a, I2 = b, HI = ε)))
            connect!(sys, Reaction(:H2 => -1, :I2 => -1, :HI => 2; law = MassAction(1.0)), :v; name = :rx)
            cs = compile(sys)
            sol = simulate(cs, (0.0, 200.0))
            gs = (gH2, gI2, gHI)
            U0 = sum(n_of(0.0)[k] * mu_(gs[k], T0) for k in 1:3)
            function Tof(ξ)
                n = n_of(ξ)
                return 298.15 + (U0 - sum(n[k] * (gs[k].hf - R̄ * 298.15) for k in 1:3)) / sum(n[k] * gs[k].cv for k in 1:3)
            end
            S(ξ) = (n = n_of(ξ); T = Tof(ξ); sum(n[k] * ms_(gs[k], T, n[k] * R̄ * T / V) for k in 1:3))
            ξeq = golden_max(S, -ε / 2 + 1e-18, min(a, b) - 1e-18)
            nf = moles(cs, sol.u[end], :v)
            @test nf.HI ≈ n_of(ξeq)[3] rtol = 1e-7
            @test temperature(cs, sol.u[end], :v) ≈ Tof(ξeq) rtol = 1e-7
            @test total_entropy(cs, sol.u[end]) - total_entropy(cs, sol.u[1]) ≈ generated_entropy(cs, sol.u[end])[:rx] rtol = 1e-9
        end
    end

    # lean H₂/O₂ in N₂ (2 % H₂, 10 % O₂), 300 K, 1 bar, 1 L; fast first-order mass action
    comb_n(ε) = (p = 1e5; V = 1e-3; N = p * V / (R̄ * 300.0);
                 (H2 = 0.02N, O2 = 0.10N, H2O = ε * N, N2 = 0.88N))
    comb_law = MassAction(1e3; orders = (H2 = 1, O2 = 1))
    gas4 = (MolarGas(H2), MolarGas(O2), MolarGas(H2O), MolarGas(N2))
    comb_U(n, T) = sum(n[k] * mu_(gas4[k], T) for k in 1:4)

    @testset "H3 adiabatic constant-volume combustion (and H10: approach to zero)" begin
        mix = IdealGasMixture(H2, O2, H2O, N2)
        n0 = comb_n(1e-20)
        sys = ThermoSystem()
        add!(sys, Body(:v, mix; T = 300.0, V = 1e-3, n = n0))
        connect!(sys, Reaction(:H2 => -2, :O2 => -1, :H2O => 2; law = comb_law), :v; name = :burn)
        cs = compile(sys)
        sol = simulate(cs, (0.0, 0.2))
        uf = sol.u[end]
        nf = moles(cs, uf, :v)
        @test nf.H2 <= 1e-12 * n0.H2                                             # H10: driven to ≤ 1e-12, no failure
        mtot = sum(charges(cs, uf, :v)[3:end])
        @test all(u -> all(x -> x >= -1e-11 * mtot, charges(cs, u, :v)[3:end]), sol.u)   # at worst a trace-floor negative
        @test length(sol.t) < 2000                                                # no stall
        @test all(u -> reaction_rates(cs, u, :burn).S_gen >= 0, sol.u)
        # complete conversion: final T from U conservation (closed form, constant c_v)
        ncomp = (H2 = 0.0, O2 = n0.O2 - n0.H2 / 2, H2O = n0.H2O + n0.H2, N2 = n0.N2)
        U0 = comb_U(n0, 300.0)
        Tf = 298.15 + (U0 - sum(ncomp[k] * (gas4[k].hf - R̄ * 298.15) for k in 1:4)) / sum(ncomp[k] * gas4[k].cv for k in 1:4)
        @test temperature(cs, uf, :v) ≈ Tf rtol = 1e-8
        @test total_entropy(cs, uf) - total_entropy(cs, sol.u[1]) ≈ generated_entropy(cs, uf)[:burn] rtol = 1e-9
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9 && cl.charge_rel <= 1e-12
    end

    @testset "H4 isothermal calorimeter: heat = formation-energy difference" begin
        mix = IdealGasMixture(H2, O2, H2O, N2)
        n0 = comb_n(1e-20)
        sys = ThermoSystem()
        add!(sys, Body(:v, mix; T = 300.0, V = 1e-3, n = n0))
        add!(sys, Reservoir(:bath; T = 300.0))
        connect!(sys, Conductor(ThermalConductance(1.0)), :v, :bath)
        connect!(sys, Reaction(:H2 => -2, :O2 => -1, :H2O => 2; law = comb_law), :v; name = :burn)
        cs = compile(sys)
        sol = simulate(cs, (0.0, 40.0))
        ncomp = (H2 = 0.0, O2 = n0.O2 - n0.H2 / 2, H2O = n0.H2O + n0.H2, N2 = n0.N2)
        @test temperature(cs, sol.u[end], :v) ≈ 300.0 rtol = 1e-9
        @test reservoir_ledger(cs, sol.u[end], :bath).U ≈ comb_U(n0, 300.0) - comb_U(ncomp, 300.0) rtol = 1e-7
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9 && cl.charge_rel <= 1e-12
    end

    @testset "H5 coupled reactions (Onsager): production ≥ 0; equilibrium independent of L" begin
        sA = gasdata(:A, (C = 4, H = 10); cp = 98.0, s0 = 300.0)
        sB = gasdata(:B, (C = 4, H = 10); cp = 98.0, s0 = 305.0)
        sC = gasdata(:C, (C = 4, H = 10); cp = 98.0, s0 = 302.0)
        mix = IdealGasMixture(sA, sB, sC)
        K1, K2 = exp(5.0 / R̄), exp(-3.0 / R̄)
        finals = []
        for L in ([2.0 0.5; 0.5 1.0], [1.0 -0.9; -0.9 1.0], [0.3 0.0; 0.0 3.0])
            sys = ThermoSystem()
            add!(sys, Body(:v, mix; T = 350.0, V = 1e-3, n = (A = 0.03, B = 1e-15, C = 1e-15)))
            connect!(sys, Reaction([[:A => -1, :B => 1], [:B => -1, :C => 1]]; law = AffinityLaw(1e-4 * L)), :v; name = :rx)
            cs = compile(sys)
            sol = simulate(cs, (0.0, 2000.0); saveat = 10.0)
            @test all(u -> entropy_production(cs, u)[:rx] >= 0, sol.u)
            n = moles(cs, sol.u[end], :v)
            @test n.B / n.A ≈ K1 rtol = 1e-7
            @test n.C / n.B ≈ K2 rtol = 1e-7
            push!(finals, n)
            cl = closure(cs, sol)
            @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9 && cl.charge_rel <= 1e-12
        end
        @test finals[1].A ≈ finals[2].A rtol = 1e-7
        @test finals[1].A ≈ finals[3].A rtol = 1e-7
    end

    @testset "08 §A on random reacting networks (A1, A3, elements, A → 0)" begin
        rng = MersenneTwister(74)
        mix = IdealGasMixture(H2, I2, HI, N2)
        for trial in 1:6
            sys = ThermoSystem()
            nb = 3
            for i in 1:nb
                add!(sys, Body(Symbol(:g, i), mix; T = 500 + 400rand(rng), V = 1e-3 * (0.5 + rand(rng)),
                               n = (H2 = 0.01rand(rng) + 1e-6, I2 = 0.01rand(rng) + 1e-6, HI = 0.01rand(rng) + 1e-6, N2 = 0.02rand(rng))))
                law = rand(rng) < 0.5 ? MassAction(0.1 + rand(rng)) : AffinityLaw(1e-5 * (0.1 + rand(rng)))
                connect!(sys, Reaction(:H2 => -1, :I2 => -1, :HI => 2; law), Symbol(:g, i); name = Symbol(:rx, i))
            end
            for (i, j) in ((1, 2), (2, 3))
                connect!(sys, Conductor(ThermalConductance(0.01 + 0.05rand(rng))), Symbol(:g, i), Symbol(:g, j))
                connect!(sys, Conductor(MassConductance(:HI, 1e-11 * (0.5 + rand(rng))); dissipation = Split(rand(rng))),
                         Symbol(:g, i), Symbol(:g, j))
                connect!(sys, Conductor(MassConductance(:H2, 1e-11 * (0.5 + rand(rng))); dissipation = Split(rand(rng))),
                         Symbol(:g, i), Symbol(:g, j))
            end
            add!(sys, Reservoir(:env; T = 600.0))
            connect!(sys, Conductor(ThermalConductance(0.02)), :g1, :env)
            cs = compile(sys)
            sol = simulate(cs, (0.0, 2e5); saveat = 2e3)
            for u in sol.u
                @test all(>=(0), values(entropy_production(cs, u)))
            end
            cl = closure(cs, sol)
            @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9 && cl.charge_rel <= 1e-11
            for i in 1:nb
                @test abs(reaction_rates(cs, sol.u[end], Symbol(:rx, i)).A[1]) <= 1e-5 * R̄ * 600.0
            end
        end
    end

    # tritium as an atomic species; 5.7 keV deposited per decay
    λT = log(2) / (12.32YEAR)
    EdepT = 5.7 * keV * N_AV                                                    # J per mol of decays

    @testset "H6 tritium decay (and H10: 60 half-lives)" begin
        Tg = gasdata(:T, (T = 1,); cp = 20.786, s0 = 153.0); He3g = gasdata(:He3, (He3 = 1,); cp = 20.786, s0 = 126.0)
        mix = IdealGasMixture(Tg, He3g)
        n0 = 1.0
        sys = ThermoSystem()
        add!(sys, Body(:tank, mix; T = 300.0, V = 1.0, n = (T = n0, He3 = 1e-20)))
        connect!(sys, NuclearSource(:T => -1, :He3 => 1; rate = FirstOrder(; half_life = 12.32YEAR), E_dep = EdepT), :tank; name = :decay)
        cs = compile(sys)
        u0 = initial_state(cs)
        du0 = rates(cs, u0)
        P0 = -boundary_ledger(cs, du0, :decay).U                                  # deposited power at t = 0
        @test P0 ≈ λT * n0 * EdepT rtol = 1e-12
        @test P0 / (n0 * molar_mass(Tg) * 1e3) ≈ 0.324 rtol = 0.01                 # ≈ 0.32 W per g of T
        @test abs(energy_residual(cs, u0).residual) <= 1e-12 * energy_residual(cs, u0).scale
        ts = [1, 5, 12.32, 30, 60 * 12.32] .* YEAR
        sol = simulate(cs, (0.0, ts[end]); saveat = ts)
        for (t, u) in zip(sol.t, sol.u)
            n = moles(cs, u, :tank)
            if t < 100YEAR
                @test n.T ≈ n0 * exp(-λT * t) rtol = 1e-6
            else                                                                     # 60 half-lives: n_T ≈ 1e-18 n0
                @test abs(n.T) <= 1e-10 * n0
            end
            @test 3n.T + 3n.He3 ≈ 3n0 rtol = 1e-12                                  # nucleons
            @test total_energy(cs, u) - total_energy(cs, u0) ≈ EdepT * (n0 - n.T) rtol = 1e-7         # dU = E_dep ξ (exact in the rates)
        end
        @test length(sol.t) == length(ts)                                             # completed, no failure
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-7 && cl.entropy_rel <= 1e-9 && cl.charge_rel <= 1e-12
    end

    @testset "H7 Bateman chain A → B → C" begin
        ad = (A100 = (M = 0.0999, A = 100), B100 = (M = 0.09989, A = 100), C100 = (M = 0.09988, A = 100))
        mk(name, atom) = liqdata(name, NamedTuple{(atom,)}((1,)); c = 25.0, atom_data = ad)
        sol_ = IdealSolution(mk(:A, :A100), mk(:B, :B100), mk(:C, :C100))
        λA, λB = 1 / 100.0, 1 / 30.0
        EA, EB = 1e8, 3e8                                                          # J/mol
        n0 = 2.0
        sys = ThermoSystem()
        add!(sys, Body(:src, sol_; T = 300.0, n = (A = n0, B = 1e-20, C = 1e-20)))
        add!(sys, Reservoir(:wall; T = 300.0))
        connect!(sys, Conductor(ThermalConductance(50.0)), :src, :wall)
        connect!(sys, NuclearSource(:A => -1, :B => 1; rate = FirstOrder(λA), E_dep = EA), :src; name = :dA)
        connect!(sys, NuclearSource(:B => -1, :C => 1; rate = FirstOrder(λB), E_dep = EB), :src; name = :dB)
        cs = compile(sys)
        ts = [10.0, 50.0, 100.0, 300.0, 1000.0]
        sol = simulate(cs, (0.0, 1000.0); saveat = ts)
        for (t, u) in zip(sol.t, sol.u)
            n = moles(cs, u, :src)
            nA = n0 * exp(-λA * t)
            nB = n0 * λA / (λB - λA) * (exp(-λA * t) - exp(-λB * t))
            @test n.A ≈ nA rtol = 1e-6
            @test n.B ≈ nB rtol = 1e-6
            @test n.C ≈ n0 - nA - nB rtol = 1e-6
            heat = -boundary_ledger(cs, u, :dA).U - boundary_ledger(cs, u, :dB).U
            @test heat ≈ EA * (n0 - nA) + EB * (n0 - nA - nB) rtol = 1e-6
        end
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9 && cl.charge_rel <= 1e-11
    end

    @testset "H8 heat source: exact exponential, boundary entropy closes" begin
        m, c, G, Q, Tr, T0 = 2.0, 500.0, 3.0, 150.0, 300.0, 280.0
        sys = ThermoSystem()
        add!(sys, Body(:slab, ConstantCapacitySolid(m, c); T = T0))
        add!(sys, Reservoir(:env; T = Tr))
        connect!(sys, Conductor(ThermalConductance(G)), :slab, :env)
        connect!(sys, HeatSource(Q), :slab; name = :heater)
        cs = compile(sys)
        τ = m * c / G
        ts = [0.1, 0.5, 1.0, 3.0] .* τ
        sol = simulate(cs, (0.0, ts[end]); saveat = ts)
        for (t, u) in zip(sol.t, sol.u)
            @test temperature(cs, u, :slab) ≈ Tr + Q / G + (T0 - Tr - Q / G) * exp(-t / τ) rtol = 1e-9
            @test -boundary_ledger(cs, u, :heater).U ≈ Q * t rtol = 1e-12
        end
        # the boundary entropy is ∫Q/T dt (negative received), not production
        @test !haskey(entropy_production(cs, sol.u[end]), :heater)
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-7 && cl.entropy_rel <= 1e-9
        # time-varying heat
        sys2 = ThermoSystem()
        add!(sys2, Body(:slab, ConstantCapacitySolid(m, c); T = T0))
        connect!(sys2, HeatSource(t -> Q * (1 + sin(t))), :slab; name = :heater)
        cs2 = compile(sys2)
        sol2 = simulate(cs2, (0.0, 5.0))
        @test total_energy(cs2, sol2.u[end]) - total_energy(cs2, sol2.u[1]) ≈ Q * (5.0 + 1 - cos(5.0)) rtol = 1e-7
    end
end
