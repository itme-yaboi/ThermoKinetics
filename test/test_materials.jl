# Materials and state functions (slice 2): ideal gas, elastic volume, composites.
# 08 §A4 (Maxwell symmetry), §A7 (uncoupled separability), §B3 (c_p − c_v).

using ForwardDiff: hessian

@testset "materials" begin
    rng = MersenneTwister(31)
    He = gas_he(); N2 = gas_n2()

    @testset "domain declarations" begin
        @test charge_symbol(Volumetric()) == :V && force_unit(Volumetric()) == "Pa"
        @test force_symbol(Volumetric()) == Symbol("-p")       # the force is −p, not p
        @test charge_symbol(Species(:He)) == :m && force_symbol(Species(:He)) == :μ
        @test Species(:He) == Species(:He) && Species(:He) != Species(:N2)
        @test !producible(Volumetric()) && !producible(Species(:He))
    end

    @testset "ideal gas state function" begin
        @test domains(He) == (Thermal(), Volumetric(), Species(:He))
        for _ in 1:200
            T = 100 + 1900rand(rng); p = 10.0^(3 + 4rand(rng)); m = 10.0^(-4 + 5rand(rng))
            E = charges_at(N2; T, p, m)
            P = forces(N2, E)
            @test P[1] ≈ T rtol = 1e-12                                 # T = ∂U/∂S
            @test -P[2] ≈ p rtol = 1e-12                                # −p = ∂U/∂V
            @test E[2] ≈ m * N2.R * T / p rtol = 1e-14                  # pV = mRT
            @test internal_energy(N2, E) ≈ m * N2.cv * T rtol = 1e-12
            # μ = g = h − Ts per unit mass, with h = c_p T and s = S/m
            @test P[3] ≈ (N2.cv + N2.R) * T - T * E[1] / m rtol = 1e-10 atol = 1e-6 * N2.cv * T
            # Euler (first-order homogeneity): U = T S − p V + μ m
            @test internal_energy(N2, E) ≈ P[1] * E[1] + P[2] * E[2] + P[3] * E[3] rtol = 1e-10
            @test valid_state(N2, E)
        end
        @test charges_at(He; T = 300.0, V = 2e-3, m = 1e-3) ≈ charges_at(He; T = 300.0, p = 1e-3 * He.R * 300 / 2e-3, m = 1e-3)
        @test !valid_state(He, [0.0, -1e-3, 1e-3])      # negative volume
        @test !valid_state(He, [0.0, 1e-3, 0.0])        # no mass: T and μ undefined
        @test_throws ArgumentError IdealGas(:X; R = -1.0, cv = 700.0)
        @test_throws ArgumentError charges_at(He; T = 300.0, m = 1.0)               # needs two of m, p, V
        @test_throws ArgumentError charges_at(He; T = 300.0, m = 1.0, p = 1e5, V = 1.0)
    end

    @testset "A4 Maxwell reciprocity: raw Hessians are symmetric" begin
        mats = (He, N2, ConstantCapacitySolid(2.0, 500.0),
                Composite(ConstantCapacitySolid(2.0, 500.0), ElasticVolume(1e-9; V0 = 1e-3, p0 = 1e5)))
        for mat in mats, _ in 1:100
            E = random_charges(rng, mat)
            H = hessian(e -> internal_energy(mat, e), collect(E))
            @test H ≈ H' rtol = 1e-12 atol = 1e-12 * maximum(abs, H)
        end
        # (∂T/∂V)_S = −(∂p/∂S)_V for the gas, i.e. ∂P_S/∂V = ∂P_V/∂S with P_V = −p
        E = charges_at(N2; T = 400.0, p = 2e5, m = 0.01)
        A = eos_matrix(N2, E)
        @test A[1, 2] == A[2, 1]
        @test A[1, 2] < 0      # compression (dV < 0) at fixed S heats the gas
    end

    @testset "homogeneity, stability, fixed-mass capacities" begin
        for _ in 1:50
            E = random_charges(rng, N2)
            A = eos_matrix(N2, E)
            # U homogeneous of degree 1 ⇒ the full Hessian is singular along E itself
            @test norm(A * E) <= 1e-8 * norm(A) * norm(E)
            @test is_stable(N2, E)                                   # fixed-mass block is PD
            @test mass_slots(N2) == (3,)
        end
        @test is_stable(ConstantCapacitySolid(1.0, 400.0), charges_at(ConstantCapacitySolid(1.0, 400.0); T = 300.0))
    end

    @testset "B3 c_p − c_v = R (Mayer)" begin
        for gas in (He, N2), _ in 1:50
            T = 100 + 1900rand(rng); p = 10.0^(3 + 4rand(rng)); m = 10.0^(-4 + 5rand(rng))
            E = charges_at(gas; T, p, m)
            Kp = capacity_matrix(gas, E; hold = :forces)     # other forces held: c_p-like
            Kv = capacity_matrix(gas, E; hold = :charges)    # other charges held: c_v-like
            @test size(Kp) == (2, 2)                          # the (S, V) block only
            @test (T / m) * (Kp[1, 1] - Kv[1, 1]) ≈ gas.R rtol = 1e-9
            @test (T / m) * Kv[1, 1] ≈ gas.cv rtol = 1e-12    # dS/dT at fixed V = m c_v / T
            @test Kp[2, 2] > Kv[2, 2] > 0                     # isothermal compressibility > isentropic
        end
    end

    @testset "A7 uncoupled separability" begin
        solid = ConstantCapacitySolid(2.0, 500.0)
        spring = ElasticVolume(2e-9; V0 = 1e-3, p0 = 1e5)
        body = Composite(solid, spring)
        @test domains(body) == (Thermal(), Volumetric())
        A = eos_matrix(body, charges_at(body, (T = 350.0,), (V = 1.1e-3,)))
        @test A[1, 2] == 0 && A[2, 1] == 0                    # block-diagonal by part
        for _ in 1:100
            E = charges_at(body, (T = 200 + 800rand(rng),), (V = 1e-3 * (0.5 + rand(rng)),))
            P = forces(body, E)
            dS, dV = randn(rng), 1e-4 * randn(rng)
            @test forces(body, E + SA[0.0, dV])[1] == P[1]    # T untouched by a volume change
            @test forces(body, E + SA[dS, 0.0])[2] == P[2]    # p untouched by an entropy change
        end
        # contrast: the ideal gas is a coupled group, so the same perturbation moves T
        E = charges_at(N2; T = 300.0, p = 1e5, m = 0.01)
        @test forces(N2, E + SA[0.0, -1e-4, 0.0])[1] > forces(N2, E)[1]

        @test_throws ArgumentError Composite(He, spring)                   # both have Volumetric
        @test_throws ArgumentError Composite(solid, ConstantCapacitySolid(1.0, 1.0))  # two thermal parts
    end

    @testset "elastic volume" begin
        e = ElasticVolume(1e-9; V0 = 1e-3, p0 = 2e5)
        E = charges_at(e; p = 3e5)
        @test -forces(e, E)[1] ≈ 3e5 rtol = 1e-12
        @test E[1] ≈ 1e-3 - 1e-9 * 1e5 rtol = 1e-12          # compressing raises p
        @test capacity_matrix(e, E)[1, 1] ≈ 1e-9 rtol = 1e-12 # dV/d(−p) = C > 0
        @test_throws ArgumentError ElasticVolume(-1.0; V0 = 1e-3)
    end
end
