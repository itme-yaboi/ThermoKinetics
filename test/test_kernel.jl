# Conductor bookkeeping identities (DESIGN §5.1, kernel A, thermal row), checked on
# the kernel itself before any bodies or integrators are involved.

@testset "kernel A, thermal row: identities" begin
    rng = MersenneTwister(1)
    for _ in 1:10_000
        G = 10.0^(4 * rand(rng) - 2)                      # 1e-2 … 1e2 W/K
        Ta = 10.0^(1 + 3 * rand(rng))                     # 10 … 1e4 K: ΔT/T up to ~1
        Tb = rand(rng) < 0.1 ? Ta : 10.0^(1 + 3 * rand(rng))
        tr = ThermoKinetics.thermal_transfer(ThermalConductance(G), Ta, Tb)
        # energy: heat deliveries sum to the released non-thermal power, which is zero here
        @test tr.Q_a + tr.Q_b == 0
        # the energy-form law is reproduced exactly at any ΔT
        @test tr.Q_b ≈ G * (Ta - Tb) rtol = 1e-13
        # entropy production is what was written to the two bodies, and is ≥ 0 exactly
        @test tr.S_gen == tr.Q_a / Ta + tr.Q_b / Tb
        @test tr.S_gen >= 0
        @test tr.S_gen ≈ G * (Ta - Tb) * (1 / Tb - 1 / Ta) rtol = 1e-10 atol = 1e-300
    end
end

@testset "kernel A through the assembled right-hand side" begin
    rng = MersenneTwister(2)
    for trial in 1:50
        sys, bodies, res, conds, _ = random_network(rng; nbodies = rand(rng, 2:8), nreservoirs = rand(rng, 0:2))
        cs = compile(sys)
        u = initial_state(cs)
        du = rates(cs, u)
        # every entropy write is accounted: Σ dS_bodies + Σ dS_reservoirs = Σ Ṡ_gen (round-off)
        dS_bodies = sum(charges(cs, du, b)[1] for b in bodies)   # the only charge is S
        dS_res = isempty(res) ? 0.0 : sum(reservoir_ledger(cs, du, r).S for r in res)
        Sgen = sum(values(entropy_production(cs, u)))
        @test dS_bodies + dS_res ≈ Sgen rtol = 1e-12 atol = 1e-15
        # energy: Σ T_i dS_i over bodies + energy into reservoirs = 0 (instantaneous, round-off)
        r = energy_residual(cs, u)
        @test abs(r.residual) <= 1e-12 * r.scale
    end
end

# Non-thermal rows with a dissipation policy (DESIGN §5.1): Split/Into (kernel A) and the
# internal thermal node (kernel B). Checked on the heat-distribution kernel directly.
@testset "dissipation policies: identities" begin
    rng = MersenneTwister(3)
    W = ThermoKinetics.WeightsPolicy
    N = ThermoKinetics.NodePolicy
    for _ in 1:10_000
        D = rand(rng) < 0.05 ? 0.0 : 10.0^(6rand(rng) - 3)         # released power ≥ 0 [W]
        Ta, Tb, Th = (10.0^(1 + 3rand(rng)) for _ in 1:3)
        w = rand(rng)
        for pol in (W(1 - w, w, 0.0), W(0.0, 0.0, 1.0),
                    N(10.0^(4rand(rng) - 2), 10.0^(4rand(rng) - 2), rand(rng) < 0.5 ? 0.0 : 10.0^(4rand(rng) - 2)))
            h = ThermoKinetics.distribute_heat(pol, D, Ta, Tb, Th)
            # energy: remainder rule, closes to round-off
            @test abs(h.Q_a + h.Q_b + h.Q_h - D) <= 4eps() * (abs(h.Q_a) + abs(h.Q_b) + abs(h.Q_h) + D)
            @test h.S_gen >= 0                                        # sign-exact ledger
            written = h.Q_a / Ta + h.Q_b / Tb + h.Q_h / Th            # entropy actually written
            @test written ≈ h.S_gen rtol = 1e-9 atol = 1e-12 * (D + abs(h.Q_a) + abs(h.Q_b)) / min(Ta, Tb, Th)
        end
        # kernel B: the node temperature balances its legs, and production splits in two parts
        pol = N(1.0 + rand(rng), 1.0 + rand(rng), 0.0)
        h = ThermoKinetics.distribute_heat(pol, D, Ta, Tb, Th)
        @test pol.G_a * (h.T_m - Ta) + pol.G_b * (h.T_m - Tb) ≈ D rtol = 1e-9 atol = 1e-9 * (Ta + Tb)
        @test h.parts.dissipation ≈ D / h.T_m rtol = 1e-12
        @test h.parts.dissipation + h.parts.conduction == h.S_gen
    end
end

@testset "volume conductor through the assembled right-hand side" begin
    rng = MersenneTwister(4)
    for _ in 1:30
        sys, bodies, res, edges = random_gas_network(rng; nbodies = rand(rng, 2:6), nreservoirs = rand(rng, 0:2))
        cs = compile(sys)
        u = initial_state(cs)
        r = energy_residual(cs, u)
        @test abs(r.residual) <= 1e-12 * r.scale
        @test all(≥(0), values(entropy_production(cs, u)))
        # volume moved out of one port arrives at the other: the volume rates sum to zero
        du = rates(cs, u)
        dV = sum(charges(cs, du, b)[2] for b in bodies) +
             sum((reservoir_charge(cs, du, rr, Volumetric()) for rr in res); init = 0.0)
        @test abs(dV) <= 1e-14 * sum(abs(flux(cs, u, e).J[1]) for (e, k) in edges if k == :piston; init = 1e-30)
    end
end
