# Conductor bookkeeping identities (DESIGN §5.1, kernel A, thermal row), checked on
# the kernel itself before any bodies or integrators are involved.

@testset "kernel A, thermal row: identities" begin
    rng = MersenneTwister(1)
    for _ in 1:10_000
        G = 10.0^(4 * rand(rng) - 2)                      # 1e-2 … 1e2 W/K
        Ta = 10.0^(1 + 3 * rand(rng))                     # 10 … 1e4 K: ΔT/T up to ~1
        Tb = rand(rng) < 0.1 ? Ta : 10.0^(1 + 3 * rand(rng))
        tr = VeinikThermo.thermal_transfer(ThermalConductance(G), Ta, Tb)
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
