# 08 §A invariants as seeded property tests on random thermal networks.
# A2 (non-thermal charge) is vacuous in slice 1: there are no non-thermal charges yet.
# A4 (Hessian symmetry) and A5 (conductance sign) are covered in test_build.jl;
# A7 (separability) needs a second domain and arrives in slice 2.

@testset "08 §A invariants" begin
    rng = MersenneTwister(20261014)

    @testset "A1 energy (isolated networks)" begin
        for trial in 1:15
            sys, bodies, _, _, p = random_network(rng; nbodies = rand(rng, 2:8))
            cs = compile(sys)
            u0 = initial_state(cs)
            # instantaneous: every connection closes exactly (round-off)
            r = energy_residual(cs, u0)
            @test abs(r.residual) <= 1e-12 * r.scale
            # integrated: drift within integrator tolerance (DESIGN §12.3), reported as closure
            sol = simulate(cs, (0.0, 5 / slowest_rate(p)))
            U0 = total_energy(cs, u0)
            @test maximum(abs(total_energy(cs, u) - U0) for u in sol.u) <= 1e-6 * U0
            cl = closure(cs, sol)
            @test cl.energy_rel <= 1e-6
        end
    end

    @testset "A1 energy closure with reservoirs" begin
        for trial in 1:10
            sys, _, res, _, p = random_network(rng; nbodies = rand(rng, 2:6), nreservoirs = rand(rng, 1:3))
            cs = compile(sys)
            sol = simulate(cs, (0.0, 2 / slowest_rate(p)))
            cl = closure(cs, sol)
            @test cl.energy_rel <= 1e-6          # bodies + energy received by reservoirs is constant
            @test cl.entropy_rel <= 1e-10        # entropy closure is a linear invariant: round-off
        end
    end

    @testset "A3 second law" begin
        for trial in 1:15
            sys, bodies, res, conds, p = random_network(rng; nbodies = rand(rng, 2:8), nreservoirs = rand(rng, 0:2))
            cs = compile(sys)
            # per-conductor production ≥ 0 at random states, including ΔT/T ~ 1
            for _ in 1:20
                u = initial_state(cs)
                for b in bodies
                    set_temperature!(cs, u, b, 10.0^(1 + 3 * rand(rng)))
                end
                @test all(≥(0), values(entropy_production(cs, u)))
            end
            # isolated part: S_total is non-decreasing along the trajectory
            if isempty(res)
                sol = simulate(cs, (0.0, 5 / slowest_rate(p)); saveat = 0.05 / slowest_rate(p))
                S = [total_entropy(cs, u) for u in sol.u]
                @test all(diff(S) .>= -1e-10 * maximum(abs, S))
                # and the increase equals the integrated, attributed production
                Sgen = sum(values(generated_entropy(cs, sol.u[end])))
                @test S[end] - S[1] ≈ Sgen rtol = 1e-9
            end
        end
    end

    @testset "A6 equilibrium (contact rule)" begin
        for trial in 1:15
            sys, bodies, _, conds, p = random_network(rng; nbodies = rand(rng, 2:8))
            cs = compile(sys)
            sol = simulate(cs, (0.0, 40 / slowest_rate(p)))
            uf = sol.u[end]
            Tf = sum(p.m .* p.c .* p.T0) / sum(p.m .* p.c)   # energy-weighted (constant c)
            for b in bodies
                @test temperature(cs, uf, b) ≈ Tf atol = 1e-6
            end
            # forces equal across every open connection
            for c in conds
                f = flux(cs, uf, c)
                @test abs(f.T_a - f.T_b) <= 1e-6
            end
        end
    end
end
