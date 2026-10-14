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

# The same invariants on networks of ideal gases joined by diathermal walls and damped
# pistons (volume conductors whose friction heat goes to an internal node).
@testset "08 §A invariants: gas networks" begin
    rng = MersenneTwister(77)
    T_END = 3000.0     # ≥ 30 of the slowest time constants for these parameter ranges

    @testset "A1, A2, A3 (isolated)" begin
        for trial in 1:12
            sys, bodies, _, edges = random_gas_network(rng; nbodies = rand(rng, 2:6))
            cs = compile(sys)
            u0 = initial_state(cs)
            sol = simulate(cs, (0.0, T_END); saveat = T_END / 200)
            cl = closure(cs, sol)
            @test cl.energy_rel <= 1e-6                                  # A1, integrated
            @test cl.entropy_rel <= 1e-10
            # A2: every non-thermal charge total is conserved to round-off
            for d in (Volumetric(), Species(:He), Species(:N2))
                Q0 = total_charge(cs, u0, d)
                @test maximum(abs(total_charge(cs, u, d) - Q0) for u in sol.u) <= 1e-13 * max(abs(Q0), 1e-30)
            end
            @test cl.charge_rel <= 1e-10
            # A3: S_total non-decreasing, per-element production ≥ 0 along the trajectory
            S = [total_entropy(cs, u) for u in sol.u]
            @test all(diff(S) .>= -1e-10 * maximum(abs, S))
            @test all(all(≥(0), values(entropy_production(cs, u))) for u in sol.u[1:20:end])
        end
    end

    @testset "A6 equilibrium: T and p equal across every open connection" begin
        for trial in 1:12
            sys, bodies, _, edges = random_gas_network(rng; nbodies = rand(rng, 2:6))
            cs = compile(sys)
            uf = simulate(cs, (0.0, T_END)).u[end]
            for (e, kind) in edges
                f = flux(cs, uf, e)
                @test abs(f.T_a - f.T_b) <= 1e-5 * f.T_a
                kind === :piston && @test abs(pressure(cs, uf, f.a) - pressure(cs, uf, f.b)) <= 1e-5 * pressure(cs, uf, f.a)
            end
        end
    end

    @testset "A1, A2 closure with reservoirs" begin
        for trial in 1:8
            sys, bodies, res, edges = random_gas_network(rng; nbodies = rand(rng, 2:5), nreservoirs = rand(rng, 1:2))
            cs = compile(sys)
            sol = simulate(cs, (0.0, T_END / 3))
            cl = closure(cs, sol)
            @test cl.energy_rel <= 1e-6
            @test cl.entropy_rel <= 1e-10
            # body volume + volume received by reservoirs: a linear invariant, so round-off, but
            # amplified by the conditioning of each implicit (Rosenbrock) solve over many steps
            @test cl.charge_rel <= 1e-10
        end
    end
end
