# Slice 6: mechanical systems (DESIGN §18). Elastic forms of motion (displacement, angle), springs
# and the pendulum potential, spring branches (converters), monotone friction laws, force and torque
# sources, gravity, gears and rack-and-pinion through transformers. 08 §G cases G3, G4, G9 here; G1,
# G2, G5–G8 are application tests (test_applications.jl).

@testset "mechanical systems" begin
    rng = MersenneTwister(20261015)

    @testset "construction and validation" begin
        @test charge_symbol(Displacement()) == :x && force_symbol(Displacement()) == :F
        @test charge_unit(Displacement()) == "m" && force_unit(Displacement()) == "N"
        @test charge_symbol(Angle()) == :θ && force_symbol(Angle()) == :M && force_unit(Angle()) == "N m"
        sp = LinearSpring(100.0)
        E = charges_at(sp; x = 0.02)
        @test forces(sp, E)[1] ≈ 2.0 && internal_energy(sp, E) ≈ 0.5 * 100.0 * 0.02^2
        @test charges_at(sp; F = 2.0)[1] ≈ 0.02
        @test forces(LinearSpring(100.0; x0 = 0.1), charges_at(sp; x = 0.1))[1] == 0.0
        ts = TorsionSpring(5.0)
        @test forces(ts, charges_at(ts; θ = 0.2))[1] ≈ 1.0 && charges_at(ts; M = 1.0)[1] ≈ 0.2
        pp = PendulumPotential(2.0, 9.81, 0.5)
        Eθ = charges_at(pp; θ = 0.3)
        @test forces(pp, Eθ)[1] ≈ 2.0 * 9.81 * 0.5 * sin(0.3)
        @test internal_energy(pp, Eθ) ≈ 2.0 * 9.81 * 0.5 * (1 - cos(0.3))
        @test eos_matrix(pp, Eθ)[1, 1] ≈ 2.0 * 9.81 * 0.5 * cos(0.3)                       # A4 (1 × 1)
        @test_throws ArgumentError LinearSpring(0.0)
        @test_throws ArgumentError TorsionSpring(-1.0)
        @test_throws ArgumentError PendulumPotential(1.0, 0.0, 1.0)
        @test_throws ArgumentError charges_at(sp; x = 0.1, F = 1.0)
        @test SpringBranch() isa Converter && TorsionBranch() isa Converter
        @test_throws ArgumentError QuadraticDrag(-1.0)
        @test_throws ArgumentError QuadraticDrag(1.0; domain = Thermal())
        @test_throws ArgumentError CoulombFriction(1.0; w_ε = 0.0)
        @test_throws ArgumentError CoulombFriction(-1.0; w_ε = 1e-3)
        @test_throws ArgumentError Conductor(QuadraticDrag(1.0))                  # dissipates: needs a policy
        @test_throws ArgumentError Source(Thermal() => 1.0)                       # heat sources: not this element
        @test_throws ArgumentError Source(Momentum() => NaN)
        @test_throws ArgumentError GravityField(Inf)
        # graph validation
        sys = ThermoSystem()
        s = add!(sys, Body(:s, LinearSpring(10.0); x = 0.0))
        r = add!(sys, Reservoir(:r, Momentum() => 0.0))
        bad1 = deepcopy(sys); connect!(bad1, Source(Momentum() => 1.0), s)     # the spring has no momentum
        @test_throws ArgumentError compile(bad1)
        bad2 = deepcopy(sys); connect!(bad2, GravityField(9.81), s)            # no inertial mass
        @test_throws ArgumentError compile(bad2)
        bad3 = deepcopy(sys); connect!(bad3, Source(Momentum() => 1.0), r)     # sources act on bodies
        @test_throws ArgumentError compile(bad3)
    end

    @testset "kernels: spring branches conserve power; friction laws dissipate (D ≥ 0) by type" begin
        for W in (SpringBranch().W, TorsionBranch().W)
            @test W == -W'
            for _ in 1:100
                P = 10 .* (rand(rng, 3) .- 0.5)
                @test abs(dot(P, W * P)) <= 1e-14 * sum(abs2, P)
            end
        end
        for law in (QuadraticDrag(0.7), CoulombFriction(3.0; w_ε = 1e-3), QuadraticDrag(2.0; domain = AngularMomentum()))
            for Δ in vcat(0.0, 1e-12, -1e-12, 10 .^ range(-8, 3; length = 30), -(10 .^ range(-8, 3; length = 30)))
                J = ThermoKinetics.fluxes(law, SVector(Δ))
                @test Δ * J[1] >= 0
                @test ThermoKinetics.fluxes(law, SVector(-Δ))[1] == -J[1]            # odd
            end
        end
    end

    @testset "sources: force and gravity ledgers close energy; momentum is accounted" begin
        sys = msd_system(; m = 2.0, k = 50.0, c = 3.0, x0 = 0.01, v0 = 0.5, F = t -> 2.0 * cos(3t), g = 9.81)
        cs = compile(sys)
        u0 = initial_state(cs)
        r = energy_residual(cs, u0)
        @test abs(r.residual) <= 1e-12 * r.scale
        du = rates(cs, u0)
        @test source_ledger(cs, du, :force).W ≈ 0.5 * 2.0 atol = 1e-12           # P·J = w F at t = 0
        @test source_ledger(cs, du, :gravity).W ≈ 0.5 * (-2.0 * 9.81) rtol = 1e-14
        @test source_ledger(cs, du, :gravity).charge ≈ 2.0 * 9.81 rtol = 1e-14    # the field receives +Mg
        sol = simulate(cs, (0.0, 5.0); reltol = 1e-10)
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-7 && cl.entropy_rel <= 1e-10
    end

    @testset "G3 two-mass chain: a normal mode released at rest stays harmonic" begin
        m, k, a = 1.0, 100.0, 0.01
        for λ in (k * (3 - sqrt(5)) / 2, k * (3 + sqrt(5)) / 2)       # eigenvalues of K = [2k −k; −k k] (M = I)
            φ = 2 - λ / k                                             # mode shape (1, φ)
            sys = ThermoSystem()
            m1 = add!(sys, Body(:m1, InertialMass(m); w = 0.0))
            m2 = add!(sys, Body(:m2, InertialMass(m); w = 0.0))
            s1 = add!(sys, Body(:s1, LinearSpring(k); x = a))
            s2 = add!(sys, Body(:s2, LinearSpring(k); x = (φ - 1) * a))   # extension x2 − x1
            frame = add!(sys, Reservoir(:frame, Momentum() => 0.0))
            connect!(sys, SpringBranch(), frame, m1, s1)
            connect!(sys, SpringBranch(), m1, m2, s2)
            cs = compile(sys)
            ω = sqrt(λ / m)
            sol = simulate(cs, (0.0, 3 * 2π / ω); saveat = 2π / ω / 16, reltol = 1e-11)
            for (t, u) in zip(sol.t, sol.u)
                x1 = charges(cs, u, s1)[1]
                x2 = x1 + charges(cs, u, s2)[1]
                @test x1 ≈ a * cos(ω * t) atol = 1e-7 * a
                @test x2 ≈ φ * a * cos(ω * t) atol = 1e-7 * a
            end
            E0 = total_energy(cs, sol.u[1])
            @test maximum(abs(total_energy(cs, u) - E0) for u in sol.u) <= 1e-8 * E0
        end
    end

    @testset "G4 spring under gravity: half the gravitational work is stored, half dissipated" begin
        m, k, g, T = 2.0, 200.0, 9.81, 300.0
        ωn = sqrt(k / m)
        xs = -m * g / k
        # undamped: oscillation about the static deflection at the unchanged frequency
        cs0 = compile(msd_system(; m, k, c = 0.0, x0 = 0.0, g))
        sol0 = simulate(cs0, (0.0, 3 * 2π / ωn); saveat = 0.05, reltol = 1e-11)
        for (t, u) in zip(sol0.t, sol0.u)
            @test charges(cs0, u, :spring)[1] ≈ xs * (1 - cos(ωn * t)) atol = 1e-7 * abs(xs)
        end
        # damped: comes to rest at x = −Mg/k; of the work M²g²/k, half is dissipated, for any c
        for ζ in (0.05, 0.5)
            c = 2ζ * sqrt(k * m)
            cs = compile(msd_system(; m, k, c, x0 = 0.0, g, T))
            sol = simulate(cs, (0.0, 40 / (ζ * ωn)); reltol = 1e-10)
            uf = sol.u[end]
            @test charges(cs, uf, :spring)[1] ≈ xs rtol = 1e-8
            Wg = source_ledger(cs, uf, :gravity).W
            @test Wg ≈ m^2 * g^2 / k rtol = 1e-7
            Q = reservoir_ledger(cs, uf, :env).U
            @test Q ≈ m^2 * g^2 / (2k) rtol = 1e-6
            @test generated_entropy(cs, uf)[:damper] ≈ Q / T rtol = 1e-6
            cl = closure(cs, sol)
            @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10
            # the steady state is found directly too
            us = steady_state(cs)
            @test charges(cs, us, :spring)[1] ≈ xs rtol = 1e-10
            @test abs(force(cs, us, :mass, Momentum())) <= 1e-12
        end
    end

    @testset "G9 nonlinear friction: quadratic drag and Coulomb friction" begin
        m, v0, T = 2.0, 10.0, 300.0
        # quadratic drag, its heat at a kernel-B node whose one leg goes to the environment
        c = 0.5
        sys = ThermoSystem()
        mass = add!(sys, Body(:mass, InertialMass(m); w = v0))
        air = add!(sys, Reservoir(:air, Momentum() => 0.0))
        env = add!(sys, Reservoir(:env; T))
        drag = connect!(sys, Conductor(QuadraticDrag(c); dissipation = InternalNode(G_h = 1.0, h = env)), mass, air)
        cs = compile(sys)
        sol = simulate(cs, (0.0, 100.0); saveat = 1.0, reltol = 1e-10)
        for (t, u) in zip(sol.t, sol.u)
            @test force(cs, u, mass, Momentum()) ≈ v0 / (1 + c * v0 * t / m) rtol = 1e-7
        end
        vf = force(cs, sol.u[end], mass, Momentum())
        @test reservoir_ledger(cs, sol.u[end], env).U ≈ 0.5m * (v0^2 - vf^2) rtol = 1e-7
        parts = entropy_production_parts(cs, sol.u[3])
        @test parts[drag].dissipation > 0 && parts[drag].conduction >= 0
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10
        # Coulomb friction (tanh-regularized): constant deceleration F_c/m until rest
        Fc, wε = 4.0, 1e-3
        sys2 = ThermoSystem()
        mass2 = add!(sys2, Body(:mass, InertialMass(m); w = v0))
        floor = add!(sys2, Reservoir(:floor, Momentum() => 0.0))
        env2 = add!(sys2, Reservoir(:env; T))
        fr = connect!(sys2, Conductor(CoulombFriction(Fc; w_ε = wε); dissipation = Into(env2)), mass2, floor)
        cs2 = compile(sys2)
        sol2 = simulate(cs2, (0.0, 10.0); saveat = [1.0, 2.0, 3.0, 4.0, 10.0], reltol = 1e-10)
        for (t, u) in zip(sol2.t[1:4], sol2.u[1:4])
            @test force(cs2, u, mass2, Momentum()) ≈ v0 - Fc * t / m rtol = 1e-8
        end
        @test abs(force(cs2, sol2.u[end], mass2, Momentum())) <= 1e-6
        @test reservoir_ledger(cs2, sol2.u[end], env2).U ≈ 0.5m * v0^2 rtol = 1e-8
        @test generated_entropy(cs2, sol2.u[end])[fr] ≈ 0.5m * v0^2 / T rtol = 1e-8
        cl2 = closure(cs2, sol2)
        @test cl2.energy_rel <= 1e-6 && cl2.entropy_rel <= 1e-10
    end

    @testset "A1, A3 on random mechanical chains (springs, dampers, drag, kernel A and B)" begin
        for trial in 1:10
            n = rand(rng, 2:5)
            sys = ThermoSystem()
            frame = add!(sys, Reservoir(:frame, Momentum() => 0.0))
            env = add!(sys, Reservoir(:env; T = 250 + 100rand(rng)))
            masses = [add!(sys, Body(Symbol(:m, i), InertialMass(0.5 + 2rand(rng)); w = randn(rng))) for i in 1:n]
            dampers = Symbol[]
            for i in 1:n
                left = i == 1 ? frame : masses[i-1]
                s = add!(sys, Body(Symbol(:s, i), LinearSpring(10 + 90rand(rng)); x = 0.02randn(rng)))
                connect!(sys, SpringBranch(), left, masses[i], s)
                law = rand(rng) < 0.7 ? LinearConductance(Momentum(), 0.1 + 2rand(rng)) : QuadraticDrag(0.1 + rand(rng))
                pol = rand(rng) < 0.5 ? Into(env) : InternalNode(G_h = 0.1 + rand(rng), h = env)
                push!(dampers, connect!(sys, Conductor(law; dissipation = pol), masses[i], left))
            end
            cs = compile(sys)
            for _ in 1:10                                # random states: energy closes, production ≥ 0
                u = initial_state(cs)
                u[1:cs.nphys] .= randn(rng, cs.nphys)
                r = energy_residual(cs, u)
                @test abs(r.residual) <= 1e-12 * r.scale
                @test all(≥(0), values(entropy_production(cs, u)))
            end
            sol = simulate(cs, (0.0, 20.0); reltol = 1e-9)
            cl = closure(cs, sol)
            @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10
        end
    end

    @testset "pendulum tolerance scaling: the potential's stiffness floor (A = 0 at θ = π/2)" begin
        m, ℓ, g = 1.0, 1.0, 9.81
        cs = compile(pendulum_system(; m, ℓ, g, θ0 = π / 2))
        tol = ThermoKinetics.abstol_vector(cs)
        @test tol[2] ≈ 1e-8 / (m * g * ℓ) rtol = 1e-12          # scaled by mgℓ, not by |A| ≈ 0
        Tp = 4 * sqrt(ℓ / g) * ellipK(sin(π / 4))
        sol = simulate(cs, (0.0, Tp); saveat = [Tp / 2, Tp], reltol = 1e-11)
        @test charges(cs, sol.u[1], :pot)[1] ≈ -π / 2 atol = 1e-7
        @test charges(cs, sol.u[2], :pot)[1] ≈ π / 2 atol = 1e-7
    end

    @testset "rack and pinion: a transformer between angular and linear momentum" begin
        J, m, R, τ = 0.02, 3.0, 0.05, 0.4
        sys = ThermoSystem()
        pinion = add!(sys, Body(:pinion, RotationalInertia(J); ω = 0.0))
        rack = add!(sys, Body(:rack, InertialMass(m); w = 0.0))
        tr = connect!(sys, Transformer(R; a = Momentum(), b = AngularMomentum()), rack, pinion)   # v = R ω
        connect!(sys, Source(AngularMomentum() => τ), pinion; name = :motor)
        cs = compile(sys)
        sol = simulate(cs, (0.0, 2.0); saveat = 0.25)
        for (t, u) in zip(sol.t, sol.u)
            ω = force(cs, u, pinion, AngularMomentum())
            @test ω ≈ τ * t / (J + m * R^2) rtol = 1e-7 atol = 1e-12
            @test force(cs, u, rack, Momentum()) ≈ R * ω rtol = 1e-7 atol = 1e-12
            @test abs(constraint_residual(cs, u, tr)) <= 1e-8
        end
        @test source_ledger(cs, sol.u[end], :motor).W ≈ total_energy(cs, sol.u[end]) rtol = 1e-6
    end

    @testset "DC motor driving a load through a gearbox" begin
        V, R, Lw, k, J, b = 12.0, 0.5, 1e-3, 0.05, 1e-4, 1e-5
        r, JL, bL = 5.0, 1e-3, 2e-4
        sys = ThermoSystem()
        src = add!(sys, Reservoir(:src, Electrical() => V))
        gnd = add!(sys, Reservoir(:gnd, Electrical() => 0.0))
        env = add!(sys, Reservoir(:env; T = 300.0))
        frame = add!(sys, Reservoir(:frame, AngularMomentum() => 0.0))
        mid = add!(sys, Junction(:mid, Electrical()))
        coil = add!(sys, Body(:coil, LinearInductor(Lw); i = 0.0))
        rotor = add!(sys, Body(:rotor, RotationalInertia(J); ω = 0.0))
        load = add!(sys, Body(:load, RotationalInertia(JL); ω = 0.0))
        connect!(sys, Conductor(Resistor(R); dissipation = Into(env)), src, mid)
        W = [0 0 -1 0; 0 0 1 0; 1 -1 0 -k; 0 0 k 0]
        connect!(sys, Converter(W, (Electrical(), Electrical(), FluxLinkage(), AngularMomentum())), mid, gnd, coil, rotor)
        connect!(sys, Conductor(LinearConductance(AngularMomentum(), b); dissipation = Into(env)), rotor, frame)
        connect!(sys, Transformer(r; a = AngularMomentum(), b = AngularMomentum()), rotor, load)    # ω_m = r ω_L
        connect!(sys, Conductor(LinearConductance(AngularMomentum(), bL); dissipation = Into(env)), load, frame)
        cs = compile(sys)
        sol = simulate(cs, (0.0, 5.0))
        uf = sol.u[end]
        ωm = V / (R * (b + bL / r^2) / k + k)                 # V = iR + kω_m, k i = (b + b_L/r²) ω_m
        @test force(cs, uf, rotor, AngularMomentum()) ≈ ωm rtol = 1e-6
        @test force(cs, uf, load, AngularMomentum()) ≈ ωm / r rtol = 1e-6
        @test force(cs, uf, coil, FluxLinkage()) ≈ (b + bL / r^2) * ωm / k rtol = 1e-6
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9
    end
end
