# Application tests (DESIGN §0.1): models a user would actually build, checked against an
# independent engineering result. One per delivered use; rows are added as slices land.

@testset "applications" begin
    @testset "transient heat transfer: composite wall" begin
        # brick + insulation, 1 m², convection on both faces (SI)
        A = 1.0
        T_in, h_in = 293.15, 8.0                          # room air
        T_out, h_out = 263.15, 25.0                       # winter outside air
        λ1, ρ1, c1, L1, n1 = 0.7, 1800.0, 840.0, 0.20, 20   # brick
        λ2, ρ2, c2, L2, n2 = 0.04, 30.0, 1400.0, 0.10, 10   # mineral-wool insulation
        T0 = 273.15

        sys = ThermoSystem()
        inside = add!(sys, Reservoir(:inside; T = T_in))
        outside = add!(sys, Reservoir(:outside; T = T_out))
        brick = rod!(sys, :brick; n = n1, L = L1, area = A, λ = λ1, ρ = ρ1, c = c1, T0, left = (inside, h_in))
        insul = rod!(sys, :insul; n = n2, L = L2, area = A, λ = λ2, ρ = ρ2, c = c2, T0, right = (outside, h_out))
        # interface: the two half cells in series
        G_if = 1 / (brick.dx / (2λ1 * A) + insul.dx / (2λ2 * A))
        iface = connect!(sys, Conductor(ThermalConductance(G_if)), brick.cells[end], insul.cells[1]; name = :interface)
        faces = vcat(brick.all_faces, [iface], insul.all_faces)
        cs = compile(sys)

        # independent result: series thermal resistances, U-value and piecewise-linear profile
        R = 1 / h_in + L1 / λ1 + L2 / λ2 + 1 / h_out      # m² K / W
        q = (T_in - T_out) / (R / A)                      # W
        Ts1 = T_in - q / (h_in * A)                       # inner surface
        Tif = Ts1 - q * L1 / (λ1 * A)                     # brick/insulation interface
        T_exact(x) = x <= L1 ? Ts1 - q * x / (λ1 * A) : Tif - q * (x - L1) / (λ2 * A)
        xs = vcat(brick.x, L1 .+ insul.x)
        cells = vcat(brick.cells, insul.cells)

        @testset "steady state" begin
            u = steady_state(cs)
            for f in faces
                @test heat_flow(cs, u, f) ≈ q rtol = 1e-9
            end
            for (cell, x) in zip(cells, xs)
                @test temperature(cs, u, cell) ≈ T_exact(x) atol = 1e-8
            end
            # loss attribution: total production is q(1/T_out − 1/T_in), and the insulation
            # (largest resistance) is where most of it happens
            Sgen = entropy_production(cs, u)
            @test sum(values(Sgen)) ≈ q * (1 / T_out - 1 / T_in) rtol = 1e-9
            S_insul = sum(Sgen[f] for f in insul.faces)
            S_brick = sum(Sgen[f] for f in brick.faces)
            @test S_insul > S_brick
        end

        @testset "warm-up transient" begin
            τ = (L2 / λ2) * ρ1 * c1 * L1                  # insulation resistance × brick capacity
            sol = simulate(cs, (0.0, 20τ); saveat = τ)
            uf = sol.u[end]
            for (cell, x) in zip(cells, xs)
                @test temperature(cs, uf, cell) ≈ T_exact(x) atol = 1e-4
            end
            # the wall starts colder than its steady profile, so over the warm-up the room
            # supplies more heat than the outside receives; the difference is the stored energy
            led_in = reservoir_ledger(cs, uf, inside)
            led_out = reservoir_ledger(cs, uf, outside)
            stored = total_energy(cs, uf) - total_energy(cs, initial_state(cs))
            @test -led_in.U - led_out.U ≈ stored rtol = 1e-6
            @test stored > 0
            cl = closure(cs, sol)
            @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10
        end
    end
end

@testset "applications: control mass" begin
    N2 = gas_n2()
    cp = N2.cv + N2.R
    T0, p_atm, V0, Th = 300.0, 101325.0, 1e-3, 400.0
    m = p_atm * V0 / (N2.R * T0)

    @testset "isochoric heating of a rigid vessel" begin
        sys = ThermoSystem()
        g = add!(sys, Body(:gas, N2; T = T0, V = V0, m))
        heater = add!(sys, Reservoir(:heater; T = Th))
        connect!(sys, Conductor(ThermalConductance(0.2)), heater, g)
        cs = compile(sys)
        sol = simulate(cs, (0.0, 200.0))
        uf = sol.u[end]
        Q = -reservoir_ledger(cs, uf, heater).U
        @test Q ≈ m * N2.cv * (Th - T0) rtol = 1e-7          # Q = m c_v ΔT
        @test pressure(cs, uf, g) ≈ p_atm * Th / T0 rtol = 1e-7
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10
    end

    @testset "isobaric heating behind a loose piston" begin
        # a gas in a cylinder, a damped piston open to the atmosphere (seal friction heats the
        # gas through an internal node), and a heater. First law for the control mass:
        # Q = ΔU + p ΔV = m c_p ΔT, regardless of friction; W = p_atm ΔV goes to the atmosphere.
        sys = ThermoSystem()
        g = add!(sys, Body(:gas, N2; T = T0, V = V0, m))
        atm = add!(sys, Reservoir(:atm; p = p_atm))
        heater = add!(sys, Reservoir(:heater; T = Th))
        piston = connect!(sys, Conductor(VolumeConductance(1e-8); dissipation = InternalNode(G_a = 1.0)), g, atm)
        connect!(sys, Conductor(ThermalConductance(0.2)), heater, g)
        cs = compile(sys)
        sol = simulate(cs, (0.0, 200.0))
        uf = sol.u[end]
        @test temperature(cs, uf, g) ≈ Th rtol = 1e-7
        @test pressure(cs, uf, g) ≈ p_atm rtol = 1e-7
        Q = -reservoir_ledger(cs, uf, heater).U
        W = reservoir_ledger(cs, uf, atm).U
        ΔV = charges(cs, uf, g)[2] - V0
        @test Q ≈ m * cp * (Th - T0) rtol = 1e-6
        @test W ≈ p_atm * ΔV rtol = 1e-6
        @test W ≈ m * N2.R * (Th - T0) rtol = 1e-6
        # second law for the whole: gas entropy rise minus what the heater gave up
        Sgen = sum(values(generated_entropy(cs, uf)))
        @test Sgen ≈ m * cp * log(Th / T0) - Q / Th rtol = 1e-6
        @test Sgen > 0
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10 && cl.charge_rel <= 1e-10
    end
end

@testset "applications: control mass, piston dynamics" begin
    @testset "piston–cylinder compression by the atmosphere, with seal friction" begin
        # a gas in a cylinder (1 bar) is compressed by a 2 kg piston that the atmosphere (3 bar)
        # pushes in; seal friction damps the oscillation and sends its heat to the wall. The gas
        # takes no heat and the piston is ideal, so the gas ends exactly on its isentrope at 3 bar,
        # whatever the dynamics. All the irreversibility is the friction heat.
        N2 = gas_n2()
        cp = N2.cv + N2.R
        A, mp, c, T0, p0, V0, patm, Twall = 0.01, 2.0, 20.0, 300.0, 1e5, 1e-3, 3e5, 300.0
        sys = ThermoSystem()
        gas = add!(sys, Body(:gas, N2; T = T0, p = p0, V = V0))
        pst = add!(sys, Body(:piston, InertialMass(mp); w = 0.0))
        atm = add!(sys, Reservoir(:atm; p = patm))
        ground = add!(sys, Reservoir(:ground, Momentum() => 0.0))
        wall = add!(sys, Reservoir(:wall; T = Twall))
        connect!(sys, Gyrator(A; a = Volumetric(), b = Momentum()), gas, pst)
        connect!(sys, Gyrator(-A; a = Volumetric(), b = Momentum()), atm, pst)
        fr = connect!(sys, Conductor(LinearConductance(Momentum(), c); dissipation = Into(wall)), pst, ground)
        cs = compile(sys)
        u0 = initial_state(cs)
        sol = simulate(cs, (0.0, 10.0))
        uf = sol.u[end]
        Tf = T0 * (patm / p0)^(N2.R / cp)
        @test pressure(cs, uf, gas) ≈ patm rtol = 1e-6
        @test temperature(cs, uf, gas) ≈ Tf rtol = 1e-6
        S0 = charges(cs, u0, gas)[1]                                       # no heat, ideal piston: S fixed
        @test abs(charges(cs, uf, gas)[1] - S0) <= 1e-12 * abs(S0)         # (to the ulps of the implicit steps)
        m = charges(cs, u0, gas)[3]
        ΔU = m * N2.cv * (Tf - T0)
        Watm = -reservoir_ledger(cs, uf, atm).U                            # work done by the atmosphere
        @test Watm ≈ patm * (V0 - charges(cs, uf, gas)[2]) rtol = 1e-6
        Qfric = reservoir_ledger(cs, uf, wall).U
        @test Qfric ≈ Watm - ΔU rtol = 1e-5                                # first law, piston at rest
        @test generated_entropy(cs, uf)[fr] ≈ Qfric / Twall rtol = 1e-6
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10
    end
end

@testset "applications: control volume" begin
    N2 = gas_n2()
    cp = N2.cv + N2.R

    @testset "adiabatic tank blowdown: the gas left behind expands isentropically" begin
        T0, p0, V, p_atm = 400.0, 5e5, 1e-2, 1e5
        sys = ThermoSystem()
        tank = add!(sys, Body(:tank, N2; T = T0, p = p0, V))
        atm = add!(sys, Reservoir(:atm; T = 300.0, p = p_atm, material = N2))
        valve = connect!(sys, Stream(HydraulicConductance(2e-8)), tank, atm)
        cs = compile(sys)
        sol = simulate(cs, (0.0, 300.0); saveat = 5.0)
        for u in sol.u
            p = pressure(cs, u, tank)
            @test temperature(cs, u, tank) ≈ T0 * (p / p0)^(N2.R / cp) rtol = 1e-7
        end
        @test pressure(cs, sol.u[end], tank) ≈ p_atm rtol = 1e-6
        # all the irreversibility (throttling + mixing) is booked where the jet ends: at the atmosphere
        @test generated_entropy(cs, sol.u[end])[valve] > 0
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9 && cl.charge_rel <= 1e-10
    end

    @testset "counter-flow heat exchanger vs ε–NTU (first-order convergence)" begin
        Th_in, Tc_in = 400.0, 300.0
        mh, mc = 1.0e-3, 1.5e-3
        Ch, Cc = mh * cp, mc * cp
        Cmin, Cr = min(Ch, Cc), min(Ch, Cc) / max(Ch, Cc)
        NTU = 2.0
        UA = NTU * Cmin
        ε_exact = (1 - exp(-NTU * (1 - Cr))) / (1 - Cr * exp(-NTU * (1 - Cr)))
        function exchanger(N)
            sys = ThermoSystem()
            hot = [add!(sys, Body(Symbol(:h, i), N2; T = Th_in, p = 1.05e5, V = 1e-4)) for i in 1:N]
            cold = [add!(sys, Body(Symbol(:c, i), N2; T = Tc_in, p = 1.05e5, V = 1e-4)) for i in 1:N]
            connect!(sys, FlowBoundary(N2; T = Th_in, p = 1.2e5, mdot = t -> mh, supply = :passive), hot[1])
            connect!(sys, FlowBoundary(N2; T = Tc_in, p = 1.2e5, mdot = t -> mc, supply = :passive), cold[1])
            for i in 1:N-1
                connect!(sys, Stream(HydraulicConductance(1e-5)), hot[i], hot[i+1])
                connect!(sys, Stream(HydraulicConductance(1e-5)), cold[i], cold[i+1])
            end
            connect!(sys, PassThroughOutlet(κ = 1e-6, p_out = 1e5), hot[N])
            connect!(sys, PassThroughOutlet(κ = 1e-6, p_out = 1e5), cold[N])
            for i in 1:N                                 # counter-flow: hot cell i faces cold cell N + 1 − i
                connect!(sys, Conductor(ThermalConductance(UA / N)), hot[i], cold[N+1-i])
            end
            cs = compile(sys)
            u = steady_state(cs; u0 = simulate(cs, (0.0, 50.0)).u[end])     # transient warm start
            Th_out, Tc_out = temperature(cs, u, hot[N]), temperature(cs, u, cold[N])
            @test Ch * (Th_in - Th_out) ≈ Cc * (Tc_out - Tc_in) rtol = 1e-8  # energy balance
            # loss attribution: total production = Σ ṁ (s_out − s_in) of the two streams
            p_h, p_c = pressure(cs, u, hot[N]), pressure(cs, u, cold[N])
            S_expected = mh * (cp * log(Th_out / Th_in) - N2.R * log(p_h / 1.2e5)) +
                         mc * (cp * log(Tc_out / Tc_in) - N2.R * log(p_c / 1.2e5))
            @test sum(values(entropy_production(cs, u))) ≈ S_expected rtol = 1e-7
            return Ch * (Th_in - Th_out) / (Cmin * (Th_in - Tc_in))
        end
        ε10, ε20 = exchanger(10), exchanger(20)
        # well-mixed cells converge at first order: the error halves when N doubles …
        @test 1.8 < (ε_exact - ε10) / (ε_exact - ε20) < 2.2
        # … so the Richardson extrapolation recovers ε–NTU (measured: 1.1e-3)
        @test abs((2ε20 - ε10) - ε_exact) < 2e-3
    end
end

@testset "applications: electrical circuits" begin
    @testset "resistor self-heating with its heat path" begin
        # a 12 V resistor whose core (internal node) is bonded to a heat-sink plate, which is
        # cooled by air at 300 K. Steady temperatures and the loss map follow from the series path.
        V, R, Gcore, Gair, Tair = 12.0, 8.0, 0.5, 0.2, 300.0
        P = V^2 / R
        sys = ThermoSystem()
        hi = add!(sys, Reservoir(:hi, Electrical() => V))
        lo = add!(sys, Reservoir(:lo, Electrical() => 0.0))
        air = add!(sys, Reservoir(:air; T = Tair))
        sink = add!(sys, Body(:sink, ConstantCapacitySolid(0.05, 900.0); T = Tair))
        res = connect!(sys, Conductor(Resistor(R); dissipation = InternalNode(G_h = Gcore, h = sink)), hi, lo)
        cool = connect!(sys, Conductor(ThermalConductance(Gair)), sink, air)
        cs = compile(sys)
        u = steady_state(cs)
        T_sink = Tair + P / Gair
        T_core = T_sink + P / Gcore
        @test temperature(cs, u, sink) ≈ T_sink rtol = 1e-10
        @test flux(cs, u, res).T_m ≈ T_core rtol = 1e-10
        # where the entropy is produced: in the core, across the bond, at the sink surface
        Sp = entropy_production_parts(cs, u)
        @test Sp[res].dissipation ≈ P / T_core rtol = 1e-10
        @test Sp[res].conduction ≈ P * (1 / T_sink - 1 / T_core) rtol = 1e-9
        @test Sp[cool].total ≈ P * (1 / Tair - 1 / T_sink) rtol = 1e-9
        @test P ≈ Tair * sum(values(entropy_production(cs, u))) rtol = 1e-9    # all of P is lost work
        # warm-up of the sink plate from ambient: first order with τ = C/G_air
        τ = 0.05 * 900.0 / Gair
        sol = simulate(cs, (0.0, 3τ); saveat = τ)
        for (t, uu) in zip(sol.t, sol.u)
            @test temperature(cs, uu, sink) ≈ Tair + (P / Gair) * (1 - exp(-t / τ)) rtol = 1e-7
        end
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10 && cl.charge_rel <= 1e-10
    end

    @testset "thermoelectric generator efficiency" begin
        α, R, K = 0.05, 1.0, 0.5                                # a module: Z = α²/(RK)
        Th, Tc = 500.0, 300.0
        Z = α^2 / (R * K)
        function generator(Rload)
            sys = ThermoSystem()
            H = add!(sys, Reservoir(:H, Thermal() => Th, Electrical() => 0.0))
            cold = add!(sys, Reservoir(:cold; T = Tc))
            sink = add!(sys, Reservoir(:sink; T = Tc))
            body = Composite(ConstantCapacitySolid(1.0, 400.0), LinearCapacitor(1e-3))
            C = add!(sys, Body(:C, body, charges_at(body, (T = Tc,), (φ = 0.0,))))
            connect!(sys, Conductor(ThermoelectricConductance(; R, K, α)), H, C)
            connect!(sys, Conductor(ThermalConductance(1e7)), C, cold)              # cold electrode on the sink plate
            connect!(sys, Conductor(Resistor(Rload); dissipation = Into(sink)), C, H)  # the load closes the circuit
            cs = compile(sys)
            u = steady_state(cs)
            du = rates(cs, u)
            Qh = -reservoir_ledger(cs, du, H).U                 # heat drawn from the hot side
            Pload = reservoir_ledger(cs, du, sink).U            # electrical power dissipated in the load
            return Qh, Pload, temperature(cs, u, C)
        end
        for m in (0.5, 1.0, 2.0)
            Qh, Pload, TC = generator(m * R)
            ΔT = Th - TC
            I = α * ΔT / (R + m * R)
            @test Pload ≈ I^2 * m * R rtol = 1e-9
            @test Qh ≈ K * ΔT + α * Th * I - 0.5 * I^2 * R rtol = 1e-9        # standard TEG heat balance
        end
        # the optimal load R_L = R√(1 + Z T̄) gives the classical maximum efficiency
        Qh, Pload, TC = generator(R * sqrt(1 + Z * (Th + Tc) / 2))
        Tbar = (Th + TC) / 2
        M = sqrt(1 + Z * Tbar)
        η_max = ((Th - TC) / Th) * (M - 1) / (M + TC / Th)
        @test Pload / Qh ≈ η_max rtol = 1e-6
        @test Pload / Qh < 1 - Tc / Th                                         # below Carnot
    end
end

@testset "applications: translational mechanics" begin
    @testset "G1 free mass–spring–damper: response, and all the energy ends up as heat" begin
        m, k, x0, T = 1.0, 100.0, 0.01, 300.0
        ωn = sqrt(k / m)
        # undamped: x0 cos ωn t, mechanical energy conserved
        cs0 = compile(msd_system(; m, k, c = 0.0, x0))
        sol0 = simulate(cs0, (0.0, 3 * 2π / ωn); saveat = 0.02, reltol = 1e-11)
        for (t, u) in zip(sol0.t, sol0.u)
            @test charges(cs0, u, :spring)[1] ≈ x0 * cos(ωn * t) atol = 1e-7 * x0
        end
        E0 = total_energy(cs0, sol0.u[1])
        @test maximum(abs(total_energy(cs0, u) - E0) for u in sol0.u) <= 1e-8 * E0
        # under-, critically and over-damped
        for ζ in (0.1, 1.0, 2.0)
            c = 2ζ * sqrt(k * m)
            σ = ζ <= 1 ? ζ * ωn : ωn * (ζ - sqrt(ζ^2 - 1))          # slowest decay rate
            cs = compile(msd_system(; m, k, c, x0, T))
            sol = simulate(cs, (0.0, 40 / σ); saveat = vcat(range(0, 8 / σ; length = 60), 40 / σ), reltol = 1e-10)
            for (t, u) in zip(sol.t, sol.u)
                @test charges(cs, u, :spring)[1] ≈ msd_x(t; ζ, ωn, x0) atol = 1e-6 * x0
            end
            uf = sol.u[end]
            Q = reservoir_ledger(cs, uf, :env).U
            @test Q ≈ 0.5k * x0^2 rtol = 1e-6                          # for any c (the mechanical D1)
            @test generated_entropy(cs, uf)[:damper] ≈ Q / T rtol = 1e-6
            cl = closure(cs, sol)
            @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10
        end
    end

    @testset "G2 forced oscillator: amplitude, phase, and power balance" begin
        m, k, ζ, F0, Ω, T = 1.0, 100.0, 0.1, 1.0, 8.0, 300.0
        c = 2ζ * sqrt(k * m)
        cs = compile(msd_system(; m, k, c, x0 = 0.0, T, F = t -> F0 * sin(Ω * t)))
        Tp = 2π / Ω
        t0 = ceil(40 / (ζ * sqrt(k / m)) / Tp) * Tp                   # transients have decayed by e⁻⁴⁰
        N = 5
        ts = collect(range(t0, t0 + N * Tp; length = 8N + 1))
        sol = simulate(cs, (0.0, ts[end]); saveat = ts, reltol = 1e-10)
        X = F0 / hypot(k - m * Ω^2, c * Ω)
        φ = atan(c * Ω, k - m * Ω^2)
        for (t, u) in zip(sol.t, sol.u)
            @test charges(cs, u, :spring)[1] ≈ X * sin(Ω * t - φ) atol = 1e-6 * X
        end
        Win = source_ledger(cs, sol.u[end], :force).W - source_ledger(cs, sol.u[1], :force).W
        Q = reservoir_ledger(cs, sol.u[end], :env).U - reservoir_ledger(cs, sol.u[1], :env).U
        @test Win ≈ 0.5c * Ω^2 * X^2 * N * Tp rtol = 1e-6               # mean input power = mean dissipation
        @test Q ≈ Win rtol = 1e-6
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10
    end

    @testset "G5 pendulum: period vs the elliptic integral, small to large amplitude; damped heat" begin
        m, ℓ, g = 1.0, 1.0, 9.81
        T0 = 2π * sqrt(ℓ / g)
        for θ0 in (0.01, 1.0, 2.5)
            Tp = 4 * sqrt(ℓ / g) * ellipK(sin(θ0 / 2))
            θ0 == 0.01 && @test Tp ≈ T0 * (1 + θ0^2 / 16) rtol = 1e-9     # the small-amplitude limit
            cs = compile(pendulum_system(; m, ℓ, g, θ0))
            sol = simulate(cs, (0.0, Tp); saveat = [Tp / 4, Tp / 2, Tp], reltol = 1e-11)
            θ = [charges(cs, u, :pot)[1] for u in sol.u]
            @test abs(θ[1]) <= 1e-7 * θ0
            @test θ[2] ≈ -θ0 atol = 1e-7 * θ0
            @test θ[3] ≈ θ0 atol = 1e-7 * θ0
            @test abs(force(cs, sol.u[3], :rod, AngularMomentum())) <= 1e-7 * θ0 * sqrt(g / ℓ)
            E0 = m * g * ℓ * (1 - cos(θ0))
            @test maximum(abs(total_energy(cs, u) - E0) for u in sol.u) <= 1e-8 * E0
        end
        θ0, b, T = 1.0, 0.5, 300.0
        cs = compile(pendulum_system(; m, ℓ, g, θ0, b, T))
        sol = simulate(cs, (0.0, 160.0); reltol = 1e-10)
        Q = reservoir_ledger(cs, sol.u[end], :env).U
        @test Q ≈ m * g * ℓ * (1 - cos(θ0)) rtol = 1e-6
        @test generated_entropy(cs, sol.u[end])[:pivot] ≈ Q / T rtol = 1e-6
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10
    end
end

@testset "applications: rotational mechanics" begin
    @testset "G6 gear train: spin-up sees the reflected inertia J₁ + J₂/r₁² + J₃/(r₁r₂)²" begin
        J1, J2, J3, r1, r2, τ = 0.01, 0.05, 0.2, 3.0, 2.0, 0.5
        sys = ThermoSystem()
        s1 = add!(sys, Body(:s1, RotationalInertia(J1); ω = 0.0))
        s2 = add!(sys, Body(:s2, RotationalInertia(J2); ω = 0.0))
        s3 = add!(sys, Body(:s3, RotationalInertia(J3); ω = 0.0))
        g1 = connect!(sys, Transformer(r1; a = AngularMomentum(), b = AngularMomentum()), s1, s2)   # ω₁ = r₁ω₂
        g2 = connect!(sys, Transformer(r2; a = AngularMomentum(), b = AngularMomentum()), s2, s3)
        connect!(sys, Source(AngularMomentum() => τ), s1; name = :motor)
        cs = compile(sys)
        Jeq = J1 + J2 / r1^2 + J3 / (r1 * r2)^2
        sol = simulate(cs, (0.0, 2.0); saveat = 0.25)
        for (t, u) in zip(sol.t, sol.u)
            ω1 = force(cs, u, s1, AngularMomentum())
            @test ω1 ≈ τ * t / Jeq rtol = 1e-7 atol = 1e-12
            @test force(cs, u, s2, AngularMomentum()) ≈ ω1 / r1 rtol = 1e-7 atol = 1e-12
            @test force(cs, u, s3, AngularMomentum()) ≈ ω1 / (r1 * r2) rtol = 1e-7 atol = 1e-12
            @test abs(constraint_residual(cs, u, g1)) <= 1e-8 && abs(constraint_residual(cs, u, g2)) <= 1e-8
        end
        uf = sol.u[end]
        @test source_ledger(cs, uf, :motor).W ≈ total_energy(cs, uf) rtol = 1e-6       # work in = kinetic energy
        @test total_energy(cs, uf) ≈ 0.5 * Jeq * (τ * 2.0 / Jeq)^2 rtol = 1e-6
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6
    end

    @testset "G7 two inertias on a torsion shaft: ω_n = √(k(J₁+J₂)/(J₁J₂))" begin
        J1, J2, k, ω0 = 0.1, 0.3, 50.0, 10.0
        sys = ThermoSystem()
        r1 = add!(sys, Body(:r1, RotationalInertia(J1); ω = ω0))
        r2 = add!(sys, Body(:r2, RotationalInertia(J2); ω = 0.0))
        shaft = add!(sys, Body(:shaft, TorsionSpring(k); θ = 0.0))
        connect!(sys, TorsionBranch(), r1, r2, shaft)              # dθ = ω₂ − ω₁
        cs = compile(sys)
        ωn = sqrt(k * (J1 + J2) / (J1 * J2))
        sol = simulate(cs, (0.0, 3 * 2π / ωn); saveat = 2π / ωn / 16, reltol = 1e-11)
        for (t, u) in zip(sol.t, sol.u)
            @test charges(cs, u, shaft)[1] ≈ -(ω0 / ωn) * sin(ωn * t) atol = 1e-7 * ω0 / ωn
            L = charges(cs, u, r1)[1] + charges(cs, u, r2)[1]
            @test L ≈ J1 * ω0 rtol = 1e-12                          # angular momentum is conserved
        end
        E0 = total_energy(cs, sol.u[1])
        @test maximum(abs(total_energy(cs, u) - E0) for u in sol.u) <= 1e-8 * E0
    end

    @testset "G8 clutch engagement: the loss ½J₁J₂Δω²/(J₁+J₂) for any clutch conductance" begin
        J1, J2, ω1, ω2, t1, T = 0.2, 0.5, 100.0, 20.0, 1.0, 300.0
        Jr = J1 * J2 / (J1 + J2)
        for G in (0.05, 5.0)
            sys = ThermoSystem()
            a = add!(sys, Body(:a, RotationalInertia(J1); ω = ω1))
            b = add!(sys, Body(:b, RotationalInertia(J2); ω = ω2))
            env = add!(sys, Reservoir(:env; T))
            cl_ = connect!(sys, Switched(Conductor(LinearConductance(AngularMomentum(), G); dissipation = Into(env));
                                         on = t -> t >= t1), a, b; name = :clutch)
            cs = compile(sys)
            tend = t1 + 40 * Jr / G
            sol = simulate(cs, (0.0, tend); saveat = [t1, tend], tstops = [t1])
            @test force(cs, sol.u[1], a, AngularMomentum()) == ω1                 # open until t₁
            ωf = (J1 * ω1 + J2 * ω2) / (J1 + J2)
            uf = sol.u[end]
            @test force(cs, uf, a, AngularMomentum()) ≈ ωf rtol = 1e-8
            @test force(cs, uf, b, AngularMomentum()) ≈ ωf rtol = 1e-8
            Q = reservoir_ledger(cs, uf, env).U
            @test Q ≈ 0.5Jr * (ω1 - ω2)^2 rtol = 1e-7
            @test generated_entropy(cs, uf)[cl_] ≈ Q / T rtol = 1e-7
            cl = closure(cs, sol)
            @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10 && cl.charge_rel <= 1e-10
        end
    end
end

@testset "applications: internal heating" begin
    @testset "slab with uniform internal heating, convective faces: T(x) = T∞ + qL/(2h) + qx(L − x)/(2λ)" begin
        L, A, λ, ρ, c, h, Tinf, q = 0.05, 1.0, 15.0, 8000.0, 500.0, 200.0, 300.0, 2e6    # q in W/m³
        T_exact(x) = Tinf + q * L / (2h) + q * x * (L - x) / (2λ)
        rise = T_exact(L / 2) - Tinf
        err = Float64[]
        for n in (10, 20, 40)
            sys = ThermoSystem()
            add!(sys, Reservoir(:coolant; T = Tinf))
            rod = rod!(sys, :slab; n, L, area = A, λ, ρ, c, T0 = Tinf, left = (:coolant, h), right = (:coolant, h))
            for cell in rod.cells
                connect!(sys, HeatSource(q * A * rod.dx), cell)
            end
            cs = compile(sys)
            u = steady_state(cs)
            # the faces carry all the heat, q L A, half to each side
            @test -heat_flow(cs, u, rod.left_face) ≈ q * L * A / 2 rtol = 1e-9
            @test heat_flow(cs, u, rod.right_face) ≈ q * L * A / 2 rtol = 1e-9
            push!(err, maximum(abs(temperature(cs, u, rod.cells[i]) - T_exact(rod.x[i])) for i in 1:n))
        end
        @test err[3] <= 1e-3 * rise
        @test 3.5 <= err[1] / err[2] <= 4.5 && 3.5 <= err[2] / err[3] <= 4.5          # second order
    end
end

@testset "applications: chemical and nuclear" begin
    @testset "water–gas shift equilibrium at 1000 K vs Gibbs minimization" begin
        CO = gasdata(:CO, (C = 1, O = 1); cp = 29.142, h_f = -110527.0, s0 = 197.660)
        CO2 = gasdata(:CO2, (C = 1, O = 2); cp = 37.135, h_f = -393522.0, s0 = 213.795)
        H2O, H2 = H2Og(), H2g()
        mix = IdealGasMixture(CO, H2O, CO2, H2)
        T, V = 1000.0, 1e-3
        n0 = (CO = 0.010, H2O = 0.015, CO2 = 1e-12, H2 = 1e-12)
        sys = ThermoSystem()
        add!(sys, Body(:reactor, mix; T, V, n = n0))
        add!(sys, Reservoir(:furnace; T))
        connect!(sys, Conductor(ThermalConductance(0.5)), :reactor, :furnace)
        connect!(sys, Reaction(:CO => -1, :H2O => -1, :CO2 => 1, :H2 => 1; law = MassAction(5.0)), :reactor; name = :shift)
        cs = compile(sys)
        sol = simulate(cs, (0.0, 300.0))
        nf = moles(cs, sol.u[end], :reactor)
        gs = MolarGas.((CO, H2O, CO2, H2))
        n_of(ξ) = (n0.CO - ξ, n0.H2O - ξ, n0.CO2 + ξ, n0.H2 + ξ)
        G(ξ) = (n = n_of(ξ); sum(n[k] * (mh_(gs[k], T) - T * ms_(gs[k], T, n[k] * R̄ * T / V)) for k in 1:4))
        ξeq = golden_max(ξ -> -G(ξ), 0.0, n0.CO - 1e-18)
        @test nf.CO2 ≈ n_of(ξeq)[3] rtol = 1e-6
        @test nf.CO ≈ n_of(ξeq)[1] rtol = 1e-6
        # the equilibrium constant from the same data, K = exp(−ΔG°/RT) (Δn = 0)
        ΔG = sum(ν * (mh_(g, T) - T * ms_(g, T, 101325.0)) for (ν, g) in zip((-1, -1, 1, 1), gs))
        @test (nf.CO2 * nf.H2) / (nf.CO * nf.H2O) ≈ exp(-ΔG / (R̄ * T)) rtol = 1e-6
        # the heat of reaction went to the furnace
        @test reservoir_ledger(cs, sol.u[end], :furnace).U ≈
              sum(n0[k] * mu_(gs[k], T) for k in 1:4) - sum(nf[k] * mu_(gs[k], T) for k in 1:4) rtol = 1e-6
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9 && cl.charge_rel <= 1e-12
    end

    λT = log(2) / (12.32YEAR)
    EdepT = 5.7 * keV * N_AV

    @testset "tritium storage bed: inventory, decay heat to the wall, bed temperature" begin
        host = liqdata(:ZrCo, (Zr = 1, Co = 1); c = 60.0)
        Tb = liqdata(:T, (T = 1,); c = 25.0, h_f = -20e3)                       # tritium bound in the bed
        He3b = liqdata(:He3, (He3 = 1,); c = 21.0)                                # helium trapped in the bed
        bed = IdealSolution(host, Tb, He3b)
        nT0 = 100e-3 / molar_mass(Tb)                                           # 100 g of tritium
        G, Tw = 2.0, 300.0
        sys = ThermoSystem()
        add!(sys, Body(:bed, bed; T = Tw, n = (ZrCo = 10.0 / molar_mass(host), T = nT0, He3 = 1e-20)))
        add!(sys, Reservoir(:wall; T = Tw))
        connect!(sys, Conductor(ThermalConductance(G)), :bed, :wall; name = :mount)
        connect!(sys, NuclearSource(:T => -1, :He3 => 1; rate = FirstOrder(λT), E_dep = EdepT), :bed; name = :decay)
        cs = compile(sys)
        ts = [0.25, 1.0, 5.0, 12.32] .* YEAR
        sol = simulate(cs, (0.0, ts[end]); saveat = ts)
        for (t, u) in zip(sol.t, sol.u)
            nT = moles(cs, u, :bed).T
            @test nT ≈ nT0 * exp(-λT * t) rtol = 1e-6
            # quasi-steady (thermal time ≪ decay time): G(T − T_w) = ξ̇ (E_dep + u_T − u_He3)
            Tbed = temperature(cs, u, :bed)
            Δu = (Tb.h_f + Tb.cv * (Tbed - 298.15)) * molar_mass(Tb) - (He3b.h_f + He3b.cv * (Tbed - 298.15)) * molar_mass(He3b)
            @test Tbed - Tw ≈ λT * nT * (EdepT + Δu) / G rtol = 1e-4
            # the heat reaching the wall is the decay energy, less what the bed stored
            @test reservoir_ledger(cs, u, :wall).U ≈ EdepT * (nT0 - nT) rtol = 1e-3
        end
        @test moles(cs, sol.u[end], :bed).T ≈ nT0 / 2 rtol = 1e-6                   # one half-life
        @test 0.32 <= λT * nT0 * EdepT / 100 <= 0.33                                 # W per g at t = 0
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9 && cl.charge_rel <= 1e-11
    end

    @testset "tritium generation in irradiated water: D(n,γ)T to saturation" begin
        wd(name, atoms) = liqdata(name, atoms; c = 75.385, h_f = -285830.0, s0 = 69.95)
        water = IdealSolution(wd(:H2O, (H = 2, O = 1)), wd(:HDO, (H = 1, D = 1, O = 1)), wd(:HTO, (H = 1, T = 1, O = 1)),
                              liqdata(:He3, (He3 = 1,); c = 20.786), liqdata(:OH, (O = 1, H = 1); c = 30.0))
        nw = 1.0 / 18.015e-3                                                    # 1 kg of water
        nD0 = 2 * 1.5576e-4 * nw                                                 # natural D/H, as HDO
        σ, φ = 0.519e-31, 1e18                                                   # m², n/(m² s)
        a = σ * φ
        sys = ThermoSystem()
        add!(sys, Body(:water, water; T = 300.0, n = (H2O = nw - nD0, HDO = nD0, HTO = 1e-20, He3 = 1e-20, OH = 1e-20)))
        add!(sys, Reservoir(:loop; T = 300.0))
        connect!(sys, Conductor(ThermalConductance(10.0)), :water, :loop)
        connect!(sys, NuclearSource(:HDO => -1, :HTO => 1; rate = FirstOrder(a), external = (nucleons = 1, charge = 0)),
                 :water; name = :capture)
        connect!(sys, NuclearSource(:HTO => -1, :He3 => 1, :OH => 1; rate = FirstOrder(λT), E_dep = EdepT), :water; name = :decay)
        cs = compile(sys)
        ts = [1.0, 10.0, 50.0] .* YEAR
        sol = simulate(cs, (0.0, ts[end]); saveat = ts)
        for (t, u) in zip(sol.t, sol.u)
            n = moles(cs, u, :water)
            nHTO = nD0 * a * (exp(-a * t) - exp(-λT * t)) / (λT - a)                  # two-member Bateman
            @test n.HTO ≈ nHTO rtol = 1e-6
            @test n.HDO ≈ nD0 * exp(-a * t) rtol = 1e-9
            @test n.He3 ≈ nD0 - n.HDO - n.HTO rtol = 1e-6
        end
        sat = a * nD0 / λT                                                        # saturation inventory
        @test moles(cs, sol.u[end], :water).HTO / sat ≈ 1 - exp(-λT * 50YEAR) rtol = 1e-4
        # after one year: inventory ≈ a nD0 t (decay still small)
        @test moles(cs, sol.u[1], :water).HTO ≈ a * nD0 * YEAR rtol = 0.03
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9 && cl.charge_rel <= 1e-12
        # per-species closure, not diluted by the water: the tritium inventory is accounted to round-off
        @test Set(keys(cl.species)) == Set([:H2O, :HDO, :HTO, :He3, :OH])
        @test cl.species[:HTO].abs <= 1e-15                                       # kg, ≈ eps × 1 kg of water
        @test cl.species[:HTO].abs == maximum(abs(total_charge(cs, u, Species(:HTO)) - total_charge(cs, sol.prob.u0, Species(:HTO)))
                                              for u in sol.u)
        @test isempty(cl.elements)                                                # no reactions here
    end
end
