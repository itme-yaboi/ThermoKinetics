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

@testset "applications: rockets" begin
    g0 = 9.80665

    "Supersonic root of the area–Mach relation, by bisection (independent of the library)."
    function mach_of(ε, k)
        f(M) = (1 / M) * ((2 / (k + 1)) * (1 + (k - 1) / 2 * M^2))^((k + 1) / (2 * (k - 1))) - ε
        lo, hi = 1.0, 100.0
        for _ in 1:200
            mid = (lo + hi) / 2
            f(mid) > 0 ? (hi = mid) : (lo = mid)
        end
        return (lo + hi) / 2
    end
    "Ideal frozen nozzle from (T_c, p_c, k, R): ṁ per throat area, v_e, p_e, vacuum c_eff."
    function ideal_nozzle(Tc, pc, k, R, At, ε)
        Γ = sqrt(k) * (2 / (k + 1))^((k + 1) / (2 * (k - 1)))
        mdot = At * pc * Γ / sqrt(R * Tc)
        Me = mach_of(ε, k)
        Te = Tc / (1 + (k - 1) / 2 * Me^2)
        pe = pc * (Te / Tc)^(k / (k - 1))
        cp = k * R / (k - 1)
        ve = sqrt(2cp * (Tc - Te))
        F = mdot * ve + pe * ε * At
        return (; mdot, Me, Te, pe, ve, F, c_eff = F / mdot, cstar = pc * At / mdot)
    end

    @testset "cold-gas thruster: isentropic tank, Δv by quadrature, thrust = ṁ c_eff" begin
        N2d = N2g()
        R = R̄ / molar_mass(N2d); cv = N2d.cv; k = (cv + R) / cv
        gas = IdealGas(:N2; R, cv)
        m_dry, T0, p0, V = 5.0, 300.0, 200e5, 0.05
        m0 = p0 * V / (R * T0)
        At, ε = 1e-5, 50.0
        sys = ThermoSystem()
        add!(sys, Body(:tank, Moving(gas; m_dry); T = T0, V, m = m0, w = 0.0))
        connect!(sys, Nozzle(At, ε), :tank; name = :jet)
        cs = compile(sys)
        c0 = ideal_nozzle(T0, p0, k, R, At, ε).c_eff
        c_eff(m) = c0 * (m / m0)^((k - 1) / 2)                     # √T on the isentrope
        τ = m0 / ideal_nozzle(T0, p0, k, R, At, ε).mdot
        ts = [0.1, 0.3, 0.6, 1.0] .* τ
        sol = simulate(cs, (0.0, ts[end]); saveat = ts)
        for u in sol.u
            m = charges(cs, u, :tank)[3]
            @test temperature(cs, u, :tank) ≈ T0 * (m / m0)^(k - 1) rtol = 1e-7      # stays on its isentrope
            ns = nozzle_state(cs, u, :jet)
            @test ns.thrust ≈ ns.mdot * c_eff(m) rtol = 1e-7
            # Δv = ∫ c_eff dm/(m_dry + m) from m to m0 (Simpson, 2000 intervals)
            n = 2000; h = (m0 - m) / n
            fq(x) = c_eff(x) / (m_dry + x)
            Δv = h / 3 * (fq(m) + fq(m0) + sum((isodd(j) ? 4 : 2) * fq(m + j * h) for j in 1:n-1))
            @test force(cs, u, :tank, Momentum()) ≈ Δv rtol = 1e-6
        end
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9 && cl.charge_rel <= 1e-10   # momentum included
    end

    # --- the bipropellant engine: H₂/O₂ gas, pump-fed, equilibrium chamber, frozen ideal nozzle ---
    H2, O2, H2O = H2g(), O2g(), H2Og()
    gH2, gO2, gH2O = MolarGas(H2), MolarGas(O2), MolarGas(H2O)
    mdot_H2, mdot_O2 = 0.05, 0.25                                    # kg/s, O/F = 5 (fuel-rich)
    T_in, p_in = 300.0, 200e5                                        # feed state (tank state in flight)
    At, ε = 1.4e-4, 40.0
    burn = MassAction(1e10)                                          # fast: an equilibrium chamber
    mixc = IdealGasMixture(H2, O2, H2O)

    "Independent hand calculation of the steady chamber: equilibrium, energy, choked flow."
    function hand_chamber()
        nH2 = mdot_H2 / molar_mass(H2); nO2 = mdot_O2 / molar_mass(O2)       # mol/s in
        gs = (gH2, gO2, gH2O)
        function resid(x)
            T, p, ξ = x[1], x[2], x[3]
            n = (nH2 - 2ξ, nO2 - ξ, 2ξ)
            N = sum(n)
            # delivered (pumped) state: isentropic from (T_in, p_in) to p
            Tdel = (T_in * (p / p_in)^(R̄ / gH2.cp), T_in * (p / p_in)^(R̄ / gO2.cp))
            Hin = nH2 * mh_(gH2, Tdel[1]) + nO2 * mh_(gO2, Tdel[2])
            Hout = sum(n[j] * mh_(gs[j], T) for j in 1:3)
            ΔG = 2mh_(gH2O, T) - 2mh_(gH2, T) - mh_(gO2, T) - T * (2ms_(gH2O, T, 101325.0) - 2ms_(gH2, T, 101325.0) - ms_(gO2, T, 101325.0))
            x_ = n ./ N
            lnQ = 2log(x_[3]) - 2log(x_[1]) - log(x_[2]) - log(p / 101325.0)
            mtot = mdot_H2 + mdot_O2
            cv_kg = sum(n[j] * gs[j].cv for j in 1:3) / mtot             # molar c_v × mol/s ÷ kg/s
            Rk = N * R̄ / mtot
            k = (cv_kg + Rk) / cv_kg
            Γ = sqrt(k) * (2 / (k + 1))^((k + 1) / (2 * (k - 1)))
            mchoke = At * p * Γ / sqrt(Rk * T)
            return [(Hin - Hout) / (abs(Hin) + 1.0), lnQ + ΔG / (R̄ * T), mchoke / mtot - 1]
        end
        x = [3000.0, 50e5, 0.95 * nO2]
        for _ in 1:100
            r = resid(x)
            J = ForwardDiff.jacobian(resid, x)
            dx = -(J \ r)
            # damped step keeping ξ inside (0, nO2)
            λ = 1.0
            while !(0 < x[3] + λ * dx[3] < min(nO2, nH2 / 2) && x[1] + λ * dx[1] > 0 && x[2] + λ * dx[2] > 0)
                λ /= 2
            end
            x .+= λ .* dx
            maximum(abs, r) < 1e-14 && break
        end
        T, p, ξ = x
        n = (nH2 - 2ξ, nO2 - ξ, 2ξ)
        mtot = mdot_H2 + mdot_O2
        cv_kg = sum(n[j] * (gH2, gO2, gH2O)[j].cv for j in 1:3) / mtot
        Rk = sum(n) * R̄ / mtot
        return (; T, p, n, x = n ./ sum(n), k = (cv_kg + Rk) / cv_kg, R = Rk)
    end

    "Static-fire model: pumped boundary feeds, the chamber, the reaction, the nozzle on a stand."
    function static_fire(; η = 1.0)
        sys = ThermoSystem()
        add!(sys, Body(:chamber, mixc; T = 3000.0, V = 2e-3, n = (H2 = 0.02, O2 = 0.002, H2O = 0.03)))
        add!(sys, Reservoir(:stand, Momentum() => 0.0))
        connect!(sys, FlowBoundary(IdealGasMixture(H2); T = T_in, p = p_in, mdot = t -> mdot_H2, supply = :pumped),
                 :chamber; name = :inj_H2)
        connect!(sys, FlowBoundary(IdealGasMixture(O2); T = T_in, p = p_in, mdot = t -> mdot_O2, supply = :pumped),
                 :chamber; name = :inj_O2)
        connect!(sys, Reaction(:H2 => -2, :O2 => -1, :H2O => 2; law = burn), :chamber; name = :combustion)
        connect!(sys, Nozzle(At, ε; η), :chamber, :stand; name = :nozzle)
        cs = compile(sys)
        u = steady_state(cs; u0 = simulate(cs, (0.0, 0.5)).u[end])
        return cs, u
    end

    hc = hand_chamber()
    ref = ideal_nozzle(hc.T, hc.p, hc.k, hc.R, At, ε)

    @testset "static fire: steady chamber and performance vs an independent hand calculation" begin
        cs, u = static_fire()
        @test temperature(cs, u, :chamber) ≈ hc.T rtol = 1e-7
        @test pressure(cs, u, :chamber) ≈ hc.p rtol = 1e-7
        n = moles(cs, u, :chamber); N = sum(values(n))
        @test n.H2O / N ≈ hc.x[3] rtol = 1e-7
        @test n.O2 / N ≈ hc.x[2] rtol = 1e-5                          # the small equilibrium residual
        ns = nozzle_state(cs, u, :nozzle)
        @test ns.mdot ≈ mdot_H2 + mdot_O2 rtol = 1e-9                   # steady: out = in
        @test ns.M_e ≈ ref.Me rtol = 1e-8
        @test ns.thrust ≈ ref.F rtol = 1e-7
        du = rates(cs, u)
        @test reservoir_charge(cs, du, :stand, Momentum()) ≈ ref.F rtol = 1e-7   # the load cell
        C_F = ref.F / (hc.p * At); Isp = ref.F / (ns.mdot * g0)
        @test ns.thrust / (pressure(cs, u, :chamber) * At) ≈ C_F rtol = 1e-7
        @test ns.thrust / (ns.mdot * g0) ≈ Isp rtol = 1e-7
        cstar = sqrt(hc.R * hc.T / hc.k) / ((2 / (hc.k + 1))^((hc.k + 1) / (2 * (hc.k - 1))))
        @test pressure(cs, u, :chamber) * At / ns.mdot ≈ cstar rtol = 1e-7   # p_c = ṁ c*/A_t
        # loss map. The chamber sits at chemical equilibrium (A ≈ 0), so the combustion
        # irreversibility is booked where the cold propellants meet it: at the injectors, whose
        # generation includes the jump from the feed's μ to the chamber's equilibrium μ (D-48).
        # The reaction element carries almost nothing; the ideal nozzle exactly nothing.
        S = entropy_production(cs, u)
        @test S[:inj_H2] > 0
        @test S[:inj_O2] > 0
        @test S[:combustion] >= 0
        @test S[:combustion] < 1e-6 * (S[:inj_H2] + S[:inj_O2])
        @test S[:nozzle] == 0
        # steady: total production = entropy out with the jet − entropy in with the feeds
        jet, iH, iO = boundary_ledger(cs, du, :nozzle), boundary_ledger(cs, du, :inj_H2), boundary_ledger(cs, du, :inj_O2)
        @test sum(values(S)) ≈ jet.S + iH.S + iO.S rtol = 1e-8           # feeds' S received is negative
        r = energy_residual(cs, u)
        @test abs(r.residual) <= 1e-10 * r.scale
        # with a lossy nozzle: same chamber, generation ṁ c_p ln(T_e/T_es) (frozen, constant c_p)
        cs2, u2 = static_fire(; η = 0.95)
        ns2 = nozzle_state(cs2, u2, :nozzle)
        @test temperature(cs2, u2, :chamber) ≈ hc.T rtol = 1e-7
        cp = hc.k * hc.R / (hc.k - 1)
        Te2 = ref.Te + 0.05 * (hc.T - ref.Te)
        @test entropy_production(cs2, u2)[:nozzle] ≈ ns2.mdot * cp * log(Te2 / ref.Te) rtol = 1e-6
        @test ns2.v_e ≈ sqrt(0.95) * ref.ve rtol = 1e-7
    end

    "Flight model: Moving tanks (one-species mixtures at the feed state) pumped into a Moving chamber
    that starts at the steady static-fire state; rigid links; gravity on every body; vacuum."
    function flight(u_static, cs_static; w0 = 0.0, gravity = true, m_dry = 50.0, mH2 = 3.0, mO2 = 15.0)
        Ec = charges(cs_static, u_static, :chamber)
        mixm = Moving(mixc)
        sys = ThermoSystem()
        add!(sys, Body(:chamber, mixm, vcat(Ec, (Ec[3] + Ec[4] + Ec[5]) * w0)))
        add!(sys, Body(:tank_H2, Moving(IdealGasMixture(H2); m_dry); T = T_in, p = p_in, m = (H2 = mH2,), w = w0))
        add!(sys, Body(:tank_O2, Moving(IdealGasMixture(O2)); T = T_in, p = p_in, m = (O2 = mO2,), w = w0))
        connect!(sys, IdealPump(t -> mdot_H2), :tank_H2, :chamber; name = :pump_H2)
        connect!(sys, IdealPump(t -> mdot_O2), :tank_O2, :chamber; name = :pump_O2)
        connect!(sys, Transformer(1.0; a = Momentum(), b = Momentum()), :tank_H2, :chamber; name = :link_H2)
        connect!(sys, Transformer(1.0; a = Momentum(), b = Momentum()), :tank_O2, :chamber; name = :link_O2)
        connect!(sys, Reaction(:H2 => -2, :O2 => -1, :H2O => 2; law = burn), :chamber; name = :combustion)
        connect!(sys, Nozzle(At, ε), :chamber; name = :nozzle)
        if gravity
            for b in (:chamber, :tank_H2, :tank_O2)
                connect!(sys, GravityField(g0), b; name = Symbol(:g_, b))
            end
        end
        return compile(sys)
    end

    @testset "flight: Beer & Johnston Sample Problem 14.9, v = u ln(m₀/(m₀ − qt)) − gt, and the altitude" begin
        cs_s, u_s = static_fire()
        u_ex = nozzle_state(cs_s, u_s, :nozzle).thrust / (mdot_H2 + mdot_O2)    # exhaust speed (c_eff, vacuum)
        cs = flight(u_s, cs_s)
        q = mdot_H2 + mdot_O2
        m0 = 50.0 + 3.0 + 15.0 + sum(charges(cs_s, u_s, :chamber)[3:5])
        ts = [2.0, 10.0, 20.0, 30.0]
        sol = simulate(cs, (0.0, ts[end]); saveat = ts)
        for (t, u) in zip(sol.t, sol.u)
            # Beer & Johnston (Vector Mechanics, Dynamics) Sample Problem 14.9
            v_exact = u_ex * log(m0 / (m0 - q * t)) - g0 * t
            x_exact = u_ex * (t - (m0 - q * t) / q * log(m0 / (m0 - q * t))) - g0 * t^2 / 2
            @test force(cs, u, :chamber, Momentum()) ≈ v_exact rtol = 1e-6
            @test force(cs, u, :tank_O2, Momentum()) ≈ v_exact rtol = 1e-6            # rigidly linked
            @test source_ledger(cs, u, :g_chamber).x ≈ x_exact rtol = 1e-6            # burnout altitude ∫ w dt
            @test temperature(cs, u, :chamber) ≈ hc.T rtol = 1e-6                      # chamber stays steady
            @test pressure(cs, u, :chamber) ≈ hc.p rtol = 1e-6
        end
        # momentum: bodies + jet + gravity fields received = 0 (the links conserve it)
        pt(u) = sum(charges(cs, u, b)[end] for b in (:chamber, :tank_H2, :tank_O2)) +
                boundary_ledger(cs, u, :nozzle).p + sum(source_ledger(cs, u, Symbol(:g_, b)).charge for b in (:chamber, :tank_H2, :tank_O2))
        @test abs(pt(sol.u[end])) <= 1e-9 * q * u_ex * ts[end]
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9 && cl.charge_rel <= 1e-10
        # Galilean invariance of the whole engine: started at 3 km/s, without gravity
        csA, csB = flight(u_s, cs_s; gravity = false), flight(u_s, cs_s; gravity = false, w0 = 3000.0)
        sA = simulate(csA, (0.0, 10.0); saveat = [10.0]); sB = simulate(csB, (0.0, 10.0); saveat = [10.0])
        @test force(csB, sB.u[end], :chamber, Momentum()) - 3000.0 ≈ force(csA, sA.u[end], :chamber, Momentum()) rtol = 1e-7
        @test temperature(csB, sB.u[end], :chamber) ≈ temperature(csA, sA.u[end], :chamber) rtol = 1e-9
        @test nozzle_state(csB, sB.u[end], :nozzle).thrust ≈ nozzle_state(csA, sA.u[end], :nozzle).thrust rtol = 1e-9
    end
end

@testset "applications: power cycles, Otto (08 E4)" begin
    # ideal Otto cycle, cold-air standard: a piston (prescribed motion, via a gyrator to a time-varying
    # momentum reservoir) compresses and expands the gas isentropically; heat is added at constant
    # volume by a prescribed HeatSource and rejected at constant volume through a switched conductor.
    # Closed form: Çengel 9e Example 9-3 (cold-air-standard result) and Eq. 9-8
    air = IdealGas(:air; R = 287.0, cv = 287.0 / 0.4)
    k = (air.cv + air.R) / air.cv
    T1, p1, m, r, qin, A = 290.0, 1e5, 1e-3, 8.0, 800e3, 0.01
    V1 = m * air.R * T1 / p1
    V2 = V1 / r
    Qin = qin * m
    tc, th, tr = 0.02, 0.01, 30.0
    t1, t2, t3, tend = tc, tc + th, 2tc + th, 2tc + th + tr
    G = 50.0                                                  # relaxation time m cv / G ≈ 0.008 s

    T2 = T1 * r^(k - 1); p2 = p1 * r^k
    T3 = T2 + qin / air.cv; p3 = p2 * T3 / T2
    T4 = T3 * r^(1 - k); p4 = p1 * T4 / T1
    qout = air.cv * (T4 - T1)
    ηref = 1 - r^(1 - k)

    w(t) = t < t1 ? -(V1 - V2) / A * (π / (2tc)) * sin(π * t / tc) :
           t < t2 ? 0.0 :
           t < t3 ? (V1 - V2) / A * (π / (2tc)) * sin(π * (t - t2) / tc) : 0.0
    sys = ThermoSystem()
    gas = add!(sys, Body(:gas, air; T = T1, p = p1, m = m))
    piston = add!(sys, Reservoir(:piston, Momentum() => w))
    sink = add!(sys, Reservoir(:sink; T = T1))
    connect!(sys, Gyrator(A; a = Volumetric(), b = Momentum()), gas, piston)             # dV/dt = A w
    connect!(sys, HeatSource(t -> t1 <= t < t2 ? Qin / th : 0.0), gas; name = :heater)
    connect!(sys, Switched(Conductor(ThermalConductance(G)); on = t -> t >= t3), gas, sink; name = :reject)
    cs = compile(sys)
    stops = [t1, t2, t3]
    sol = simulate(cs, (0.0, tend); saveat = [0.0, t1, t2, t3, tend], tstops = stops)
    u0, u1, u2, u3, u4 = sol.u

    @testset "states at the end of each phase match the closed form" begin
        for (u, T, p) in ((u1, T2, p2), (u2, T3, p3), (u3, T4, p4), (u4, T1, p1))
            @test temperature(cs, u, gas) ≈ T rtol = 1e-6
            @test pressure(cs, u, gas) ≈ p rtol = 1e-6
        end
        @test temperature(cs, u4, gas) ≈ T1 rtol = 1e-8
        @test charges(cs, u1, gas)[2] ≈ V2 rtol = 1e-8
        @test charges(cs, u2, gas)[2] ≈ V2 rtol = 1e-8
        @test charges(cs, u4, gas)[2] ≈ V1 rtol = 1e-8
    end

    @testset "compression and expansion are isentropic" begin
        S(u) = charges(cs, u, gas)[1]
        @test abs(S(u1) - S(u0)) <= 1e-8 * m * air.cv
        @test abs(S(u3) - S(u2)) <= 1e-8 * m * air.cv
    end

    @testset "net work and thermal efficiency" begin
        Wnet = reservoir_ledger(cs, u4, piston).U - reservoir_ledger(cs, u0, piston).U
        Qin_l = -boundary_ledger(cs, u4, :heater).U
        Qout = reservoir_ledger(cs, u4, sink).U
        @test Qin_l ≈ Qin rtol = 1e-6
        @test Qout ≈ qout * m rtol = 1e-6
        @test Wnet ≈ Qin_l - Qout rtol = 1e-6
        @test Wnet / Qin_l ≈ ηref rtol = 1e-6
        @test Wnet / Qin_l ≈ 0.565 rtol = 2e-3                 # the textbook's rounding (56.5 %)
    end

    @testset "loss map: only the rejection conductor produces entropy" begin
        gen = generated_entropy(cs, u4)
        @test haskey(gen, :reject)
        @test !haskey(gen, :heater)
        @test !haskey(entropy_production(cs, u4), :heater)
        Sexp = m * air.cv * log(T1 / T4) + qout * m / T1
        @test gen[:reject] > 0
        @test gen[:reject] ≈ Sexp rtol = 1e-6
        @test sum(values(gen)) ≈ gen[:reject] rtol = 1e-12
    end

    cl = closure(cs, sol)
    @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9
end

@testset "applications: real-fluid Brayton cycles and a steam turbine (slice 10a)" begin
    # Independent reference: Clapeyron's own pressure-based property functions, with an isentropic
    # temperature found here by bisection on Clapeyron.entropy(model, p, T, z). The library goes
    # through the volume-based Helmholtz energy and its own inversion (D-50).
    clz(f, y) = (n = collect(y) ./ collect(Clapeyron.mw(f.model) ./ 1000); n ./ sum(n))
    Mbar(f, z) = sum(z .* collect(Clapeyron.mw(f.model) ./ 1000))
    h_pT(f, p, T, z) = Clapeyron.enthalpy(f.model, p, T, z) / Mbar(f, z)
    s_pT(f, p, T, z) = Clapeyron.entropy(f.model, p, T, z) / Mbar(f, z)
    function T_isentropic(f, s, p, z; lo = 100.0, hi = 3000.0)
        for _ in 1:200
            mid = (lo + hi) / 2
            s_pT(f, p, mid, z) > s ? (hi = mid) : (lo = mid)
        end
        return (lo + hi) / 2
    end

    "Closed Brayton loop: a cold body (compressor inlet) and a hot body (turbine inlet), an ideal
    compressor cold → hot and an ideal turbine hot → cold at the same ṁ, and the two bodies on
    reservoirs at T1 and T3 through stiff conductances. Returns rates at the relaxed state."
    function brayton(f, T1, p1, T3, rp, y; mdot = 0.01, V = 0.05, G = 1e5)
        m = y isa Real ? 1.0 : NamedTuple{f.names}(Tuple(y))
        Ec = charges_at(f; T = T1, p = p1, m = m)
        Eh = charges_at(f; T = T3, p = rp * p1, m = m)
        sc, sh = V / Ec[2], V / Eh[2]                            # scale both to volume V
        sys = ThermoSystem()
        add!(sys, Body(:cold, f, Ec .* sc))
        add!(sys, Body(:hot, f, Eh .* sh))
        add!(sys, Reservoir(:sink; T = T1))
        add!(sys, Reservoir(:source; T = T3))
        connect!(sys, Conductor(ThermalConductance(G)), :sink, :cold; name = :cooler)
        connect!(sys, Conductor(ThermalConductance(G)), :source, :hot; name = :heater)
        connect!(sys, IdealPump(t -> mdot), :cold, :hot; name = :comp)
        connect!(sys, IdealPump(t -> mdot), :hot, :cold; name = :turb)
        cs = compile(sys)
        sol = simulate(cs, (0.0, 20.0); saveat = [20.0])
        u = sol.u[end]
        du = rates(cs, u)
        W = pump_work(cs, du)
        return (; cs, u, du, Wc = W[:comp], Wt = W[:turb], Qin = -reservoir_ledger(cs, du, :source).U,
                Qout = reservoir_ledger(cs, du, :sink).U, mdot, cl = closure(cs, sol),
                T2 = ThermoKinetics._eval(cs, u, :comp).T_star, T4 = ThermoKinetics._eval(cs, u, :turb).T_star)
    end

    "The same cycle by hand, from the bodies' actual states, with Clapeyron's (p, T) functions only."
    function by_hand(f, r, y)
        z = clz(f, y)
        # qualified: Clapeyron also exports `temperature` and `pressure`
        T1, p1 = ThermoKinetics.temperature(r.cs, r.u, :cold), ThermoKinetics.pressure(r.cs, r.u, :cold)
        T3, p3 = ThermoKinetics.temperature(r.cs, r.u, :hot), ThermoKinetics.pressure(r.cs, r.u, :hot)
        h1, s1 = h_pT(f, p1, T1, z), s_pT(f, p1, T1, z)
        h3, s3 = h_pT(f, p3, T3, z), s_pT(f, p3, T3, z)
        T2 = T_isentropic(f, s1, p3, z); T4 = T_isentropic(f, s3, p1, z)
        h2, h4 = h_pT(f, p3, T2, z), h_pT(f, p1, T4, z)
        return (; T2, T4, wc = h2 - h1, wt = h3 - h4, qin = h3 - h2, η = ((h3 - h4) - (h2 - h1)) / (h3 - h2))
    end

    function check_cycle(f, r, hd; Tsink, rtol = 1e-6)
        @test r.Wc ≈ r.mdot * hd.wc rtol = rtol                  # compressor work, isentropic
        @test -r.Wt ≈ r.mdot * hd.wt rtol = rtol                 # turbine work (negative: an ideal turbine)
        @test r.T2 ≈ hd.T2 rtol = 1e-8
        @test r.T4 ≈ hd.T4 rtol = 1e-8
        @test r.Qin ≈ r.mdot * hd.qin rtol = rtol                # steady: the heater supplies h3 − h2
        Wnet = -(r.Wc + r.Wt)
        @test Wnet ≈ r.Qin - r.Qout rtol = 1e-6                  # first law over the steady loop
        @test Wnet / r.Qin ≈ hd.η rtol = rtol
        # loss map: Gouy–Stodola between the two reservoirs, W_Carnot − W = T_sink Σ Ṡ_gen
        Tsrc = ThermoKinetics.reservoir_force(r.cs.reservoirs[findfirst(==(:source), r.cs.reservoir_names)], Thermal())
        Sgen = sum(values(entropy_production(r.cs, r.u)))
        @test r.Qin * (1 - Tsink / Tsrc) - Wnet ≈ Tsink * Sgen rtol = 1e-6
        @test r.cl.energy_rel <= 1e-6 && r.cl.entropy_rel <= 1e-9
    end

    @testset "Çengel 9e Example 9-6: ideal Brayton, air, rp = 8, 300 K → 1300 K" begin
        air = RealFluid(Clapeyron.GERG2008(["nitrogen", "oxygen", "argon"]); species = (:N2, :O2, :Ar))
        y = (0.7552, 0.2314, 0.0134)                             # dry air, mass fractions
        r = brayton(air, 300.0, 100e3, 1300.0, 8.0, y)
        hd = by_hand(air, r, y)
        check_cycle(air, r, hd; Tsink = 300.0)
        # vs the textbook (air tables, variable specific heats): T2 = 540 K, T4 = 770 K, rbw = 0.403,
        # η = 0.426; the bodies sit within ~1 K of 300 K and 1300 K (stiff, not infinite, conductances)
        @test r.T2 ≈ 540.0 rtol = 1e-2
        @test r.T4 ≈ 770.0 rtol = 1e-2
        @test r.Wc / -r.Wt ≈ 0.403 rtol = 1e-2
        @test -(r.Wc + r.Wt) / r.Qin ≈ 0.426 rtol = 1e-2
    end

    @testset "closed helium Brayton loop (fusion balance of plant)" begin
        he = RealFluid(Clapeyron.GERG2008(["helium"]); species = :He)
        r = brayton(he, 300.0, 3e6, 800.0, 2.5, 1.0)
        hd = by_hand(he, r, (1.0,))
        check_cycle(he, r, hd; Tsink = 300.0)
        k = 5 / 3                                                # helium is nearly ideal: η ≈ 1 − rp^((1−k)/k)
        @test -(r.Wc + r.Wt) / r.Qin ≈ 1 - 2.5^((1 - k) / k) rtol = 1e-2
    end

    @testset "simple supercritical CO₂ Brayton cycle near the critical point" begin
        co2 = RealFluid(Clapeyron.SingleFluid("carbon dioxide"); species = :CO2)
        r = brayton(co2, 305.0, 7.8e6, 823.0, 2.6, 1.0; V = 0.01)
        hd = by_hand(co2, r, (1.0,))
        check_cycle(co2, r, hd; Tsink = 305.0)
        @test r.Wc / -r.Wt < 0.3                                 # dense near-critical compression: low back work
    end

    @testset "Çengel 9e Example 7-5: ideal steam turbine, 5 MPa 450 °C → 1.4 MPa" begin
        water = RealFluid(Clapeyron.IAPWS95(); species = :H2O)
        sys = ThermoSystem()
        add!(sys, Body(:inlet, water; T = 723.15, p = 5e6, m = 1.0))
        add!(sys, Body(:exhaust, water; T = 500.0, p = 1.4e6, m = 1.0))
        connect!(sys, IdealPump(t -> 1.0), :inlet, :exhaust; name = :turbine)     # 1 kg/s
        cs = compile(sys)
        u = initial_state(cs)
        E = charges(cs, u, :inlet)
        h1 = (ThermoKinetics.internal_energy(water, E) + 5e6 * E[2]) / E[3]
        @test h1 ≈ 3317.2e3 rtol = 1e-3                           # Table A-6
        @test E[1] / E[3] ≈ 6.8210e3 rtol = 1e-3
        # w = h1 − h2s. The textbook finds h2s by linear interpolation in s between table rows; h(s) is
        # convex along an isobar, so that overestimates h2s by ~1 kJ/kg (IAPWS-95 gives w = 350.9 kJ/kg)
        @test -pump_work(cs, rates(cs, u))[:turbine] ≈ 349.8e3 rtol = 5e-3
        st = ThermoKinetics.isentropic_state(water, E[1] / E[3], 1.4e6)
        @test 520 < st.T < 550                                    # superheated at 1.4 MPa (T_sat = 468 K): ≈ 262 °C
    end
end
