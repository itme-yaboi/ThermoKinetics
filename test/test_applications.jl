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
