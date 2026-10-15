# 08 §E: Carnot cycles. An ideal-gas body whose volume is driven through a closed path by a
# prescribed-velocity drive (a time-varying momentum reservoir, joined by a gyrator), with
# switched thermal contacts to a hot and a cold reservoir.

@testset "08 §E cycles" begin
    N2 = gas_n2()
    Th, Tc = 500.0, 300.0
    V1, p1 = 1e-3, 5e5
    m = p1 * V1 / (N2.R * Th)
    V2 = 2V1
    V3 = V2 * (Th / Tc)^(N2.cv / N2.R)                    # reversible adiabats of the ideal gas
    V4 = V1 * (Th / Tc)^(N2.cv / N2.R)
    tl = 100.0                                            # duration of each leg [s]
    period = 4tl

    function carnot(G; ncycles = 2)
        rates = ((V2 - V1) / tl, (V3 - V2) / tl, (V4 - V3) / tl, (V1 - V4) / tl)
        leg(t) = mod(floor(Int, t / tl), 4) + 1
        sys = ThermoSystem()
        gas = add!(sys, Body(:gas, N2; T = Th, V = V1, m))
        hot = add!(sys, Reservoir(:hot; T = Th))
        cold = add!(sys, Reservoir(:cold; T = Tc))
        drive = add!(sys, Reservoir(:drive, Momentum() => t -> rates[leg(t)]))
        connect!(sys, Gyrator(1.0; a = Volumetric(), b = Momentum()), gas, drive)     # dV/dt = w(t)
        connect!(sys, Switched(Conductor(ThermalConductance(G)); on = t -> leg(t) == 1), hot, gas; name = :hotc)
        connect!(sys, Switched(Conductor(ThermalConductance(G)); on = t -> leg(t) == 3), cold, gas; name = :coldc)
        cs = compile(sys)
        stops = collect(tl:tl:ncycles*period)
        # ledger control (the default) matters here: the drive makes V linear in t and the adiabats keep
        # S fixed, so without it a whole leg is one step and ∫p dV a one-step quadrature (D-29, D-31)
        sol = simulate(cs, (0.0, ncycles * period); saveat = vcat(0.0, stops), tstops = stops)
        # the last cycle, from ledger differences (periodic by then: each isotherm resets the gas)
        i0 = length(sol.t) - 4
        u(k) = sol.u[i0+k]
        Δ(f) = f(u(4)) - f(u(0))
        W = Δ(v -> reservoir_ledger(cs, v, drive).U)                 # work delivered to the drive
        Qh = -Δ(v -> reservoir_ledger(cs, v, hot).U)
        Qc = Δ(v -> reservoir_ledger(cs, v, cold).U)
        Sgen = Δ(v -> sum(values(generated_entropy(cs, v))))
        S(v) = charges(cs, v, gas)[1]
        return (; W, Qh, Qc, Sgen, ΔU = Δ(v -> total_energy(cs, v)),
                ΔS1 = S(u(1)) - S(u(0)), ΔS3 = S(u(3)) - S(u(2)), cl = closure(cs, sol))
    end

    ηC = 1 - Tc / Th
    r10, r100 = carnot(10.0), carnot(100.0)

    @testset "E1 Carnot: η → 1 − T″/T′ and Ṡ_gen → 0 as the bridges become reversible" begin
        for r in (r10, r100)
            @test abs(r.ΔU) <= 1e-6 * r.Qh                          # a closed cycle
            @test r.W ≈ r.Qh - r.Qc rtol = 1e-6                     # first law over the cycle
            @test r.W / r.Qh < ηC                                   # never above Carnot
            @test r.Sgen > 0
            @test r.cl.energy_rel <= 1e-6 && r.cl.entropy_rel <= 1e-9
        end
        e10, e100 = ηC - r10.W / r10.Qh, ηC - r100.W / r100.Qh
        @test e10 < 5e-3
        @test 5 < e10 / e100 < 20                                   # the deficit falls like 1/G
        @test 5 < r10.Sgen / r100.Sgen < 20
    end

    @testset "E2 Carnot with bridges: η = 1 − T₂/T₁ at the internal temperatures (490)" begin
        for r in (r10, r100)
            T1 = r.Qh / r.ΔS1                                       # mean T at which the gas took in heat
            T2 = r.Qc / (-r.ΔS3)                                    # … and gave it out
            @test Tc < T2 && T1 < Th
            @test r.W / r.Qh ≈ 1 - T2 / T1 rtol = 1e-6
            # the lost work is T_c times the entropy produced in the bridges (Gouy–Stodola); exact for a
            # periodic cycle, and this one repeats to ~4e-8 of Q_h (lost work: 0.3 J and 0.03 J of 138 J)
            @test abs((r.Qh * ηC - r.W) - Tc * r.Sgen) <= 1e-7 * r.Qh
            @test r.Qh * ηC - r.W > 0
        end
    end
end

@testset "08 E3 rectangular cycle: Veinik's switched cycle converter (§62, Figs. 69–70)" begin
    # A body with P = B e T (`EntropicCharge`, here on the electrical form of motion: "a piece of
    # metal with interrupter wires") is switched between charge reservoirs at P′ and P″ and kept
    # in thermal contact with a reservoir. Legs: 1→2 at P′, 2→3 at constant e, 3→4 at P″, 4→1 at
    # constant e. Net: charge falls from P′ to P″ and lifts entropy (a heat pump driven by charge).
    C, B, P1, P2, e1, e2 = 10.0, 1.0, 400.0, 300.0, 1.0, 1.2       # J/K, V/(C K), V, V, C, C
    mat = EntropicCharge(Electrical(); C, B, T_ref = 300.0)
    Ts = (P1 / (B * e1), P1 / (B * e2), P2 / (B * e2), P2 / (B * e1))   # T at states 1, 2, 3, 4
    tl = 1.0; period = 4tl; ncycles = 4
    leg(t) = mod(floor(Int, t / tl), 4) + 1

    "Run the cycle. tracking: the reservoir temperature follows the reversible path (matched
    switching); otherwise a hot and a cold reservoir at fixed temperatures. reversed: run 1→4→3→2→1."
    function e3(G_th, G_e; tracking = true, reversed = false)
        path = reversed ? (Ts[1], Ts[4], Ts[3], Ts[2], Ts[1]) : (Ts[1], Ts[2], Ts[3], Ts[4], Ts[1])
        Tpath(t) = (k = leg(t); f = t / tl - floor(t / tl); path[k] + f * (path[k+1] - path[k]))
        # forward: leg 1 at P′, leg 3 at P″; reversed: leg 2 at P″ (charge in), leg 4 at P′ (charge out)
        onP1 = reversed ? (t -> leg(t) == 4) : (t -> leg(t) == 1)
        onP2 = reversed ? (t -> leg(t) == 2) : (t -> leg(t) == 3)
        sys = ThermoSystem()
        add!(sys, Body(:b, mat, charges_at(mat; T = Ts[1], e = e1)))
        add!(sys, Reservoir(:R1, Electrical() => P1))
        add!(sys, Reservoir(:R2, Electrical() => P2))
        if tracking
            add!(sys, Reservoir(:heat, Thermal() => Tpath))
            connect!(sys, Conductor(ThermalConductance(G_th)), :heat, :b; name = :thermal)
            sink = :heat
        else
            # the forward run rejects entropy on legs 1, 2 (to a cold reservoir below the path) and
            # absorbs it on legs 3, 4 (from a hot one above it)
            add!(sys, Reservoir(:hot; T = 410.0)); add!(sys, Reservoir(:cold; T = 240.0))
            rej(t) = leg(t) in (1, 2)
            connect!(sys, Switched(Conductor(ThermalConductance(G_th)); on = t -> !rej(t)), :hot, :b; name = :thot)
            connect!(sys, Switched(Conductor(ThermalConductance(G_th)); on = rej), :cold, :b; name = :tcold)
            sink = :cold
        end
        connect!(sys, Switched(Conductor(LinearConductance(Electrical(), G_e); dissipation = Into(sink)); on = onP1),
                 :R1, :b; name = :tap1)
        connect!(sys, Switched(Conductor(LinearConductance(Electrical(), G_e); dissipation = Into(sink)); on = onP2),
                 :R2, :b; name = :tap2)
        cs = compile(sys)
        stops = collect(tl:tl:ncycles*period)
        sol = simulate(cs, (0.0, ncycles * period); saveat = vcat(0.0, stops), tstops = stops)
        i0 = length(sol.t) - 4                                     # the last cycle
        Δ(f) = f(sol.u[i0+4]) - f(sol.u[i0])
        e_out1 = -Δ(v -> reservoir_charge(cs, v, :R1, Electrical()))   # charge supplied by R1
        e_in2 = Δ(v -> reservoir_charge(cs, v, :R2, Electrical()))      # charge received by R2
        Qres = sum(Δ(v -> reservoir_ledger(cs, v, r).U) for r in (tracking ? (:heat,) : (:hot, :cold)))
        Sgen = Δ(v -> sum(values(generated_entropy(cs, v))))
        Sb = [charges(cs, sol.u[i0+k], :b)[1] for k in 0:4]
        return (; cs, sol, e_out1, e_in2, Qres, Sgen, ΔSthrough = sum(abs, diff(Sb)), Sb,
                eb = [charges(cs, sol.u[i0+k], :b)[2] for k in 0:4], cl = closure(cs, sol))
    end

    r1, r10 = e3(1e3, 10.0), e3(1e4, 100.0)
    @testset "forward: charge falls P′ → P″ and lifts entropy" begin
        for r in (r1, r10)
            q1, q2 = P1 * r.e_out1, P2 * r.e_in2                  # charge work taken at P′, given at P″
            @test r.e_out1 > 0
            @test r.e_in2 > 0
            @test r.e_out1 ≈ r.e_in2 rtol = 1e-6                  # periodic: the body returns to its charge
            @test (q1 - q2) / q1 ≈ 1 - P2 / P1 rtol = 1e-6        # η = 1 − P″/P′ on the charge leg (494), (498)
            @test r.Qres ≈ q1 - q2 rtol = 1e-6                    # ∮T dS = −∮P de: the cycle rejects q as heat
            @test abs(r.Sb[end] - r.Sb[1]) <= 1e-6 * r.ΔSthrough   # entropy returns: a closed cycle
            @test r.Sgen > 0
            @test r.cl.energy_rel <= 1e-6 && r.cl.entropy_rel <= 1e-9 && r.cl.charge_rel <= 1e-10
            @test all(u -> all(>=(0), values(entropy_production(r.cs, u))), r.sol.u)
        end
        # the body's path is the rectangle: e1 → e2 at P′, e2 → e1 at P″ (quasi-static for the stiffer coupling;
        # at G_th = 1e3 W/K the thermal lag shortens the charge swing by ~2 %)
        @test r10.e_out1 ≈ e2 - e1 rtol = 1e-2
        @test r10.eb[3] ≈ e2 rtol = 1e-2
        @test r10.eb[5] ≈ e1 rtol = 1e-2
        # reversible limit: matched-force switching and tracking heat exchange; production falls like 1/G
        @test r10.Sgen / r10.ΔSthrough < 1e-3
        @test 5 < r1.Sgen / r10.Sgen < 20
    end

    @testset "fixed-temperature reservoirs: production does not vanish" begin
        f1, f10 = e3(1e3, 10.0; tracking = false), e3(1e4, 100.0; tracking = false)
        for r in (f1, f10)
            q1, q2 = P1 * r.e_out1, P2 * r.e_in2
            @test (q1 - q2) / q1 ≈ 1 - P2 / P1 rtol = 1e-6        # the force-ratio efficiency still holds (05)
            @test r.Qres ≈ q1 - q2 rtol = 1e-6
            @test r.cl.energy_rel <= 1e-6 && r.cl.entropy_rel <= 1e-9
        end
        @test f10.Sgen > 0.5 * f1.Sgen                            # the ΔT to fixed reservoirs stays
        @test f10.Sgen > 100 * r10.Sgen
    end

    @testset "reversed: entropy falls and lifts charge P″ → P′" begin
        r = e3(1e4, 100.0; reversed = true)
        @test r.e_out1 < 0                                       # R1 receives …
        @test r.e_in2 < 0                                        # … R2 supplies
        @test r.e_out1 ≈ r.e_in2 rtol = 1e-6
        q1, q2 = -P1 * r.e_out1, -P2 * r.e_in2                    # charge work given at P′, taken at P″
        @test r.Qres ≈ -(q1 - q2) rtol = 1e-6                     # the reservoir supplies the lifted charge's work
        @test r.Sgen > 0
        @test r.Sgen / r.ΔSthrough < 1e-3
        @test r.cl.energy_rel <= 1e-6 && r.cl.entropy_rel <= 1e-9
    end

    @testset "EntropicCharge state function" begin
        rng = MersenneTwister(91)
        for _ in 1:100
            T = 100 + 900rand(rng); e = 4 * (rand(rng) - 0.5)
            E = charges_at(mat; T, e)
            F = forces(mat, E)
            @test F[1] ≈ T rtol = 1e-12
            @test F[2] ≈ B * e * T rtol = 1e-12 atol = 1e-12        # P = B e T
            @test internal_energy(mat, E) ≈ C * T rtol = 1e-12
            @test is_stable(mat, E)
            A = eos_matrix(mat, E)
            @test det(A) ≈ B * T^2 / C rtol = 1e-9                  # convex
        end
        @test_throws ArgumentError EntropicCharge(Thermal(); C = 1.0, B = 1.0)
        @test_throws ArgumentError EntropicCharge(Electrical(); C = 1.0, B = -1.0)
    end
end
