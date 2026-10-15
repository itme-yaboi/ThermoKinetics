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
        # ledger_control: the drive makes V linear in t and the adiabats keep S fixed, so without
        # error control on the ledgers a whole leg is one step and ∫p dV is a one-step quadrature (D-29)
        sol = simulate(cs, (0.0, ncycles * period); saveat = vcat(0.0, stops), tstops = stops, ledger_control = true)
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
