# Slice 5: kinetic forms of motion (momentum, angular momentum, flux linkage), ideal converters
# (F = W P, W antisymmetric: gyrators, pistons, motors, inductor branches), zero-capacity
# junctions and transformers (mass-matrix DAE), switched connections, time-varying reservoirs.

@testset "kinetic domains, converters, junctions, transformers" begin
    N2 = gas_n2()

    @testset "construction and validation" begin
        @test charge_symbol(Momentum()) == :p && force_symbol(Momentum()) == :w && force_unit(Momentum()) == "m/s"
        @test charge_symbol(AngularMomentum()) == :L && force_unit(AngularMomentum()) == "rad/s"
        @test charge_symbol(FluxLinkage()) == :Λ && force_symbol(FluxLinkage()) == :i && force_unit(FluxLinkage()) == "A"
        m = InertialMass(2.0)
        E = charges_at(m; w = 3.0)
        @test E[1] ≈ 6.0 && forces(m, E)[1] ≈ 3.0 && internal_energy(m, E) ≈ 9.0     # ½ m w²
        @test forces(RotationalInertia(0.5), charges_at(RotationalInertia(0.5); ω = 4.0))[1] ≈ 4.0
        @test forces(LinearInductor(1e-3), charges_at(LinearInductor(1e-3); i = 2.0))[1] ≈ 2.0
        @test_throws ArgumentError InertialMass(0.0)
        @test_throws ArgumentError Converter([0.0 1.0; 1.0 0.0], (Volumetric(), Momentum()))      # not antisymmetric
        @test_throws ArgumentError Converter([0.0 1.0; -1.0 0.0], (Volumetric(),))                # wrong size
        @test_throws ArgumentError Converter([0.0 1.0; -1.0 0.0], (Thermal(), Momentum()))        # entropy is never converted
        @test Gyrator(0.01; a = Volumetric(), b = Momentum()) isa Converter
        @test_throws ArgumentError Transformer(0.0; a = Electrical(), b = Electrical())
        @test_throws ArgumentError Junction(:j, Thermal())                                        # no zero-capacity thermal node
        @test_throws ArgumentError Junction(:j)
    end

    @testset "LC oscillator: an inductor branch is a 3-port converter; energy is conserved exactly" begin
        C, L, V0 = 1e-3, 1e-2, 5.0
        sys = ThermoSystem()
        cap = add!(sys, Body(:cap, LinearCapacitor(C); φ = V0))
        ind = add!(sys, Body(:ind, LinearInductor(L); i = 0.0))
        gnd = add!(sys, Reservoir(:gnd, Electrical() => 0.0))
        br = connect!(sys, InductorBranch(), cap, gnd, ind)
        cs = compile(sys)
        u0 = initial_state(cs)
        r = energy_residual(cs, u0)
        @test abs(r.residual) <= 1e-12 * max(r.scale, 1e-30)
        ω0 = 1 / sqrt(L * C)
        T = 2π / ω0
        sol = simulate(cs, (0.0, 3T); saveat = T / 8, reltol = 1e-10)
        for (t, u) in zip(sol.t, sol.u)
            @test force(cs, u, cap, Electrical()) ≈ V0 * cos(ω0 * t) atol = 1e-6 * V0
            @test force(cs, u, ind, FluxLinkage()) ≈ V0 * C * ω0 * sin(ω0 * t) atol = 1e-6 * V0 * C * ω0
        end
        E0 = total_energy(cs, u0)
        @test maximum(abs(total_energy(cs, u) - E0) for u in sol.u) <= 1e-7 * E0
        @test isempty(entropy_production(cs, u0))                  # an ideal converter produces nothing
    end

    @testset "RLC: all the stored energy ends up as heat" begin
        C, L, R, V0, T = 1e-3, 1e-2, 1.0, 5.0, 300.0
        sys = ThermoSystem()
        cap = add!(sys, Body(:cap, LinearCapacitor(C); φ = V0))
        ind = add!(sys, Body(:ind, LinearInductor(L); i = 0.0))
        mid = add!(sys, Junction(:mid, Electrical()))           # zero-capacity node between R and L: a DAE
        env = add!(sys, Reservoir(:env; T))
        gnd = add!(sys, Reservoir(:gnd, Electrical() => 0.0))
        res = connect!(sys, Conductor(Resistor(R); dissipation = Into(env)), cap, mid)
        connect!(sys, InductorBranch(), mid, gnd, ind)
        cs = compile(sys)
        τ = 2L / R
        sol = simulate(cs, (0.0, 30τ))
        uf = sol.u[end]
        @test reservoir_ledger(cs, uf, env).U ≈ 0.5 * C * V0^2 rtol = 1e-6
        @test generated_entropy(cs, uf)[res] ≈ 0.5 * C * V0^2 / T rtol = 1e-6
        @test abs(rates(cs, uf)[algebraic_index(cs, mid, Electrical())]) <= 1e-12        # KCL at the junction
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9
    end

    @testset "gas spring: a frictionless piston on a gas never changes the gas entropy" begin
        A, mp = 0.01, 2.0
        sys = ThermoSystem()
        gas = add!(sys, Body(:gas, N2; T = 300.0, p = 2e5, V = 1e-3))
        pst = add!(sys, Body(:piston, InertialMass(mp); w = 0.0))
        atm = add!(sys, Reservoir(:atm; p = 1e5))
        connect!(sys, Gyrator(A; a = Volumetric(), b = Momentum()), gas, pst)       # dV = A w, dp = A p_gas
        connect!(sys, Gyrator(-A; a = Volumetric(), b = Momentum()), atm, pst)      # the atmosphere pushes back
        cs = compile(sys)
        u0 = initial_state(cs)
        sol = simulate(cs, (0.0, 1.0); saveat = 0.01)
        S0 = charges(cs, u0, gas)[1]
        # isentropic: its entropy rate is exactly zero; the dense LU of the implicit steps leaves ~1 ulp
        @test all(abs(charges(cs, u, gas)[1] - S0) <= 1e-13 * abs(S0) for u in sol.u)
        @test maximum(force(cs, u, pst, Momentum()) for u in sol.u) > 0.1         # it really oscillates
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10
    end

    @testset "DC motor: inductor, resistor, gyrator to a rotor, viscous load (junction DAE)" begin
        V, R, Lw, k, J, b = 12.0, 0.5, 1e-3, 0.05, 1e-4, 1e-5
        sys = ThermoSystem()
        src = add!(sys, Reservoir(:src, Electrical() => V))
        gnd = add!(sys, Reservoir(:gnd, Electrical() => 0.0))
        env = add!(sys, Reservoir(:env; T = 300.0))
        frame = add!(sys, Reservoir(:frame, AngularMomentum() => 0.0))
        mid = add!(sys, Junction(:mid, Electrical()))
        coil = add!(sys, Body(:coil, LinearInductor(Lw); i = 0.0))
        rotor = add!(sys, Body(:rotor, RotationalInertia(J); ω = 0.0))
        connect!(sys, Conductor(Resistor(R); dissipation = Into(env)), src, mid)
        # winding: (φ_mid, φ_gnd, i, ω): dq_mid = −i, dq_gnd = i, dΛ = φ_mid − φ_gnd − kω, dL = k i
        W = [0 0 -1 0; 0 0 1 0; 1 -1 0 -k; 0 0 k 0]
        connect!(sys, Converter(W, (Electrical(), Electrical(), FluxLinkage(), AngularMomentum())), mid, gnd, coil, rotor)
        connect!(sys, Conductor(LinearConductance(AngularMomentum(), b); dissipation = Into(env)), rotor, frame)
        cs = compile(sys)
        sol = simulate(cs, (0.0, 20.0))
        uf = sol.u[end]
        ω_ss = V / (R * b / k + k)                                    # V = iR + kω, k i = b ω
        @test force(cs, uf, rotor, AngularMomentum()) ≈ ω_ss rtol = 1e-6
        @test force(cs, uf, coil, FluxLinkage()) ≈ b * ω_ss / k rtol = 1e-6
        @test force(cs, uf, mid, Electrical()) ≈ V - R * b * ω_ss / k rtol = 1e-6
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9
    end

    @testset "ideal transformer (index-reduced DAE): RC charging sees C₁ + C₂/r²" begin
        C1, C2, r, R, V = 1e-3, 4e-3, 2.0, 10.0, 6.0
        sys = ThermoSystem()
        src = add!(sys, Reservoir(:src, Electrical() => V))
        env = add!(sys, Reservoir(:env; T = 300.0))
        c1 = add!(sys, Body(:c1, LinearCapacitor(C1); φ = 0.0))
        c2 = add!(sys, Body(:c2, LinearCapacitor(C2); φ = 0.0))
        connect!(sys, Conductor(Resistor(R); dissipation = Into(env)), src, c1)
        tr = connect!(sys, Transformer(r; a = Electrical(), b = Electrical()), c1, c2)     # φ₁ = r φ₂
        cs = compile(sys)
        τ = R * (C1 + C2 / r^2)
        sol = simulate(cs, (0.0, 10τ); saveat = τ / 2)
        for (t, u) in zip(sol.t, sol.u)
            φ1, φ2 = force(cs, u, c1, Electrical()), force(cs, u, c2, Electrical())
            @test φ1 ≈ V * (1 - exp(-t / τ)) atol = 1e-6 * V
            @test abs(φ1 - r * φ2) <= 1e-8 * V                              # the constraint holds
            @test abs(constraint_residual(cs, u, tr)) <= 1e-8 * V
        end
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-9
        @test generated_entropy(cs, sol.u[end]) |> values |> sum ≈ (0.5 * (C1 + C2 / r^2) * V^2) / 300.0 rtol = 1e-5
        # an initial state violating the constraint is rejected
        sys2 = ThermoSystem()
        d1 = add!(sys2, Body(:d1, LinearCapacitor(C1); φ = 1.0))
        d2 = add!(sys2, Body(:d2, LinearCapacitor(C2); φ = 0.0))
        connect!(sys2, Transformer(r; a = Electrical(), b = Electrical()), d1, d2)
        @test_throws ArgumentError simulate(compile(sys2), (0.0, 1.0))
    end

    @testset "switched connections and time-varying reservoirs" begin
        m, c, G, T0 = 1.0, 500.0, 10.0, 300.0
        sys = ThermoSystem()
        b = add!(sys, Body(:b, ConstantCapacitySolid(m, c); T = T0))
        hot = add!(sys, Reservoir(:hot, Thermal() => t -> 400.0 + 10.0 * (t >= 100.0)))   # a step at t = 100
        k = connect!(sys, Switched(Conductor(ThermalConductance(G)); on = t -> 20.0 <= t < 200.0), hot, b)
        cs = compile(sys)
        sol = simulate(cs, (0.0, 300.0); saveat = [10.0, 20.0, 100.0, 200.0, 300.0], tstops = [20.0, 100.0, 200.0])
        τ = m * c / G
        T(t) = t < 20 ? T0 : t < 100 ? 400 - 100exp(-(t - 20) / τ) :
               (T100 = 400 - 100exp(-80 / τ); t < 200 ? 410 - (410 - T100) * exp(-(t - 100) / τ) :
                                              410 - (410 - T100) * exp(-100 / τ))
        for (t, u) in zip(sol.t, sol.u)
            @test temperature(cs, u, b) ≈ T(t) rtol = 1e-7
        end
        @test temperature(cs, sol.u[end], b) == temperature(cs, sol.u[4], b)       # off after t = 200: frozen
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10
    end
end
