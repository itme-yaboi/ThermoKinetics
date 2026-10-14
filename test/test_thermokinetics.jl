@testset "08 §C thermokinetics" begin
    # a steel rod between two reservoirs (SI)
    λ, ρ, c = 45.0, 7850.0, 470.0
    L, F, n = 0.5, 1e-3, 25
    T1, T2 = 500.0, 300.0

    function rod_between_reservoirs()
        sys = ThermoSystem()
        r1 = add!(sys, Reservoir(:hot; T = T1))
        r2 = add!(sys, Reservoir(:cold; T = T2))
        rd = rod!(sys, :rod; n, L, area = F, λ, ρ, c, T0 = 400.0, left = r1, right = r2)
        return sys, r1, r2, rd
    end

    @testset "C1 steady conduction" begin
        sys, r1, r2, rd = rod_between_reservoirs()
        cs = compile(sys)
        u = steady_state(cs)
        Q = λ * F * (T1 - T2) / L
        # linear T(x) (exact for the cell-centred finite-volume chain)
        for (cell, x) in zip(rd.cells, rd.x)
            @test temperature(cs, u, cell) ≈ T1 + (T2 - T1) * x / L atol = 1e-8
        end
        # Q̇ = λFΔT/Δx through every face, including the boundary half-cells
        for face in rd.all_faces
            @test flux(cs, u, face).Q ≈ Q rtol = 1e-9
        end
        # total production Q̇(1/T₂ − 1/T₁), attributed over the faces
        Sgen = sum(values(entropy_production(cs, u)))
        @test Sgen ≈ Q * (1 / T2 - 1 / T1) rtol = 1e-9
        # the transient run reaches the same state
        τ = ρ * c * L^2 / λ
        sol = simulate(cs, (0.0, 3τ))
        for cell in rd.cells
            @test temperature(cs, sol.u[end], cell) ≈ temperature(cs, u, cell) atol = 1e-5
        end
        cl = closure(cs, sol)
        @test cl.energy_rel <= 1e-6 && cl.entropy_rel <= 1e-10
    end

    @testset "C2 entropy flux growth" begin
        sys, r1, r2, rd = rod_between_reservoirs()
        cs = compile(sys)
        u = steady_state(cs)
        Q = λ * F * (T1 - T2) / L
        Js_prev = nothing
        for face in rd.all_faces
            f = flux(cs, u, face)
            # entropy leaves the upstream node at its T and arrives at the downstream node at its T
            @test f.S_out_of_a * f.T_a ≈ Q rtol = 1e-9
            @test f.S_into_b * f.T_b ≈ Q rtol = 1e-9
            @test f.S_into_b / f.S_out_of_a ≈ f.T_a / f.T_b rtol = 1e-12
            # J_s grows across every face (the generated entropy joins the flux) …
            @test f.S_into_b > f.S_out_of_a
            # … and is continuous through every cell at steady state (what enters a cell leaves it)
            Js_prev === nothing || @test f.S_out_of_a ≈ Js_prev rtol = 1e-9
            Js_prev = f.S_into_b
        end
        # whole rod: J_s,out / J_s,in = T₁/T₂
        rates_ = rates(cs, u)
        Js_in = -reservoir_ledger(cs, rates_, r1).S
        Js_out = reservoir_ledger(cs, rates_, r2).S
        @test Js_out / Js_in ≈ T1 / T2 rtol = 1e-9
        # energy flux is conserved: the hot reservoir gives what the cold one gets
        @test -reservoir_ledger(cs, rates_, r1).U ≈ reservoir_ledger(cs, rates_, r2).U rtol = 1e-9
    end

    @testset "C8 irreversibility efficiency" begin
        Ta, Tb, G = 900.0, 300.0, 2.0
        sys = ThermoSystem()
        a = add!(sys, Reservoir(:a; T = Ta))
        b = add!(sys, Reservoir(:b; T = Tb))
        k = connect!(sys, Conductor(ThermalConductance(G)), a, b)
        cs = compile(sys)
        du = rates(cs, initial_state(cs))
        S1 = -reservoir_ledger(cs, du, a).S          # entropy leaving the hot reservoir
        Q1 = -reservoir_ledger(cs, du, a).U          # P₁ΔE
        Q2 = Tb * S1                                 # P₂ΔE: the transported charge at the lower force
        @test Q1 ≈ Ta * S1 rtol = 1e-14
        @test Q2 / Q1 ≈ Tb / Ta rtol = 1e-14         # (216)
        Qd = Q1 - Q2                                 # ΔP ΔE, dissipated at Tb
        @test entropy_production(cs, initial_state(cs))[k] ≈ Qd / Tb rtol = 1e-12
        # a system with no bodies still integrates its ledgers
        sol = simulate(cs, (0.0, 10.0))
        @test generated_entropy(cs, sol.u[end])[k] ≈ 10 * Qd / Tb rtol = 1e-10
    end
end
