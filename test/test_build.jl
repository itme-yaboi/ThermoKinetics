@testset "build and validation" begin
    @testset "domain declaration" begin
        d = Thermal()
        @test charge_symbol(d) == :S && charge_unit(d) == "J/K"
        @test force_symbol(d) == :T && force_unit(d) == "K"
        @test producible(d)
    end

    @testset "ConstantCapacitySolid state function" begin
        mat = ConstantCapacitySolid(2.0, 450.0)            # 2 kg, 450 J/(kg K)
        for T in (50.0, 300.0, 1234.5)
            E = charges_at(mat; T)
            @test temperature(mat, E) ≈ T rtol = 1e-14       # P = ∂U/∂E inverts charges_at
            @test internal_energy(mat, E) ≈ 2.0 * 450.0 * T rtol = 1e-14
            A = eos_matrix(mat, E)
            @test A == A'                                     # A4 (trivially 1×1 here)
            @test 1 / A[1, 1] ≈ 2.0 * 450.0 / T rtol = 1e-12  # capacity dS/dT = m c / T
        end
        @test_throws ArgumentError ConstantCapacitySolid(-1.0, 450.0)
        @test_throws ArgumentError ConstantCapacitySolid(1.0, 0.0)
        @test_throws ArgumentError charges_at(mat; T = -5.0)
    end

    @testset "graph construction errors" begin
        sys = ThermoSystem()
        a = add!(sys, Body(:a, ConstantCapacitySolid(1.0, 500.0); T = 300.0))
        r = add!(sys, Reservoir(:r; T = 400.0))
        @test_throws ArgumentError add!(sys, Body(:a, ConstantCapacitySolid(1.0, 500.0); T = 300.0))
        @test_throws ArgumentError connect!(sys, Conductor(ThermalConductance(1.0)), a, :nope)
        @test_throws ArgumentError connect!(sys, Conductor(ThermalConductance(1.0)), a, a)
        c = connect!(sys, Conductor(ThermalConductance(1.0)), a, r; name = :c1)
        @test c == :c1
        @test_throws ArgumentError connect!(sys, Conductor(ThermalConductance(1.0)), c, a)  # element is not a node
        @test_throws ArgumentError Reservoir(:bad; T = 0.0)
    end

    @testset "A5: conductance sign is a type invariant" begin
        @test_throws ArgumentError ThermalConductance(-1e-3)
        @test_throws ArgumentError ThermalConductance(NaN)
        @test ThermalConductance(0.0) isa ThermalConductance                  # an open contact is allowed
    end

    @testset "parallel conductors are allowed" begin
        sys = ThermoSystem()
        a = add!(sys, Body(:a, ConstantCapacitySolid(1.0, 500.0); T = 300.0))
        b = add!(sys, Body(:b, ConstantCapacitySolid(1.0, 500.0); T = 400.0))
        connect!(sys, Conductor(ThermalConductance(1.0)), a, b)
        connect!(sys, Conductor(ThermalConductance(2.0)), a, b)
        cs = compile(sys)
        @test length(conductor_names(cs)) == 2
        @test heat_flow(cs, initial_state(cs), conductor_names(cs)[2]) ≈ 2.0 * (300.0 - 400.0)
    end

    @testset "SI only: no technical-system constants in src" begin
        srcdir = joinpath(pkgdir(ThermoKinetics), "src")
        for (root, _, files) in walkdir(srcdir), f in files
            txt = read(joinpath(root, f), String)
            @test !occursin(r"kcal|kgf|\b427\b"i, txt)
        end
    end
end

@testset "slice 2 construction and validation" begin
    N2 = gas_n2()
    solid = ConstantCapacitySolid(1.0, 500.0)

    @testset "reservoirs" begin
        r = Reservoir(:atm; T = 300.0, p = 101325.0)
        @test domains(r) == (Thermal(), Volumetric())
        @test ThermoKinetics.reservoir_force(r, Volumetric()) == -101325.0     # stored as the force −p
        @test domains(Reservoir(:hp; p = 2e5)) == (Volumetric(),)
        @test_throws ArgumentError Reservoir(:none)
        @test_throws ArgumentError Reservoir(:bad; p = -1.0)
        raw = Reservoir(:raw, Thermal() => 300.0, Volumetric() => -1e5)
        @test domains(raw) == (Thermal(), Volumetric())
    end

    @testset "dissipation policy is required for non-thermal conductors" begin
        @test_throws ArgumentError Conductor(VolumeConductance(1e-9))                 # no policy given
        @test_throws ArgumentError Conductor(ThermalConductance(1.0); dissipation = Split(0.5))
        @test_throws ArgumentError VolumeConductance(-1e-9)
        @test_throws ArgumentError Split(1.5)
        @test_throws ArgumentError InternalNode()                                     # no legs
        @test_throws ArgumentError InternalNode(G_h = 1.0)                            # leg to h without h
        @test_throws ArgumentError InternalNode(G_a = -1.0)
    end

    @testset "graph validation of non-thermal conductors" begin
        sys = ThermoSystem()
        g1 = add!(sys, Body(:g1, N2; T = 300.0, p = 1e5, m = 1e-3))
        g2 = add!(sys, Body(:g2, N2; T = 300.0, p = 2e5, m = 1e-3))
        s = add!(sys, Body(:s, solid; T = 300.0))
        hp = add!(sys, Reservoir(:hp; p = 3e5))                                       # no thermal DOF
        sink = add!(sys, Reservoir(:sink; T = 300.0))
        # Into(h) adds a third port; h must be a distinct node
        c = connect!(sys, Conductor(VolumeConductance(1e-9); dissipation = Into(sink)), g1, g2)
        @test ThermoKinetics.port_node(sys, c, :h) == sink
        @test_throws ArgumentError connect!(sys, Conductor(VolumeConductance(1e-9); dissipation = Into(g1)), g1, g2)
        @test_throws ArgumentError connect!(sys, Conductor(VolumeConductance(1e-9); dissipation = Into(:nope)), g1, g2)
        @test compile(sys) isa CompiledSystem

        # internal isolation (n ≤ m): a solid has no volumetric DOF
        sys2 = deepcopy(sys)
        connect!(sys2, Conductor(VolumeConductance(1e-9); dissipation = Split(0.5)), g1, s)
        @test_throws ArgumentError compile(sys2)
        # dissipation may not be routed to a node without a thermal DOF
        sys3 = deepcopy(sys)
        connect!(sys3, Conductor(VolumeConductance(1e-9); dissipation = Split(0.5)), g1, hp)
        @test_throws ArgumentError compile(sys3)
        sys4 = deepcopy(sys)
        connect!(sys4, Conductor(VolumeConductance(1e-9); dissipation = InternalNode(G_a = 1.0, G_b = 1.0)), g1, hp)
        @test_throws ArgumentError compile(sys4)
        # … but zero weight / no leg on it is fine
        sys5 = deepcopy(sys)
        connect!(sys5, Conductor(VolumeConductance(1e-9); dissipation = InternalNode(G_a = 1.0)), g1, hp)
        connect!(sys5, Conductor(VolumeConductance(1e-9); dissipation = Split(0.0)), g2, hp)
        @test compile(sys5) isa CompiledSystem
    end

    @testset "every exported name has a docstring (the HTML reference lists them all)" begin
        undocumented = [n for n in names(ThermoKinetics) if n != :ThermoKinetics && !Base.Docs.hasdoc(ThermoKinetics, n)]
        @test isempty(undocumented)
    end
end
