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
        srcdir = joinpath(pkgdir(VeinikThermo), "src")
        for (root, _, files) in walkdir(srcdir), f in files
            txt = read(joinpath(root, f), String)
            @test !occursin(r"kcal|kgf|\b427\b"i, txt)
        end
    end
end
