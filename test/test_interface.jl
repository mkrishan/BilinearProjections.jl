@testset "Interface" begin
    C = CrossConstraint()
    x = [1.0, 2.0]
    y = [3.0, 4.0]

    xout = similar(x)
    yout = similar(y)
    info = project!((xout, yout), C, (x, y))
    @test abs(dot(xout, yout)) <= 1e-10
    @test info.distance2 ≈ distance2(C, (x, y))
    @test isfeasible(C, (xout, yout); atol=1e-10)

    H = HyperbolicParaboloid(1.0)
    zout = Ref(0.0)
    infoH = project!((xout, yout, zout), H, (x, y, 1.0))
    @test dot(xout, yout) ≈ zout[] atol=1e-8
    @test infoH.converged

    @test_throws DimensionMismatch project(C, ([1.0], [1.0, 2.0]))
    @test_throws ArgumentError HyperbolicParaboloid(0.0)
    @test_throws ArgumentError HyperbolicParaboloid(1.0; beta=0.0)
    @test_throws ArgumentError project(C, (ComplexF64[1], ComplexF64[1]))
end
