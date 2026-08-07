@testset "CanonicalHyperbolicParaboloid" begin
    H = CanonicalHyperbolicParaboloid(5.0; beta=1.0)

    # Numerical examples from the paper.
    (u, v, z), info = project(H, ([2.0], [-3.0], 4.0); return_info=true)
    @test u[1] ≈ 4.20311 atol=2e-5
    @test v[1] ≈ -1.96830 atol=2e-5
    @test z ≈ 1.37919 atol=2e-5
    @test norm(u)^2 - norm(v)^2 ≈ 10z atol=1e-7
    @test info.converged

    u, v, z = project(H, ([0.0], [-3.0], 3.0))
    @test u[1] ≈ 0.0 atol=1e-12
    @test v[1] ≈ -1.80187 atol=2e-5
    @test z ≈ -0.32467 atol=2e-5

    P = projection_set(H, ([0.0], [sqrt(32.0)], 6.0))
    @test P isa SphereProjection
    pplus = P([1.0])
    pminus = P([-1.0])
    @test pplus[1][1] ≈ sqrt(18.0) atol=1e-10
    @test pminus[1][1] ≈ -sqrt(18.0) atol=1e-10
    @test pplus[2][1] ≈ sqrt(8.0) atol=1e-10
    @test pplus[3] ≈ 1.0 atol=1e-12

    P0 = projection_set(H, ([0.0], [0.0], 6.0))
    @test P0 isa SphereProjection
    p = P0([1.0])
    @test p[1][1] ≈ sqrt(10.0) atol=1e-10
    @test p[2][1] == 0.0
    @test p[3] == 1.0

    pzero = project(H, ([0.0], [0.0], 4.0))
    @test pzero[1] == [0.0]
    @test pzero[2] == [0.0]
    @test pzero[3] == 0.0
end

@testset "HyperbolicParaboloid" begin
    H = HyperbolicParaboloid(2.0; beta=0.75)
    x = [1.0, -0.5, 2.0]
    y = [-1.0, 1.5, 0.25]
    z = 0.8
    (xh, yh, zh), info = project(H, (x, y, z); return_info=true)
    @test dot(xh, yh) ≈ 2zh atol=1e-8
    @test info.feasibility_residual <= 1e-8
    @test info.stationarity_residual <= 1e-7

    # Exceptional x = -y branch.
    a = [1.0, 0.0]
    P = projection_set(H, (a, -a, 4.0))
    @test P isa SphereProjection
    q = P([0.0, 1.0])
    @test dot(q[1], q[2]) ≈ 2q[3] atol=1e-9
end
