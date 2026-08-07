@testset "CrossConstraint" begin
    C = CrossConstraint()

    x = [1.0, 2.0, -1.0]
    y = [3.0, -1.0, 4.0]
    (xh, yh), info = project(C, (x, y); return_info=true)
    @test abs(dot(xh, yh)) <= 1e-10
    @test info.multiplicity == :singleton
    @test info.feasibility_residual <= 1e-10
    @test info.stationarity_residual <= 1e-8
    @test info.distance2 ≈ norm(xh - x)^2 + norm(yh - y)^2

    # Idempotence on a feasible pair.
    xf = [1.0, 0.0]
    yf = [0.0, 2.0]
    xfp, yfp = project(C, (xf, yf))
    @test xfp ≈ xf
    @test yfp ≈ yf

    # Nonunique branch: x = y.
    z = [1.0, 0.0]
    P = projection_set(C, (z, z))
    @test P isa SphereProjection
    @test P.info.multiplicity == :continuum
    p1 = P([1.0, 0.0])
    p2 = P([0.0, 1.0])
    @test abs(dot(p1[1], p1[2])) <= 1e-12
    @test abs(dot(p2[1], p2[2])) <= 1e-12
    @test norm(p1[1] - z)^2 + norm(p1[2] - z)^2 ≈ P.info.distance2
    @test norm(p2[1] - z)^2 + norm(p2[2] - z)^2 ≈ P.info.distance2

    # Orthogonal equivariance under a fixed planar rotation.
    theta = 0.37
    Q = [cos(theta) -sin(theta); sin(theta) cos(theta)]
    a = [1.2, -0.7]
    b = [0.4, 2.1]
    ah, bh = project(C, (a, b))
    qah, qbh = project(C, (Q * a, Q * b))
    @test qah ≈ Q * ah atol=1e-9
    @test qbh ≈ Q * bh atol=1e-9

    # Deterministic property sweep without a random-number dependency.
    for k in 1:20
        a = [sin(k + j / 3) for j in 1:7]
        b = [cos(k / 2 + 2j / 5) for j in 1:7]
        ah, bh = project(C, (a, b))
        @test abs(dot(ah, bh)) <= 1e-8 * (1 + norm(ah) * norm(bh))
    end
end
