@testset "BilinearLevelSet" begin
    for gamma in (2.0, -1.5)
        C = BilinearLevelSet(gamma)
        x = [1.0, -2.0, 0.5]
        y = [0.3, 1.5, -1.0]
        (xh, yh), info = project(C, (x, y); return_info=true)
        @test dot(xh, yh) ≈ gamma atol=1e-8
        @test info.converged
        @test info.multiplicity == :singleton
        @test info.stationarity_residual <= 1e-7
    end

    # Positive level and x = -y: complete nonunique family.
    Cp = BilinearLevelSet(1.0)
    x = [1.0, 0.0]
    Pp = projection_set(Cp, (x, -x))
    @test Pp isa SphereProjection
    for d in ([1.0, 0.0], [0.0, 1.0])
        u, v = Pp(d)
        @test dot(u, v) ≈ 1.0 atol=1e-10
        @test norm(u - x)^2 + norm(v + x)^2 ≈ Pp.info.distance2 atol=1e-10
    end

    # Negative level and x = y: symmetric nonunique family.
    Cn = BilinearLevelSet(-1.0)
    Pn = projection_set(Cn, (x, x))
    @test Pn isa SphereProjection
    u, v = Pn([0.0, 1.0])
    @test dot(u, v) ≈ -1.0 atol=1e-10

    # Level zero delegates to the cross.
    Cz = BilinearLevelSet(0.0)
    u, v = project(Cz, ([1.0, 2.0], [3.0, 4.0]))
    @test abs(dot(u, v)) <= 1e-10

    # Matrices use the Frobenius inner product.
    A = [1.0 2.0; 3.0 4.0]
    B = [0.5 -1.0; 2.0 0.25]
    Ah, Bh = project(BilinearLevelSet(3.0), (A, B))
    @test dot(Ah, Bh) ≈ 3.0 atol=1e-8
end
