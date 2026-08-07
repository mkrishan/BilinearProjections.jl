@testset "Independent low-dimensional checks" begin
    # In X = R the cross is the union of the coordinate axes.
    xh, yh = project(CrossConstraint(), ([3.0], [1.0]))
    @test xh ≈ [3.0]
    @test yh ≈ [0.0]

    # Compare the hyperbola formula against a direct one-dimensional grid over
    # y = gamma/x, including both signs of x.
    gamma = 2.0
    x0 = [3.0]
    y0 = [1.0]
    xh, yh = project(BilinearLevelSet(gamma), (x0, y0))
    objective = (xh[1] - x0[1])^2 + (yh[1] - y0[1])^2
    grid_best = Inf
    for sign in (-1.0, 1.0), t in range(0.05, 8.0; length=20_000)
        x = sign * t
        y = gamma / x
        grid_best = min(grid_best, (x - x0[1])^2 + (y - y0[1])^2)
    end
    @test objective <= grid_best + 1e-6

    # Direct grid check for the canonical paraboloid z=(u^2-v^2)/(2alpha).
    alpha = 1.5
    beta = 0.8
    H = CanonicalHyperbolicParaboloid(alpha; beta)
    input = ([1.2], [-0.7], 0.4)
    projected = project(H, input)
    objectiveH = (projected[1][1] - input[1][1])^2 +
                 (projected[2][1] - input[2][1])^2 +
                 beta^2 * (projected[3] - input[3])^2
    grid_bestH = Inf
    for u in range(-2.0, 3.0; length=500), v in range(-2.5, 2.0; length=500)
        z = (u^2 - v^2) / (2alpha)
        value = (u - input[1][1])^2 + (v - input[2][1])^2 + beta^2 * (z - input[3])^2
        grid_bestH = min(grid_bestH, value)
    end
    @test objectiveH <= grid_bestH + 2e-4
end
