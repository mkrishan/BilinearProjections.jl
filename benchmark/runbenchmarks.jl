using BilinearProjections
using LinearAlgebra

function benchmark_projector(name, C, point; samples=20, repeats=100)
    # Warm-up
    project(C, point)
    times = Float64[]
    for _ in 1:samples
        elapsed = @elapsed begin
            for _ in 1:repeats
                project(C, point)
            end
        end
        push!(times, elapsed / repeats)
    end
    sort!(times)
    median = times[(length(times) + 1) ÷ 2]
    println(rpad(name, 28), " median seconds/call = ", median)
end

n = 10_000
x = [sin(i / 17) for i in 1:n]
y = [cos(i / 23) for i in 1:n]

benchmark_projector("CrossConstraint", CrossConstraint(), (x, y))
benchmark_projector("BilinearLevelSet", BilinearLevelSet(2.0), (x, y))
benchmark_projector("HyperbolicParaboloid", HyperbolicParaboloid(2.0), (x, y, 0.5))
