using Pkg

root = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(root)
Pkg.instantiate()
Pkg.test()

for example in ("cross.jl", "hyperbola.jl", "paraboloid.jl", "nonunique.jl")
    println("\nRunning example: ", example)
    include(joinpath(root, "examples", example))
end

println("\nBuilding documentation")
Pkg.activate(joinpath(root, "docs"))
Pkg.instantiate()
include(joinpath(root, "docs", "make.jl"))
