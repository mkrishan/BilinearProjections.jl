using Test
using Aqua
using BilinearProjections
using LinearAlgebra

include("test_cross.jl")
include("test_level_set.jl")
include("test_paraboloid.jl")
include("test_interface.jl")
include("test_global_checks.jl")

@testset "Aqua" begin
    Aqua.test_all(BilinearProjections)
end
