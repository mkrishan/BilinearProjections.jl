using BilinearProjections
using LinearAlgebra

x = [1.0, 0.0]
P = projection_set(CrossConstraint(), (x, x))

println(P)
println("default projection = ", P.default)
println("direction e1 = ", P([1.0, 0.0]))
println("direction e2 = ", P([0.0, 1.0]))
