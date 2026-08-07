using BilinearProjections
using LinearAlgebra

C = CrossConstraint()
x = [1.0, 2.0]
y = [3.0, 1.0]

(xhat, yhat), info = project(C, (x, y); return_info=true)

println("dot(xhat, yhat) = ", dot(xhat, yhat))
println("distance squared = ", info.distance2)
println("branch = ", info.branch)
