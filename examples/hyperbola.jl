using BilinearProjections
using LinearAlgebra

C = BilinearLevelSet(2.0)
x = [1.0, -2.0, 0.5]
y = [0.3, 1.5, -1.0]

(xhat, yhat), info = project(C, (x, y); return_info=true)

println("dot(xhat, yhat) = ", dot(xhat, yhat))
println("root iterations = ", info.iterations)
println("converged = ", info.converged)
