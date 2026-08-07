using BilinearProjections
using LinearAlgebra

H = CanonicalHyperbolicParaboloid(5.0; beta=1.0)
(u, v, z), info = project(H, ([2.0], [-3.0], 4.0); return_info=true)

println("projection = ", (u, v, z))
println("constraint residual = ", norm(u)^2 - norm(v)^2 - 10z)
println("branch = ", info.branch)
