# BilinearProjections.jl

BilinearProjections.jl provides exact projection operators for several
nonconvex bilinear constraint sets in real Hilbert spaces.

The package is designed around two principles:

1. implement the complete mathematical branch structure, including exceptional
   cases;
2. preserve genuinely set-valued projections through `SphereProjection`.

```@example intro
using BilinearProjections, LinearAlgebra

C = CrossConstraint()
x = [1.0, 2.0]
y = [3.0, 1.0]
(xh, yh), info = project(C, (x, y); return_info=true)
(dot(xh, yh), info.distance2, info.branch)
```

## Implemented sets

```math
\begin{aligned}
C_0 &= \{(x,y):\langle x,y\rangle=0\},\\
C_\gamma &= \{(x,y):\langle x,y\rangle=\gamma\},\\
\widetilde C_\alpha
&= \{(u,v,z):\|u\|^2-\|v\|^2=2\alpha z\},\\
C_\alpha
&= \{(x,y,z):\langle x,y\rangle=\alpha z\}.
\end{aligned}
```

All paraboloid projections use the weighted metric

```math
\|(x,y,z)\|_\beta^2=\|x\|^2+\|y\|^2+\beta^2|z|^2.
```
