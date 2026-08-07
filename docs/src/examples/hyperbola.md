# Hyperbola

```@example hyperbola
using BilinearProjections, LinearAlgebra

C = BilinearLevelSet(2.0)
x = [1.0, -2.0, 0.5]
y = [0.3, 1.5, -1.0]

(xh, yh), info = project(C, (x, y); return_info=true)
(dot(xh, yh), info.iterations, info.converged)
```

Negative levels use the same public interface:

```@example hyperbola
Cn = BilinearLevelSet(-1.0)
project(Cn, (x, y))
```
