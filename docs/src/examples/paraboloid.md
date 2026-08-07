# Hyperbolic paraboloid

## Canonical coordinates

```@example paraboloid
using BilinearProjections, LinearAlgebra

H = CanonicalHyperbolicParaboloid(5.0; beta=1.0)
(u, v, z), info = project(H, ([2.0], [-3.0], 4.0); return_info=true)
(norm(u)^2 - norm(v)^2, 10z, info.branch)
```

## Bilinear coordinates

```@example paraboloid
B = HyperbolicParaboloid(2.0; beta=0.75)
x = [1.0, -0.5, 2.0]
y = [-1.0, 1.5, 0.25]
(xh, yh, zh) = project(B, (x, y, 0.8))
(dot(xh, yh), 2zh)
```
