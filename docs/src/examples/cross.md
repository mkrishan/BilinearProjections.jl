# Cross

```@example cross
using BilinearProjections, LinearAlgebra

x = [1.0, 2.0, -1.0]
y = [3.0, -1.0, 4.0]

(xh, yh), info = project(CrossConstraint(), (x, y); return_info=true)

(dot(xh, yh), info.distance2, info.stationarity_residual)
```

The operation is symmetric: both vectors move by the least total squared
amount required to become orthogonal.
