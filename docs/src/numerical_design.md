# Numerical design

## Cross

After the 45-degree rotation, projecting onto the cross becomes

```math
\min_{u,v}\ \|u-u_0\|^2+\|v-v_0\|^2
\quad\text{subject to}\quad \|u\|=\|v\|.
```

When both canonical inputs are nonzero, the optimal common radius is

```math
r=\frac{\|u_0\|+\|v_0\|}{2}.
```

This avoids cancellation in a direct multiplier formula.

## Hyperbola and paraboloid

The generic cases reduce to continuous strictly monotone scalar equations on
`(-1,1)`. The implementation uses bracketed bisection instead of unrestricted
quartic or quintic root extraction. This enforces the mathematically relevant
multiplier interval and gives predictable convergence.

## Diagnostics

`ProjectionInfo` records:

- the recovered Lagrange multiplier;
- squared projection distance;
- feasibility residual;
- stationarity residual;
- branch and multiplicity;
- root-solver iterations and convergence status.
