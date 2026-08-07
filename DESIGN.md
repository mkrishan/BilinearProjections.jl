# Package design

## Public layers

1. **Constraint objects** define the mathematical set and metric.
2. **`projection_set`** returns the complete metric projection.
3. **`project`** applies a deterministic or user-specified selector.
4. **`ProjectionInfo`** exposes numerical and variational diagnostics.

This separation is intentional. A conventional single-valued function cannot
faithfully represent the exceptional cases of these nonconvex projectors.

## Adding another constraint family

Implement the following methods:

```julia
residual(C::NewConstraint, point)
BilinearProjections._residual_scale(C::NewConstraint, point)
projection_set(C::NewConstraint, point; kwargs...)
```

Reuse the internal helpers when the set is isometrically equivalent to a
canonical difference-of-squared-norms constraint. Return an internal raw
singleton or raw sphere family, map it to the public coordinates, and construct
`ProjectionInfo` from the final point.

## Intended future extensions

- weighted bilinear level sets;
- generalized constraints `dot(A*x, B*y) = gamma` under suitable structure;
- interval constraints `lower <= dot(x,y) <= upper`;
- matrix product constraints;
- batched CPU and GPU kernels;
- ChainRules support on smooth singleton strata;
- integration adapters for proximal and feasibility-algorithm packages.

Automatic differentiation must not assign an ordinary derivative across a
multivalued or branch-singular stratum without an explicit selection rule.
