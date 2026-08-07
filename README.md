# BilinearProjections.jl

[![CI](https://github.com/mkrishan/BilinearProjections.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/mkrishan/BilinearProjections.jl/actions/workflows/CI.yml)
[![Docs](https://github.com/mkrishan/BilinearProjections.jl/actions/workflows/Docs.yml/badge.svg)](https://github.com/mkrishan/BilinearProjections.jl/actions/workflows/Docs.yml)
[![codecov](https://codecov.io/gh/mkrishan/BilinearProjections.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mkrishan/BilinearProjections.jl)

**Exact, set-valued-aware projections onto nonconvex bilinear constraints in Julia.**

The initial release implements projections onto

- the cross `dot(x,y) = 0`;
- bilinear level sets `dot(x,y) = gamma`;
- canonical rectangular hyperbolic paraboloids
  `norm(u)^2 - norm(v)^2 = 2alpha*z`;
- bilinear-coordinate hyperbolic paraboloids `dot(x,y) = alpha*z`.

The implementation follows the complete Hilbert-space projection formulas of
Bauschke, Krishan Lal, and Wang. Nonunique projections are represented as
parameterized projection sets rather than silently collapsed to one point.

## Installation

Until the package is registered, install it directly from GitHub:

```julia
pkg> add https://github.com/mkrishan/BilinearProjections.jl
```

For a local checkout:

```julia
pkg> develop path/to/BilinearProjections.jl
```

## Quick start

### Cross: symmetric minimal-change orthogonalization

```julia
using BilinearProjections, LinearAlgebra

C = CrossConstraint()
x = [1.0, 2.0]
y = [3.0, 1.0]

(xhat, yhat), info = project(C, (x, y); return_info=true)

@assert isapprox(dot(xhat, yhat), 0.0; atol=1e-10)
info.distance2
info.multiplier
```

This solves

```math
\min_{u,v}\; \|u-x\|^2+\|v-y\|^2
\quad\text{subject to}\quad \langle u,v\rangle=0.
```

### Hyperbola / bilinear level set

```julia
Cgamma = BilinearLevelSet(2.0)
(xhat, yhat), info = project(Cgamma, (x, y); return_info=true)

@assert isapprox(dot(xhat, yhat), 2.0; atol=1e-9)
```

### Hyperbolic paraboloid

```julia
H = HyperbolicParaboloid(5.0; beta=1.0)
z = 4.0

(xhat, yhat, zhat), info = project(H, (x, y, z); return_info=true)
@assert isapprox(dot(xhat, yhat), 5.0 * zhat; atol=1e-9)
```

The canonical coordinate representation is also available:

```julia
Hcanonical = CanonicalHyperbolicParaboloid(5.0; beta=1.0)
(uhat, vhat, zhat) = project(Hcanonical, ([2.0], [-3.0], 4.0))
@assert isapprox(norm(uhat)^2 - norm(vhat)^2, 10.0 * zhat; atol=1e-9)
```

## Nonunique projections

Some degenerate inputs have infinitely many nearest points. Use
`projection_set` to retain the full solution family:

```julia
P = projection_set(CrossConstraint(), ([1.0, 0.0], [1.0, 0.0]))

P isa SphereProjection
P.default
P([0.0, 1.0])                  # select a unit-sphere parameter
project(CrossConstraint(),
        ([1.0, 0.0], [1.0, 0.0]);
        selector=UnitDirection([0.0, 1.0]))
```

`project` returns a deterministic default selection. `projection_set` exposes
both the default and the complete unit-sphere parameterization.

## API

```julia
projection_set(C, point; atol=nothing, rtol=nothing, maxiters=200)
project(C, point; selector=DefaultSelector(), return_info=false, ...)
project!(output, C, point; ...)
residual(C, point)
isfeasible(C, point; atol=0, rtol=sqrt(eps()))
distance2(C, point; ...)
```

For pair constraints, `project!` accepts `(xout, yout)`. For paraboloids it
accepts `(xout, yout, Ref(zout))`.

## Numerical design

- Cross projections use a stable 45-degree coordinate rotation and a closed
  form in the canonical variables.
- Hyperbola and paraboloid projections solve the mathematically distinguished
  monotone scalar equation on `(-1,1)` by bracketed bisection.
- Exceptional branches are implemented directly from the complete formulas.
- Inputs may be vectors, matrices, views, or other real `AbstractArray`s. The
  Frobenius inner product is used for matrix inputs.
- `Float32`, `Float64`, and `BigFloat` are supported.

## Development

```julia
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

See `CONTRIBUTING.md` for testing and style expectations.

## License

MIT. See `LICENSE`.
