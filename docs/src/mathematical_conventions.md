# Mathematical conventions

The package works with nonempty real `AbstractArray`s. Julia's `dot` and
`norm` define the Hilbert-space inner product and norm. Therefore matrix inputs
are interpreted using the Frobenius inner product.

The 45-degree change of variables

```math
u=\frac{x+y}{\sqrt2},\qquad
v=\frac{-x+y}{\sqrt2}
```

turns a bilinear constraint into a difference-of-squared-norms constraint:

```math
\langle x,y\rangle=\frac12(\|u\|^2-\|v\|^2).
```

This isometry is used internally for stable closed forms and for the complete
exceptional-case analysis.

## Output precision

Integer inputs are promoted to floating point. `Float32`, `Float64`, and
`BigFloat` inputs preserve their floating precision after promotion with scalar
parameters.
