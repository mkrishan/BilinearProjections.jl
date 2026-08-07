# Nonunique projections

A projection onto a closed nonconvex set can contain more than one nearest
point. BilinearProjections.jl represents this explicitly.

```@example families
using BilinearProjections, LinearAlgebra

C = CrossConstraint()
x = [1.0, 0.0]
P = projection_set(C, (x, x))

P isa SphereProjection
p1 = P([1.0, 0.0])
p2 = P([0.0, 1.0])
(dot(p1[1], p1[2]), dot(p2[1], p2[2]))
```

A `SphereProjection` contains:

- `generator`: the full unit-sphere parameterization;
- `default`: a deterministic selection;
- `extras`: any separately represented points, when needed;
- `description`: a human-readable description;
- `info`: common projection metadata.

Use `UnitDirection(direction)` with `project` to request a particular member:

```@example families
project(C, (x, x); selector=UnitDirection([0.0, 1.0]))
```

Directions are normalized internally. A zero direction is rejected.
