# Contributing

Thank you for contributing to BilinearProjections.jl.

## Development setup

```julia
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

## Requirements for a new projector

A pull request adding a constraint family must include:

1. a precise mathematical definition of the set and metric;
2. a complete treatment of singleton and non-singleton branches;
3. feasibility, KKT/stationarity, idempotence, and symmetry tests;
4. an independent low-dimensional numerical check where practical;
5. API documentation and one executable example;
6. a primary mathematical reference or a self-contained derivation.

## Numerical standards

- Prefer bracketed scalar root solvers when the theory provides a monotone root.
- Do not select an arbitrary polynomial root without proving it is the relevant multiplier.
- Avoid exact equality tests except where the formula genuinely has an exact algebraic branch.
- Report convergence and residual metadata.
- Preserve full set-valued answers when projections are nonunique.

## Style

- Use descriptive mathematical names.
- Keep exported APIs small.
- Add docstrings to exported types and methods.
- Run `Aqua.test_all(BilinearProjections)` through the standard test suite.

## Changes

Add user-visible changes to `CHANGELOG.md` under an `Unreleased` heading.
