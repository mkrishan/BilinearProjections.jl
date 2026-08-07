# Benchmarks

Run from the repository root:

```bash
julia --project benchmark/runbenchmarks.jl
```

The benchmark intentionally uses only Julia's standard timing facilities, so
it introduces no package dependency. For publication-grade benchmarking,
repeat the measurements in a controlled environment and consider wrapping the
same calls with BenchmarkTools.jl.
