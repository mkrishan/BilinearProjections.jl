using Pkg
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))

using Documenter
using BilinearProjections

makedocs(
    sitename = "BilinearProjections.jl",
    modules = [BilinearProjections],
    format = Documenter.HTML(
        canonical = "https://mkrishan.github.io/BilinearProjections.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Mathematical conventions" => "mathematical_conventions.md",
        "Nonunique projections" => "nonunique_projections.md",
        "Numerical design" => "numerical_design.md",
        "Examples" => [
            "Cross" => "examples/cross.md",
            "Hyperbola" => "examples/hyperbola.md",
            "Hyperbolic paraboloid" => "examples/paraboloid.md",
        ],
        "API" => "api.md",
        "References" => "references.md",
    ],
    warnonly = [:missing_docs],
)

deploydocs(
    repo = "github.com/mkrishan/BilinearProjections.jl.git",
    devbranch = "main",
)
