using Documenter, DocStringExtensions, DynamicEnergyBudgets

makedocs(
    modules = [DynamicEnergyBudgets],
    doctest = false,
    clean = false,
    sitename = "DynamicEnergyBudgets.jl",
    format = :html,
    pages = Any[
        "Home" => "index.md",
    ]
)

deploydocs(
    repo = "github.com/rafaqz/DynamicEnergyBudgets.jl.git",
    osname = "linux",
    julia = "0.6",
    target = "build",
    deps = nothing,
    make = nothing
)
