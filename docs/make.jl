using Documenter, DocStringExtensions, DynamicEnergyBudgets

makedocs(
    modules = [DynamicEnergyBudgets],
    sitename = "DynamicEnergyBudgets.jl",
    pages = Any[
        "Home" => "index.md",
    ],
    clean = false,
)

deploydocs(
    repo = "github.com/rafaqz/DynamicEnergyBudgets.jl.git",
)
