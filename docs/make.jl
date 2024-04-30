using Revise
using Documenter
using Radiant

makedocs(
    sitename = "Radiant.jl",
    pages = [
        "User Interface" => [
            "Computation_Unit" => "computation_unit.md",
            "Cross_Sections" => "cross_sections.md",
            "Solvers" => [
                "Discrete_Ordinates" => "discrete_ordinates.md"
                "Solvers" => "solvers.md"
            ],
            "Geometry" => "geometry.md",
            "Interactions" => [
                "Annihilation" => "annihilation.md",
                "Auger " => "auger.md",
                "Bremsstrahlung" => "bremsstrahlung.md",
                "Compton" => "compton.md",
                "Elastic_Leptons" => "elastic_leptons.md",
                "Fluorescence" => "fluorescence.md",
                "Inelastic_Leptons" => "inelastic_leptons.md",
                "Pair_Production" => "pair_production.md",
                "Photoelectric" => "photoelectric.md",
                "Rayleigh" => "rayleigh.md",
            ],
            "Material" => "material.md",
            "Sources" => [
                "Surface_Source" => "surface_source.md",
                "Volume_Source" => "volume_source.md",
                "Fixed_Sources" => "fixed_source.md"
            ]
        ]
    ]
    )

deploydocs(
    repo = "github.com/CBienvenue/Radiant.jl"
)