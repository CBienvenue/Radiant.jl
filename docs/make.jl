using Revise
using Documenter
using Radiant
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

makedocs(
    format = Documenter.HTML(assets=String["assets/citations.css"],),
    plugins=[bib],
    sitename = "Radiant",
    pages = [
        "Installation Guide" => ["Quick Start" => "quick_start.md"],

        "Theory and Methodology" => [
            "1 Radiation Transport" => [
                "1.1 The Transport Equation" => "theory_transport_equation.md",
                "1.2 Energy Discretization"  => "theory_energy_discretization.md",
                "1.3 Angular Discretization" => "theory_angular_discretization.md",
                "1.4 Spatial Discretization" => "theory_spatial_discretization.md",
                "1.5 Closure Relations"      => "theory_closure_relations.md",
                "1.6 Sweeping Process"       => "theory_sweeping_process.md"
            ],
            "2 Cross-Sections" => [
                "2.1 Photons Interactions" => [
                    "2.1.1 Rayleigh Scattering"  => "theory_rayleigh.md",  
                    "2.1.2 Compton Scattering"   => "theory_compton.md",
                    "2.1.3 Photoelectric Effect" => "theory_photoelectric.md",
                    "2.1.4 Pair Production"      => "theory_pair_production.md",
                ],
                "2.2 Electrons and Positrons Interactions" => [
                    "2.2.1 Impact Ionization"  => "theory_impact_ionization.md",
                    "2.2.2 Elastic Scattering" => "theory_elastic_scattering.md",
                    "2.2.3 Bremsstrahlung"     => "theory_bremsstrahlung.md",
                    "2.2.4 Annihilation"       => "theory_annihilation.md",
                ],
                "2.3 Atomic Relaxation" => "theory_atomic_relaxation.md",
            ],
            "3 Mathematical Tools" => [
                "3.1 Polynomials"    => "theory_polynomials.md",
                "3.2 Quadratures"    => "theory_quadrature.md",
                "3.3 Interpolations" => "theory_interpolation.md"
            ],
            "References" => "theory_references.md"
        ],

        "User's Guide" => [
            "1 Getting Started with Radiant" => "user_guide_getting_started.md",
            "2 Materials"                    => "user_guide_materials.md",
            "3 Particles"                    => "user_guide_particles.md",
            "4 Cross-Sections"               => "user_guide_cross_sections.md",
            "5 Geometry"                     => "user_guide_geometry.md",
            "6 Solvers"                      => "user_guide_solvers.md",
            "7 Fixed External Sources"       => "user_guide_fixed_external_sources.md",
            "8 Transport Calculations"       => "user_guide_transport_calculations.md",
        ],

        #"Examples" => [],

        "API" => [
            "Computation_Unit" => "api_computation_unit.md",
            "Cross_Sections"   => "api_cross_sections.md",
            "Solvers" => [
                "Discrete_Ordinates" => "api_discrete_ordinates.md"
                "Solvers"            => "api_solvers.md"
            ],
            "Geometry"     => "api_geometry.md",
            "Interactions" => [
                "Annihilation"      => "api_annihilation.md",
                "Auger "            => "api_auger.md",
                "Bremsstrahlung"    => "api_bremsstrahlung.md",
                "Compton"           => "api_compton.md",
                "Elastic_Leptons"   => "api_elastic_leptons.md",
                "Fluorescence"      => "api_fluorescence.md",
                "Inelastic_Leptons" => "api_inelastic_leptons.md",
                "Pair_Production"   => "api_pair_production.md",
                "Photoelectric"     => "api_photoelectric.md",
                "Rayleigh"          => "api_rayleigh.md",
            ],
            "Material" => "api_material.md",
            "Sources" => [
                "Surface_Source" => "api_surface_source.md",
                "Volume_Source"  => "api_volume_source.md",
                "Fixed_Sources"  => "api_fixed_source.md"
            ]
        ]
    ]
    )

deploydocs(
    repo = "github.com/CBienvenue/Radiant.jl",
    devurl = "dev",
    versions = ["stable" => devurl, devurl, devurl => devurl]
)