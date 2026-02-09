module Radiant

    #----
    # External packages
    #----
    import Base.println
    using Printf: @sprintf
    using LinearAlgebra
    using JLD2
    using SpecialFunctions

    #----
    # Include files
    #----
    radiant_src = Dict{String,Vector{String}}()
    radiant_src["cross_sections/"] = [
        "read_fmac_m.jl",
        "write_fmac_m.jl",
        "generate_cross_sections.jl",
        "custom_cross_sections.jl",
        "energy_group_structure.jl",
        "atomic_electron_cascades.jl",
        "multigroup.jl",
        "feed.jl",
        "mean_excitation_energy.jl",
        "atomic_weight.jl",
        "density.jl",
        "nuclei_density.jl",
        "plasma_energy.jl",
        "electron_subshells.jl",
        "orbital_compton_profiles.jl",
        "resonance_energy.jl",
        "fermi_density_effect.jl",
        "atomic_number.jl",
        "transport_correction.jl",
        "angular_fokker_planck_decomposition.jl",
        "bethe.jl",
        "moller.jl",
        "bhabha.jl",
        "kawrakow_correction.jl",
        "moliere_screening.jl",
        "mott.jl",
        "sommerfield.jl",
        "poskus.jl",
        "ratio_positron_electron_bremsstrahlung.jl",
        "seltzer_berger.jl",
        "klein_nishina.jl",
        "waller_hartree.jl",
        "impulse_approximation.jl",
        "biggs_lighthill.jl",
        "photoelectric_per_subshell.jl",
        "sauter.jl",
        "rayleigh.jl",
        "high_energy_coulomb_correction.jl",
        "baro.jl",
        "heitler.jl",
        "relaxation.jl",
        "interaction_interdependances.jl",
        "inelastic_collision_heavy_particle.jl",
        "soft_catastrophic_cutoff.jl"
    ]
    radiant_src["particle_transport/"] = [
        "geometry.jl",
        "volume_source.jl",
        "surface_source.jl",
        "scattering_source.jl",
        "particle_source.jl",
        "fokker_planck_source.jl",
        "fokker_planck_scattering_matrix.jl",
        "fokker_planck_finite_difference.jl",
        "fokker_planck_differential_quadrature.jl",
        "fokker_planck_galerkin.jl",
        "electromagnetic_scattering_matrix.jl",
        "transport.jl",
        "sn_flux.jl",
        "pn_flux.jl",
        "sn_one_speed.jl",
        "pn_one_speed.jl",
        "sn_sweep_1D.jl",
        "sn_sweep_2D.jl",
        "sn_sweep_3D.jl",
        "sn_1D_bte.jl",
        "sn_1D_bfp.jl",
        "sn_2D_bte.jl",
        "sn_2D_bfp.jl",
        "sn_3D_bte.jl",
        "sn_3D_bfp.jl",
        "adaptive.jl",
        "scheme_weights.jl",
        "map_moments.jl",
        "constant_linear.jl",
        "angular_polynomial_basis.jl",
        "livolant.jl",
        "energy_deposition.jl",
        "charge_deposition.jl",
        "flux.jl",
        "pn_1D_bte.jl",
        "pn_1D_bfp.jl",
        "pn_2D_bte.jl",
        "pn_3D_bte.jl",
        "pn_sweep_1D.jl",
        "pn_sweep_2D.jl",
        "pn_sweep_3D.jl",
        "restricted_to_full_domain_matrix.jl",
        "pn_weights.jl",
        "gn_flux.jl",
        "gn_one_speed.jl",
        "gn_sweep_1D.jl",
        "gn_sweep_2D.jl",
        "gn_sweep_3D.jl",
        "gn_1D_bte.jl",
        "gn_1D_bfp.jl",
        "gn_2D_bte.jl",
        "gn_2D_bfp.jl",
        "gn_3D_bte.jl",
        "gn_3D_bfp.jl"
    ]
    radiant_src["structures/"] = [
        "Particle.jl",
        "Interaction.jl",
        "Elastic_Collision.jl",
        "Inelastic_Collision.jl",
        "Bremsstrahlung.jl",
        "Compton.jl",
        "Pair_Production.jl",
        "Photoelectric.jl",
        "Annihilation.jl",
        "Rayleigh.jl",
        "Relaxation.jl",
        "Multigroup_Cross_Sections.jl",
        "Material.jl",
        "Cross_Sections.jl",
        "Geometry.jl",
        "Discrete_Ordinates.jl",
        "Spherical_Harmonics.jl",
        "Galerkin.jl",
        "Solvers.jl",
        "Surface_Source.jl",
        "Volume_Source.jl",
        "Source.jl",
        "Sources.jl",
        "Fixed_Sources.jl",
        "Flux_Per_Particle.jl",
        "Flux.jl",
        "Computation_Unit.jl"
    ]
    radiant_src["tools/"] = [
        "quadrature.jl",
        "gauss_legendre.jl",
        "gauss_lobatto.jl",
        "gauss_legendre_chebychev.jl",
        "lebedev.jl",
        "carlson.jl",
        "legendre_polynomials.jl",
        "jacobi_polynomials.jl",
        "ferrer_associated_legendre.jl",
        "real_spherical_harmonics.jl",
        "heaviside.jl",
        "interpolation.jl",
        "integrals.jl",
        "spline.jl",
        "voronoi.jl",
        "newton_bisection.jl",
        "one_space.jl",
        "cache.jl",
        "find_package_root.jl",
        "python_method_notation.jl",
        "factorial_factor.jl",
        "double_factorial.jl"
    ]
    for folder in ["structures/","tools/","cross_sections/","particle_transport/"]
        for file in radiant_src[folder] include(string(folder,file)) end
    end

    #----
    # Export objects
    #----
    export Particle, Photon, Electron, Positron, Proton, Antiproton, Alpha, Muon, Antimuon
    export Elastic_Collision,Inelastic_Collision,Bremsstrahlung,Compton,Pair_Production,Photoelectric,Annihilation,Rayleigh,Relaxation,Fluorescence,Auger
    export Material,Cross_Sections,Geometry,Discrete_Ordinates,Solvers,Surface_Source,Volume_Source,Fixed_Sources,Computation_Unit,Spherical_Harmonics,Galerkin

    #----
    # Execution of Radiant script files
    #----
    export @radiant_input
    macro radiant_input()
        return quote
            if abspath(PROGRAM_FILE) == @__FILE__
                using Radiant
                Radiant.run_script(@__FILE__)
                exit()
            end
        end
    end
    function run_script(script::AbstractString)
        isfile(script) || error("Input script not found: $script")
        include(script)
    end
end
