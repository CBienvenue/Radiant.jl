module Radiant

    #----
    # External packages
    #----
    import Base.println
    using Printf: @sprintf
    using LinearAlgebra
    using JLD2

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
        "interaction_interdependances.jl"
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
        "compute_flux.jl",
        "compute_one_speed.jl",
        "compute_sweep_1D.jl",
        "compute_sweep_2D.jl",
        "compute_sweep_3D.jl",
        "flux_1D_BTE.jl",
        "flux_1D_BFP.jl",
        "flux_2D_BTE.jl",
        "flux_2D_BFP.jl",
        "flux_3D_BTE.jl",
        "flux_3D_BFP.jl",
        "adaptive.jl",
        "scheme_weights.jl",
        "map_moments.jl",
        "constant_linear.jl",
        "angular_polynomial_basis.jl",
        "livolant.jl",
        "energy_deposition.jl",
        "charge_deposition.jl",
        "flux.jl"
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
        "ferrer_associated_legendre.jl",
        "real_spherical_harmonics.jl",
        "heaviside.jl",
        "interpolation.jl",
        "integrals.jl",
        "spline.jl",
        "voronoi.jl",
        "newton_bissection.jl",
        "one_space.jl",
        "cache.jl",
        "find_package_root.jl",
        "python_method_notation.jl"
    ]
    for folder in ["structures/","tools/","cross_sections/","particle_transport/"]
        for file in radiant_src[folder] include(string(folder,file)) end
    end

    #----
    # Export objects
    #----
    export Particle, Photon, Electron, Positron, Proton, Antiproton, Alpha, Muon, Antimuon
    export Elastic_Collision,Inelastic_Collision,Bremsstrahlung,Compton,Pair_Production,Photoelectric,Annihilation,Rayleigh,Relaxation,Fluorescence,Auger
    export Material,Cross_Sections,Geometry,Discrete_Ordinates,Solvers,Surface_Source,Volume_Source,Fixed_Sources,Computation_Unit

end
