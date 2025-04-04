module Radiant

    #----
    # External packages
    #----
    import Base.println
    using Printf: @sprintf
    using LinearAlgebra
    using JLD2

    #----
    # Find root of Radiant
    #----
    function find_package_root()
        current_dir = @__DIR__
        while current_dir != "/"
            if isfile(joinpath(current_dir, "Project.toml"))
                return current_dir
            end
            current_dir = dirname(current_dir)
        end
        error("Package root not found.")
    end

    #----
    # Global unique ID generator
    #----
    const unique_id_counter = Ref(0)
    function generate_unique_id()
        unique_id_counter[] += 1
        return unique_id_counter[]
    end

    #----
    # Cache
    #----
    const cache_files = Ref{Dict{String,Any}}(Dict())

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
        "baro_coefficients.jl",
        "transport_correction.jl",
        "angular_fokker_planck_decomposition.jl",
        "bethe.jl"
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
        "Elastic_Leptons.jl",
        "Inelastic_Leptons.jl",
        "Bremsstrahlung.jl",
        "Compton.jl",
        "Pair_Production.jl",
        "Photoelectric.jl",
        "Annihilation.jl",
        "Rayleigh.jl",
        "Fluorescence.jl",
        "Auger.jl",
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
        "one_space.jl"
    ]
    for folder in ["structures/","tools/","cross_sections/","particle_transport/"]
        for file in radiant_src[folder] include(string(folder,file)) end
    end

    #----
    # Export objects
    #----
    export Particle, Photon, Electron, Positron, Proton, Antiproton, Alpha, Muon, Antimuon
    export Elastic_Leptons, Inelastic_Leptons, Bremsstrahlung,Compton,Pair_Production,Photoelectric,Annihilation,Rayleigh,Fluorescence,Auger
    export Material,Cross_Sections,Geometry,Discrete_Ordinates,Solvers,Surface_Source,Volume_Source,Fixed_Sources,Computation_Unit

    #----
    # To use Python-like notation for methods
    #----
    RadiantObject = Union{
        Material,
        Cross_Sections,
        Multigroup_Cross_Sections,
        Geometry,
        Discrete_Ordinates,
        Solvers,
        Surface_Source,
        Volume_Source,
        Fixed_Sources,
        Source,
        Sources,
        Computation_Unit,
        Flux_Per_Particle,
        Flux,
        Interaction,
        Particle
    }

    """
        Base.getproperty(object::RadiantObject, s::Symbol)

    Special function to enable calling a method "m(this::T,...)" on a struct "S" of type "T"
    using the notation "S.m(...)" for a subset of struct.

    """
    function Base.getproperty(object::RadiantObject, s::Symbol)
        if hasfield(typeof(object), s)
            return getfield(object, s)
        else
            # Check if the function is defined in the Radiant module
            if isdefined(Radiant, s) && getproperty(Radiant, s) isa Function
                try
                    return function(x...)
                        func = getproperty(Radiant, s)
                        return func(object, x...)
                    end
                catch e
                    # Log the original error and rethrow it
                    type = typeof(object)
                    @error "Failed to invoke function '$s' on object of type '$type'. Error: $e"
                    rethrow()
                end
            else
                type = typeof(object)
                throw(error("The object of type '$type' has no function '$s'."))
            end
        end
    end

end
