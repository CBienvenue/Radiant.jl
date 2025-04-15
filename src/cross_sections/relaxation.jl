"""
    relaxation(Z::Int64,Ei::Float64,Ecutoff::Float64,Ec::Float64,Ef⁻::Float64,Ef⁺::Float64,
    δi::Int64,incoming_particle::Particle,produced_particle::Particle,ηmin::Float64=0.01)

Gives the relaxation (fluorescence and Auger electron production) cross-section.

# Input Argument(s)
- `Z::Int64` : atomic number.
- `Ei::Float64` : incoming particle energy.
- `Ecutoff::Float64` : cutoff energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic events.
- `Ef⁻::Float64` : higher energy bound.
- `Ef⁺::Float64` : lower energy bound. 
- `δi::Int64` : subshell index.
- `incoming_particle::Particle` : incoming particle.
- `produced_particle::Particle` : particle produced following relaxation cascades.
- `ηmin::Float64 = 0.001` : minimum probability of the production of specific relaxation 
  particle production following electron cascades.

# Output Argument(s)
- `σs::Float64` : relaxation cross-sections.

# Reference(s)
- Perkins et al. (1991), Tables and Graphs of Atomic Subshell and Relaxation Data Derived
  from the LLNL Evaluated Atomic Data Library (EADL), Z = 1 - 100.
- Bienvenue et al. (2025), Toward highly accurate multigroup coupled 
  photon-electron-positron cross-sections for the Boltzmann Fokker-Planck equation.

"""
function relaxation(Z::Int64,Ei::Float64,Ecutoff::Float64,Ec::Float64,Ef⁻::Float64,Ef⁺::Float64,δi::Int64,incoming_particle::Particle,produced_particle::Particle,ηmin::Float64=0.01)

    #----
    # Compute and save in cache, if not already in cache
    #----
    if ~haskey(cache_radiant[],"relaxation_data") || ~haskey(cache_radiant[]["relaxation_data"],Z)
        ΔE_auger, η_auger, ΔE_florescence, η_fluoresence = atomic_electron_cascades(Z,Ecutoff,ηmin)
        if ~haskey(cache_radiant[],"relaxation_data")
            data = Dict{Int64,Dict}()
            cache_radiant[]["relaxation_data"] = data
        end
        data = cache_radiant[]["relaxation_data"]
        data[Z] = Dict{String,Vector{Vector{Float64}}}()
        data[Z]["ΔE_auger"] = ΔE_auger
        data[Z]["η_auger"] = η_auger
        data[Z]["ΔE_florescence"] = ΔE_florescence
        data[Z]["η_fluoresence"] = η_fluoresence
        cache_radiant[]["relaxation_data"] = data
    end

    #----
    # Extract data from cache
    #----
    data = cache_radiant[]["relaxation_data"]

    #----
    # Extract data
    #----
    _,Zi,Ui,_,_,_ = electron_subshells(Z)
    σ_per_subshell = 0

    # Inelastic collisionnal electron scattering
    if is_electron(incoming_particle)
        σ_per_subshell += Zi[δi] * integrate_moller_per_subshell(Z,Ei,0,Ui[δi],Ei-Ec,Ei)

    # Inelastic collisionnal positron scattering
    elseif is_positron(incoming_particle)
        σ_per_subshell += Zi[δi] * integrate_bhabha_per_subshell(Z,Ei,0,Ui[δi],Ei-Ec,Ei)

    # Photoelectric
    elseif is_photon(incoming_particle)
        σ_per_subshell += photoelectric_per_subshell(Z,Ei,δi)

    else
        error("Unknown incoming particle.")
    end

    #----
    # Fluorescence production cross section
    #----
    σs = 0
    if is_photon(produced_particle)
        ΔE_florescence = data[Z]["ΔE_florescence"][δi]
        η_fluoresence = data[Z]["η_fluoresence"][δi]
        Nt = length(ΔE_florescence)
        for δj in range(1,Nt)
            if Ef⁺ ≤ ΔE_florescence[δj] ≤ Ef⁻
                σs += η_fluoresence[δj] * σ_per_subshell
            end
        end

    #----
    # Auger electron production cross-section
    #----
    elseif is_electron(produced_particle)
        ΔE_auger = data[Z]["ΔE_auger"][δi]
        η_auger = data[Z]["η_auger"][δi]
        Nt = length(ΔE_auger)
        for δj in range(1,Nt)
            if Ef⁺ ≤ ΔE_auger[δj] ≤ Ef⁻
                σs += η_auger[δj] * σ_per_subshell
            end
        end
    else
        error("Relaxation produce either photons or electrons.")
    end
    return σs
end
