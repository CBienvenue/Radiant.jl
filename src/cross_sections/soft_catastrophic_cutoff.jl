"""
    soft_catastrophic_cutoff(Ei::Float64,Ei⁻::Float64,Ei⁺::Float64,Ei²⁺::Float64,
    scattering_model::String)

Compute the cutoff energy between soft and catastrophic events.

# Input Argument(s)
- `Ei::Float64` : energy of the incoming particle.
- `Ei⁻::Float64` : upper bound of the energy group g containing the particle energy.
- `Ei⁺::Float64` : lower bound of the energy group g containing the particle energy.
- `Ei²⁺::Float64` : lower bound of the energy group g+1.
- `scattering_model::String` : scattering model.

# Output Argument(s)
- `Ec::Float64`: cutoff energy between soft and catastrophic events.

# Reference(s)
- Bienvenue et al. (2025), Toward highly accurate multigroup coupled photon-electron-positron
  cross-sections for the Boltzmann Fokker-Planck equation.

"""
function soft_catastrophic_cutoff(Ei::Float64,Ei⁻::Float64,Ei⁺::Float64,Ei²⁺::Float64,scattering_model::String)
    #----
    # Boltzmann Fokker-Planck solver (mixed soft and catastrphic events)
    #----
    if scattering_model == "BFP"
        Ec = Ei*(Ei⁺-Ei²⁺)/(Ei⁻-Ei⁺) - (Ei⁺^2-Ei⁻*Ei²⁺)/(Ei⁻-Ei⁺)
    #----
    # Fokker-Planck solver (soft events only)
    #----
    elseif scattering_model == "FP"
        Ec = 0.0
    #----
    # Boltzmann solver (catastrophic events only)
    #----
    elseif scattering_model == "BTE"
        Ec = Ei
    else
        error("Unknown scattering model $scattering_model.")
    end
    return Ec
end