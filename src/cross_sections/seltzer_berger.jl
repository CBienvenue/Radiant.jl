"""
    seltzer_berger_cross_section(Z::Int64,Ei::Float64,Eγ::Float64,particle::Particles)

Gives the bremsstrahlung scattering cross-sections based on Seltzer and Berger tables.

# Input Argument(s)
- `Z::Int64` : atomic number.
- `Ei::Float64` : incoming particle energy.
- `Eγ::Float64` : Bremsstrahlung photon energy.
- `particle::Particle` : incoming particle.

# Output Argument(s)
- `σs::Float64` : Bremsstrahlung scattering cross-sections

# Reference(s)
- Seltzer and Berger (1986), Bremsstrahlung energy spectra from electrons with kinetic
  energy 1 keV - 10 GeV incident on screened nuclei and orbital electrons of neutral atoms
  with Z = 1-100.

"""
function seltzer_berger_cross_section(Z::Int64,Ei::Float64,Eγ::Float64,particle::Particle)

    #----
    # Extract data
    #----
    data = fast_load("bremsstrahlung_data_seltzer_berger_1986.jld2")
    E = data["incident_electron_energy"] # (in mₑc²)
    if ~(E[1] ≤ Ei ≤ E[end]) error("Bremsstrahlung tables from Seltzer and Berger are defined between 1 keV and 10 GeV.") end
    r = data["radiative_energy_fraction"]
    ri = Eγ/Ei
    if ~(r[1] ≤ ri ≤ r[end]) error("Ratio between the emitted bremsstrahlung photon and the incoming particle should be between 0 and 1.") end

    #----
    # Extract scaled cross-sections
    #----
    if ~haskey(data["scaled_cross_sections"],Z) error("Undefined scaled Bremsstrahlung cross-sections for Z = $Z") end
    χ = data["scaled_cross_sections"][Z]

    #----
    # Correction for positrons
    #----
    if is_positron(particle) 
        ratio = ratio_positron_electron_bremsstrahlung(Z,Ei) 
    elseif is_electron(particle) 
        ratio = 1 
    else 
        error("Unknown particle.") 
    end

    #----
    # Compute the scattering cross-section
    #----
    β² = Ei*(Ei+2)/(Ei+1)^2
    i = searchsortedfirst(E,Ei)
    j = searchsortedfirst(r,ri)
    if (i == 1 && E[1] ≈ Ei) i=2; Ei=E[1] end
    if (j == 1 && r[1] ≈ ri) j=2; ri=r[1] end
    if (i == length(E)+1 && E[end] ≈ Ei) i=length(E); Ei=E[end] end
    if (j == length(r)+1 && r[end] ≈ ri) j=length(r); ri=r[end] end
    if (i == 1 || i > length(E)) error("Interpolation value is outside the interpolation vector.") end
    if (j == 1 || j > length(r)) error("Interpolation value is outside the interpolation vector.") end
    χij⁻ = log_interpolation(Ei,E[i-1:i],χ[i-1:i,j-1])
    χij⁺ = log_interpolation(Ei,E[i-1:i],χ[i-1:i,j])
    χij = linear_interpolation(ri,r[j-1:j],[χij⁻,χij⁺])
    σs = ratio * χij/1e27/Eγ * Z^2/β²
    return σs
end

"""
    seltzer_berger_stopping_power(Z::Int64,Ei::Float64,particle::Particles)

Gives the bremsstrahlung stopping powers based on Seltzer and Berger tables.

# Input Argument(s)
- `Z::Int64` : atomic number.
- `Ei::Float64` : incoming particle energy.
- `particle::Particle` : incoming particle.

# Output Argument(s)
- `S::Float64` : Bremsstrahlung stopping powers

# Reference(s)
- Seltzer and Berger (1986), Bremsstrahlung energy spectra from electrons with kinetic
  energy 1 keV - 10 GeV incident on screened nuclei and orbital electrons of neutral atoms
  with Z = 1-100.

"""
function seltzer_berger_stopping_power(Z::Int64,Ei::Float64,particle::Particle)

    #----
    #Extract data
    #----
    data = fast_load("bremsstrahlung_data_seltzer_berger_1986.jld2")
    E = data["incident_electron_energy"] # (in mₑc²)
    if ~(E[1] ≤ Ei ≤ E[end]) error("Bremsstrahlung tables from Seltzer and Berger are defined between 1 keV and 10 GeV.") end
    if ~haskey(data["scaled_stopping_powers"],Z) error("Undefined scaled Bremsstrahlung cross-sections for Z = $Z") end
    ϕ = data["scaled_stopping_powers"][Z]

    #----
    # Compute spline to interpolate the scaled stopping power
    #----
    spline_ssp = fast_spline(E,ϕ,Z,"scaled_stopping_powers")

    #----
    # Correction for positrons
    #----
    if is_positron(particle) 
        ratio = ratio_positron_electron_bremsstrahlung(Z,Ei) 
    elseif is_electron(particle) 
        ratio = 1 
    else 
        error("Unknown particle.") 
    end

    #----
    # Compute the stopping power
    #----
    ϕ_spline = spline_ssp(Ei)
    α = 1/137
    rₑ = 2.81794092e-13 # (in cm)
    S = ratio * ϕ_spline * α * rₑ^2 * Z^2 * (Ei+1)
    return S
end