"""
    heitler(Z::Int64,Ei::Float64,Ef::Float64)

Gives the annihilation scattering cross-section based on Heitler model.

# Input Argument(s)
- `Z::Int64` : atomic numbers of the element.
- `Ei::Float64` : incoming positron energy.
- `Ef::Float64` : outgoing lowest-energy photon energy.

# Output Argument(s)
- `σs::Float64` : annihilation scattering cross-section.

# Reference(s)
- Heitler (1954), The Quantum Theory of Radiation.
- Nelson et al. (1985), The EGS4 Code System, Technical Report SLAC-265.
- Salvat (2019), PENELOPE-2018: A Code System for Monte-Carlo Simulation of Electron and
  Photon Transport.

"""
function heitler(Z::Int64,Ei::Float64,Ef::Float64)

    # Initialization
    rₑ = 2.81794092e-13 # (in cm)
    γ = Ei+1
    ζmin = 1/(γ+1+sqrt(γ^2-1))

    # Compute the direction cosine
    if ζmin ≤ Ef/(Ei+2) ≤ 1/2  # Lowest energy photon
        ζ = Ef/(Ei+2)
    elseif 1/2 ≤ Ef/(Ei+2) ≤ 1-ζmin # Highest energy photon
        ζ = (Ei+2-Ef)/(Ei+2)
    else
        error("Value of ζ is out of bounds. $Ef $Ei")
    end

    # Compute annihilation scattering cross-section
    S_heitler(ζ) = -(γ+1)^2 + (γ^2+4*γ+1)/ζ - 1/ζ^2
    σs = Z*π*rₑ^2/((γ+1)^2*(γ^2-1)) * (S_heitler(ζ)+S_heitler(1-ζ))
    return σs
end

"""
    integrated_heitler(Z::Int64,Ei::Float64,Ef::Float64)

Gives the absorption cross-section based on Heitler model.

# Input Argument(s)
- `Z::Int64` : atomic numbers of the element.
- `Ei::Float64` : incoming positron energy.

# Output Argument(s)
- `σt::Float64` : annihilation absorption cross-section.

# Reference(s)
- Heitler (1954), The Quantum Theory of Radiation.
- Nelson et al. (1985), The EGS4 Code System, Technical Report SLAC-265.
- Salvat (2019), PENELOPE-2018: A Code System for Monte-Carlo Simulation of Electron and
  Photon Transport.

"""
function integrated_heitler(Z::Int64,Ei::Float64)

    # Initialization
    γ = Ei+1
    rₑ = 2.81794092e-13 # (in cm)

    # Compute the absorption annihilation cross-section
    σa = Z*π*rₑ^2/((γ+1)*(γ^2-1)) * ((γ^2+4*γ+1)*log(γ+sqrt(γ^2-1))-(3+γ)*sqrt(γ^2-1))
    return σa
end

"""
    angular_heitler(Ei::Float64,Ef::Float64)

Gives the Legendre moments of the Heitler angular distribution.

# Input Argument(s)
- `Ei::Float64` : incoming positron energy.
- `Ef::Float64` : outgoing lowest-energy photon energy.
- `L::Int64` : Legendre truncation order.

# Output Argument(s)
- `Wl::Float64` : Legendre moments of the Heitler angular distribution.

# Reference(s)
- Heitler (1954), The Quantum Theory of Radiation.
- Nelson et al. (1985), The EGS4 Code System, Technical Report SLAC-265.
- Salvat (2019), PENELOPE-2018: A Code System for Monte-Carlo Simulation of Electron and
  Photon Transport.

"""
function angular_heitler(Ei::Float64,Ef::Float64,L::Int64)

    # Initialization
    Wl = zeros(L+1)
    γ = Ei+1
    ζmin = 1/(γ+1+sqrt(γ^2-1))

    # Compute the direction cosine
    if ζmin ≤ Ef/(Ei+2) ≤ 1/2  # Lowest energy photon
        ζ = Ef/(Ei+2)  
        μ = (γ+1-1/ζ)/sqrt(γ^2-1)
    elseif 1/2 ≤ Ef/(Ei+2) ≤ 1-ζmin # Highest energy photon
        ζ = (Ei+2-Ef)/(Ei+2)
        μ = (γ+1-1/(1-ζ))/sqrt(γ^2-1)
    else
        error("Value of ζ is out of bounds.")
    end

    # Compute the Legendre moments of the angular distribution
    Plμ = legendre_polynomials_up_to_L(L,μ)
    for l in range(0,L) Wl[l+1] += Plμ[l+1] end
    return Wl
end
