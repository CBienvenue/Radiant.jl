"""
    impulse_approximation(Z::Int64,δi::Int64,Ei::Float64,Ef::Float64)

Gives the Legendre moments of the Compton differential cross-section based on the
relativistic impulse approximation.

# Input Argument(s)
- `Z::Int64` : atomic number.
- `δi::Int64` : subshell index
- `Ei::Float64` : incoming photon energy.
- `Ef::Float64` : outgoing photon energy.

# Output Argument(s)
- `σℓ::Float64` : Legendre moments of the Compton differential cross-section.

# Reference(s)
- Brusa et al. (1996), Fast sampling algorithm for the simulation of photon Compton
  scattering.

"""
function impulse_approximation(Z::Int64,L::Int64,δi::Int64,Ei::Float64,Ef::Float64)

    #----
    # Initialization
    #----
    Nμ = 80
    μ,w = quadrature(Nμ,"gauss-lobatto")
    mₑc² = 0.510999
    J₀i = orbital_compton_profiles(Z)[δi]
    _,Zi,Ui,_,_,_ = electron_subshells(Z)

    #----
    # Change of units
    #----
    mₑ = 1                        # (a₀)
    c = 137.03599908388762        # (a₀×Eₕ/ħ)
    rₑ = 5.325135459237564e-5     # (mₑ)
    a₀ = 5.29177210903e-11        # (a₀)
    # Conversion to atomic units (MeV) -> (Eₕ)
    Ei = Ei * 1e6 / 27.211386245988 * mₑc²
    Ef = Ef * 1e6 / 27.211386245988 * mₑc²
    Ui .= Ui * 1e6 ./ 27.211386245988 * mₑc²

    #----
    # Compute the scattering cross-section
    #----
    σℓ = zeros(L+1)
    if Ei - Ef - Ui[δi] ≥ 0
        for n in range(1,Nμ)
            Pℓμ = legendre_polynomials(L,μ[n])
            Ec = Ei*mₑ*c^2/(mₑ*c^2+Ei*(1-μ[n]))
            pz = (Ei*Ef*(1-μ[n])-mₑ*c^2*(Ei-Ef))/(c*sqrt(Ei^2+Ef^2-2*Ei*Ef*μ[n]))
            Ji = J₀i*(1+2*J₀i*abs(pz))*exp(1/2-1/2*(1+2*J₀i*abs(pz))^2)
            σs = w[n] * π * rₑ^2 * Ef/Ei * (Ec/Ei+Ei/Ec+μ[n]^2-1) * 1/sqrt(Ei^2+Ef^2-2*Ei*Ef*μ[n]+Ei^2*(Ef-Ec)^2/Ec^2) * Zi[δi]*Ji *mₑ*c
            σs *= (a₀ ^ 2) / 27.211386245988 * 100^2 * 1e6 * mₑc² # Change of unit
            for ℓ in range(0,L) σℓ[ℓ+1] += Pℓμ[ℓ+1] * σs end
        end
    end
    return σℓ
end