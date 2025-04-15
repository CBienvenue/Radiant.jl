"""
    klein_nishina(Ei::Float64,Ef::Float64)

Gives the Klein-Nishina differential cross-section.

# Input Argument(s)
- `Ei::Float64` : incoming photon energy.
- `Ef::Float64` : outgoing photon energy.

# Output Argument(s)
- `σs::Float64` : Klein-Nishina differential cross-section.

# Reference(s)
- Lorence et al. (1989), Physics guide to CEPXS: a multigroup coupled electron-photon
  cross-section generating code.

"""
function klein_nishina(Ei::Float64,Ef::Float64)
    if Ei/(1+2*Ei) ≤ Ef ≤ Ei
        rₑ = 2.81794092e-13 # (in cm)
        σs = π * rₑ^2 / Ei^2 * (Ei/Ef + Ef/Ei - 2*(1/Ef-1/Ei) + (1/Ef-1/Ei)^2)
    end
    return σs
end

"""
    klein_nishina(Ei::Float64,Ef::Float64)

Gives the Legendre moments of the Klein-Nishina angular distribution.

# Input Argument(s)
- `Ei::Float64` : incoming photon energy.
- `Ef::Float64` : outgoing particle energy.
- `L::Int64` : Legendre truncation order.
- `particle::Particle` : outgoing particle.

# Output Argument(s)
- `Wℓ::Float64` : Legendre moments of the Klein-Nishina angular distribution.

# Reference(s)
- Lorence et al. (1989), Physics guide to CEPXS: a multigroup coupled electron-photon
  cross-section generating code.

"""
function angular_klein_nishina(Ei::Float64,Ef::Float64,L::Int64,particle::Particle)

    # Initialization
    Wℓ = zeros(L+1)
    if Ei < Ef return Wℓ end

    # Compute the direction cosine
    if is_photon(particle)
        if ~(Ei/(1+2*Ei) ≤ Ef ≤ Ei) return Wℓ end
        μ = 1 + 1/Ei - 1/Ef
    elseif is_electron(particle)
        if ~(0 ≤ Ef ≤ 2*Ei^2/(1+2*Ei)) return Wℓ end
        μ = (1 + Ei)/Ei * 1/sqrt(2/Ef+1)
    end

    # Compute the Legendre moments angular distribution
    Pℓμ = legendre_polynomials(L,μ)
    for ℓ in range(0,L) Wℓ[ℓ+1] += Pℓμ[ℓ+1] end
    return Wℓ
end