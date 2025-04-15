"""
    sommerfield(Ei::Float64,L::Int64)

Gives the Legendre moments of the Sommerfield angular distribution.

# Input Argument(s)
- `Ei::Float64` : incoming particle energy.
- `L::Int64` : Legendre truncation order.

# Output Argument(s)
- `Wℓ::Vector{Wℓ}` : Sommerfield angular distribution.

# Reference(s)
- Lorence et al. (1989), Physics guide to CEPXS: a multigroup coupled electron-photon

"""
function sommerfield(Ei::Float64,L::Int64)
    β = sqrt(Ei*(Ei+2)/(Ei+1)^2)
    Wℓ = zeros(L+1)
    for ℓ in range(0,L)
        if ℓ == 0
            Wℓ[ℓ+1] = 1
        elseif ℓ == 1
            Wℓ[ℓ+1] = (2*β + (1-β^2)*log((1-β)/(1+β)))/(2*β^2)
        else
            Wℓ[ℓ+1] = ((2*ℓ-1)*Wℓ[ℓ] - ℓ*β*Wℓ[ℓ-1])/((ℓ-1)*β)
        end
    end
    return Wℓ
end