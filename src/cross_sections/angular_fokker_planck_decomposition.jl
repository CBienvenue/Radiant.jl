"""
    angular_fokker_planck_decomposition(interaction::Interaction,L::Int64,Σt::Float64,
    Σsℓ::Vector{Float64},α::Float64)

Decompose the elastic cross-sections between soft and catastrophic, such as the
catastrophic part is dealt with the Boltzmann operator and the soft part, with the angular
Fokker-Planck operator.

# Input Argument(s)
- `interaction::Interaction` : type of interaction.
- `L::Int64` : Legendre truncation order.
- `Σt::Float64` : total cross-section.
- `Σsℓ::Vector{Float64}` : Legendre moments of the scattering cross-section.
- `α::Float64` : momentum transfer.
- `scattering_model::String` : scattering_model type.
- `is_AFP_decomposition::Bool` : boolean indicating if angular Fokker-Planck decomposition
  method is employed.

# Output Argument(s)
- `Σt::Float64` : total cross-section.
- `Σsℓ::Vector{Float64}` : Legendre moments of the scattering cross-section.
- `α::Float64` : momentum transfer.

# Reference(s)
- Landesman and Morel (1989), Angular Fokker-Planck Decomposition and Representation
  Techniques.

"""
function angular_fokker_planck_decomposition(interaction::Interaction,L::Int64,Σt::Float64,Σsℓ::Vector{Float64},T::Float64)

# Find the last non-zero Legendre moment
ℓmax = L
for ℓ in range(L,1,step=-1)
    ℓmax = ℓ
    is_momentum_transfer_ok = true
    for l in range(0,ℓ-1)
        if Σsℓ[l+1] - Σsℓ[ℓ+1] - (Σsℓ[ℓ]-Σsℓ[ℓ+1])/(2*ℓ)*(ℓ*(ℓ+1)-l*(l+1)) < 0
            is_momentum_transfer_ok = false
            break
        end
    end
    if Σsℓ[ℓ+1] != 0.0 && is_momentum_transfer_ok
        Σsℓ[ℓ+2:end] .= 0
        break 
    end
end

# Momentum transfer, total cross-sections and elastic scattering cross-sections modification
if ℓmax > 1 && interaction.is_AFP_decomposition
    Tt = (Σsℓ[ℓmax]-Σsℓ[ℓmax+1])/(2*ℓmax)
    T += Tt
    Σt -= Tt*ℓmax*(ℓmax+1) + Σsℓ[ℓmax+1]
    for ℓ in range(0,ℓmax)
        Σsℓ[ℓ+1] -= Tt*(ℓmax*(ℓmax+1)-ℓ*(ℓ+1)) + Σsℓ[ℓmax+1]
    end 
end

return Σt, Σsℓ, T
end