"""
    transport_correction(interaction::Interaction,L::Int64,Σt::Float64,
    Σsℓ::Vector{Float64},α::Float64,scattering_model::String)

Compute the transport correction and/or the decomposition of elastic scattering between
soft and catastrophic components.

# Input Argument(s)
- `interaction::Interaction` : type of interaction.
- `L::Int64` : Legendre truncation order.
- `Σt::Float64` : total cross-section.
- `Σsℓ::Vector{Float64}` : Legendre moments of the scattering cross-section.
- `α::Float64` : momentum transfer.
- `scattering_model::String` : scattering_model type.

# Output Argument(s)
- `Σt::Float64` : total cross-section.
- `Σsℓ::Vector{Float64}` : Legendre moments of the scattering cross-section.
- `α::Float64` : momentum transfer.

# Reference(s)
- Hébert (2016), Applied Reactor Physics.

"""
function transport_correction(interaction::Interaction,L::Int64,Σt::Float64,Σsℓ::Vector{Float64})

# Find the last non-zero Legendre moment
ℓmax = L
for ℓ in range(L,1,step=-1)
    ℓmax = ℓ
    if Σsℓ[ℓ+1] != 0.0 && Σsℓ[ℓ+1] < minimum(Σsℓ[1:ℓ])
        Σsℓ[ℓ+2:end] .= 0
        break
    end
end

# Extended transport correction (ETC) for elastic scattering
if interaction.is_ETC
    Σt -= Σsℓ[ℓmax+1]
    Σsℓ[1:ℓmax+1] .-= Σsℓ[ℓmax+1]
end

return Σt, Σsℓ
end