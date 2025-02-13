"""
    transport_correction(interaction::Interaction,L::Int64,Σt::Float64,
    Σsℓ::Vector{Float64},α::Float64,solver::String)

Compute the transport correction and/or the decomposition of elastic scattering between
soft and catastrophic components.

# Input Argument(s)
- `interaction::Interaction` : type of interaction.
- `L::Int64` : Legendre truncation order.
- `Σt::Float64` : total cross-section.
- `Σsℓ::Vector{Float64}` : Legendre moments of the scattering cross-section.
- `α::Float64` : momentum transfer.
- `solver::String` : solver type.

# Output Argument(s)
- `Σt::Float64` : total cross-section.
- `Σsℓ::Vector{Float64}` : Legendre moments of the scattering cross-section.
- `α::Float64` : momentum transfer.

# Reference(s)
- Hébert (2016), Applied Reactor Physics.
- Landesman and Morel (1989), Angular Fokker-Planck Decomposition and Representation
  Techniques.

"""
function transport_correction(interaction::Interaction,L::Int64,Σt::Float64,Σsℓ::Vector{Float64},T::Float64,solver::String)

if interaction.is_elastic && solver ∈ ["BFP","BTE"]

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

    # Momentum transfer
    if ℓmax > 1 && solver == "BFP" && interaction.is_AFP
        Tt = (Σsℓ[ℓmax]-Σsℓ[ℓmax+1])/(2*ℓmax)
        T += Tt
        Σt -= Tt*ℓmax*(ℓmax+1) + Σsℓ[ℓmax+1]
        for ℓ in range(0,ℓmax)
            Σsℓ[ℓ+1] -= Tt*(ℓmax*(ℓmax+1)-ℓ*(ℓ+1)) + Σsℓ[ℓmax+1]
        end 
    end

end

return Σt, Σsℓ, T
end