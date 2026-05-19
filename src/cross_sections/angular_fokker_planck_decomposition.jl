"""
    angular_fokker_planck_decomposition(interaction::Interaction,L::Int64,Σt::Float64,
    Σsl::Vector{Float64},α::Float64)

Decompose the elastic cross-sections between soft and catastrophic, such as the
catastrophic part is dealt with the Boltzmann operator and the soft part, with the angular
Fokker-Planck operator.

# Input Argument(s)
- `interaction::Interaction` : type of interaction.
- `L::Int64` : Legendre truncation order.
- `Σt::Float64` : total cross-section.
- `Σsl::Vector{Float64}` : Legendre moments of the scattering cross-section.
- `α::Float64` : momentum transfer.
- `scattering_model::String` : scattering_model type.
- `is_AFP_decomposition::Bool` : boolean indicating if angular Fokker-Planck decomposition
  method is employed.

# Output Argument(s)
- `Σt::Float64` : total cross-section.
- `Σsl::Vector{Float64}` : Legendre moments of the scattering cross-section.
- `α::Float64` : momentum transfer.

# Reference(s)
- Landesman and Morel (1989), Angular Fokker-Planck Decomposition and Representation
  Techniques.

"""
function angular_fokker_planck_decomposition(interaction::Interaction,L::Int64,Σt::Float64,Σsl::Vector{Float64},T::Float64)

# Find the last non-zero Legendre moment
lmax = L
for l in range(L,1,step=-1)
    lmax = l
    is_momentum_transfer_ok = true
    for li in range(0,l-1)
        if Σsl[li+1] - Σsl[l+1] - (Σsl[l]-Σsl[l+1])/(2*l)*(l*(l+1)-li*(li+1)) < 0
            is_momentum_transfer_ok = false
            break
        end
    end
    if Σsl[l+1] != 0.0 && is_momentum_transfer_ok
        Σsl[l+2:end] .= 0
        break 
    end
end

# Momentum transfer, total cross-sections and elastic scattering cross-sections modification
if lmax > 1 && interaction.is_AFP_decomposition
    Tt = (Σsl[lmax]-Σsl[lmax+1])/(2*lmax)
    T += Tt
    Σt -= Tt*lmax*(lmax+1) + Σsl[lmax+1]
    for l in range(0,lmax)
        Σsl[l+1] -= Tt*(lmax*(lmax+1)-l*(l+1)) + Σsl[lmax+1]
    end 
end

return Σt, Σsl, T
end