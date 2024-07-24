function transport_correction(interaction::Interaction,L::Int64,Σt::Float64,Σsℓ::Vector{Float64},α::Float64,solver)

if interaction.is_elastic && solver ∈ ["BFP","BTE"]

    # Find the last non-zero Legendre moment
    ℓmax = L
    for ℓ in range(L,1,step=-1)
        ℓmax = ℓ
        if Σsℓ[ℓ+1] != 0.0 break end
    end

    # Mott momentum transfer transpot correction
    if ℓmax > 1 && solver == "BFP" && interaction.is_AFP
        αt = (Σsℓ[ℓmax]-Σsℓ[ℓmax+1])/ℓmax
        α += αt
        Σt -= αt/2*ℓmax*(ℓmax+1)
        for ℓ in range(0,ℓmax)
            Σsℓ[ℓ+1] -= αt/2*(ℓmax*(ℓmax+1)-ℓ*(ℓ+1))
        end 
    end

    # Extended transport correction (ETC) for elastic scattering
    if interaction.is_ETC
        Σt -= Σsℓ[ℓmax+1]
        Σsℓ[1:ℓmax+1] .-= Σsℓ[ℓmax+1]
    end

end

return Σt, Σsℓ, α
end