"""
    multigroup(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,state_of_matter::String,
    Eiᵇ::Vector{Float64},Efᵇ::Vector{Float64},L::Int64,interaction::Interaction,
    full_type::String,incoming_particle::Particle,scattered_particle::Particle,
    particles::Vector{Particle},Npts::Int64,isStandard::Bool,
    interactions::Vector{Interaction})

Produce the multigroup macroscopic cross sections.

# Input Argument(s)
- 'Z::Int64' : atomic number of the element.
- 'ωz::Vector{Float64}' : weight fraction of the element(s) composing the material.
- 'ρ::Float64' : density of the material [in g/cm³].
- 'state_of_matter::String' : state of matter (solid, liquid or gaz).
- 'Eiᵇ::Vector{Float64}': energy boundaries of the incoming particle [in MeV].
- 'Efᵇ::Vector{Float64}': energy boundaries of the outgoing particle [in MeV].
- 'L::Int64': Legendre truncation order.
- 'interaction::Interaction': structure containing information about the interaction.
- 'full_type::String' : type of interaction (scattering, production or absorption).
- 'incoming_particle::Particle' : incoming particle in the interaction.
- 'scattered_particle::Particle' : scattered particle in the interaction.
- 'particles::Vector{Particle}' : list of particles involved in the interaction.
- 'Npts::Int64' : number of points in the quadrature.
- 'isStandard::Bool' : indicate if a new interaction is defined or not when calling
  multigroup function.
- 'interactions::Vector{Interaction}' : list of all interactions that are taken into
  account for the cross-sections library.

# Output Argument(s)
- 'Σsℓ::Array{Float64,3}': Legendre moments of the differential cross section [in cm⁻¹].
- 'Σt::Vector{Float64}': total cross sections [in cm⁻¹].
- 'Σa::Vector{Float64}': absorption cross sections [in cm⁻¹].
- 'Σs::Vector{Float64}': secondary production cross sections [in cm⁻¹].
- 'Σe::Vector{Float64}': energy deposition cross sections [in MeV × cm⁻¹].
- 'Σc::Vector{Float64}': charge deposition cross sections [in cm⁻¹].
- 'S::Vector{Float64}': stopping power [MeV × cm⁻¹].
- 'T::Vector{Float64}': momentum transfer [in cm⁻¹].

# Reference(s)
- Lorence (1989), Physics Guide to CEPXS: A Multigroup Coupled Electron-Photon
  Cross-Section Generating Code.

"""
function multigroup(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,state_of_matter::String,Eiᵇ::Vector{Float64},Efᵇ::Vector{Float64},L::Int64,interaction::Interaction,full_type::String,incoming_particle::Particle,scattered_particle::Particle,particles::Vector{Particle},Npts::Int64,isStandard::Bool,interactions::Vector{Interaction})

if isStandard
    println("Start of ",interaction.name," calculations.") 
end

# Initialization
mₑc² = 0.510999
Ngi = length(Eiᵇ)-1; Ngf = length(Efᵇ)-1
Σt = zeros(Ngi); Σtₑ = zeros(Ngi); Σa = zeros(Ngi); Σs = zeros(Ngi); Σe = zeros(Ngi+1); Σc = zeros(Ngi+1); S = zeros(Ngi+1); Sm = zeros(Ngi); T = zeros(Ngi)
Σsℓ = zeros(Ngi,Ngf,L+1); Σsₑ = zeros(Ngi,Ngf)
𝓕 = zeros(Ngf+1,L+1); 𝓕ₑ = zeros(Ngf+1)
charge_in = incoming_particle.get_charge()
charge_out = scattered_particle.get_charge()
type = string(full_type[1])


# Multigroup cross sections preparation

# Change of units (MeV → mₑc²)
E_in = Eiᵇ./mₑc²; E_out = Efᵇ./mₑc²;

# Incoming particle energy spectrum
is_dirac, Np, q_type = in_distribution_dispatch(interaction)
if is_dirac Np = 1; u = [0]; w = [2] else u,w = quadrature(Np,q_type) end

# Preloading data for calculations
if (interaction.is_preload_data) preload_data_dispatch(interaction,Z,E_in[1],E_in[end],L,ωz,ρ,E_out,incoming_particle,full_type,interactions) end

# Compute cross sections in each energy group
@inbounds for gi in range(1,Ngi)

    # Initial energy group
    Ei⁻ = E_in[gi]; Ei⁺ = E_in[gi+1]
    if (gi < Ngi) Ei²⁺ = E_in[gi+2] else Ei²⁺ = 0.0 end
    ΔEi = Ei⁻ - Ei⁺

    # Total and scattering differential cross section
    for ni in range(1,Np)

        # Incoming energy
        Ei = (u[ni]*ΔEi + (Ei⁻+Ei⁺))/2

        # Boundary between catastrophic and soft interactions
        Ec = Ei * (Ei⁺-Ei²⁺)/ (Ei⁻-Ei⁺) - (Ei⁺^2-Ei⁻*Ei²⁺)/(Ei⁻-Ei⁺)
        if (interaction.scattering_model == "FP") Ec = 0.0 end

        # Total cross sections
        if type != "P" && ~(interaction.scattering_model == "FP" && type == "S")
            Nz = length(Z)
            Σtᵢ = 0.0
            for i in range(1,Nz)
                Σtᵢ += 1/2 * w[ni] * tcs_dispatch(interaction,Ei,Z[i],Ec,i,incoming_particle,E_in[end],E_out,full_type) * nuclei_density(Z[i],ρ) * ωz[i]
            end
            if is_dirac Σtᵢ /= ΔEi end
            Σt[gi] += Σtᵢ
            if (~interaction.is_elastic)
                if typeof(interaction) == Pair_Production
                    Σtₑ[gi] += Σtᵢ * (Ei-2)
                elseif typeof(interaction) == Annihilation
                    Σtₑ[gi] += Σtᵢ * (Ei+2)
                else
                    Σtₑ[gi] += Σtᵢ * Ei
                end
            end
        end

        # Scattering cross sections
        if type == "A" continue end # No scattering for absorption interaction
        if ~(interaction.scattering_model == "FP" && type == "S")
            𝓕, 𝓕ₑ = feed(Z,ωz,ρ,L,Ei,E_out,Ngf,interaction,gi,Ngi,particles,Npts,full_type,incoming_particle,scattered_particle,E_in,Ec)
            if is_dirac 𝓕 ./= ΔEi; 𝓕ₑ ./= ΔEi end
            for gf in range(1,Ngf)
                Σsℓ[gi,gf,:] += 1/2 * w[ni] * 𝓕[gf,:]
                if (~interaction.is_elastic) Σsₑ[gi,gf] += 1/2 * w[ni] * 𝓕ₑ[gf] end
            end
        end

        # Momentum transfer
        if  (interaction.name == "mott" && interaction.scattering_model == "FP") || (interaction.is_CSD && type != "P")
            Nz = length(Z)
            T[gi] = 0.0
            for i in range(1,Nz)
                T[gi] += 1/2 * w[ni] * mt_dispatch(interaction) * nuclei_density(Z[i],ρ) * ωz[i]
            end
            if is_dirac T[gi] ./= ΔEi end
        end

        # Stopping power
        if (interaction.is_CSD) && type != "P"
            Sm[gi] += 1/2 * w[ni] * sp_dispatch(interaction,Z,ωz,ρ,state_of_matter,Ei,Ec,incoming_particle,E_out)
            if is_dirac Sm[gi] ./= ΔEi end
        end

    end

    # Stopping power at boundaries
    if (interaction.is_CSD) && type != "P" && full_type != "Aₐ"
        S[gi] = sp_dispatch(interaction,Z,ωz,ρ,state_of_matter,Ei⁻,Ei⁺,incoming_particle,E_out)
        if (gi == Ngi) S[gi+1] += sp_dispatch(interaction,Z,ωz,ρ,state_of_matter,Ei⁺,0.0,incoming_particle,E_out) end
    end

    # Elastic transport corrections
    Σt[gi],Σsℓ[gi,gi,:],T[gi] = transport_correction(interaction,L,Σt[gi],Σsℓ[gi,gi,:],T[gi],interaction.scattering_model)

end

if typeof(interaction) == Annihilation
    if full_type == "P_pp"
        α_p = 0
        β_p = 1/2
    else
        α_p = 1
        β_p = 0
    end
elseif typeof(interaction) == Pair_Production
    α_p = 0
    β_p = 0 # no electron extracted from medium in pair production
else
    α_p = -charge_in
    β_p = -charge_out
end

@inbounds for gi in range(1,Ngi)

    if type == "A"

        # Absorption cross section
        Σa[gi] = Σt[gi]

        # Energy deposition cross sections
        Σe[gi] = Σtₑ[gi]

        # No secondary production cross section
        # ∅

        # Charge deposition cross sections
        Σc[gi] = Σt[gi] * α_p

    elseif type == "S"

        # Absorption cross section
        Σa[gi] = Σt[gi] - sum(Σsℓ[gi,:,1])

        # Energy deposition cross sections
        Σe[gi] = Σtₑ[gi] - sum(Σsₑ[gi,:]) + Sm[gi]

        # No secondary production cross section
        # ∅

        # Charge deposition cross sections
        Σc[gi] = Σt[gi] * α_p - sum(Σsℓ[gi,:,1]) * β_p

    elseif type == "P"

        # Total cross-section
        Σt[gi] = 0.0

        # No absorption cross section
        # ∅

        # Energy deposition cross sections
        Σe[gi] = -sum(Σsₑ[gi,:])

        # Secondary production cross section
        Σs[gi] = sum(Σsℓ[gi,:,1])

        # Charge deposition cross sections
        Σc[gi] = -sum(Σsℓ[gi,:,1]) * β_p

    end

    # Particle conservation
    if Σa[gi] != Σt[gi] + Σs[gi] - sum(Σsℓ[gi,:,1])
        error("Particle conservation is not satisfied: ",[Σa[gi],Σt[gi],Σs[gi],sum(Σsℓ[gi,:,1])])
    end

end

# Slowing-down under the cutoff contributions
if type != "P"
    Σe[Ngi+1] = S[Ngi+1] * E_in[end]/(E_in[end-1]-E_in[end])
    Σc[Ngi+1] = S[Ngi+1] * (-charge_in)/(E_in[end-1]-E_in[end])
end

# Change of units (mₑc² → MeV)
Σe *= mₑc²; S *= mₑc²;

if isStandard println("End of ",interaction.name," calculations."); println() end

return Σsℓ, Σt, Σa, Σs, Σe, Σc, S, T
end