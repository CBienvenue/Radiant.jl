"""
    multigroup(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,state_of_matter::String,
    Eiᵇ::Vector{Float64},Efᵇ::Vector{Float64},L::Int64,interaction::Interaction,
    full_type::String,incoming_particle::Particle,scattered_particle::Particle,
    particles::Vector{Particle},interactions::Vector{Interaction})

Produce the multigroup macroscopic cross sections.

# Input Argument(s)
- `Z::Int64` : atomic number of the element.
- `ωz::Vector{Float64}` : weight fraction of the element(s) composing the material.
- `ρ::Float64` : density of the material [in g/cm³].
- `state_of_matter::String` : state of matter (solid, liquid or gaz).
- `Eiᵇ::Vector{Float64}`: energy boundaries of the incoming particle [in MeV].
- `Efᵇ::Vector{Float64}`: energy boundaries of the outgoing particle [in MeV].
- `L::Int64`: Legendre truncation order.
- `interaction::Interaction`: structure containing information about the interaction.
- `full_type::String` : type of interaction (scattering, production or absorption).
- `incoming_particle::Particle` : incoming particle in the interaction.
- `scattered_particle::Particle` : scattered particle in the interaction.
- `particles::Vector{Particle}` : list of particles involved in the interaction.
- `interactions::Vector{Interaction}` : list of all interactions that are taken into
  account for the cross-sections library.

# Output Argument(s)
- `Σsl::Array{Float64,3}`: Legendre moments of the differential cross section [in cm⁻¹].
- `Σt::Vector{Float64}`: total cross sections [in cm⁻¹].
- `Σa::Vector{Float64}`: absorption cross sections [in cm⁻¹].
- `Σe::Vector{Float64}`: energy deposition cross sections [in MeV × cm⁻¹].
- `Σc::Vector{Float64}`: charge deposition cross sections [in cm⁻¹].
- `Sb::Vector{Float64}`: stopping power at boundaries [MeV × cm⁻¹].
- `S::Vector{Float64}`: stopping power [MeV × cm⁻¹].
- `T::Vector{Float64}`: momentum transfer [in cm⁻¹].

# Reference(s)
- Lorence et al. (1989), Physics Guide to CEPXS: A Multigroup Coupled Electron-Photon
  Cross-Section Generating Code.
- Bienvenue et al. (2025), Toward highly accurate multigroup coupled photon-electron-positron
  cross-sections for the Boltzmann Fokker-Planck equation.

"""
function multigroup(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,state_of_matter::String,Eiᵇ::Vector{Float64},Efᵇ::Vector{Float64},L::Int64,interaction::Interaction,full_type::String,incoming_particle::Particle,scattered_particle::Particle,particles::Vector{Particle},interactions::Vector{Interaction},I_eff::Float64=NaN)

#----
# Initialization
#----
mₑc² = 0.510999
Nz = length(Z)
Ngi = length(Eiᵇ)-1; Ngf = length(Efᵇ)-1
Σt = zeros(Ngi); Σtₑ = zeros(Ngi); Σa = zeros(Ngi); Σe = zeros(Ngi+1); Σc = zeros(Ngi+1); Sb = zeros(Ngi+1); S = zeros(Ngi); T = zeros(Ngi)
Σsl = zeros(Ngi,Ngf,L+1); Σsₑ = zeros(Ngi,Ngf)
𝓕 = zeros(Ngf+1,L+1); 𝓕ₑ = zeros(Ngf+1)
charge_in = incoming_particle.get_charge()
charge_out = scattered_particle.get_charge()
q_deposited, q_extracted = get_charge_variation(interaction,full_type,charge_in,charge_out)
type = string(full_type[1])
scattering_model, is_CSD, is_AFP, is_AFP_decomposition, is_ETC = interaction.get_scattering_model()
is_elastic = is_elastic_scattering(interaction)
ΔQ = get_mass_energy_variation(interaction,type,false)
is_subshells = interaction.get_is_subshells_dependant()
E_in = Eiᵇ./mₑc²; E_out = Efᵇ./mₑc²; # Change of units (MeV → mₑc²)

#----
# Multigroup cross sections preparation
#----

# Incoming particle energy discretization
is_dirac, Np, q_type = in_distribution_dispatch(interaction)
if is_dirac Np = 1; u = [0]; w = [2] else u,w = quadrature(Np,q_type) end

# Compute cross sections, stopping powers and momentum transfers in each energy group
for gi in range(1,Ngi)

    # Incoming particle energy group
    Ei⁻ = E_in[gi]; Ei⁺ = E_in[gi+1]
    if (gi < Ngi) Ei²⁺ = E_in[gi+2] else Ei²⁺ = 0.0 end
    ΔEi = Ei⁻ - Ei⁺

    # Quadrature over energy group
    for ni in range(1,Np)

        # Incoming energy
        Ei = (u[ni]*ΔEi + (Ei⁻+Ei⁺))/2

        # Boundary between catastrophic and soft interactions
        Ec = soft_catastrophic_cutoff(Ei,Ei⁻,Ei⁺,Ei²⁺,scattering_model)

        # Scattering cross sections
        if type ∈ ["S","P"] && scattering_model != "FP"
            𝓕, 𝓕ₑ = feed(Z,ωz,ρ,L,Ei,E_out,Ngf,interaction,gi,Ngi,particles,full_type,incoming_particle,scattered_particle,E_in,Ec,is_elastic,is_subshells)
            if is_dirac 𝓕 ./= ΔEi; 𝓕ₑ ./= ΔEi end
            for gf in range(1,Ngf)
                Σsl[gi,gf,1:L+1] += w[ni]/2 * 𝓕[gf,1:L+1]
                Σsₑ[gi,gf] += w[ni]/2 * 𝓕ₑ[gf] 
            end
        end

        # Absorption and total cross sections
        if type ∈ ["S","A"] && scattering_model != "FP"
            Σtᵢ = 0.0
            for i in range(1,Nz)
                Σtᵢ += w[ni]/2 * tcs_dispatch(interaction,Ei,Z[i],Ec,i,incoming_particle,E_in[end],E_out) * nuclei_density(Z[i],ρ) * ωz[i]
            end
            if is_dirac Σtᵢ /= ΔEi end
            Σt[gi] += Σtᵢ
            Σtₑ[gi] += Σtᵢ * (Ei-ΔQ)
            Σa[gi] = Σtᵢ - w[ni]/2 * sum(𝓕[:,1])
        end

        # Momentum transfer
        if is_AFP && type == "S" && scattering_model != "BTE"
            T[gi] = 0.0
            for i in range(1,Nz)
                T[gi] += w[ni]/2 * mt_dispatch(interaction) * nuclei_density(Z[i],ρ) * ωz[i]
            end
            if is_dirac T[gi] ./= ΔEi end
        end

        # Stopping power
        if is_CSD && type == "S" && scattering_model != "BTE"
            S[gi] += w[ni]/2 * sp_dispatch(interaction,Z,ωz,ρ,state_of_matter,Ei,Ec,incoming_particle,E_out,I_eff)
            if is_dirac S[gi] ./= ΔEi end
        end
    end

    # Stopping power at energy group boundaries
    if is_CSD && type == "S" && scattering_model != "BTE"
        Ec = soft_catastrophic_cutoff(Ei⁻,Ei⁻,Ei⁺,Ei²⁺,scattering_model)
        Sb[gi] = sp_dispatch(interaction,Z,ωz,ρ,state_of_matter,Ei⁻,Ec,incoming_particle,E_out,I_eff)
        if (gi == Ngi)
            Ec = soft_catastrophic_cutoff(Ei⁺,Ei⁻,Ei⁺,Ei²⁺,scattering_model)
            Sb[gi+1] = sp_dispatch(interaction,Z,ωz,ρ,state_of_matter,Ei⁺,Ec,incoming_particle,E_out,I_eff)
        end
    end

    # Energy deposition cross sections
    Σe[gi] = Σtₑ[gi] - sum(Σsₑ[gi,:]) + S[gi]

    # Charge deposition cross sections
    Σc[gi] = Σt[gi] * q_deposited - sum(Σsl[gi,:,1]) * q_extracted

    # Extended transport corrections
    if scattering_model ∈ ["BFP","BTE"] && is_ETC
        Σt[gi],Σsl[gi,gi,:] = transport_correction(interaction,L,Σt[gi],Σsl[gi,gi,:])
    end

    # Elastic decomposition in soft and catastrophic components
    if scattering_model == "BFP" && is_AFP_decomposition
        Σt[gi],Σsl[gi,gi,:],T[gi] = angular_fokker_planck_decomposition(interaction,L,Σt[gi],Σsl[gi,gi,:],T[gi])
    end
end

# Contributions of particle slowing-down to the cutoff energy to the energy and charge deposition
if type ∈ ["S","A"]
    Σe[Ngi+1] = Sb[Ngi+1] * E_in[end]/(E_in[end-1]-E_in[end])
    Σc[Ngi+1] = Sb[Ngi+1] * (-charge_in)/(E_in[end-1]-E_in[end])
end

Σe *= mₑc²; Sb *= mₑc²; S *= mₑc² ; # Change of units (mₑc² → MeV)

return Σsl, Σt, Σa, Σe, Σc, Sb, S, T
end