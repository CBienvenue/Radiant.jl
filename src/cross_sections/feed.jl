"""
    feed(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,L::Int64,Ei::Float64,
    Eout::Vector{Float64},Ng::Int64,interaction::Interaction,gi::Int64,Ngi::Int64,
    particles::Vector{Particle},Npts::Int64,type::String,incoming_particle::Particle,
    scattered_particle::Particle,Ein::Vector{Float64},Ec::Float64)

Calculate the feed function 𝓕 (normalized probability of scattering from Ei into each
group gf) for each Legendremoments up to order L. Also calculate the energy weighted
feed function 𝓕ₑ for energy-deposition cross section.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic number of the element(s) composing the material.
- `ωz::Vector{Float64}` : weight fraction of the element(s) composing the material.
- `ρ::Float64` : density of the material [in g/cm³].
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : energy of the incoming particle [in mₑc²].
- `Eout::Vector{Float64}` : energy group boundaries [in mₑc²].
- `Ng::Int64` : number of groups.
- `interaction::Interaction` : interaction informations.
- `gi::Int64` : incoming particle group index.
- `Ngi::Int64` :  number of groups for the incoming particle.
- `particles::Vector{Particle}` : list of the particles imply in the interaction.
- `Npts::Int64` : number of points in the quadrature.
- `type::String` : type of interaction (scattering or production).
- `incoming_particle::Particle` : incoming particle.
- `scattered_particle::Particle` : scattered particle.
- `Ein::Vector{Float64}` : energy group boundaries corresponding to the incoming
  particle [in mₑc²].
- `Ec::Float64` : cutoff energy between soft and catastrophic interaction.

# Output Argument(s)
- `𝓕::Array{Float64}` : feed function.
- `𝓕ₑ::Vector{Float64}` : energy weighted feed function.

# Reference(s)
- MacFarlane et al. (2021) : The NJOY Nuclear Data Processing System, Version 2012.

"""
function feed(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,L::Int64,Ei::Float64,Eout::Vector{Float64},Ng::Int64,interaction::Interaction,gi::Int64,Ngi::Int64,particles::Vector{Particle},Npts::Int64,type::String,incoming_particle::Particle,scattered_particle::Particle,Ein::Vector{Float64},Ec::Float64)

#----
# Initialization
#----

𝓕 = zeros(Ng+1,L+1)
𝓕ₑ = zeros(Ng+1)

# Outgoing particle energy spectrum
is_dirac, Np, q_type = out_distribution_dispatch(interaction,type)
if is_dirac Np = 1; u = [0]; w = [2] else u,w = quadrature(Np,q_type) end

#----
# Feed function over all groups and under the cutoff energy
#----

# Loop over the coumpound elements
Nz = length(Z)
@inbounds for i in range(1,Nz)

    # Loop over subshells
    if interaction.is_subshells_dependant
        Nshells,Zi,Ui,Ti,ri,subshells = electron_subshells(Z[i])
    else
        Ui = [0.0]; Nshells = 1; Zi = Z[i]; Ti = [0.0]; ri = [0.0]; subshells = [""]; # Free atomic electron
    end

    for gf in range(1,Ng), δi in range(1,Nshells)
        
        # Final energy group
        Ef⁻ = Eout[gf]; Ef⁺ = Eout[gf+1]
        Ef⁻,Ef⁺,isSkip = bounds_dispatch(interaction,Ef⁻,Ef⁺,Ei,gi,gf,type,Ui[δi],Ec,incoming_particle)
        if isSkip continue end
        ΔEf = Ef⁻ - Ef⁺
        
        # Integration over the energy group
        𝓕i = zeros(L+1)
        𝓕iₑ = 0
        for n in range(1,Np)
            Ef = (u[n]*ΔEf + (Ef⁻+Ef⁺))/2
            if ΔEf > 0

                # Compute Legendre angular flux moments
                Σsᵢ = w[n]/2 .* dcs_dispatch(interaction,L,Ei,Ef,Z[i],scattered_particle,type,i,particles,Ein,Z,Ef⁻,Ef⁺,δi,Ui[δi],Zi[δi],Ti[δi],Ec,incoming_particle) * nuclei_density(Z[i],ρ) * ωz[i]

                if ~is_dirac Σsᵢ *= ΔEf  end
                𝓕i .+= Σsᵢ
                if typeof(interaction) == Annihilation && type ∈ ["P_inel","P_brems"]
                    𝓕iₑ += Σsᵢ[1] * (Ef-1)
                else
                    𝓕iₑ += Σsᵢ[1] * Ef
                end
            end
        end
        𝓕[gf,:] .+= 𝓕i
        𝓕ₑ[gf] += 𝓕iₑ
    end
end

return 𝓕, 𝓕ₑ
end