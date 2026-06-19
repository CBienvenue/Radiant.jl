"""
    feed(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,L::Int64,Ei::Float64,
    Eout::Vector{Float64},Ng::Int64,interaction::Interaction,gi::Int64,Ngi::Int64,
    particles::Vector{Particle},type::String,incoming_particle::Particle,
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
- `type::String` : type of interaction (scattering or production).
- `incoming_particle::Particle` : incoming particle.
- `scattered_particle::Particle` : scattered particle.
- `Ein::Vector{Float64}` : energy group boundaries corresponding to the incoming
  particle [in mₑc²].
- `Ec::Float64` : cutoff energy between soft and catastrophic interaction.
- `is_elastic_scattering::Bool` : boolean indicating if the scattering is elastic.
- `is_subshells::Bool` : boolean indicating if the cross-sections are subshells dependant.

# Output Argument(s)
- `𝓕::Array{Float64}` : feed function.
- `𝓕ₑ::Vector{Float64}` : energy weighted feed function.

# Reference(s)
- MacFarlane et al. (2021) : The NJOY Nuclear Data Processing System, Version 2012.

"""
function feed(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,L::Int64,Ei::Float64,Eout::Vector{Float64},Ng::Int64,interaction::Interaction,gi::Int64,Ngi::Int64,particles::Vector{Particle},type::String,incoming_particle::Particle,scattered_particle::Particle,Ein::Vector{Float64},Ec::Float64,is_elastic::Bool,is_subshells::Bool)

#----
# Initialization
#----
𝓕 = zeros(Ng+1,L+1)
𝓕ₑ = zeros(Ng+1)
ΔQ = get_mass_energy_variation(interaction,type,true)

# Outgoing particle energy spectrum
is_dirac, Np, q_type = out_distribution_dispatch(interaction,type)
if is_dirac Np = 1; u = [0]; w = [2] else u,w = quadrature(Np,q_type) end

#----
# Feed function over all groups and under the cutoff energy
#----

# Heavy inelastic S (same scattered particle)
is_heavy_inelastic_S = (is_proton(incoming_particle) || is_alpha(incoming_particle)) && (incoming_particle == scattered_particle)

# Loop over the coumpound elements
Nz = length(Z)
for i in range(1,Nz)

    nd = nuclei_density(Z[i],ρ) * ωz[i]

    # Loop over subshells and outgoing groups
    Nshells,Zi,Ui,Ti,ri,_ = electron_subshells(Z[i],~is_subshells)
    for gf in range(1,Ng), δi in range(1,Nshells)

        # Final energy group
        Ef⁻ = Eout[gf]; Ef⁺ = Eout[gf+1]
        Ef⁻,Ef⁺,isSkip = bounds_dispatch(interaction,Ef⁻,Ef⁺,Ei,gi,gf,type,Ui[δi],Ec,incoming_particle)
        if isSkip continue end
        ΔEf = Ef⁻ - Ef⁺

        # Integration over the energy group
        𝓕i = zeros(L+1)
        𝓕iₑ = 0

        # For heavy particles, compute the analytic singular contribution once per group
        analytic_A = 0.0
        M₁ = 0.0
        cache = nothing
        if is_heavy_inelastic_S
            cache = HeavyInelasticCache(Zi[δi], Ei, incoming_particle)
            analytic_A = integrate_A_over_W2_per_subshell(cache, Ef⁻, Ef⁺) * nd
            M₁ = feed_first_moment_heavy_particle(cache, Ef⁻, Ef⁺) * nd
        end

        # Quadrature integration (regular remainder for heavy particles, full for others)
        for n in range(1,Np)

            # Outgoing particle energy group
            if (is_elastic) Ef = Ei else Ef = (u[n]*ΔEf + (Ef⁻+Ef⁺))/2 end

            # Compute Legendre angular flux moments
            Σsᵢ = ΔEf .* w[n]/2 .* dcs_dispatch(interaction,L,Ei,Ef,Z[i],scattered_particle,type,i,particles,Ein,Ef⁻,Ef⁺,δi,Ui[δi],Zi[δi],Ti[δi],ri[δi],Ec,incoming_particle) * nd
            if is_dirac Σsᵢ /= ΔEf  end
            𝓕i .+= Σsᵢ
            𝓕iₑ += Σsᵢ[1] * (Ef + ΔQ)
        end

        # Add analytic singular contribution for heavy particles
        if is_heavy_inelastic_S
            𝓕i .+= analytic_A .* ones(L+1)
            for l in range(0,L)
                𝓕i[l+1] += integrate_leading_1overW_per_subshell(cache, Ef⁻, Ef⁺, l) * nd
            end
            σ_analytic = feed_analytical_heavy_particle(cache, Ef⁻, Ef⁺) * nd
            𝓕iₑ = Ei * σ_analytic - M₁

            # Consistency check: ensure analytic σ equals numeric l=0 moment
            if ~isapprox(𝓕i[1], σ_analytic; rtol=1e-3, atol=1e-12)
                rel = abs(𝓕i[1] - σ_analytic) / max(abs(σ_analytic), 1e-300)
                @warn "Heavy inelastic feed: analytic/numeric mismatch" Z=Z[i] δi=δi gf=gf σ_analytic=σ_analytic numeric_l0=𝓕i[1] rel_diff=rel
                𝓕i[1] = σ_analytic
            end
        end
        𝓕[gf,:] .+= 𝓕i
        𝓕ₑ[gf] += 𝓕iₑ
    end
end
return 𝓕, 𝓕ₑ
end