"""
    mott(Z::Int64,Ei::Float64,particle::Particle,L::Int64,
    Ecutoff::Union{Missing,Float64}=missing,is_seltzer_correction::Bool=true,
    is_kawrakow_correction::Bool=true,is_subshell_inelastic::Bool=true,
    model::String="boschini")

Gives the Legendre moments of the Mott cross-sections.

# Input Argument(s)
- `Z::Int64` : atomic number.
- `Ei::Float64` : incoming particle energy.
- `particle::Particle` : incoming particle.
- `L::Int64` : Legendre truncation order.
- `Ecutoff::Union{Missing,Float64}` : cutoff energy (lowest energy for transport).
- `is_seltzer_correction::Bool` : boolean to enable or not the Seltzer correction to
  Molière screening factor.
- `is_kawrakow_correction::Bool` : boolean to enable or not the Kawrakow correction.
- `is_subshell_inelastic::Bool` : boolean to enable or not subshell-dependant Kawrakow
  correction.
- `model::String` : model, which can have the following values:
    - `"boschini"` : screened Mott based on Boschini parameters.
    - `"rutherford"` : screened Rutherford.

# Output Argument(s)
- `σℓ::Vector{Float64}` : Legendre moments of the Mott cross-sections.

# Reference(s)
- Boschini et al. (2013), An expression for the Mott cross section of electrons and
  positrons on nuclei with Z up to 118.
- Bienvenue et al. (2025), Toward highly accurate multigroup coupled 
  photon-electron-positron cross-sections for the Boltzmann Fokker-Planck equation.

"""
function mott(Z::Int64,Ei::Float64,particle::Particle,L::Int64,Ecutoff::Union{Missing,Float64}=missing,is_seltzer_correction::Bool=true,is_kawrakow_correction::Bool=true,is_subshell_inelastic::Bool=true,model::String="boschini")

    #----
    # Initialization
    #----
    rₑ = 2.81794092E-13 # (in cm)
    β² = Ei*(Ei+2)/(Ei+1)^2
    β = sqrt(β²)
    β₀ = 0.7181287
    η = moliere_screening(Z,Ei,is_seltzer_correction)
    σℓ = zeros(L+1)

    #----
    # Import and compute data
    #----

    # Compute and save in cache, if not already in cache
    if ~haskey(cache_radiant[],"Cℓk") || length(cache_radiant[]["Cℓk"][:,1]) < L+1 
        Cℓk = zeros(L+1,div(L,2)+1)
        for ℓ in range(0,L), k in range(0,div(L,2))
            Cℓk[ℓ+1,k+1] = (-1)^k * exp( sum(log.(1:2*ℓ-2*k)) - sum(log.(1:k)) - sum(log.(1:ℓ-k)) - sum(log.(1:ℓ-2*k)) )
        end
        cache_radiant[]["Cℓk"] = Cℓk
    end
    if ~haskey(cache_radiant[],"Cℓki") || length(cache_radiant[]["Cℓki"][:,1,1,1]) < L+1
        Cℓki = zeros(L+1,div(L,2)+1,2,L+2)
        for ℓ in range(0,L), k in range(0,div(L,2)), i in range(0,1), g in range(0,ℓ-2*k+i)
            Cℓki[ℓ+1,k+1,i+1,g+1] = 2 * exp( sum(log.(1:ℓ-2*k+i)) - sum(log.(1:g)) - sum(log.(1:ℓ-2*k+i-g)) ) * (-1)^g
        end
        cache_radiant[]["Cℓki"] = Cℓki
    end

    # Extract data
    Cℓk = cache_radiant[]["Cℓk"]
    Cℓki = cache_radiant[]["Cℓki"]
    if model == "boschini"
        data = fast_load("mott_data_boschini_2013.jld2")
        if is_electron(particle)
            particle_name = "electrons"
        elseif is_positron(particle)
            particle_name = "positrons"
        else
            error("Unknown particle")
        end
        if ~haskey(data,particle_name) error("Mott cross-sections are only available for electrons and positrons.") end
        bjk = data[particle_name][Z]
        ai = zeros(5)
        for j in range(0,4), k in range(1,6)
            ai[j+1] += bjk[j+1,k] * (β-β₀)^(k-1)
        end
    elseif model == "rutherford"
        ai = zeros(5)
        ai[1] = 1
    else
        error("Unknown elastic model $model.")
    end
    αi = [ai[1]+ai[3]+ai[5],-(ai[3]+2*ai[5]),ai[5],ai[2]+ai[4],-ai[4]]

    #----
    # Kawrakow correction
    #----
    ξ = 0
    if is_kawrakow_correction
        if ismissing(Ecutoff) error("The cutoff energy is required for Kawrakow correction.") end
        ξ = kawrakow_correction(Z,Ei,Ecutoff,η,particle,ai,model,is_subshell_inelastic)
    end

    #----
    # Legendre moments of the scattering cross-sections.
    #----
    μmax = 1; μmin = -1
    𝒢₁⁺ = 𝒢₁(μmax,2+L,1+2*η,-1)
    𝒢₁⁻ = 𝒢₁(μmin,2+L,1+2*η,-1)
    𝒢₂⁺ = 𝒢₂(sqrt(1-μmin),2*(2+L),2*η,1)
    𝒢₂⁻ = 𝒢₂(sqrt(1-μmax),2*(2+L),2*η,1)
    for ℓ in range(0,L)
        for k in range(0,div(ℓ,2))
            σℓk = 0.0
            # Compute I₁ -----
            I₁ = zeros(3)
            for i in range(0,2)
                I₁[i+1] += 𝒢₁⁺[ℓ-2*k+i+1] - 𝒢₁⁻[ℓ-2*k+i+1]
                I₁[i+1] *= αi[i+1]
            end 
            σℓk += sum(I₁)
            # Compute I₂ -----
            I₂ = zeros(2)
            for i in range(0,1)
                for g in range(0,ℓ-2*k+i)
                    I₂[i+1] += Cℓki[ℓ+1,k+1,i+1,g+1] * ( 𝒢₂⁺[2*(1+g)] - 𝒢₂⁻[2*(1+g)])
                end
                I₂[i+1] *= αi[i+4]
            end
            σℓk += sum(I₂)
            σℓ[ℓ+1] += Cℓk[ℓ+1,k+1] * σℓk
        end
        σℓ[ℓ+1] *= 1/(2^ℓ)
    end
    Γ = 2*π*rₑ^2*Z*(Z+ξ)/(β²*Ei*(Ei+2))
    σℓ .*= Γ

    #----
    # Correction to deal with high-order Legendre moments instabilities
    #----
    for ℓ in range(1,L)
        if abs(σℓ[1]) < abs(σℓ[ℓ+1])
            σℓ[ℓ+1:end] .= 0.0
            break
        end
    end

    return σℓ
end