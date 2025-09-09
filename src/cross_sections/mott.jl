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
- `σl::Vector{Float64}` : Legendre moments of the Mott cross-sections.

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
    σl = zeros(L+1)

    #----
    # Import and compute data
    #----

    # Compute and save in cache, if not already in cache
    if ~haskey(cache_radiant[],"Clk") || length(cache_radiant[]["Clk"][:,1]) < L+1 
        Clk = zeros(L+1,div(L,2)+1)
        for l in range(0,L), k in range(0,div(L,2))
            Clk[l+1,k+1] = (-1)^k * exp( sum(log.(1:2*l-2*k)) - sum(log.(1:k)) - sum(log.(1:l-k)) - sum(log.(1:l-2*k)) )
        end
        cache_radiant[]["Clk"] = Clk
    end
    if ~haskey(cache_radiant[],"Clki") || length(cache_radiant[]["Clki"][:,1,1,1]) < L+1
        Clki = zeros(L+1,div(L,2)+1,2,L+2)
        for l in range(0,L), k in range(0,div(L,2)), i in range(0,1), g in range(0,l-2*k+i)
            Clki[l+1,k+1,i+1,g+1] = 2 * exp( sum(log.(1:l-2*k+i)) - sum(log.(1:g)) - sum(log.(1:l-2*k+i-g)) ) * (-1)^g
        end
        cache_radiant[]["Clki"] = Clki
    end

    # Extract data
    Clk = cache_radiant[]["Clk"]
    Clki = cache_radiant[]["Clki"]
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
    for l in range(0,L)
        for k in range(0,div(l,2))
            σlk = 0.0
            # Compute I₁ -----
            I₁ = zeros(3)
            for i in range(0,2)
                I₁[i+1] += 𝒢₁⁺[l-2*k+i+1] - 𝒢₁⁻[l-2*k+i+1]
                I₁[i+1] *= αi[i+1]
            end 
            σlk += sum(I₁)
            # Compute I₂ -----
            I₂ = zeros(2)
            for i in range(0,1)
                for g in range(0,l-2*k+i)
                    I₂[i+1] += Clki[l+1,k+1,i+1,g+1] * ( 𝒢₂⁺[2*(1+g)] - 𝒢₂⁻[2*(1+g)])
                end
                I₂[i+1] *= αi[i+4]
            end
            σlk += sum(I₂)
            σl[l+1] += Clk[l+1,k+1] * σlk
        end
        σl[l+1] *= 1/(2^l)
    end
    Γ = 2*π*rₑ^2*Z*(Z+ξ)/(β²*Ei*(Ei+2))
    σl .*= Γ

    #----
    # Correction to deal with high-order Legendre moments instabilities
    #----
    for l in range(1,L)
        if abs(σl[1]) < abs(σl[l+1])
            σl[l+1:end] .= 0.0
            break
        end
    end

    return σl
end