
mutable struct Compton <: Interaction

    # Variable(s)
    name::String
    incoming_particle::Vector{String}
    interaction_particles::Vector{String}
    interaction_types::Dict{Tuple{String,String},Vector{String}}
    energy_integration_method::String
    Cℓk::Array{Float64}
    Cℓki::Array{Float64}
    is_CSD::Bool
    is_AFP::Bool
    is_elastic::Bool
    is_preload_data::Bool
    is_subshells_dependant::Bool
    incoherent_scattering_factor::Function

    # Constructor(s)
    function Compton()
        this = new()
        this.name = "compton"
        this.interaction_types = Dict(("photons","photons") => ["S"],("photons","electrons") => ["P"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.energy_integration_method = "quadrature"
        this.is_CSD = false
        this.is_AFP = false
        this.is_elastic = false
        this.is_preload_data = true
        this.is_subshells_dependant = false
        return this
    end
end

# Method(s)
function in_distribution(this::Compton)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

function out_distribution(this::Compton)
    if this.energy_integration_method == "analytic"
        is_dirac = true
        N = 1
        quadrature = "dirac"
    elseif this.energy_integration_method == "quadrature"
        is_dirac = false
        N = 8
        quadrature = "gauss-legendre"
    else
        error("Unknown energy integration method.")
    end
    return is_dirac, N, quadrature
end

function bounds(this::Compton,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String)
    # Scattered photon
    if type == "S" 
        Ef⁻ = min(Ei,Ef⁻)
        Ef⁺ = max(Ei/(1+2*Ei),Ef⁺)
        if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    # Produced electron
    elseif type == "P" 
        Ef⁻ = min(2*Ei^2/(1+2*Ei),Ef⁻)
        if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    else
        error("Unknown type of method for Klein-Nishina scattering.")
    end
    return Ef⁻,Ef⁺,isSkip
end

function dcs(this::Compton,L::Int64,Ei::Float64,Ef::Float64,type::String,Efmax::Float64,Efmin::Float64,Z::Int64,iz::Int64)

    # Initialization
    rₑ = 2.81794092e-13 # (in cm)
    σs = 0
    σℓ = zeros(L+1)

    # Scattered photon
    if type == "S"

        # Exact integration
        if this.energy_integration_method == "analytic"
            Γ = Z * π*rₑ^2/Ei^2
            αi = [1,Ei-2-2/Ei,2/Ei+1/Ei^2,1/Ei]
            𝒢γ = zeros(L+4)
            @inbounds for i in range(-(L+2),1)
                if i == -1
                    𝒢γ[(L+2)+i+1] = log(Efmax) - log(Efmin)
                else
                    𝒢γ[(L+2)+i+1] = ((Efmax)^(i+1) - (Efmin)^(i+1))/(i+1)
                end
            end
            @inbounds for ℓ in range(0,L)
                for k in range(0,div(ℓ,2))
                    σℓk = 0.0
                    for i in range(0,ℓ-2k)
                        σℓki = 0.0
                        for j in range(0,3)
                            σℓki += αi[j+1] * 𝒢γ[(L+2)+(j-i-2)+1]
                        end
                        σℓk += (1+1/Ei)^(ℓ-2*k-i) * this.Cℓki[ℓ+1,k+1,i+1] * σℓki
                    end 
                    σℓ[ℓ+1] += this.Cℓk[ℓ+1,k+1] * σℓk
                end
                σℓ[ℓ+1] *= Γ/(2^ℓ)
            end
            # Correction to deal with high-order Legendre moments
            for ℓ in range(1,L)
                if abs(σℓ[1]) < abs(σℓ[ℓ+1])
                    σℓ[ℓ+1:end] .= 0.0
                    break
                end
            end
        # Quadrature-based integration
        elseif this.energy_integration_method == "quadrature"
            # Compute the differential scattering cross section
            if Ei/(1+2*Ei) ≤ Ef ≤ Ei
                σs = π * rₑ^2 / Ei^2 * (Ei/Ef + Ef/Ei - 2*(1/Ef-1/Ei) + (1/Ef-1/Ei)^2)
            end
            # Compute the Legendre moments of the flux
            if σs != 0
                μ = max(min(1 + 1/Ei - 1/Ef,1),-1)
                S = this.incoherent_scattering_factor(iz,Ei,μ)
                Pℓμ = legendre_polynomials(L,μ)
                for ℓ in range(0,L) σℓ[ℓ+1] += Pℓμ[ℓ+1] * σs * S end
            end
        else
            error("Unknown energy integration method.")
        end

    # Produced electron
    elseif type == "P"
        # Exact integration
        if this.energy_integration_method == "analytic"
            Γ = Z*π*rₑ^2/Ei^2
            αi = [1,Ei^2-2*Ei-2,-1/Ei,1+2/Ei+1/Ei^2,-(Ei-2-2/Ei)]
            @inbounds for ℓ in range(0,L)
                for k in range(0,div(ℓ,2))
                    σℓk = 0.0
                    if mod(ℓ-2k,2) == 0
                        m = div(ℓ-2*k,2)
                        for i in range(0,1)
                            σℓk += αi[i+1] * (𝒢₃(i-2,-m,1,2,-1,Ei,1/Efmin)-𝒢₃(i-2,-m,1,2,-1,Ei,1/Efmax))
                        end
                        for i in range(0,2)
                            σℓk += αi[i+3] * (𝒢₃(i-3,-m,1,2,0,1,1/Efmin)-𝒢₃(i-3,-m,1,2,0,1,1/Efmax))
                        end
                    else
                        m = div(ℓ-2*k-1,2)
                        for i in range(0,1)
                            σℓk += αi[i+1] * (𝒢₄(2-i,m,1,2,-1,Ei,1/Efmin)-𝒢₄(2-i,m,1,2,-1,Ei,1/Efmax))
                        end
                        for i in range(0,2)
                            σℓk += αi[i+3] * (𝒢₄(3-i,m,1,2,0,1,1/Efmin)-𝒢₄(3-i,m,1,2,0,1,1/Efmax))
                        end
                    end
                    σℓ[ℓ+1] += (1+1/Ei)^(ℓ-2*k) * this.Cℓk[ℓ+1,k+1] * σℓk
                end
                σℓ[ℓ+1] *= Γ/(2^ℓ)
            end 
            # Correction to deal with high-order Legendre moments
            for ℓ in range(1,L)
                if abs(σℓ[1]) < abs(σℓ[ℓ+1])
                    σℓ[ℓ+1:end] .= 0.0
                    break
                end
            end
        # Quadrature-based integration
        elseif this.energy_integration_method == "quadrature"
            # Change of variable
            Eₑ = copy(Ef)
            Ef = Ei - Eₑ
            # Compute the differential scattering cross section
            if 0 ≤ Eₑ ≤ 2*Ei^2/(1+2*Ei)
                σs = π * rₑ^2 / Ei^2 * (Ei/Ef + Ef/Ei - 2*(1/Ef-1/Ei) + (1/Ef-1/Ei)^2)
            end
            # Compute the Legendre moments of the flux
            if σs != 0
                μ = max(min((1 + Ei)/Ei * 1/sqrt(2/Eₑ+1),1),-1)
                μγ = max(min(1 + 1/Ei - 1/Ef,1),-1)
                S = this.incoherent_scattering_factor(iz,Ei,μγ)
                Pℓμ = legendre_polynomials(L,μ)
                for ℓ in range(0,L) σℓ[ℓ+1] += Pℓμ[ℓ+1] * σs * S end
            end
        else
            error("Unknown energy integration method.")
        end
    else
        error("Unknown interaction.")
    end
    return σℓ
end

function tcs(this::Compton,Ei::Float64,Z::Int64,Eout::Vector{Float64},iz::Int64)

    #=
    rₑ = 2.81794092E-13 # (in cm)
    f_klein_nishina(x) = (Ei-2-2/Ei)*log(x) + (1/(2*Ei))*x^2 + (2/Ei+1/Ei^2)*x - 1/x
    σt = Z*π*rₑ^2/Ei^2*( f_klein_nishina(Ei) - f_klein_nishina(Ei/(1+2*Ei)) )
    return σt
    =#

    σt = 0.0
    rₑ = 2.81794092E-13 # (in cm)
    Ngf = length(Eout)-1
    is_dirac, Np, q_type = out_distribution(this)
    if is_dirac Np = 1; u = [0]; w = [2] else u,w = quadrature(Np,q_type) end
    for gf in range(1,Ngf+1)
        Ef⁻ = Eout[gf]
        if (gf != Ngf+1) Ef⁺ = Eout[gf+1] else Ef⁺ = 0.0 end
        Ef⁻,Ef⁺,isSkip = bounds(this,Ef⁻,Ef⁺,Ei,"S")
        if isSkip continue end
        ΔEf = Ef⁻ - Ef⁺
        for n in range(1,Np)
            Ef = (u[n]*ΔEf + (Ef⁻+Ef⁺))/2

            # Compute the differential scattering cross section
            if Ei/(1+2*Ei) ≤ Ef ≤ Ei
                μ = max(min(1 + 1/Ei - 1/Ef,1),-1)
                S = this.incoherent_scattering_factor(iz,Ei,μ)
                σs = S * π * rₑ^2 / Ei^2 * (Ei/Ef + Ef/Ei - 2*(1/Ef-1/Ei) + (1/Ef-1/Ei)^2)
            end
            
            # Cross-sections
            σt += ΔEf/2 * w[n] * σs
        end
    end

    return σt
    
end

function preload_data(this::Compton,L::Int64,Z::Vector{Int64})

    # Precompute angular integration factors
    if this.energy_integration_method == "analytic"
        this.Cℓk = zeros(L+1,div(L,2)+1)
        this.Cℓki = zeros(L+1,div(L,2)+1,L+1)
        for ℓ in range(0,L), k in range(0,div(L,2))
            this.Cℓk[ℓ+1,k+1] = (-1)^k * exp( sum(log.(1:2*ℓ-2*k)) - sum(log.(1:k)) - sum(log.(1:ℓ-k)) - sum(log.(1:ℓ-2*k)) )
            for i in range(0,L)
                this.Cℓki[ℓ+1,k+1,i+1] = (-1)^i * exp( sum(log.(1:ℓ-2k)) -  sum(log.(1:i)) -  sum(log.(1:ℓ-2k-i)) )
            end
        end
    end
    
    # Incoherent scattering factor
    Nz = length(Z)
    x = Vector{Vector{Float64}}(undef,Nz)
    S = Vector{Vector{Float64}}(undef,Nz)
    S_spline = Vector{Function}(undef,Nz)
    for iz in range(1,Nz)
        data = load("./data/compton_factors_JENDL5.jld2")
        x[iz] = data["x"][Z[iz]]
        S[iz] = data["F"][Z[iz]]
        S_spline[iz] = cubic_hermite_spline(x[iz],S[iz])
    end
    this.incoherent_scattering_factor = function incoherent_scattering_factor(iz::Int64,Ei::Float64,μ::Float64)
        hc = 1/20.60744 # (hc in mₑc² × Å)
        xi = 2*Ei/(hc)*sqrt((1-μ)/2) * sqrt((1+(Ei^2+2*Ei)*(1-μ)/2))/(1+Ei*(1-μ))
        Si = S_spline[iz](xi)
        return Si
    end
    
end
