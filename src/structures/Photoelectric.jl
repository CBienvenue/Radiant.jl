"""
    Photoelectric

Structure used to define parameters for production of multigroup photoelectric cross-sections.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{String,String},Vector{String}} = Dict(("photons","photons") => ["A"],("photons","electrons") => ["P"])`: Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `("photons","photons") => ["A"]`: absorption of incoming photon.
    - `("photons","electrons") => ["P"]`: produced photo-electron.

"""
mutable struct Photoelectric <: Interaction

    # Variable(s)
    name::String
    incoming_particle::Vector{String}
    interaction_particles::Vector{String}
    interaction_types::Dict{Tuple{String,String},Vector{String}}
    is_CSD::Bool
    is_AFP::Bool
    is_elastic::Bool
    is_preload_data::Bool
    is_subshells_dependant::Bool
    photoelectric_cross_sections::Function
    model::String
    Cℓk::Array{Float64}

    # Constructor(s)
    function Photoelectric()
        this = new()
        this.name = "photoelectric"
        this.interaction_types = Dict(("photons","photons") => ["A"],("photons","electrons") => ["P"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.is_CSD = false
        this.is_AFP = false
        this.is_elastic = false
        this.is_preload_data = true
        this.model = "jendl5"
        if this.model == "jendl5"
            this.is_subshells_dependant = true
        elseif this.model == "biggs_lighthill"
            this.is_subshells_dependant = false
        else    
            error("Unknown photoelectric model.")
        end
        return this
    end

end

# Method(s)
"""
    set_interaction_types(this::Photoelectric,interaction_types::Dict{Tuple{String,String},Vector{String}})

To define the interaction types for photoelectric processes.

# Input Argument(s)
- `this::Photoelectric`: photoelectric structure.
- `interaction_types::Dict{Tuple{String,String},Vector{String}}`: Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `("photons","photons") => ["A"]`: absorption of incoming photon.
    - `("photons","electrons") => ["P"]`: produced photo-electron.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> photoelectric = Photoelectric()
julia> photoelectric.set_interaction_types( Dict(("photons","photons") => ["A"],("photons","electrons") => ["P"]) ) # Full photoelectric phenomenon (default case).
```
"""
function set_interaction_types(this::Photoelectric,interaction_types::Dict{Tuple{String,String},Vector{String}})
    this.interaction_types = interaction_types
end

function in_distribution(this::Photoelectric)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

function out_distribution(this::Photoelectric,type::String)
    if type == "A"
        is_dirac = false
        N = 8
        quadrature = "gauss-legendre"
    elseif type == "P"
        is_dirac = true
        N = 1
        quadrature = "dirac"
    else
        error("Unknown type")
    end
    return is_dirac, N, quadrature
end

function bounds(this::Photoelectric,Ef⁻::Float64,Ef⁺::Float64,gi::Int64,type::String,Ui::Float64,E_in::Vector{Float64})
    if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    return Ef⁻,Ef⁺,isSkip
end

function dcs(this::Photoelectric,L::Int64,Ei::Float64,Z::Int64,iz::Int64,δi::Int64,Ef⁻::Float64,Ef⁺::Float64)

    # Absorption cross section
    Nshells,Zi,Ui,Ti,ri,subshells = electron_subshells(Z)
    if Ef⁺ < Ei-Ui[δi] ≤ Ef⁻
        if this.model == "jendl5"
            σa = this.photoelectric_cross_sections(iz,Ei,subshells[δi])
        elseif this.model == "biggs_lighthill"
            σa = this.photoelectric_cross_sections(iz,Ei)
        else
            error("Unknown photoelectric model.")
        end
    else
        σa = 0.0
    end

    # Angular distribution
    σℓ = zeros(L+1)
    γ = Ei+1
    β = sqrt(Ei*(Ei+2)/(Ei+1)^2)
    Γ = 1/(4/(3*(1-β^2)^2)+γ*(γ-1)*(γ-2)/(2*β^3)*(2*β/(1-β^2)-log((1+β)/(1-β))))
    α = [1,γ*(γ-1)*(γ-2)/2]
    ΔG3 = zeros(L+3,2)
    @inbounds for i in range(0,L+2), j in range(0,1)
        ΔG3[i+1,j+1] =  𝒢₃(i,j-4,1,-β,0,1,1)-𝒢₃(i,j-4,1,-β,0,1,-1)
    end
    @inbounds for ℓ in range(0,L)
        for k in range(0,div(ℓ,2))
            σℓk = 0.0
            for i in range(0,1), j in range(0,1)
                σℓk += α[i+1] * (-1)^j * ΔG3[ℓ-2*k+2*j+1,i+1]
            end
            σℓ[ℓ+1] += this.Cℓk[ℓ+1,k+1] * σℓk
        end
        σℓ[ℓ+1] *= Γ/(2^ℓ) * σa
    end

    # Correction to deal with high-order Legendre moments
    for ℓ in range(1,L)
        if abs(σℓ[1]) < abs(σℓ[ℓ+1])
            σℓ[ℓ+1:end] .= 0.0
            break
        end
    end

    return σℓ
end

function tcs(this::Photoelectric,Ei::Float64,Z::Int64,iz::Int64)

    if this.model == "jendl5"
        Nshells,Zi,Ui,Ti,ri,subshells = electron_subshells(Z)
        σt = 0.0
        for δi in range(1,Nshells)
            σt += this.photoelectric_cross_sections(iz,Ei,subshells[δi])
        end
    elseif this.model == "biggs_lighthill"
        σt = this.photoelectric_cross_sections(iz,Ei)
    else
        error("Unknown photoelectric model.")
    end

    return σt
end

function preload_data(this::Photoelectric,Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,L::Int64)
    this.preload_photoelectric_cross_sections(Z,ρ)
    # Precompute angular integration factors
    this.Cℓk = zeros(L+1,div(L,2)+1)
    for ℓ in range(0,L), k in range(0,div(L,2))
        this.Cℓk[ℓ+1,k+1] = (-1)^k * exp( sum(log.(1:2*ℓ-2*k)) - sum(log.(1:k)) - sum(log.(1:ℓ-k)) - sum(log.(1:ℓ-2*k)) )
    end
end

function preload_photoelectric_cross_sections(this::Photoelectric,Z::Vector{Int64},ρ::Float64)

    if this.model == "jendl5"

        data = load("./data/photoelectric_JENDL5.jld2")
        Nz = length(Z)
        E = Vector{Dict{String,Vector{Float64}}}(undef,Nz)
        σ = Vector{Dict{String,Vector{Float64}}}(undef,Nz)
        photoelectric_spline = Vector{Dict{String,Function}}(undef,Nz)
        for iz in range(1,Nz)
            E[iz] = data["E"][Z[iz]]
            σ[iz] = data["σ"][Z[iz]]
            photoelectric_spline[iz] = Dict()
            for subshells in keys(E[iz])
                photoelectric_spline[iz][subshells] = cubic_hermite_spline(E[iz][subshells],σ[iz][subshells])
            end
        end

        # Return the interpolation function
        this.photoelectric_cross_sections = function photoelectric_cross_sections_per_subshell(iz::Int64,Ei::Float64,subshells::String)
            if Ei > E[iz][subshells][1]
                return photoelectric_spline[iz][subshells](Ei)
            else
                return 0.0
            end
        end

    elseif this.model == "biggs_lighthill"

        data = load("./data/photoelectric_biggs_lighthill_1988.jld2")
        Nz = length(Z)
        E⁻ = Vector{Vector{Float64}}(undef,Nz)
        M = Vector{Array{Float64}}(undef,Nz)
        for iz in range(1,Nz)
            E⁻[iz] = data["E"][Z[iz]]
            M[iz] = data["M"][Z[iz]]
        end

        # Return the interpolation function
        this.photoelectric_cross_sections = function biggs_lighthill_cross_sections(iz::Int64,Ei::Float64)
            
            # Extract the 4 parameters from M matrix corresponding to the input energy
            A = Vector{Float64}(undef,4)
            N_interval = length(E⁻[iz])
            for i in range(1,N_interval)
                if i == N_interval && Ei >= E⁻[iz][i] || Ei >= E⁻[iz][i] && Ei < E⁻[iz][i+1] 
                    A = M[iz][i,:]
                    break
                end
            end

            # Absorption cross-sections
            σ = 0.0
            for i in range(1,4)
                σ += A[i]/Ei^i
            end
            return σ * ρ / nuclei_density(Z[iz],ρ)
        end
    else
        error("Unknown photoelectric model.")
    end
end