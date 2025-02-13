"""
    Photoelectric

Structure used to define parameters for production of multigroup photoelectric cross-sections.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}} = Dict((Photon,Photon) => ["A"],(Photon,Electron) => ["P"])` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `(Photon,Photon) => ["A"]` : absorption of incoming photon.
    - `(Photon,Electron) => ["P"]` : produced photo-electron.

"""
mutable struct Photoelectric <: Interaction

    # Variable(s)
    name::String
    incoming_particle::Vector{Type}
    interaction_particles::Vector{Type}
    interaction_types::Dict{Tuple{Type,Type},Vector{String}}
    is_CSD::Bool
    is_AFP::Bool
    is_elastic::Bool
    is_preload_data::Bool
    is_subshells_dependant::Bool
    photoelectric_cross_sections::Function
    model::String
    Cℓk::Array{Float64}
    scattering_model::String

    # Constructor(s)
    function Photoelectric()
        this = new()
        this.name = "photoelectric"
        this.interaction_types = Dict((Photon,Photon) => ["A"],(Photon,Electron) => ["P"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.is_CSD = false
        this.is_AFP = false
        this.is_elastic = false
        this.is_preload_data = true
        this.set_model("jendl5")
        this.scattering_model = "BTE"
        return this
    end

end

# Method(s)
"""
    set_interaction_types(this::Photoelectric,interaction_types::Dict{Tuple{DataType,DataType},Vector{String}})

To define the interaction types for photoelectric processes.

# Input Argument(s)
- `this::Photoelectric` : photoelectric structure.
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}}` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `(Photon,Photon) => ["A"]` : absorption of incoming photon.
    - `(Photon,Electron) => ["P"]` : produced photo-electron.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> photoelectric = Photoelectric()
julia> photoelectric.set_interaction_types( Dict((Photon,Photon) => ["A"],(Photon,Electron) => ["P"]) ) # Full photoelectric phenomenon (default case).
```
"""
function set_interaction_types(this::Photoelectric,interaction_types)
    this.interaction_types = interaction_types
end

"""
    set_model(this::Photoelectric,model::String)

To define the photoelectric model.

# Input Argument(s)
- `this::Photoelectric` : photoelectric structure.
- `model::String` : cross-section model:
    - `jendl5` : evaluated subshell-dependent cross-sections.
    - `biggs_lighthill` : Biggs and Lighthill cross-sections.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> photoelectric = Photoelectric()
julia> photoelectric.set_model("biggs_lighthill")
```
"""
function set_model(this::Photoelectric,model::String)
    if lowercase(model) ∉ ["jendl5","biggs_lighthill"] error("Unkown photoelectric model: $model.") end
    this.model = lowercase(model)
    if this.model == "jendl5"
        this.is_subshells_dependant = true
    elseif this.model == "biggs_lighthill"
        this.is_subshells_dependant = false
    else    
        error("Unknown photoelectric model.")
    end
end

"""
    in_distribution(this::Photoelectric)

Describe the energy discretization method for the incoming particle in the photoelectric
interaction.

# Input Argument(s)
- `this::Photoelectric` : photoelectric structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function in_distribution(this::Photoelectric)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    out_distribution(this::Photoelectric,type::String)

Describe the energy discretization method for the outgoing particle in the photoelectric
interaction.

# Input Argument(s)
- `this::Photoelectric` : photoelectric structure.
- `type::String` : type of interaction.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
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

"""
    bounds(this::Photoelectric,Ef⁻::Float64,Ef⁺::Float64)

Gives the integration energy bounds for the outgoing particle for photoelectric
interaction. 

# Input Argument(s)
- `this::Photoelectric` : photoelectric structure.
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.

# Output Argument(s)
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `isSkip::Bool` : define if the integration is skipped or not.

"""
function bounds(this::Photoelectric,Ef⁻::Float64,Ef⁺::Float64)
    if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    return Ef⁻,Ef⁺,isSkip
end

"""
    dcs(this::Photoelectric,L::Int64,Ei::Float64,Z::Int64,iz::Int64,δi::Int64,Ef⁻::Float64,
    Ef⁺::Float64)

Gives the Legendre moments of the scattering cross-sections for photoelectric interaction. 

# Input Argument(s)
- `this::Photoelectric` : photoelectric structure.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.
- `iz::Int64` : index of the element in the material.
- `δi::Int64` : subshell index.
- `Ef⁻::Float64` : upper bounds associated with the outgoing particle.
- `Ef⁺::Float64` : lower bounds associated with the outgoing particle.

# Output Argument(s)
- `σℓ::Vector{Float64}` : Legendre moments of the scattering cross-sections.

"""
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

"""
    tcs(this::Photoelectric,Ei::Float64,Z::Int64,iz::Int64)

Gives the total cross-section for photoelectric interaction. 

# Input Argument(s)
- `this::Photoelectric` : photoelectric structure. 
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.
- `iz::Int64` : index of the element in the material.

# Output Argument(s)
- `σt::Float64` : total cross-section.

"""
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

"""
   preload_data(this::Photoelectric,Z::Vector{Int64},ρ::Float64,L::Int64)

Preload data for multigroup photoelectric calculations.

# Input Argument(s)
- `this::Photoelectric` : photoelectric structure. 
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `ρ::Float64` : material density.
- `L::Int64` : Legendre truncation order.

# Output Argument(s)
N/A

"""
function preload_data(this::Photoelectric,Z::Vector{Int64},ρ::Float64,L::Int64)
    this.preload_photoelectric_cross_sections(Z,ρ)
    # Precompute angular integration factors
    this.Cℓk = zeros(L+1,div(L,2)+1)
    for ℓ in range(0,L), k in range(0,div(L,2))
        this.Cℓk[ℓ+1,k+1] = (-1)^k * exp( sum(log.(1:2*ℓ-2*k)) - sum(log.(1:k)) - sum(log.(1:ℓ-k)) - sum(log.(1:ℓ-2*k)) )
    end
end

"""
    preload_photoelectric_cross_sections(this::Photoelectric,Z::Vector{Int64},ρ::Float64)

Preload data for photoelectric cross-sections per subshells.

# Input Argument(s)
- `this::Photoelectric` : photoelectric structure. 
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `ρ::Float64` : material density.

# Output Argument(s)
N/A

"""
function preload_photoelectric_cross_sections(this::Photoelectric,Z::Vector{Int64},ρ::Float64)

    if this.model == "jendl5"

        path = joinpath(find_package_root(), "data", "photoelectric_JENDL5.jld2")
        data = load(path)
        Nz = length(Z)
        E = Vector{Dict{String,Vector{Float64}}}(undef,Nz)
        σ = Vector{Dict{String,Vector{Float64}}}(undef,Nz)
        photoelectric_spline = Vector{Dict{String,Function}}(undef,Nz)
        for iz in range(1,Nz)
            E[iz] = data["E"][Z[iz]]
            σ[iz] = data["σ"][Z[iz]]

            # Temporary fix - Delete additionnal data in photoelectric_JENDL5...
            for subshells in keys(E[iz])
                index = length(E[iz][subshells])
                for i in range(2,length(E[iz][subshells]))
                    if E[iz][subshells][i-1] > E[iz][subshells][i] index = i-1; break end
                end
                E[iz][subshells] = E[iz][subshells][1:index]
                σ[iz][subshells] = σ[iz][subshells][1:index] 
            end

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

        path = joinpath(find_package_root(), "data", "photoelectric_biggs_lighthill_1988.jld2")
        data = load(path)
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