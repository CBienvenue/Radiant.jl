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
    is_AFP_decomposition::Bool
    is_elastic::Bool
    is_subshells_dependant::Bool
    model::String
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
        this.is_AFP_decomposition = false
        this.is_elastic = false
        this.set_model("epdl97")
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
    - `epdl97` : evaluated subshell-dependent cross-sections.
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
    if lowercase(model) ∉ ["epdl97","biggs_lighthill"] error("Unkown photoelectric model: $model.") end
    this.model = lowercase(model)
    if this.model == "epdl97"
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
- `δi::Int64` : subshell index.
- `Ef⁻::Float64` : upper bounds associated with the outgoing particle.
- `Ef⁺::Float64` : lower bounds associated with the outgoing particle.

# Output Argument(s)
- `σℓ::Vector{Float64}` : Legendre moments of the scattering cross-sections.

"""
function dcs(this::Photoelectric,L::Int64,Ei::Float64,Z::Int64,δi::Int64,Ef⁻::Float64,Ef⁺::Float64)

    # Absorption cross section
    _,_,Ui,_,_,_ = electron_subshells(Z)
    if Ef⁺ < Ei-Ui[δi] ≤ Ef⁻
        if this.model == "epdl97"
            σa = photoelectric_per_subshell(Z,Ei,δi)
        elseif this.model == "biggs_lighthill"
            σa = biggs_lighthill(Z,Ei)
        else
            error("Unknown photoelectric model.")
        end
    else
        σa = 0.0
    end

    # Angular distribution
    Wℓ = sauter(Ei,L)

    # Legendre moments of the scattering cross-section
    σℓ = zeros(L+1)
    for ℓ in range(0,L) σℓ[ℓ+1] = σa * Wℓ[ℓ+1] end
    return σℓ
end

"""
    tcs(this::Photoelectric,Ei::Float64,Z::Int64)

Gives the total cross-section for photoelectric interaction. 

# Input Argument(s)
- `this::Photoelectric` : photoelectric structure. 
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.

# Output Argument(s)
- `σt::Float64` : total cross-section.

"""
function tcs(this::Photoelectric,Ei::Float64,Z::Int64)
    if this.model == "epdl97"
        Nshells,_,_,_,_,_ = electron_subshells(Z)
        σt = 0.0
        for δi in range(1,Nshells)
            σt += photoelectric_per_subshell(Z,Ei,δi)
        end
    elseif this.model == "biggs_lighthill"
        σt = biggs_lighthill(Z,Ei)
    else
        error("Unknown photoelectric model.")
    end
    return σt
end

"""
    acs(this::Photoelectric,Ei::Float64,Z::Int64)

Gives the absorption cross-section for photoelectric interaction. 

# Input Argument(s)
- `this::Photoelectric` : photoelectric structure. 
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.

# Output Argument(s)
- `σa::Float64` : total cross-section.

"""
function acs(this::Photoelectric,Ei::Float64,Z::Int64)
    return tcs(this,Ei,Z)
end