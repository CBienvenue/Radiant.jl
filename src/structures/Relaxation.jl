"""
    Relaxation

Structure used to define parameters for relaxation cascades (florescence and Auger electron
production)

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `ηmin::Float64 = 0.001` : minimum probability of the production of specific relaxation 
  particle production following electron cascades.
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}} = Dict((Photon,Electron) => ["P"],(Electron,Electron) => ["P"],(Positron,Electron) => ["P"],(Photon,Photon) => ["P"],(Electron,Photon) => ["P"],(Positron,Photon) => ["P"])` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `(Photon,Electron) => ["P"]` : production of Auger electron following incident photon ionization of subshells (by photoelectric effect).
    - `(Electron,Electron) => ["P"]` : production of Auger electron following incident electrons ionization of subshells (by Møller interaction).
    - `(Positron,Electron) => ["P"]` : production of Auger electron following incident positrons ionization of subshells (by Bhabha interaction).
    - `(Photon,Photon) => ["P"]` : production of fluorescence following incident photon ionization of subshells (by photoelectric effect).
    - `(Electron,Photon) => ["P"]` : production of fluorescence following incident electrons ionization of subshells (by Møller interaction).
    - `(Positron,Photon) => ["P"]` : production of fluorescence following incident positrons ionization of subshells (by Bhabha interaction).

"""
mutable struct Relaxation <: Interaction

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
    photoelectric_cross_sections::Function
    inelastic::Vector{Vector{Function}}
    η::Vector{Vector{Vector{Float64}}}
    ΔE::Vector{Vector{Vector{Float64}}}
    ηmin::Float64
    scattering_model::String
    inelastic_collision_model::Union{Inelastic_Collision,Missing}
    photoelectric_model::Union{Photoelectric,Missing}

    # Constructor(s)
    function Relaxation(interaction_types = Dict((Photon,Electron) => ["P"],(Electron,Electron) => ["P"],(Positron,Electron) => ["P"],(Photon,Photon) => ["P"],(Electron,Photon) => ["P"],(Positron,Photon) => ["P"]))
        this = new()
        this.name = "relaxation"
        this.interaction_types = interaction_types
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.is_CSD = false
        this.is_AFP = false
        this.is_AFP_decomposition = false
        this.is_elastic = false
        this.is_subshells_dependant = true
        this.set_minimum_probability(0.001)
        this.scattering_model = "BTE"
        this.inelastic_collision_model = missing
        this.photoelectric_model = missing
        return this
    end
end

# Method(s)
function Fluorescence() Relaxation(Dict((Photon,Photon) => ["P"],(Electron,Photon) => ["P"],(Positron,Photon) => ["P"])) end
function Auger() Relaxation(Dict((Photon,Electron) => ["P"],(Electron,Electron) => ["P"],(Positron,Electron) => ["P"])) end

"""
    set_interaction_types(this::Relaxation,interaction_types::Dict{Tuple{DataType,DataType},Vector{String}})

To define the interaction types for Relaxation processes.

# Input Argument(s)
- `this::Relaxation` : Relaxation structure.
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}} = Dict((Photon,Electron) => ["P"],(Electron,Electron) => ["P"],(Positron,Electron) => ["P"],(Photon,Photon) => ["P"],(Electron,Photon) => ["P"],(Positron,Photon) => ["P"])` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `(Photon,Electron) => ["P"]` : production of Auger electron following incident photon ionization of subshells (by photoelectric effect).
    - `(Electron,Electron) => ["P"]` : production of Auger electron following incident electrons ionization of subshells (by Møller interaction).
    - `(Positron,Electron) => ["P"]` : production of Auger electron following incident positrons ionization of subshells (by Bhabha interaction).
    - `(Photon,Photon) => ["P"]` : production of fluorescence following incident photon ionization of subshells (by photoelectric effect).
    - `(Electron,Photon) => ["P"]` : production of fluorescence following incident electrons ionization of subshells (by Møller interaction).
    - `(Positron,Photon) => ["P"]` : production of fluorescence following incident positrons ionization of subshells (by Bhabha interaction).

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> relaxation = Relaxation()
julia> relaxation.set_interaction_types( Dict((Electron,Electron) => ["P"]) ) # Only cascades following Møller interactions.
```
"""
function set_interaction_types(this::Relaxation,interaction_types)
    this.interaction_types = interaction_types
end

"""
    set_minimum_probability(this::Relaxation,ηmin::Real)

To define the minimum probability of a specific Relaxation electron production.

# Input Argument(s)
- `this::Relaxation` : Relaxation structure.
- `ηmin::Real` : minimum probability of a specific Relaxation electron production.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> relaxation = Relaxation()
julia> relaxation.set_minimum_probability(0.1) # Only Relaxation electron with probability greater than 10%.
```
"""
function set_minimum_probability(this::Relaxation,ηmin::Real)
    if ~(0 ≤ ηmin ≤ 1) error("Probability should be between 0 and 1.") end
    this.ηmin = ηmin
end

"""
    in_distribution(this::Relaxation)

Describe the energy discretization method for the incoming particle in the Relaxation
electron production interaction.

# Input Argument(s)
- `this::Relaxation` : Relaxation structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function in_distribution(this::Relaxation)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    out_distribution(this::Relaxation)

Describe the energy discretization method for the outgoing particle in the Relaxation
electron production interaction.

# Input Argument(s)
- `this::Relaxation` : Relaxation structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function out_distribution(this::Relaxation)
    is_dirac = true
    N = 1
    quadrature = "dirac"
    return is_dirac, N, quadrature
end

"""
    bounds(this::Relaxation,Ef⁻::Float64,Ef⁺::Float64)

Gives the integration energy bounds for the outgoing particle for Relaxation electron
production. 

# Input Argument(s)
- `this::Relaxation` : Relaxation structure.
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.

# Output Argument(s)
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `isSkip::Bool` : define if the integration is skipped or not.

"""
function bounds(this::Relaxation,Ef⁻::Float64,Ef⁺::Float64)
    if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    return Ef⁻,Ef⁺,isSkip
end

"""
    dcs(this::Relaxation,L::Int64,Ei::Float64,Ecutoff::Float64,Z::Int64,δi::Int64,
    Ef⁻::Float64,Ef⁺::Float64,incoming_particle::Particle,produced_particle::Particle,
    ηmin::Float64)

Gives the Legendre moments of the scattering cross-sections for relaxation production. 

# Input Argument(s)
- `this::Relaxation` : Relaxation structure.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : incoming particle energy.
- `Ecutoff::Float64` : cutoff energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `Z::Int64` : atomic number.
- `δi::Int64` : index of the subshell.
- `Ef⁻::Float64` : upper bound of the outgoing energy group.
- `Ef⁺::Float64` : lower bound of the outgoing energy group.
- `incoming_particle::Particle` : incoming particle.
- `produced_particle::Particle` : produced particle.

# Output Argument(s)
- `σℓ::Vector{Float64}` : Legendre moments of the scattering cross-sections.

"""
function dcs(this::Relaxation,L::Int64,Ei::Float64,Ecutoff::Float64,Ec::Float64,Z::Int64,δi::Int64,Ef⁻::Float64,Ef⁺::Float64,incoming_particle::Particle,produced_particle::Particle)

    # Initialization
    σℓ = zeros(L+1)
    is_relaxation = false
    if is_positron(incoming_particle) || is_electron(incoming_particle)
        if ~ismissing(this.inelastic_collision_model) && this.inelastic_collision_model.is_subshells_dependant
            is_relaxation = true
        end
    elseif is_photon(incoming_particle)
        if ~ismissing(this.photoelectric_model) && this.photoelectric_model.is_subshells_dependant
            is_relaxation = true
        end
    else
        error("Unknown interaction")
    end

    if is_relaxation

        # Compute absorption cross-section
        _,Zi,Ui,Ti,_,_ = electron_subshells(Z)
        σ_per_subshell = 0

        # Inelastic collisionnal electron scattering
        if is_electron(incoming_particle)
            σ_per_subshell += Zi[δi] * integrate_moller_per_subshell(Ei,0,Ui[δi],Ti[δi],Ei-Ec,Ei,this.inelastic_collision_model.is_focusing_møller)

        # Inelastic collisionnal positron scattering
        elseif is_positron(incoming_particle)
            σ_per_subshell += Zi[δi] * integrate_bhabha_per_subshell(Z,Ei,0,Ui[δi],Ei-Ec,Ei)

        # Photoelectric
        elseif is_photon(incoming_particle)
            σ_per_subshell += photoelectric_per_subshell(Z,Ei,δi)

        else
            error("Unknown incoming particle.")
        end

        # Relaxation cross-section
        σs = relaxation(Z,Ecutoff,Ef⁻,Ef⁺,δi,produced_particle,σ_per_subshell,this.ηmin)

        # Legendre moments of the cross-section (isotropic)
        σℓ[1] = σs

    end

    return σℓ
end