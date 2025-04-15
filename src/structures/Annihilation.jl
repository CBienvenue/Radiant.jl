"""
    Annihilation

Structure used to define parameters for production of multigroup annihilation cross-sections.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}} = Dict((Positron,Positron) => ["A"],(Positron,Photon) => ["P","P_inel","P_brems"],(Photon,Photon) => ["P_pp"])` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `(Positron,Positron) => ["A"]` : absorption of the incoming positron.
    - `(Positron,Photon) => ["P"]` : production of the two photons following annihilation.
    - `(Positron,Photon) => ["P_inel"]` : production of annihilation photons from inelastic collisional positrons absorption following scattering under the cutoff energy.
    - `(Positron,Photon) => ["P_brems"]` : production of annihilation photons from Bremsstrahlung positrons absorption following scattering under the cutoff energy.
    - `(Photon,Photon) => ["P_pp"]` : production of annihilation photons from absorption of positrons following their production under the cutoff energy.

"""
mutable struct Annihilation <: Interaction

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
    scattering_model::String
    inelastic_collision_model::Inelastic_Collision
    bremsstrahlung_model::Bremsstrahlung
    pair_production_model::Pair_Production

    # Constructor(s)
    function Annihilation()
        this = new()
        this.name = "annihilation"
        this.interaction_types = Dict((Positron,Positron) => ["A"],(Positron,Photon) => ["P","P_inel","P_brems"],(Photon,Photon) => ["P_pp"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.is_CSD = false
        this.is_AFP = false
        this.is_AFP_decomposition = false
        this.is_elastic = false
        this.is_subshells_dependant = false
        this.scattering_model = "BTE"
        this.inelastic_collision_model = Inelastic_Collision()
        this.bremsstrahlung_model = Bremsstrahlung()
        this.pair_production_model = Pair_Production()
        return this
    end
end

# Method(s)
"""
    set_interaction_types(this::Annihilation,interaction_types::Dict{Tuple{DataType,DataType},Vector{String}})

To define the interaction types for annihilation processes.

# Input Argument(s)
- `this::Annihilation` : annihilation structure.
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}}` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `(Positron,Positron) => ["A"]` : absorption of the incoming positron.
    - `(Positron,Photon) => ["P"]` : production of the two photons following annihilation.
    - `(Positron,Photon) => ["P_inel"]` : production of annihilation photons from inelastic collisional positrons absorption following scattering under the cutoff energy.
    - `(Positron,Photon) => ["P_brems"]` : production of annihilation photons from Bremsstrahlung positrons absorption following scattering under the cutoff energy.
    - `(Photon,Photon) => ["P_pp"]` : production of annihilation photons from absorption of positrons following their production under the cutoff energy.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> annihilation = Annihilation()
julia> annihilation.set_interaction_types( Dict((Positron,Positron) => ["A"]) ) # Annihilation is set to be only absorption of positrons without any production of photons
```
"""
function set_interaction_types(this::Annihilation,interaction_types)
    this.interaction_types = interaction_types
end

"""
    in_distribution(this::Annihilation)

Describe the energy discretization method for the incoming particle in the annihilation
interaction.

# Input Argument(s)
- `this::Annihilation` : annihilation structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function in_distribution(this::Annihilation)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    out_distribution(this::Annihilation,type::String)

Describe the energy discretization method for the outgoing particle in the annihilation
interaction.

# Input Argument(s)
- `this::Annihilation` : annihilation structure.
- `type::String` : type of interaction.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function out_distribution(this::Annihilation,type::String)
    if type ∈ ["A","P"]
        is_dirac = false
        N = 8
        quadrature = "gauss-legendre"
    else
        is_dirac = true
        N = 1
        quadrature = "dirac"
    end
    return is_dirac, N, quadrature
end

"""
    bounds(this::Annihilation,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String)

Gives the integration energy bounds for the outgoing particle for annihilation. 

# Input Argument(s)
- `this::Annihilation` : annihilation structure.
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `Ei::Float64` : energy of the incoming particle.
- `type::String` : type of interaction.

# Output Argument(s)
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `isSkip::Bool` : define if the integration is skipped or not.

"""
function bounds(this::Annihilation,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String)
    γ = Ei+1
    ζmin = 1/(γ+1+sqrt(γ^2-1))
    if type == "P"
        Ef⁻ = min(Ef⁻,(γ+1)*(1-ζmin))
        Ef⁺ = max(Ef⁺,(γ+1)*ζmin)
        if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    elseif type ∈ ["P_inel","P_brems","P_pp"]
        if (Ef⁻-Ef⁺ < 0 || ~(Ef⁺ < 1 < Ef⁻) ) isSkip = true else isSkip = false end
    else
        error("Unknown type of method for annihilation.")
    end
    return Ef⁻,Ef⁺,isSkip
end

"""
    dcs(this::Annihilation,L::Int64,Ei::Float64,Ef::Float64,type::String,Z::Vector{Int64},
    iz::Int64,Ein::Vector{Float64},Ec::Float64)

Gives the Legendre moments of the scattering cross-sections for annihilation. 

# Input Argument(s)
- `this::Annihilation` : annihilation structure.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : incoming particle energy.
- `Ef::Float64` : outgoing particle energy.
- `type::String` : type of interaction.
- `Z::Int64` : atomic numbers of the element.
- `Ecutoff::Float64` : cutoff energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.

# Output Argument(s)
- `σℓ::Vector{Float64}` : Legendre moments of the scattering cross-sections.

"""
function dcs(this::Annihilation,L::Int64,Ei::Float64,Ef::Float64,type::String,Z::Int64,Ein::Vector{Float64},Ec::Float64)

    #----
    # Initialization
    #----
    σℓ = zeros(L+1)

    #----
    # Two annihilation photons
    #----
    if type == "P"

        # Scattering cross-section
        σs = heitler(Z,Ei,Ef)

        # Angular distribution
        Wℓ = angular_heitler(Ei,Ef,L)

        # Compute the Legendre moments of the cross-section
        for ℓ in range(0,L) σℓ[ℓ+1] = σs * Wℓ[ℓ+1] end

    #----
    # Annihilation of positrons scattered under the cutoff from inelastic collisionnal interaction
    #----
    elseif type == "P_inel"
        σa = this.inelastic_collision_model.acs(Ei,Ec,Positron(),Z,Ein[end])
        σℓ[1] += 2 * σa

    #----
    # Annihilation of positrons scattered under the cutoff from Bremsstrahlung interaction
    #----
    elseif type == "P_brems"
        σa = this.bremsstrahlung_model.acs(Ei,Z,Ec,Positron(),Ein[end])
        σℓ[1] += 2 * σa

    #----
    # Annihilation of positrons produced under the cutoff following pair production event
    #----
    elseif type == "P_pp"
        σa = this.pair_production_model.acs(Ei,Z,[Ein[end]])
        σℓ[1] += 2 * σa

    else
        error("Unknown interaction.")
    end
    return σℓ
end

"""
    tcs(this::Annihilation,Ei::Float64,Z::Int64)

Gives the total cross-section for annihilation. 

# Input Argument(s)
- `this::Annihilation` : annihilation structure.
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.

# Output Argument(s)
- `σt::Float64` : total cross-section.

"""
function tcs(this::Annihilation,Ei::Float64,Z::Int64)
    return integrated_heitler(Z,Ei)
end

"""
    acs(this::Annihilation,Ei::Float64,Z::Int64)

Gives the absorption cross-section for annihilation. 

# Input Argument(s)
- `this::Annihilation` : annihilation structure.
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.

# Output Argument(s)
- `σa::Float64` : absorption cross-section.

"""
function acs(this::Annihilation,Ei::Float64,Z::Int64)
    return tcs(this,Ei,Z)
end