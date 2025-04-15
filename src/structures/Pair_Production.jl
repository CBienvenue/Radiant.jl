"""
    Pair_Production

Structure used to define parameters for production of multigroup pair production cross-sections.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}} = Dict((Photon,Photon) => ["A"],(Photon,Electron) => ["P"],(Photon,Positron) => ["P"])` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `(Photon,Photon) => ["A"]` : absorption of incoming photon.
    - `(Photon,Electron) => ["P"]` : produced electron.
    - `(Photon,Positron) => ["P"]` : produced positron.
- `angular_scattering_type::String=modified_dipole` : type of angular scattering, which can takes the following values:
    - `angular_scattering_type = modified_dipole` : modified dipôle distribution, based on Poskus (2019) shape functions.
    - `angular_scattering_type = sommerfield` : Sommerfield distribution.

"""
mutable struct Pair_Production <: Interaction

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
    angular_scattering_type::String
    is_triplet_contribution::Bool
    scattering_model::String
    is_positron::Bool

    # Constructor(s)
    function Pair_Production()
        this = new()
        this.name = "pair_production"
        this.interaction_types = Dict((Photon,Photon) => ["A"],(Photon,Electron) => ["P"],(Photon,Positron) => ["P"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.is_CSD = false
        this.is_AFP = false
        this.is_AFP_decomposition = false
        this.is_elastic = false
        this.is_subshells_dependant = false
        this.set_angular_scattering_type("sommerfield")
        this.scattering_model = "BTE"
        this.is_positron = true
        return this
    end
end

# Method(s)
"""
    set_interaction_types(this::Pair_Production,interaction_types::Dict{Tuple{DataType,DataType},Vector{String}})

To define the interaction types for pair production processes.

# Input Argument(s)
- `this::Pair_Production` : pair production structure.
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}}` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `(Photon,Photon) => ["A"]` : absorption of incoming photon.
    - `(Photon,Electron) => ["P"]` : produced electron.
    - `(Photon,Positron) => ["P"]` : produced positron.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> pair_production = Pair_Production()
julia> pair_production.set_interaction_types( Dict((Electron,Electron) => ["S"]) ) # Only electron scattering, with photon absorption.
```
"""
function set_interaction_types(this::Pair_Production,interaction_types)
    this.interaction_types = interaction_types
end

"""
    set_angular_scattering_type(this::Pair_Production,angular_scattering_type::String)

To define the pair_production photons angular distribution.

# Input Argument(s)
- `this::Pair_Production` : pair_production structure.
- `angular_scattering_type::String` : angular scattering type.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> pair_production = Pair_Production()
julia> pair_production.set_angular_scattering_type("sommerfield")
```
"""
function set_angular_scattering_type(this::Pair_Production,angular_scattering_type::String)
    if lowercase(angular_scattering_type) ∉ ["modified_dipole","sommerfield"] error("Undefined angular distribution.") end
    this.angular_scattering_type = lowercase(angular_scattering_type)
end

"""
    in_distribution(this::Pair_Production)

Describe the energy discretization method for the incoming particle in the pair production
interaction.

# Input Argument(s)
- `this::Pair_Production` : pair production structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function in_distribution(this::Pair_Production)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    out_distribution(this::Pair_Production)

Describe the energy discretization method for the outgoing particle in the pair production
interaction.

# Input Argument(s)
- `this::Pair_Production` : pair production structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function out_distribution(this::Pair_Production)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    bounds(this::Pair_Production,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String)

Gives the integration energy bounds for the outgoing particle for pair production. 

# Input Argument(s)
- `this::Pair_Production` : pair production structure.
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `Ei::Float64` : energy of the incoming particle.
- `type::String` : type of interaction.

# Output Argument(s)
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `isSkip::Bool` : define if the integration is skipped or not.

"""
function bounds(this::Pair_Production,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String)
    # Electron/positron production
    if type == "P" || type == "A"
        Ef⁻ = min(Ef⁻,Ei-2)
        if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    else
        error("Unknown type of method for pair production.")
    end
    return Ef⁻,Ef⁺,isSkip
end

"""
    dcs(this::Pair_Production,L::Int64,Ei::Float64,Ef::Float64,Z::Int64,particle::Particle,
    type::String,iz::Int64)

Gives the Legendre moments of the scattering cross-sections for pair production. 

# Input Argument(s)
- `this::Pair_Production` : pair production structure.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : incoming particle energy.
- `Ef::Float64` : outgoing particle energy.
- `Z::Int64` : atomic number.
- `type::String` : type of interaction.
- `iz::Int64` : index of the element in the material.
- `particles::Vector{Particle}` : list of particles.

# Output Argument(s)
- `σℓ::Vector{Float64}` : Legendre moments of the scattering cross-sections.

"""
function dcs(this::Pair_Production,L::Int64,Ei::Float64,Ef::Float64,Z::Int64,type::String,iz::Int64,particles::Vector{Particle})

    # Validation
    if type != "P" error("Unknown interaction.") end
    σℓ = zeros(L+1)
    if Ei-2 < 0 return σℓ end

    # Electron/Positron production
    σs = baro(Z,Ei,Ef)

    # If no explicit positron transport, generate two electrons
    if ~this.is_positron σs *= 2 end

    # Compute angular distribution
    if this.angular_scattering_type == "sommerfield"
        Wℓ = sommerfield(Ei,L)
    elseif this.angular_scattering_type == "modified_dipole"
        Wℓ = poskus(Z,Ei,Ef,L)
    else
        error("Unkown angular distribution.")
    end

    # Compute Legendre moments of the scattering cross-section
    for ℓ in range(0,L) σℓ[ℓ+1] += σs * Wℓ[ℓ+1] end
    return σℓ
end

"""
    tcs(this::Pair_Production,Ei::Float64,Z::Int64,Ecutoff::Float64)

Gives the absorption cross-section for pair production. 

# Input Argument(s)
- `this::Pair_Production` : pair production structure. 
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.
- `Eout::Vector{Float64}` : outgoing energy boundaries.

# Output Argument(s)
- `σt::Float64` : absorption cross-section.

"""
function tcs(this::Pair_Production,Ei::Float64,Z::Int64,Eout::Vector{Float64})

    # Initialization
    σt = 0.0
    if Ei-2 < 0 return σt end

    # Compute total cross-sections
    Ngf = length(Eout)-1
    is_dirac, Np, q_type = out_distribution(this)
    if is_dirac Np = 1; u = [0]; w = [2] else u,w = quadrature(Np,q_type) end
    for gf in range(1,Ngf+1)
        Ef⁻ = Eout[gf]
        if (gf != Ngf+1) Ef⁺ = Eout[gf+1] else Ef⁺ = 0.0 end
        Ef⁻,Ef⁺,isSkip = bounds(this,Ef⁻,Ef⁺,Ei,"A")
        if isSkip continue end
        ΔEf = Ef⁻ - Ef⁺
        for n in range(1,Np)
            Ef = (u[n]*ΔEf + (Ef⁻+Ef⁺))/2
            σs = baro(Z,Ei,Ef)
            σt += ΔEf/2 * w[n] * σs
        end
     end
    return σt
end

"""
    acs(this::Pair_Production,Ei::Float64,Z::Int64,Ecutoff::Float64)

Gives the absorption cross-section for pair production. 

# Input Argument(s)
- `this::Pair_Production` : pair production structure. 
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.
- `Eout::Vector{Float64}` : outgoing energy boundaries.

# Output Argument(s)
- `σa::Float64` : absorption cross-section.

"""
function acs(this::Pair_Production,Ei::Float64,Z::Int64,Eout::Vector{Float64})
    return tcs(this,Ei,Z,Eout)
end