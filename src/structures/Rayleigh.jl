"""
    Rayleigh

Structure used to define parameters for production of multigroup Rayleigh cross-sections.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}} = Dict((Photon,Photon) => ["S"])` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `(Photon,Photon) => ["S"]` : elastic scattering of photons

"""
mutable struct Rayleigh <: Interaction

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
    is_ETC::Bool
    scattering_model::String

    # Constructor(s)
    function Rayleigh()
        this = new()
        this.name = "rayleigh"
        this.interaction_types = Dict((Photon,Photon) => ["S"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.is_CSD = false
        this.is_AFP = false
        this.is_AFP_decomposition = false
        this.is_elastic = true
        this.is_subshells_dependant = false
        this.is_ETC = true
        this.scattering_model = "BTE"
        return this
    end

end

# Method(s)
"""
    in_distribution(this::Rayleigh)

Describe the energy discretization method for the incoming particle in the Rayleigh
interaction.

# Input Argument(s)
- `this::Rayleigh` : Rayleigh structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function in_distribution(this::Rayleigh)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    out_distribution(this::Rayleigh)

Describe the energy discretization method for the outgoing particle in the Rayleigh
interaction.

# Input Argument(s)
- `this::Rayleigh` : Rayleigh structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function out_distribution(this::Rayleigh)
    is_dirac = true
    N = 1
    quadrature = "dirac"
    return is_dirac, N, quadrature
end

"""
    bounds(this::Rayleigh,Ef⁻::Float64,Ef⁺::Float64,gi::Int64,gf::Int64)

Gives the integration energy bounds for the outgoing particle for Rayleigh. 

# Input Argument(s)
- `this::Rayleigh` : Rayleigh structure.
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `gi::Int64` : group of the incoming particle.
- `gf::Int64` : group of the outgoing particle.

# Output Argument(s)
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `isSkip::Bool` : define if the integration is skipped or not.

"""
function bounds(this::Rayleigh,Ef⁻::Float64,Ef⁺::Float64,gi::Int64,gf::Int64)
    if (gf != gi) isSkip = true else isSkip = false end
    return Ef⁻,Ef⁺,isSkip
end

"""
    dcs(this::Rayleigh,L::Int64,Ei::Float64,Ef::Float64,Z::Int64,particle::Particle,
    type::String,iz::Int64)

Gives the Legendre moments of the scattering cross-sections for Rayleigh. 

# Input Argument(s)
- `this::Rayleigh` : Rayleigh structure.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.
- `iz::Int64` : index of the element in the material.

# Output Argument(s)
- `σl::Vector{Float64}` : Legendre moments of the scattering cross-sections.

"""
function dcs(this::Rayleigh,L::Int64,Ei::Float64,Z::Int64,iz::Int64)

    # Initialization
    rₑ = 2.81794092E-13 # (in cm)
    σl = zeros(L+1)

    # Angular quadrature
    μ,w = quadrature(32,"gauss-legendre")
    for n in range(1,32)
        σs = rayleigh(Z,Ei,μ[n])
        Plμ = legendre_polynomials_up_to_L(L,μ[n])
        for l in range(0,L) σl[l+1] += w[n] * σs * Plμ[l+1] end
    end
    return σl
end

"""
    tcs(this::Rayleigh,Ei::Float64,Z::Int64,iz::Int64)

Gives the total cross-section for Rayleigh. 

# Input Argument(s)
- `this::Rayleigh` : Rayleigh structure. 
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.
- `iz::Int64` : index of the element in the material.

# Output Argument(s)
- `σt::Float64` : total cross-section.

"""
function tcs(this::Rayleigh,Ei::Float64,Z::Int64,iz::Int64)
    σt = dcs(this,0,Ei,Z,iz)[1]
    return σt
end

"""
    acs(this::Rayleigh)

Gives the absorption cross-section for Rayleigh. 

# Input Argument(s)
- `this::Rayleigh` : Rayleigh structure.

# Output Argument(s)
- `σa::Float64` : absorption cross-section.

"""
function acs(this::Rayleigh)
    return 0
end