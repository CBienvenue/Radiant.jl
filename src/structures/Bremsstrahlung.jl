"""
    Bremsstrahlung

Structure used to define parameters for production of multigroup bremsstrahlung cross-sections.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}} = Dict((Electron,Electron) => ["S"],(Electron,Photon) => ["P"],(Positron,Positron) => ["S"],(Positron,Photon) => ["P"])` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `(Electron,Electron) => ["S"]` : scattering of incident electron following Bremsstrahlung interaction.
    - `(Electron,Photon) => ["P"]` : produced photon following Bremsstrahlung interaction by incident electron.
    - `(Positron,Positron) => ["S"]` : scattering of incident positron following Bremsstrahlung interaction.
    - `(Positron,Photon) => ["P"]` : produced photon following Bremsstrahlung interaction by incident positron.
- `angular_scattering_type::String=modified_dipole` : type of angular scattering, which can takes the following values:
    - `angular_scattering_type = modified_dipole` : modified dipôle distribution, based on Poskus (2019) shape functions.
    - `angular_scattering_type = sommerfield` : Sommerfield distribution.

"""
mutable struct Bremsstrahlung <: Interaction

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
    scattering_model::String

    # Constructor(s)
    function Bremsstrahlung()
        this = new()
        this.name = "bremsstrahlung"
        this.interaction_types = Dict((Electron,Electron) => ["S"],(Electron,Photon) => ["P"],(Positron,Positron) => ["S"],(Positron,Photon) => ["P"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.is_CSD = true
        this.is_AFP = true
        this.is_AFP_decomposition = false
        this.is_elastic = false
        this.is_subshells_dependant = false
        this.set_angular_scattering_type("sommerfield")
        this.set_scattering_model("BFP")
        return this
    end

end

# Method(s)
"""
    set_interaction_types(this::Bremsstrahlung,interaction_types::Dict{Tuple{DataType,DataType},Vector{String}})

To define the interaction types for bremsstrahlung processes.

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure.
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}}` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `(Electron,Electron) => ["S"]` : scattering of incident electron following Bremsstrahlung interaction.
    - `(Electron,Photon) => ["P"]` : produced photon following Bremsstrahlung interaction by incident electron.
    - `(Positron,Positron) => ["S"]` : scattering of incident positron following Bremsstrahlung interaction.
    - `(Positron,Photon) => ["P"]` : produced photon following Bremsstrahlung interaction by incident positron.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> bremsstrahlung = Bremsstrahlung()
julia> bremsstrahlung.set_interaction_types( Dict((Electron,Electron) => ["S"]) ) # Only electron scattering, with photon absorption.
```
"""
function set_interaction_types(this::Bremsstrahlung,interaction_types)
    this.interaction_types = interaction_types
end

"""
    set_angular_scattering_type(this::Bremsstrahlung,angular_scattering_type::String)

To define the bremsstrahlung photons angular distribution.

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure.
- `angular_scattering_type::String` : angular scattering type.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> bremsstrahlung = Bremsstrahlung()
julia> bremsstrahlung.set_angular_scattering_type("sommerfield")
```
"""
function set_angular_scattering_type(this::Bremsstrahlung,angular_scattering_type::String)
    if lowercase(angular_scattering_type) ∉ ["modified_dipole","sommerfield"] error("Undefined angular distribution.") end
    this.angular_scattering_type = lowercase(angular_scattering_type)
end

"""
    scattering_model(this::Bremsstrahlung,scattering_model::String)

To define the solver for bremsstrahlung scattering.

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure.
- `solver::String` : solver for bremsstrahlung scattering, which can be:
    - `BFP` : Boltzmann Fokker-Planck solver.
    - `FP` : Fokker-Planck solver.
    - `BTE` : Boltzmann solver.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> bremsstrahlung = Bremsstrahlung()
julia> bremsstrahlung.set_scattering_model("FP")
```
"""
function set_scattering_model(this::Bremsstrahlung,scattering_model::String)
    if uppercase(scattering_model) ∉ ["BFP","FP","BTE"] error("Unknown scattering model (should be BFP, FP or BTE).") end
    this.scattering_model = uppercase(scattering_model)
end

"""
    in_distribution(this::Bremsstrahlung)

Describe the energy discretization method for the incoming particle in the bremsstrahlung
interaction.

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function in_distribution(this::Bremsstrahlung)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    out_distribution(this::Bremsstrahlung)

Describe the energy discretization method for the outgoing particle in the bremsstrahlung
interaction.

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function out_distribution(this::Bremsstrahlung)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    bounds(this::Bremsstrahlung,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String,
    Ec::Float64)

Gives the integration energy bounds for the outgoing particle for bremsstrahlung. 

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure.
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `Ei::Float64` : energy of the incoming particle.
- `type::String` : type of interaction.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.

# Output Argument(s)
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `isSkip::Bool` : define if the integration is skipped or not.

"""
function bounds(this::Bremsstrahlung,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String,Ec::Float64)
    # Scattered electron or positron
    if type == "S" 
        Ef⁻ = min(Ef⁻,Ec)
        if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    # Produced photon
    elseif type == "P" 
        Ef⁻ = min(Ef⁻,Ei)
        if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    else
        error("Unknown type of method for Bremsstrahlung scattering.")
    end
    return Ef⁻,Ef⁺,isSkip
end

"""
    dcs(this::Bremsstrahlung,L::Int64,Ei::Float64,Ef::Float64,Z::Int64,particle::Particle,
    type::String,iz::Int64)

Gives the Legendre moments of the scattering cross-sections for bremsstrahlung. 

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : incoming particle energy.
- `Ef::Float64` : outgoing particle energy.
- `Z::Int64` : atomic number.
- `particle::Particle` : incoming particle.
- `type::String` : type of interaction.
- `iz::Int64` : index of the element in the material.

# Output Argument(s)
- `σl::Vector{Float64}` : Legendre moments of the scattering cross-sections.

"""
function dcs(this::Bremsstrahlung,L::Int64,Ei::Float64,Ef::Float64,Z::Int64,particle::Particle,type::String,iz::Int64)

    # Inititalisation
    σl = zeros(L+1)

    # Scattered electron
    if type == "S" 

        # Compute the differential scattering cross section
        σs = seltzer_berger_cross_section(Z,Ei,Ei-Ef,particle)

        # Compute the angular distribution, μ = 1.0 (Forward peaked)
        Wl = ones(L+1)

        # Compute the Legendre moments of the cross-section
        for l in range(0,L) σl[l+1] += σs * Wl[l+1] end

    # Produced photons
    elseif type == "P"

        # Compute the differential scattering cross section
        σs = seltzer_berger_cross_section(Z,Ei,Ef,particle)

        # Compute the angular distribution
        if this.angular_scattering_type == "sommerfield"
            Wl = sommerfield(Ei,L)
        elseif this.angular_scattering_type == "modified_dipole"
            Wl = poskus(Z,Ei,Ef,L)
        else
            error("Unkown angular distribution.")
        end

        # Compute the Legendre moments of the cross-section
        for l in range(0,L) σl[l+1] += σs * Wl[l+1] end
        
    else
        error("Unknown interaction.")
    end
    return σl
end

"""
    tcs(this::Bremsstrahlung,Ei::Float64,Z::Int64,Ec::Float64,iz::Int64,particle::Particle,
    type::String,Eout::Vector{Float64})

Gives the total cross-section for bremsstrahlung. 

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure. 
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `iz::Int64` : index of the element in the material.
- `particle::Particle` : incoming particle.
- `type::String` : type of interaction.
- `Eout::Vector{Float64}` : outgoing energy boundaries.

# Output Argument(s)
- `σt::Float64` : total cross-section.

"""
function tcs(this::Bremsstrahlung,Ei::Float64,Z::Int64,Ec::Float64,particle::Particle,Eout::Vector{Float64})

    # Compute total cross section
    σt = 0.0
    Ngf = length(Eout)-1
    is_dirac, Np, q_type = out_distribution(this)
    if is_dirac Np = 1; u = [0]; w = [2] else u,w = quadrature(Np,q_type) end
    for gf in range(1,Ngf+1)
        Ef⁻ = Eout[gf]
        if (gf != Ngf+1) Ef⁺ = Eout[gf+1] else Ef⁺ = 0.0 end
        Ef⁻,Ef⁺,isSkip = bounds(this,Ef⁻,Ef⁺,Ei,"S",Ec)
        if isSkip continue end
        ΔEf = Ef⁻ - Ef⁺
        for n in range(1,Np)
            Ef = (u[n]*ΔEf + (Ef⁻+Ef⁺))/2
            Eγ = Ei-Ef
            if Ei ≥ Ef && ΔEf ≥ 0
                σt += ΔEf/2 * w[n] * seltzer_berger_cross_section(Z,Ei,Eγ,particle)
            end
        end
    end
    return σt
end

"""
    acs(this::Bremsstrahlung,Ei::Float64,Z::Int64,Ec::Float64,particle::Particle,
    type::String,Eout::Vector{Float64})

Gives the absorption cross-section for bremsstrahlung. 

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure. 
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.
- `type::String` : type of interaction.
- `Ecutoff::Float64` : cutoff energy,

# Output Argument(s)
- `σa::Float64` : absorption cross-section.

"""
function acs(this::Bremsstrahlung,Ei::Float64,Z::Int64,Ec::Float64,particle::Particle,Ecutoff::Float64)
    return tcs(this,Ei,Z,Ec,particle,[Ecutoff])
end

"""
    sp(this::Bremsstrahlung,Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,Ei::Float64,
    Ec::Float64,Eout::Vector{Float64},particle::Particle)

Gives the stopping power for bremsstrahlung.

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure. 
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `ωz::Vector{Float64}` : weight fraction of the elements composing the material.
- `ρ::Float64` : material density.
- `Ei::Float64` : incoming particle energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `Eout::Vector{Float64}` : energy boundaries associated with outgoing particles.
- `particle::Particle` : incoming particle.

# Output Argument(s)
- `S::Float64` : stopping power.

"""
function sp(this::Bremsstrahlung,Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,Ei::Float64,Ec::Float64,Eout::Vector{Float64},particle::Particle)

    # Initialization
    𝒩ₙ = nuclei_density.(Z,ρ)
    Nz = length(Z)

    # Compute the total stopping power 
    St = 0.0
    for iz in range(1,Nz)
        St += ωz[iz] * 𝒩ₙ[iz] * seltzer_berger_stopping_power(Z[iz],Ei,particle)
    end

    # Compute the catastrophic stopping power
    Sc = 0.0
    Ngf = length(Eout)-1
    is_dirac, Np, q_type = out_distribution(this)
    if is_dirac Np = 1; u = [0]; w = [2] else u,w = quadrature(Np,q_type) end
    for iz in range(1,Nz)
        for gf in range(1,Ngf+1)
            Ef⁻ = Eout[gf]
            if (gf != Ngf+1) Ef⁺ = Eout[gf+1] else Ef⁺ = 0.0 end
            Ef⁻,Ef⁺,isSkip = bounds(this,Ef⁻,Ef⁺,Ei,"S",Ec)
            if isSkip continue end
            ΔEf = Ef⁻ - Ef⁺
            for n in range(1,Np)
                Ef = (u[n]*ΔEf + (Ef⁻+Ef⁺))/2
                Eγ = Ei-Ef
                if Ei ≥ Ef && ΔEf ≥ 0
                    Sc += ωz[iz] * ΔEf/2 * w[n] * 𝒩ₙ[iz] * Eγ * seltzer_berger_cross_section(Z[iz],Ei,Eγ,particle)
                end
            end
        end
    end

    # Compute the soft stopping powers
    S = St-Sc
    return S
end

"""
    mt(this::Bremsstrahlung)

Gives the momentum transfer for bremsstrahlung.

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure. 

# Output Argument(s)
- `T::Float64` : momentum transfer.

"""
function mt(this::Bremsstrahlung)
    T = 0.0 # Because μ = 1.0
    return T
end