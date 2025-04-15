"""
    Elastic_Leptons

Structure used to define parameters for production of multigroup elastic cross-sections for leptons.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}} = Dict((Electron,Electron) => ["S"],(Positron,Positron) => ["S"])` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `(Electron,Electron) => ["S"]` : elastic interaction of electrons.
    - `(Positron,Positron) => ["S"]` : elastic interaction of positrons.

"""
mutable struct Elastic_Leptons <: Interaction

    # Variable(s)
    name::String
    interaction_types::Dict{Tuple{Type,Type},Vector{String}}
    incoming_particle::Vector{Type}
    interaction_particles::Vector{Type}
    is_CSD::Bool
    is_elastic::Bool
    is_ETC::Bool
    is_AFP::Bool
    is_AFP_decomposition::Bool
    is_subshells_dependant::Bool
    is_kawrakow_correction::Bool
    is_seltzer_correction::Bool
    is_subshell_inelastic::Bool
    model::String
    solver::String
    scattering_model::String

    # Constructor(s)
    function Elastic_Leptons()
        this = new()
        this.name = "elastic_leptons"
        this.interaction_types = Dict((Electron,Electron) => ["S"],(Positron,Positron) => ["S"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.is_CSD = false
        this.is_elastic = true
        this.is_ETC = true
        this.is_AFP = false
        this.is_AFP_decomposition = true
        this.is_subshells_dependant = false
        this.model = "boschini"
        this.is_kawrakow_correction = true
        this.is_seltzer_correction = true
        this.solver = "BFP"
        this.scattering_model="BFP"
        return this
    end

end

# Method(s)
"""
    set_interaction_types(this::Elastic_Leptons,interaction_types::Dict{Tuple{DataType,DataType},Vector{String}})

To define the interaction types for Elastic_Leptons processes.

# Input Argument(s)
- `this::Elastic_Leptons` : elastic leptons structure.
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}}` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `(Electron,Electron) => ["S"]` : elastic interaction of electrons.
    - `(Positron,Positron) => ["S"]` : elastic interaction of positrons.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_leptons = Elastic_Leptons()
julia> elastic_leptons.set_interaction_types( Dict((Positron,Positron) => ["S"]) ) # Elastic only for positrons
```
"""
function set_interaction_types(this::Elastic_Leptons,interaction_types)
    this.interaction_types = interaction_types
end

"""
    set_model(this::Elastic_Leptons,model::String,is_KC::Bool=true)

To define the elastic scattering model.

# Input Argument(s)
- `this::Elastic_Leptons` : elastic leptons structure.
- `model::String` : model of elastic scattering:
    - `rutherford` : screened Rutherford cross-sections.
    - `boschini` : screened Mott cross-sections.
- `is_KC::Bool` : Apply Karakow correction (true) or not (false).
- `is_SC::Bool` : Apply Seltzer correction (true) or not (false).

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_leptons = Elastic_Leptons()
julia> elastic_leptons.set_model("rutherford")
```
"""
function set_model(this::Elastic_Leptons,model::String,is_KC::Bool=false,is_SC::Bool=false)
    if lowercase(model) ∉ ["rutherford","boschini"] error("Unkown elastic model: $model.") end
    this.model = lowercase(model)
    this.is_kawrakow_correction = is_KC
    this.is_seltzer_correction = is_SC
end

"""
    set_transport_correction(this::Elastic_Leptons,is_ETC::Bool)

Enable or not extended transport correcton.

# Input Argument(s)
- `this::Elastic_Leptons` : elastic leptons structure.
- `is_ETC::Bool` : Enable (true) or not (false) extended transport correcton.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_leptons = Elastic_Leptons()
julia> elastic_leptons.is_ETC(false)
```
"""
function set_transport_correction(this::Elastic_Leptons,is_ETC::Bool)
    this.is_ETC = is_ETC
end

"""
    set_solver(this::Elastic_Leptons,solver::String)

Dictate how the cross-sections is distributed to the Boltzmann and the Fokker-Planck operators.

# Input Argument(s)
- `this::Elastic_Leptons` : elastic leptons structure.
- `solver::String` :  model of elastic scattering:
    - `BTE` : the cross-sections are made to be used with the Boltzmann operator.
    - `FP` : the cross-sections are made to be used with the Fokker-Planck operator.
    - `BFP` : the cross-sections are distributed to the Boltzmann and Fokker-Planck operator.
# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_leptons = Elastic_Leptons()
julia> elastic_leptons.is_AFP(false)
```
"""
function set_solver(this::Elastic_Leptons,solver::String)
    if uppercase(solver) ∉ ["BTE","FP","BFP"] error("Unkown elastic model: $model.") end
    if uppercase(solver) ∈ ["FP","BFP"] this.is_AFP = true else this.is_AFP = false end
    this.solver = uppercase(solver)
    this.scattering_model = uppercase(solver)
end

"""
    in_distribution(this::Elastic_Leptons)

Describe the energy discretization method for the incoming particle in the elastic lepton
interaction.

# Input Argument(s)
- `this::Elastic_Leptons` : elastic lepton structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function in_distribution(this::Elastic_Leptons)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    out_distribution(this::Elastic_Leptons)

Describe the energy discretization method for the outgoing particle in the elastic lepton
interaction.

# Input Argument(s)
- `this::Elastic_Leptons` : elastic lepton structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function out_distribution(this::Elastic_Leptons)
    is_dirac = true
    N = 1
    quadrature = "dirac"
    return is_dirac, N, quadrature
end

"""
    bounds(this::Elastic_Leptons,Ef⁻::Float64,Ef⁺::Float64,gi::Int64,gf::Int64)

Gives the integration energy bounds for the outgoing particle for elastic lepton
interaction. 

# Input Argument(s)
- `this::Elastic_Leptons` : elastic lepton structure.
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `gi::Int64` : group of the incoming particle.
- `gf::Int64` : group of the outgoing particle.

# Output Argument(s)
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `isSkip::Bool` : define if the integration is skipped or not.

"""
function bounds(this::Elastic_Leptons,Ef⁻::Float64,Ef⁺::Float64,gi::Int64,gf::Int64)
    if (gf != gi) isSkip = true else isSkip = false end # Elastic scattering only
    return Ef⁻,Ef⁺,isSkip
end

"""
    dcs(this::Elastic_Leptons,L::Int64,Ei::Float64,Z::Int64,particle::Particle,
    Ecutoff::Float64,iz::Int64)

Gives the Legendre moments of the scattering cross-sections for elastic lepton
interaction. 

# Input Argument(s)
- `this::Annihilation` : elastic lepton structure.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.
- `particle::Particle` : incoming particle.
- `Ecutoff::Float64` : cutoff energy.
- `iz::Int64` : index of the element in the material.

# Output Argument(s)
- `σℓ::Vector{Float64}` : Legendre moments of the scattering cross-sections.

"""
function dcs(this::Elastic_Leptons,L::Int64,Ei::Float64,Z::Int64,particle::Particle,Ecutoff::Float64)
    σℓ = mott(Z,Ei,particle,L,Ecutoff,this.is_seltzer_correction,this.is_kawrakow_correction,this.is_subshell_inelastic,this.model)
    return σℓ
end

"""
    tcs(this::Elastic_Leptons,Ei::Float64,Z::Int64,particle::Particle,Ecutoff::Float64,iz::Int64)

Gives the total cross-section for elastic lepton interaction. 

# Input Argument(s)
- `this::Elastic_Leptons` : elastic lepton structure.
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.
- `particle::Particle` : incoming particle.
- `Ecutoff::Float64` : cutoff energy.
- `iz::Int64` : element index in the material.

# Output Argument(s)
- `σt::Float64` : total cross-section.

"""
function tcs(this::Elastic_Leptons,Ei::Float64,Z::Int64,particle::Particle,Ecutoff::Float64)
    σt = mott(Z,Ei,particle,0,Ecutoff,this.is_seltzer_correction,this.is_kawrakow_correction,this.is_subshell_inelastic,this.model)[1]
    return σt
end

"""
    acs(this::Elastic_Leptons)

Gives the absorption cross-section for elastic lepton interaction. 

# Input Argument(s)
- `this::Elastic_Leptons` : elastic lepton structure.

# Output Argument(s)
- `σt::Float64` : absorption cross-section.

"""
function acs(this::Elastic_Leptons)
    return 0
end