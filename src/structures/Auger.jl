"""
    Auger

Structure used to define parameters for production of multigroup annihilation cross-sections.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `ηmin::Float64 = 0.001` : minimum probability of the production of specific Auger electrons following electron cascades.
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}} = Dict((Photon,Electron) => ["P"],(Electron,Electron) => ["P"],(Positron,Electron) => ["P"])` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `(Photon,Electron) => ["P"]` : production of Auger electron following incident photon ionization of subshells (by photoelectric effect).
    - `(Electron,Electron) => ["P"]` : production of Auger electron following incident electrons ionization of subshells (by Møller interaction).
    - `(Positron,Electron) => ["P"]` : production of Auger electron following incident positrons ionization of subshells (by Bhabha interaction).

"""
mutable struct Auger <: Interaction

    # Variable(s)
    name::String
    incoming_particle::Vector{Type}
    interaction_particles::Vector{Type}
    interaction_types::Dict{Tuple{Type,Type},Vector{String}}
    is_CSD::Bool
    is_AFP::Bool
    is_AFP_decomposition::Bool
    is_elastic::Bool
    is_preload_data::Bool
    is_subshells_dependant::Bool
    photoelectric_cross_sections::Function
    inelastic::Vector{Vector{Function}}
    η::Vector{Vector{Vector{Float64}}}
    ΔE::Vector{Vector{Vector{Float64}}}
    ηmin::Float64
    scattering_model::String

    # Constructor(s)
    function Auger()
        this = new()
        this.name = "auger"
        this.interaction_types = Dict((Photon,Electron) => ["P"],(Electron,Electron) => ["P"],(Positron,Electron) => ["P"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.is_CSD = false
        this.is_AFP = false
        this.is_AFP_decomposition = false
        this.is_elastic = false
        this.is_preload_data = true
        this.is_subshells_dependant = true
        this.set_minimum_probability(0.001)
        this.scattering_model = "BTE"
        return this
    end

end

# Method(s)
"""
    set_interaction_types(this::Auger,interaction_types::Dict{Tuple{DataType,DataType},Vector{String}})

To define the interaction types for Auger processes.

# Input Argument(s)
- `this::Auger` : Auger structure.
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}}` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `(Photon,Electron) => ["P"]` : production of Auger electron following incident photon ionization of subshells (by photoelectric effect).
    - `(Electron,Electron) => ["P"]` : production of Auger electron following incident electrons ionization of subshells (by Møller interaction).
    - `(Positron,Electron) => ["P"]` : production of Auger electron following incident positrons ionization of subshells (by Bhabha interaction).

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> auger = Auger()
julia> auger.set_interaction_types( Dict((Electron,Electron) => ["P"]) ) # Only cascades following Møller interactions.
```
"""
function set_interaction_types(this::Auger,interaction_types)
    this.interaction_types = interaction_types
end

"""
    set_minimum_probability(this::Auger,ηmin::Real)

To define the minimum probability of a specific Auger electron production.

# Input Argument(s)
- `this::Auger` : Auger structure.
- `ηmin::Real` : minimum probability of a specific Auger electron production.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> auger = Auger()
julia> auger.set_minimum_probability(0.1) # Only Auger electron with probability greater than 10%.
```
"""
function set_minimum_probability(this::Auger,ηmin::Real)
    if ~(0 ≤ ηmin ≤ 1) error("Probability should be between 0 and 1.") end
    this.ηmin = ηmin
end

"""
    in_distribution(this::Auger)

Describe the energy discretization method for the incoming particle in the Auger
electron production interaction.

# Input Argument(s)
- `this::Auger` : Auger structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function in_distribution(this::Auger)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    out_distribution(this::Auger)

Describe the energy discretization method for the outgoing particle in the Auger
electron production interaction.

# Input Argument(s)
- `this::Auger` : Auger structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function out_distribution(this::Auger)
    is_dirac = true
    N = 1
    quadrature = "dirac"
    return is_dirac, N, quadrature
end

"""
    bounds(this::Auger,Ef⁻::Float64,Ef⁺::Float64)

Gives the integration energy bounds for the outgoing particle for Auger electron
production. 

# Input Argument(s)
- `this::Auger` : Auger structure.
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.

# Output Argument(s)
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `isSkip::Bool` : define if the integration is skipped or not.

"""
function bounds(this::Auger,Ef⁻::Float64,Ef⁺::Float64)
    if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    return Ef⁻,Ef⁺,isSkip
end

"""
    dcs(this::Auger,L::Int64,Ei::Float64,Z::Int64,iz::Int64,δi::Int64,Ef⁻::Float64,
    Ef⁺::Float64,particle::Particle)

Gives the Legendre moments of the scattering cross-sections for Auger electron production. 

# Input Argument(s)
- `this::Auger` : Auger structure.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.
- `iz::Int64` : index of the element in the material.
- `δi::Int64` : index of the subshell.
- `Ef⁻::Float64` : upper bound of the outgoing energy group.
- `Ef⁺::Float64` : lower bound of the outgoing energy group.
- `particle::Particle` : incoming particle.

# Output Argument(s)
- `σℓ::Vector{Float64}` : Legendre moments of the scattering cross-sections.

"""
function dcs(this::Auger,L::Int64,Ei::Float64,Z::Int64,iz::Int64,δi::Int64,Ef⁻::Float64,Ef⁺::Float64,particle::Particle)

    # Photoelectric cross section
    if is_photon(particle)
        _,_,_,_,_,subshells = electron_subshells(Z)
        σa = this.photoelectric_cross_sections(iz,Ei,subshells[δi])
    # Load election impact ionization
    elseif get_type(particle) ∈ [Electron,Positron]
        σa = this.inelastic[iz][δi](Ei)
    end

    # Fluorescence cross section
    σs = 0
    Nt = length(this.ΔE[iz][δi])
    for δj in range(1,Nt)
        if Ef⁺ ≤ this.ΔE[iz][δi][δj] ≤ Ef⁻
            σs += this.η[iz][δi][δj] * σa
        end
    end

    # Angular distribution (isotropic)
    σℓ = zeros(L+1)
    σℓ[1] = σs

    return σℓ
end

"""
    preload_data(this::Auger,Z::Vector{Int64},Emax::Float64,Emin::Float64,L::Int64,
    type::String,Eout::Vector{Float64},interactions::Vector{Interaction})

Preload data for multigroup Auger electron production calculations. 

# Input Argument(s)
- `this::Auger` : Auger structure. 
- `Z::Vector{Int64}` : atomic numbers of the elements in the material. 
- `ρ::Float64` : material density.
- `particle::Particle` : incoming particle.
- `Ecutoff::Float64` : cutoff energy.

# Output Argument(s)
N/A

"""
function preload_data(this::Auger,Z::Vector{Int64},ρ::Float64,particle::Particle,Ecutoff::Float64)
    
    if is_photon(particle)
        # Load photoelectric cross-sections
        photoelectric = Photoelectric()
        photoelectric.model = "jendl5"
        photoelectric.is_subshells_dependant = true
        photoelectric.preload_photoelectric_cross_sections(Z,ρ)
        this.photoelectric_cross_sections = photoelectric.photoelectric_cross_sections
    elseif is_electron(particle)
        # Load election impact ionization
        rₑ = 2.81794092E-13       # (in cm)
        Nz = length(Z)
        this.inelastic = Vector{Vector{Function}}(undef,Nz)
        for iz in range(1,Nz)
            Nshells,Zi,Ui,Ti,ri,subshells = electron_subshells(Z[iz])
            this.inelastic[iz] = Vector{Function}(undef,Nshells)
            for δi in range(1,Nshells)
                this.inelastic[iz][δi] = function inelastic_for_auger(Ei::Float64)
                    Wmax = (Ei-Ui[δi])/2; Wmin = 0.0
                    if Wmax > Wmin
                        β² = Ei*(Ei+2)/(Ei+1)^2
                        Jauger(x) = -1/(x+Ui[δi]) + 1/(Ei-x) + x/(Ei+1)^2 - (2*Ei+1)/(Ei+1)^2 * (log(x+Ui[δi])-log(Ei-x))/(Ei+Ui[δi])
                        return 2*π*rₑ^2/β² * Zi[δi] * (Jauger(Wmax)-Jauger(Wmin))
                    else
                        return 0.0
                    end
                end
            end
        end
    elseif is_positron(particle)
        # Load positron impact ionization
        rₑ = 2.81794092E-13       # (in cm)
        Nz = length(Z)
        this.inelastic = Vector{Vector{Function}}(undef,Nz)
        for iz in range(1,Nz)
            Nshells,Zi,Ui,Ti,ri,subshells = electron_subshells(Z[iz])
            this.inelastic[iz] = Vector{Function}(undef,Nshells)
            for δi in range(1,Nshells)
                this.inelastic[iz][δi] = function inelastic_for_auger2(Ei::Float64)
                    Wmax = Ei-Ui[δi]; Wmin = Ui[δi]
                    if Wmax > Wmin
                        γ = Ei+1
                        β² = Ei*(Ei+2)/(Ei+1)^2
                        b = ((γ-1)/γ)^2
                        b1 = b * (2*(γ+1)^2-1)/(γ^2-1)
                        b2 = b * (3*(γ+1)^2+1)/(γ+1)^2
                        b3 = b * (2*(γ-1)*γ)/(γ+1)^2
                        b4 = b * (γ-1)^2/(γ+1)^2
                        Jauger2(x) = -1/x - b1*log(x)/Ei + b2*x/Ei^2 - b3*x^2/(2*Ei^3) + b4*x^3/(3*Ei^4)
                        return 2*π*rₑ^2/β² * Zi[δi] * (Jauger2(Wmax)-Jauger2(Wmin))
                    else
                        return 0.0
                    end
                end
            end
        end
    else
        error("Unknown particles.")
    end

    # Load fluorescence data
    this.ΔE, this.η = atomic_electron_cascades("auger",Z,Ecutoff,this.ηmin)
    
end