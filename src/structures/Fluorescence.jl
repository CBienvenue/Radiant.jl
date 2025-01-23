"""
    Fluorescence

Structure used to define parameters for production of multigroup fluorescence cross-sections.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}} = Dict((Photon,Photon) => ["P"],(Electron,Photon) => ["P"],(Positron,Photon) => ["P"])` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `(Photon,Photon) => ["P"]` : production of fluorescence following incident photon ionization of subshells (by photoelectric effect).
    - `(Electron,Photon) => ["P"]` : production of fluorescence following incident electrons ionization of subshells (by Møller interaction).
    - `(Positron,Photon) => ["P"]` : production of fluorescence following incident positrons ionization of subshells (by Bhabha interaction).
- `ηmin::Float64=0.001` : minimum probability of the production of specific fluorescence photon following electron cascades.

"""
mutable struct Fluorescence <: Interaction

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
    inelastic::Vector{Vector{Function}}
    η::Vector{Vector{Vector{Float64}}}
    ΔE::Vector{Vector{Vector{Float64}}}
    ηmin::Float64
    scattering_model::String

    # Constructor(s)
    function Fluorescence()
        this = new()
        this.name = "fluorescence"
        this.interaction_types = Dict((Photon,Photon) => ["P"],(Electron,Photon) => ["P"],(Positron,Photon) => ["P"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.is_CSD = false
        this.is_AFP = false
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
    set_interaction_types(this::Fluorescence,interaction_types::Dict{Tuple{DataType,DataType},Vector{String}})

To define the interaction types for fluorescence processes.

# Input Argument(s)
- `this::Fluorescence` : fluorescence structure.
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}}` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `(Photon,Photon) => ["P"]` : production of fluorescence following incident photon ionization of subshells (by photoelectric effect).
    - `(Electron,Photon) => ["P"]` : production of fluorescence following incident electrons ionization of subshells (by Møller interaction).
    - `(Positron,Photon) => ["P"]` : production of fluorescence following incident positrons ionization of subshells (by Bhabha interaction).

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> fluorescence = Fluorescence()
julia> fluorescence.set_interaction_types( Dict((Electron,Photon) => ["P"]) ) # Only cascades following Møller interactions.
```
"""
function set_interaction_types(this::Fluorescence,interaction_types)
    this.interaction_types = interaction_types
end

"""
    set_minimum_probability(this::Fluorescence,ηmin::Real)

To define the minimum probability of a specific fluorescence production.

# Input Argument(s)
- `this::Fluorescence` : fluorescence structure.
- `ηmin::Real` : minimum probability of a specific fluorescence production.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> fluorescence = Fluorescence()
julia> fluorescence.set_minimum_probability(0.1) # Only fluorescence photons with probability greater than 10%.
```
"""
function set_minimum_probability(this::Fluorescence,ηmin::Real)
    if ~(0 ≤ ηmin ≤ 1) error("Probability should be between 0 and 1.") end
    this.ηmin = ηmin
end

function in_distribution(this::Fluorescence)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

function out_distribution(this::Fluorescence)
    is_dirac = true
    N = 1
    quadrature = "dirac"
    return is_dirac, N, quadrature
end

function bounds(this::Fluorescence,Ef⁻::Float64,Ef⁺::Float64,gi::Int64,type::String,Ui::Float64,E_in::Vector{Float64})
    if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    return Ef⁻,Ef⁺,isSkip
end

function dcs(this::Fluorescence,L::Int64,Ei::Float64,Z::Int64,iz::Int64,δi::Int64,Ef⁻::Float64,Ef⁺::Float64,particle::Particle)

    # Photoelectric cross section
    if is_photon(particle)
        Nshells,Zi,Ui,Ti,ri,subshells = electron_subshells(Z)
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

function preload_data(this::Fluorescence,Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,L::Int64,particle::Particle,Ecutoff::Float64)
    
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
                this.inelastic[iz][δi] = function inelastic_for_fluorescence(Ei::Float64)
                    Wmax = (Ei-Ui[δi])/2; Wmin = 0.0
                    if Wmax > Wmin
                        β² = Ei*(Ei+2)/(Ei+1)^2
                        Jfluorescence(x) = -1/(x+Ui[δi]) + 1/(Ei-x) + x/(Ei+1)^2 - (2*Ei+1)/(Ei+1)^2 * (log(x+Ui[δi])-log(Ei-x))/(Ei+Ui[δi])
                        return 2*π*rₑ^2/β² * Zi[δi] * (Jfluorescence(Wmax)-Jfluorescence(Wmin))
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
                this.inelastic[iz][δi] = function inelastic_for_fluorescence2(Ei::Float64)
                    Wmax = Ei-Ui[δi]; Wmin = Ui[δi]
                    if Wmax > Wmin
                        γ = Ei+1
                        β² = Ei*(Ei+2)/(Ei+1)^2
                        b = ((γ-1)/γ)^2
                        b1 = b * (2*(γ+1)^2-1)/(γ^2-1)
                        b2 = b * (3*(γ+1)^2+1)/(γ+1)^2
                        b3 = b * (2*(γ-1)*γ)/(γ+1)^2
                        b4 = b * (γ-1)^2/(γ+1)^2
                        Jfluorescence2(x) = -1/x - b1*log(x)/Ei + b2*x/Ei^2 - b3*x^2/(2*Ei^3) + b4*x^3/(3*Ei^4)
                        return 2*π*rₑ^2/β² * Zi[δi] * (Jfluorescence2(Wmax)-Jfluorescence2(Wmin))
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
    this.ΔE, this.η = atomic_electron_cascades("fluorescence",Z,Ecutoff,this.ηmin)
    
end