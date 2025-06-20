"""
Inelastic_Collision

Structure used to define parameters for production of multigroup inelastic collisional cross-sections.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}} = Dict((Positron,Positron) => ["S"],(Positron,Electron) => ["P"],(Electron,Electron) => ["S","P"])` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `(Positron,Positron) => ["S"]` : scattering of incident positrons
    - `(Electron,Electron) => ["S"]` : scattering of incident electrons
    - `(Positron,Electron) => ["P"]` : production of electrons by incident positrons
    - `(Electron,Electron) => ["P"]` : production of electrons by incident electrons

"""
mutable struct Inelastic_Collision <: Interaction

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
    is_shell_correction::Bool
    is_focusing_møller::Bool
    is_hydrogenic_distribution_term::Bool
    is_distant_collision::Bool
    density_correction::String
    scattering_model::String

    # Constructor(s)
    function Inelastic_Collision()
        this = new()
        this.name = "Inelastic_Collision"
        this.interaction_types = Dict((Positron,Positron) => ["S"],(Positron,Electron) => ["P"],(Electron,Electron) => ["S","P"],(Proton,Proton) => ["S"],(Proton,Electron) => ["P"],(Alpha,Alpha) => ["S"],(Alpha,Electron) => ["S"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.set_density_correction("fano")
        this.is_CSD = true
        this.is_AFP = true
        this.is_AFP_decomposition = false
        this.is_elastic = false
        this.is_focusing_møller = false
        this.is_hydrogenic_distribution_term = false
        this.set_is_subshells_dependant(true)
        this.set_is_shell_correction(true)
        this.set_scattering_model("BFP")
        return this
    end
end

# Method(s)
"""
    set_interaction_types(this::Inelastic_Collision,interaction_types::Dict{Tuple{DataType,DataType},Vector{String}})

To define the interaction types for inelastic collision processes.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}}` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `(Positron,Positron) => ["S"]` : scattering of incident positrons
    - `(Electron,Electron) => ["S"]` : scattering of incident electrons
    - `(Positron,Electron) => ["P"]` : production of electrons by incident positrons
    - `(Electron,Electron) => ["P"]` : production of electrons by incident electrons

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_collision = Inelastic_Collision()
julia> elastic_collision.set_interaction_types( Dict((Electron,Electron) => ["S"]) ) # No knock-on electrons.
```
"""
function set_interaction_types(this::Inelastic_Collision,interaction_types)
    this.interaction_types = interaction_types
end

"""
    set_density_correction(this::Inelastic_Collision,density_correction::String)

Set the Fermi density correction.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `density_correction::String` : type of density effect:
    - `fano` : Fano density effect.
    - `sternheimer` : Sternheimer semi-empirical density effect.
    - `none` : no density effect.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_collision = Inelastic_Collision()
julia> elastic_collision.set_density_correction("sternheimer")
```
"""
function set_density_correction(this::Inelastic_Collision,density_correction::String)
    if lowercase(density_correction) ∉ ["fano","sternheimer","none"] error("Unknown '$density_correction' density correction.") end
    this.density_correction = lowercase(density_correction)
end

"""
    set_is_shell_correction(this::Inelastic_Collision,is_shell_correction::Bool)

Activate or desactivate shell correction for stopping powers.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `is_shell_correction::Bool` : activate (true) or desactivate (false) the shell correction.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_collision = Inelastic_Collision()
julia> elastic_collision.set_is_shell_correction(false)
```
"""
function set_is_shell_correction(this::Inelastic_Collision,is_shell_correction::Bool)
    this.is_shell_correction = is_shell_correction
end

"""
    set_is_subshells_dependant(this::Inelastic_Collision,is_subshells_dependant::Bool)

Compute the inelastic cross-sections assuming bounded or unbounded electrons.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `is_subshells_dependant::Bool` : bounded (true) or unbounded (false) electrons.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_collision = Inelastic_Collision()
julia> elastic_collision.set_is_subshells_dependant(false)
```
"""
function set_is_subshells_dependant(this::Inelastic_Collision,is_subshells_dependant::Bool)
    this.is_subshells_dependant = is_subshells_dependant
end

"""
    scattering_model(this::Inelastic_Collision,scattering_model::String)

To define the solver for inelastic scattering.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic structure.
- `solver::String` : solver for inelastic scattering, which can be:
    - `BFP` : Boltzmann Fokker-Planck solver.
    - `FP` : Fokker-Planck solver.
    - `BTE` : Boltzmann solver.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_collision = Inelastic_Collision()
julia> elastic_collision.set_scattering_model("FP")
```
"""
function set_scattering_model(this::Inelastic_Collision,scattering_model::String)
    if uppercase(scattering_model) ∉ ["BFP","FP","BTE"] error("Unknown scattering model (should be BFP, FP or BTE).") end
    this.scattering_model = uppercase(scattering_model)
end

"""
    in_distribution(this::Inelastic_Collision)

Describe the energy discretization method for the incoming particle in the inelastic collision
interaction.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function in_distribution(this::Inelastic_Collision)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    out_distribution(this::Inelastic_Collision)

Describe the energy discretization method for the outgoing particle in the inelastic collision
interaction.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function out_distribution(this::Inelastic_Collision)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    bounds(this::Inelastic_Collision,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String,
    Ec::Float64,Ui::Float64,particle::Particle)

Gives the integration energy bounds for the outgoing particle for inelastic collision
interaction. 

# Input Argument(s)
- `this::Annihilation` : annihilation structure.
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `Ei::Float64` : energy of the incoming particle.
- `type::String` : type of interaction.
- `Ec::Float64` : energy cutoff between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.

# Output Argument(s)
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `isSkip::Bool` : define if the integration is skipped or not.

"""
function bounds(this::Inelastic_Collision,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String,Ec::Float64,Ui::Float64,particle::Particle)
    if is_electron(particle)
        # Scattered electron
        if type == "S"
            Ef⁻ = min(Ef⁻,Ec,Ei-Ui)
            Ef⁺ = max(Ef⁺,(Ei-Ui)/2)
        # Knock-on electron
        elseif type == "P"
            Ef⁻ = min(Ef⁻,(Ei-Ui)/2)
            Ef⁺ = max(Ef⁺,0.0)
        else
            error("Unknown type.")
        end
    elseif is_positron(particle)
        # Scattered positron/proton/alpha
        if type == "S"
            Ef⁻ = min(Ef⁻,Ec,Ei-Ui)
            Ef⁺ = max(Ef⁺,0.0)
        # Knock-on electron
        elseif type == "P"
            Ef⁻ = min(Ef⁻,Ei-Ui)
            Ef⁺ = max(Ef⁺,0.0)
        else
            error("Unknown type.")
        end
    elseif is_proton(particle) || is_alpha(particle)
        mₑ = 0.51099895069
        ratio_mass = get_mass(particle)/mₑ
        γ = (Ei + ratio_mass)/ratio_mass
        β² = (γ^2-1)/γ^2
        R = 1/(1+(1/ratio_mass)^2+2*γ/ratio_mass)
        W_ridge = 2*β²*γ^2*R
        # Scattered heavy particle
        if type == "S"
            Ef⁻ = min(Ef⁻,Ec,Ei-Ui)
            Ef⁺ = max(Ef⁺,0.0,Ei-W_ridge)
        # Knock-on electron
        elseif type == "P"
            Ef⁻ = min(Ef⁻,Ei-Ui,W_ridge-Ui)
            Ef⁺ = max(Ef⁺,0.0)
        else
            error("Unknown type.")
        end
    else
        error("Unknown particle")
    end
    if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    return Ef⁻,Ef⁺,isSkip
end

"""
    dcs(this::Inelastic_Collision,L::Int64,Ei::Float64,Ef::Float64,type::String,
    particle::Particle,Ui::Float64,Zi::Real,Ti::Float64)

Gives the Legendre moments of the scattering cross-sections for inelastic collision
interaction. 

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `Z::Int64` : atomic number.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : incoming particle energy.
- `Ef::Float64` : outgoing particle energy.
- `type::String` : type of interaction.
- `particle::Particle` : incoming particle.
- `δi::Int64` : subshell index.
- `Ui::Float64` : binding energy of the subshell.
- `Zi::Int64` : number of electrons in the subshell.
- `Ti::Float64` : mean kinetic energy of electrons in the subshell.
- `ri::Float64` : average radius per subshell.

# Output Argument(s)
- `σℓ::Vector{Float64}` : Legendre moments of the scattering cross-sections.

"""
function dcs(this::Inelastic_Collision,Z::Int64,L::Int64,Ei::Float64,Ef::Float64,type::String,particle::Particle,δi::Int64,Ui::Float64,Zi::Real,Ti::Float64,ri::Float64)

    #----
    # Initialization
    #----
    σs = 0.0
    σℓ = zeros(L+1)
    if type == "S"
        E⁺ = Ef
        E⁻ = Ei-Ui-E⁺
        W = E⁻+Ui
    elseif type == "P"
        E⁻ = Ef
        E⁺ = Ei-Ui-E⁻
        W = E⁻+Ui
    else
        error("Unknown type")
    end
    if W < 0 return σℓ end

    #----
    # Close collisions
    #----
    if is_electron(particle)
        σs += moller(Zi,Ei,W,Ui,Ti,this.is_focusing_møller,this.is_hydrogenic_distribution_term)
    elseif is_positron(particle)
        σs += bhabha(Zi,Ei,W)
    elseif is_proton(particle) || is_alpha(particle)
        σs += inelastic_collision_heavy_particle(Zi,Ei,W,particle)
    else
        error("Unknown particle")
    end

    #----
    # Compute the angular distribution
    #----
    Wℓ = angular_moller(Ei,Ef,L)

    #----
    # Compute the Legendre moments of the cross-section
    #----
    for ℓ in range(0,L) σℓ[ℓ+1] += Wℓ[ℓ+1] * σs end
    return σℓ
end

"""
    tcs(this::Inelastic_Collision,Ei::Float64,Ec::Float64,particle::Particle,Z::Int64)

Gives the total cross-section for inelastic collision interaction. 

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure. 
- `Ei::Float64` : incoming particle energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.
- `Z::Int64` : atomic number.
- `Eout::Vector{Float64}` : energy boundaries associated with outgoing particles.

# Output Argument(s)
- `σt::Float64` : total cross-section.

"""
function tcs(this::Inelastic_Collision,Ei::Float64,Ec::Float64,particle::Particle,Z::Int64)

    #----
    # Initialization
    #----
    σt = 0

    #----
    # Close collisions
    #----
    if is_electron(particle)
        σt += integrate_moller(Z,Ei,0,Ei-Ec,this.is_focusing_møller,this.is_hydrogenic_distribution_term)
    elseif is_positron(particle)
        σt += integrate_bhabha(Z,Ei,0,Ei-Ec)
    elseif is_proton(particle) || is_alpha(particle)
        σt += integrate_inelastic_collision_heavy_particle(Z,Ei,0,particle,Ei-Ec)
    else
        error("Unknown particle")
    end

    return σt
end

"""
    acs(this::Inelastic_Collision,Ei::Float64,Ec::Float64,particle::Particle,Z::Int64)

Gives the absorption cross-section for inelastic collision interaction. 

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure. 
- `Ei::Float64` : incoming particle energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.
- `Z::Int64` : atomic number.
- `Ecutoff::Float64` : cutoff energy.

# Output Argument(s)
- `σa::Float64` : absorption cross-section.

"""
function acs(this::Inelastic_Collision,Ei::Float64,Ec::Float64,particle::Particle,Z::Int64,Ecutoff::Float64)

    #----
    # Close collisions
    #----
    if is_electron(particle)
        σa = integrate_moller(Z,Ei,0,Ei-min(Ec,Ecutoff),this.is_focusing_møller,this.is_hydrogenic_distribution_term)
    elseif is_positron(particle)
        σa = integrate_bhabha(Z,Ei,0,Ei-min(Ec,Ecutoff))
    else
        error("Unknown particle")
    end
    return σa
    
end

"""
    sp(this::Inelastic_Collision,Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,
    state_of_matter::String,Ei::Float64,Ec::Float64,particle::Particle)

Gives the stopping power for inelastic collision interaction.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure. 
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `ωz::Vector{Float64}` : weight fraction of the elements composing the material.
- `ρ::Float64` : material density.
- `state_of_matter::String` : state of matter.
- `Ei::Float64` : incoming particle energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.
- `Eout::Vector{Float64}` : energy boundaries associated with outgoing particles.

# Output Argument(s)
- `S::Float64` : stopping power.

"""
function sp(this::Inelastic_Collision,Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,state_of_matter::String,Ei::Float64,Ec::Float64,particle::Particle)

    # Compute the total cross-section
    Stot = bethe(Z,ωz,ρ,Ei,particle,this.density_correction,state_of_matter)
    
    # Compute the catastrophic Møller- or Bhabha- derived stopping power
    Sc = 0
    Nz = length(Z)
    for i in range(1,Nz)

        # Close collision
        if is_electron(particle)
            Sc += ωz[i] * nuclei_density(Z[i],ρ) * integrate_moller(Z[i],Ei,1,Ei-Ec,this.is_focusing_møller,this.is_hydrogenic_distribution_term)
        elseif is_positron(particle)
            Sc += ωz[i] * nuclei_density(Z[i],ρ) * integrate_bhabha(Z[i],Ei,1,Ei-Ec)
        elseif is_proton(particle) || is_alpha(particle)
            Sc += ωz[i] * nuclei_density(Z[i],ρ) * integrate_inelastic_collision_heavy_particle(Z[i],Ei,1,particle,Ei-Ec)
        end
    end

    # Compute the soft stopping power
    Sr = Stot - Sc

    return Sr
end

"""
    mt(this::Inelastic_Collision)

Gives the momentum transfer for inelastic collision interaction.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure. 

# Output Argument(s)
- `T::Float64` : momentum transfer.

"""
function mt(this::Inelastic_Collision)
    T = 0.0
    return T
end

