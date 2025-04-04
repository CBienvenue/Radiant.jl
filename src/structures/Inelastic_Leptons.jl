"""
Inelastic_Leptons

Structure used to define parameters for production of multigroup inelastic collisional cross-sections for leptons.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}} = Dict((Positron,Positron) => ["S"],(Positron,Electron) => ["P"],(Electron,Electron) => ["S","P"])` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `(Positron,Positron) => ["S"]` : scattering of incident positrons
    - `(Electron,Electron) => ["S"]` : scattering of incident electrons
    - `(Positron,Electron) => ["P"]` : production of electrons by incident positrons
    - `(Electron,Electron) => ["P"]` : production of electrons by incident electrons

"""
mutable struct Inelastic_Leptons <: Interaction

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
    is_shell_correction::Bool
    density_correction::String
    scattering_model::String

    # Constructor(s)
    function Inelastic_Leptons()
        this = new()
        this.name = "inelastic_leptons"
        this.interaction_types = Dict((Positron,Positron) => ["S"],(Positron,Electron) => ["P"],(Electron,Electron) => ["S","P"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.set_density_correction("fano")
        this.is_CSD = true
        this.is_AFP = true
        this.is_AFP_decomposition = false
        this.is_elastic = false
        this.is_preload_data = false
        this.set_is_subshells_dependant(true)
        this.set_is_shell_correction(true)
        this.set_scattering_model("BFP")
        return this
    end
end

# Method(s)
"""
    set_interaction_types(this::Inelastic_Leptons,interaction_types::Dict{Tuple{DataType,DataType},Vector{String}})

To define the interaction types for inelastic leptons processes.

# Input Argument(s)
- `this::Inelastic_Leptons` : inelastic leptons structure.
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}}` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `(Positron,Positron) => ["S"]` : scattering of incident positrons
    - `(Electron,Electron) => ["S"]` : scattering of incident electrons
    - `(Positron,Electron) => ["P"]` : production of electrons by incident positrons
    - `(Electron,Electron) => ["P"]` : production of electrons by incident electrons

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_leptons = Inelastic_Leptons()
julia> elastic_leptons.set_interaction_types( Dict((Electron,Electron) => ["S"]) ) # No knock-on electrons.
```
"""
function set_interaction_types(this::Inelastic_Leptons,interaction_types)
    this.interaction_types = interaction_types
end

"""
    set_density_correction(this::Inelastic_Leptons,density_correction::String)

Set the Fermi density correction.

# Input Argument(s)
- `this::Inelastic_Leptons` : inelastic leptons structure.
- `density_correction::String` : type of density effect:
    - `fano` : Fano density effect.
    - `sternheimer` : Sternheimer semi-empirical density effect.
    - `none` : no density effect.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_leptons = Inelastic_Leptons()
julia> elastic_leptons.set_density_correction("sternheimer")
```
"""
function set_density_correction(this::Inelastic_Leptons,density_correction::String)
    if lowercase(density_correction) ∉ ["fano","sternheimer","none"] error("Unknown '$density_correction' density correction.") end
    this.density_correction = lowercase(density_correction)
end

"""
    set_is_shell_correction(this::Inelastic_Leptons,is_shell_correction::Bool)

Activate or desactivate shell correction for stopping powers.

# Input Argument(s)
- `this::Inelastic_Leptons` : inelastic leptons structure.
- `is_shell_correction::Bool` : activate (true) or desactivate (false) the shell correction.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_leptons = Inelastic_Leptons()
julia> elastic_leptons.set_is_shell_correction(false)
```
"""
function set_is_shell_correction(this::Inelastic_Leptons,is_shell_correction::Bool)
    this.is_shell_correction = is_shell_correction
end

"""
    set_is_subshells_dependant(this::Inelastic_Leptons,is_subshells_dependant::Bool)

Compute the inelastic cross-sections assuming bounded or unbounded electrons.

# Input Argument(s)
- `this::Inelastic_Leptons` : inelastic leptons structure.
- `is_subshells_dependant::Bool` : bounded (true) or unbounded (false) electrons.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_leptons = Inelastic_Leptons()
julia> elastic_leptons.set_is_subshells_dependant(false)
```
"""
function set_is_subshells_dependant(this::Inelastic_Leptons,is_subshells_dependant::Bool)
    this.is_subshells_dependant = is_subshells_dependant
end

"""
    scattering_model(this::Inelastic_Leptons,scattering_model::String)

To define the solver for inelastic scattering.

# Input Argument(s)
- `this::Inelastic_Leptons` : inelastic structure.
- `solver::String` : solver for inelastic scattering, which can be:
    - `BFP` : Boltzmann Fokker-Planck solver.
    - `FP` : Fokker-Planck solver.
    - `BTE` : Boltzmann solver.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_leptons = Inelastic_Leptons()
julia> elastic_leptons.set_scattering_model("FP")
```
"""
function set_scattering_model(this::Inelastic_Leptons,scattering_model::String)
    if uppercase(scattering_model) ∉ ["BFP","FP","BTE"] error("Unknown scattering model (should be BFP, FP or BTE).") end
    this.scattering_model = uppercase(scattering_model)
end

"""
    in_distribution(this::Inelastic_Leptons)

Describe the energy discretization method for the incoming particle in the inelastic lepton
interaction.

# Input Argument(s)
- `this::Inelastic_Leptons` : inelastic lepton structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function in_distribution(this::Inelastic_Leptons)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    out_distribution(this::Inelastic_Leptons)

Describe the energy discretization method for the outgoing particle in the inelastic lepton
interaction.

# Input Argument(s)
- `this::Inelastic_Leptons` : inelastic lepton structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function out_distribution(this::Inelastic_Leptons)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    bounds(this::Inelastic_Leptons,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String,
    Ec::Float64,Ui::Float64,particle::Particle)

Gives the integration energy bounds for the outgoing particle for inelastic lepton
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
function bounds(this::Inelastic_Leptons,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String,Ec::Float64,Ui::Float64,particle::Particle)
    if is_electron(particle)
        # Scattered electron
        if type == "S"
            Ef⁻ = min(Ef⁻,Ec,Ei-Ui)
            Ef⁺ = max(Ef⁺,(Ei-Ui)/2)
            if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
        # Knock-on electron
        elseif type == "P"
            Ef⁻ = min(Ef⁻,(Ei-Ui)/2)
            Ef⁺ = max(Ef⁺,0.0)
            if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
        else
            error("Unknown type of method for Møller scattering.")
        end
    elseif is_positron(particle)
        # Scattered positron
        if type == "S"
            Ef⁻ = min(Ef⁻,Ec,Ei-Ui)
            Ef⁺ = max(Ef⁺,0.0)
            if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
        # Knock-on electron
        elseif type == "P"
            Ef⁻ = min(Ef⁻,Ei-Ui)
            Ef⁺ = max(Ef⁺,0.0)
            if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
        else
            error("Unknown type of method for Bhabha scattering.")
        end
    else
        error("Unknown particle")
    end
    return Ef⁻,Ef⁺,isSkip
end

"""
    dcs(this::Inelastic_Leptons,L::Int64,Ei::Float64,Ef::Float64,type::String,
    particle::Particle,Ui::Float64,Zi::Real,Ti::Float64)

Gives the Legendre moments of the scattering cross-sections for inelastic lepton
interaction. 

# Input Argument(s)
- `this::Inelastic_Leptons` : inelastic lepton structure.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : incoming particle energy.
- `Ef::Float64` : outgoing particle energy.
- `type::String` : type of interaction.
- `particle::Particle` : incoming particle.
- `Ui::Float64` : binding energy of the subshell.
- `Zi::Int64` : number of electrons in the subshell.
- `Ti::Float64` : mean kinetic energy of electrons in the subshell.

# Output Argument(s)
- `σℓ::Vector{Float64}` : Legendre moments of the scattering cross-sections.

"""
function dcs(this::Inelastic_Leptons,L::Int64,Ei::Float64,Ef::Float64,type::String,particle::Particle,Ui::Float64,Zi::Real,Ti::Float64)

    # Initialization
    rₑ = 2.81794092e-13 # (in cm)
    γ = Ei+1
    β² = Ei*(Ei+2)/(Ei+1)^2
    σs = 0.0
    σℓ = zeros(L+1)
    if type == "S"
        Ep = Ef
        W = (Ei-Ui)-Ep
    elseif type == "P"
        W = Ef
        Ep = (Ei-Ui)-W
    else
        error("Unknown type")
    end

    # Close collisions
    if W ≥ 0
        if is_electron(particle)
           F = 1/(W+Ui)^2 + 1/(Ei-W)^2 + 1/(Ei+1)^2 - (2*Ei+1)/(Ei+1)^2 * 1/((Ei-W)*(W+Ui))
        elseif is_positron(particle)
            b = ((γ-1)/γ)^2
            b1 = b * (2*(γ+1)^2-1)/(γ^2-1)
            b2 = b * (3*(γ+1)^2+1)/(γ+1)^2
            b3 = b * (2*(γ-1)*γ)/(γ+1)^2
            b4 = b * (γ-1)^2/(γ+1)^2
            F = (1 - b1*((W+Ui)/Ei) + b2*((W+Ui)/Ei)^2 - b3*((W+Ui)/Ei)^3 + b4*((W+Ui)/Ei)^4)/(W+Ui)^2
        else
            error("Unknown particle")
        end
        σs += 2*π*rₑ^2/β² * Zi * F
    end

    # Compute the Legendre moments of the flux
    if type == "S" 
        μ = sqrt((Ep*(Ei+2))/(Ei*(Ep+2)))
        Pℓμ = legendre_polynomials(L,μ)
        for ℓ in range(0,L) σℓ[ℓ+1] += Pℓμ[ℓ+1] * σs end
    elseif type == "P"
        μ = sqrt((W*(Ei+2))/(Ei*(W+2)))
        Pℓμ = legendre_polynomials(L,μ)
        for ℓ in range(0,L) σℓ[ℓ+1] += Pℓμ[ℓ+1] * σs end
    else
        error("Unknown type")
    end

    return σℓ
end

"""
    tcs(this::Inelastic_Leptons,Ei::Float64,Ec::Float64,particle::Particle,Z::Int64)

Gives the total cross-section for inelastic lepton interaction. 

# Input Argument(s)
- `this::Inelastic_Leptons` : inelastic lepton structure. 
- `Ei::Float64` : incoming particle energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.
- `Z::Int64` : atomic number.

# Output Argument(s)
- `σt::Float64` : total cross-section.

"""
function tcs(this::Inelastic_Leptons,Ei::Float64,Ec::Float64,particle::Particle,Z::Int64)

    # Inititalisation
    rₑ = 2.81794092E-13       # (in cm)
    γ = Ei+1
    β² = Ei*(Ei+2)/(Ei+1)^2
    Nshells,Zi,Ui,Ti,ri,subshells = electron_subshells(Z)
    σt = 0.0

    # Close collisions
    if is_electron(particle)
        for δi in range(1,Nshells)
            Wmax = (Ei-Ui[δi])/2
            Wmin = Ei-Ec
            if Wmax > Wmin
                J₀⁻(x) = -1/(x+Ui[δi]) + 1/(Ei-x) + x/(Ei+1)^2 - (2*Ei+1)/(Ei+1)^2 * (log(x+Ui[δi])-log(Ei-x))/(Ei+Ui[δi])
                σt += 2*π*rₑ^2/β² * Zi[δi] * (J₀⁻(Wmax)-J₀⁻(Wmin))
            end
        end
    elseif is_positron(particle)
        b = ((γ-1)/γ)^2
        b1 = b * (2*(γ+1)^2-1)/(γ^2-1)
        b2 = b * (3*(γ+1)^2+1)/(γ+1)^2
        b3 = b * (2*(γ-1)*γ)/(γ+1)^2
        b4 = b * (γ-1)^2/(γ+1)^2
        for δi in range(1,Nshells)
            Wmax = Ei-Ui[δi]
            Wmin = Ei-Ec
            if Wmax > Wmin
                J₀⁺(x) = -1/(x+Ui[δi]) - b1*log(x+Ui[δi])/Ei + b2*x/Ei^2 - b3*(x^2/2+Ui[δi]*x)/Ei^3 + b4*(x^3/3+Ui[δi]*x^2+Ui[δi]^2*x)/Ei^4
                σt += 2*π*rₑ^2/β² * Zi[δi] * (J₀⁺(Wmax)-J₀⁺(Wmin))
            end
        end
    else
        error("Unknown particle")
    end
    return σt
end

"""
    sp(this::Inelastic_Leptons,Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,
    state_of_matter::String,Ei::Float64,Ec::Float64,particle::Particle)

Gives the stopping power for inelastic lepton interaction.

# Input Argument(s)
- `this::Inelastic_Leptons` : inelastic lepton structure. 
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `ωz::Vector{Float64}` : weight fraction of the elements composing the material.
- `ρ::Float64` : material density.
- `state_of_matter::String` : state of matter.
- `Ei::Float64` : incoming particle energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.

# Output Argument(s)
- `S::Float64` : stopping power.

"""
function sp(this::Inelastic_Leptons,Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,state_of_matter::String,Ei::Float64,Ec::Float64,particle::Particle)

    # Initialization
    rₑ = 2.81794092e-13 # (in cm)
    γ = Ei+1
    β² = Ei*(Ei+2)/(Ei+1)^2

    # Compute the total cross-section
    Stot = bethe(Z,ωz,ρ,Ei,particle)
    
    # Compute the catastrophic Møller- or Bhabha- derived stopping power
    Sc = 0
    Nz = length(Z)
    for i in range(1,Nz)
        Nshells,Zi,Ui,Ti,ri,subshells = electron_subshells(Z[i])
        if is_electron(particle)
            for δi in range(1,Nshells)
                Wmax = (Ei-Ui[δi])/2
                Wmin = Ei-Ec
                if Wmax > Wmin
                    J₁⁻(x) = log(x+Ui[δi]) + log(Ei-x) + (Ei+Ui[δi])/(Ei-x) + x*(x+2*Ui[δi])/(2*(Ei+1)^2) + (2*Ei+1)/(Ei+1)^2 * log(Ei-x)
                    Sc += 2*π*rₑ^2/β² * ωz[i] * nuclei_density(Z[i],ρ) * Zi[δi] * (J₁⁻(Wmax)-J₁⁻(Wmin))
                end
            end
        elseif is_positron(particle)
            b = ((γ-1)/γ)^2
            b1 = b * (2*(γ+1)^2-1)/(γ^2-1)
            b2 = b * (3*(γ+1)^2+1)/(γ+1)^2
            b3 = b * (2*(γ-1)*γ)/(γ+1)^2
            b4 = b * (γ-1)^2/(γ+1)^2
            for δi in range(1,Nshells)
                Wmax = Ei-Ui[δi]
                Wmin = Ei-Ec
                if Wmax > Wmin
                    J₁⁺(x) = log(x+Ui[δi]) - b1*x/Ei + b2*(x^2/2+Ui[δi]*x)/Ei^2 - b3*(x^3/3+Ui[δi]*x^2+Ui[δi]^2*x)/Ei^3 + b4*(x^4/4+Ui[δi]*x^3+3*Ui[δi]^2*x^2/2+Ui[δi]^3*x)/Ei^4
                    Sc += 2*π*rₑ^2/β² * ωz[i] * nuclei_density(Z[i],ρ) * Zi[δi] * (J₁⁺(Wmax)-J₁⁺(Wmin))
                end
            end
        end
    end

    # Compute the soft stopping power
    Sr = Stot - Sc

    return Sr
end

"""
    mt(this::Inelastic_Leptons)

Gives the momentum transfer for inelastic lepton interaction.

# Input Argument(s)
- `this::Inelastic_Leptons` : inelastic lepton structure. 

# Output Argument(s)
- `T::Float64` : momentum transfer.

"""
function mt(this::Inelastic_Leptons)
    T = 0.0
    return T
end

