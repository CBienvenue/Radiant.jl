"""
    Compton

Structure used to define parameters for production of multigroup Compton cross-sections.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}} = Dict((Photon,Photon) => ["S"],(Photon,Electron) => ["P"])` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `(Photon,Photon) => ["S"]` : scattering of incident photon following Compton interaction.
    - `(Photon,Electron) => ["P"]` : produced electron following Compton interaction.

"""
mutable struct Compton <: Interaction

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
    is_waller_hartree_factor::Bool
    scattering_model::String
    model::String

    # Constructor(s)
    function Compton()
        this = new()
        this.name = "compton"
        this.interaction_types = Dict((Photon,Photon) => ["S"],(Photon,Electron) => ["P"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.is_CSD = false
        this.is_AFP = false
        this.is_AFP_decomposition = false
        this.is_elastic = false
        this.is_subshells_dependant = false
        this.model = "waller-hartree"
        this.scattering_model = "BTE"
        return this
    end
end

# Method(s)
"""
    set_interaction_types(this::Compton,interaction_types::Dict{Tuple{DataType,DataType},Vector{String}})

To define the interaction types for Compton processes.

# Input Argument(s)
- `this::Compton` : compton structure.
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}}` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `(Photon,Photon)   => ["S"]` : scattering of incident photon following Compton interaction.
    - `(Photon,Electron) => ["P"]` : produced electron following Compton interaction.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> compton = Compton()
julia> compton.set_interaction_types( Dict((Photon,Photon) => ["S"]) ) # Electron are absorbed following Compton interaction.
```
"""
function set_interaction_types(this::Compton,interaction_types)
    this.interaction_types = interaction_types
end

"""
    set_model(this::Compton,model::String)

To define the Compton physics model.

# Input Argument(s)
- `this::Compton` : compton structure.
- `model::String` : Compton physics interaction model, which can be:
    - `klein-nishina`  : Klein-Nishina model.
    - `waller-hartree` : Klein-Nishina model, with Waller-Hartree incoherent scattering function.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> compton = Compton()
julia> compton.set_model("klein-nishina")
```
"""
function set_model(this::Compton,model::String)
    if lowercase(model) ∉ ["klein-nishina","waller-hartree"] error("Unknown Compton model.") end
    this.model = lowercase(model)
end

"""
    in_distribution(this::Compton)

Describe the energy discretization method for the incoming particle in the Compton
interaction.

# Input Argument(s)
- `this::Compton` : Compton structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function in_distribution(this::Compton)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    out_distribution(this::Compton)

Describe the energy discretization method for the outgoing particle in the Compton
interaction.

# Input Argument(s)
- `this::Compton` : Compton structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function out_distribution(this::Compton)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    bounds(this::Compton,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String)

Gives the integration energy bounds for the outgoing particle for Compton interaction. 

# Input Argument(s)
- `this::Compton` : Compton structure.
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `Ei::Float64` : energy of the incoming particle.
- `type::String` : type of interaction.

# Output Argument(s)
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `isSkip::Bool` : define if the integration is skipped or not.

"""
function bounds(this::Compton,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String)
    # Scattered photon
    if type == "S" 
        Ef⁻ = min(Ei,Ef⁻)
        if this.model != "impulse_approximation"
            Ef⁺ = max(Ei/(1+2*Ei),Ef⁺)
        end
        if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    # Produced electron
    elseif type == "P" 
        if this.model != "impulse_approximation"
            Ef⁻ = min(2*Ei^2/(1+2*Ei),Ei,Ef⁻)
        else
            Ef⁻ = min(Ei,Ef⁻)
        end
        if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    else
        error("Unknown type of method for Klein-Nishina scattering.")
    end
    return Ef⁻,Ef⁺,isSkip
end

"""
    dcs(this::Compton,L::Int64,Ei::Float64,Ef::Float64,type::String,Z::Int64,iz::Int64,
    δi::Int64)

Gives the Legendre moments of the scattering cross-sections for Compton interaction. 

# Input Argument(s)
- `this::Compton` : Compton structure.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : incoming particle energy.
- `Ef::Float64` : outgoing particle energy.
- `type::String` : type of interaction.
- `Z::Int64` : atomic number.
- `iz::Int64` : index of the element.
- `δi::Int64` : subshell index.
- `particle::Particle` : outgoing particle.

# Output Argument(s)
- `σℓ::Vector{Float64}` : Legendre moments of the scattering cross-sections.

"""
function dcs(this::Compton,L::Int64,Ei::Float64,Ef::Float64,type::String,Z::Int64,iz::Int64,δi::Int64,particle::Particle)

    if this.model ∈ ["klein-nishina","waller-hartree"]

        #----
        # Initialization
        #----
        σℓ = zeros(L+1)

        #----
        # Scattered photon
        #----
        if type == "S"

            # Compute the differential scattering cross section
            if this.model == "klein-nishina"
                σs = Z * klein_nishina(Ei,Ef)
            elseif this.model == "waller-hartree"
                σs = waller_hartree(Z,Ei,Ef)
            end

            # Compute the Legendre moments of the flux
            Wℓ = angular_klein_nishina(Ei,Ef,L,particle)

            # Compute the Legendre moments of the cross-section
            for ℓ in range(0,L) σℓ[ℓ+1] = σs * Wℓ[ℓ+1] end

        #----
        # Produced electron
        #----
        elseif type == "P"
            
            # Compute the differential scattering cross section
            if this.model == "klein-nishina"
                σs = Z * klein_nishina(Ei,Ei-Ef)
            elseif this.model == "waller-hartree"
                σs = waller_hartree(Z,Ei,Ei-Ef)
            end

            # Compute the Legendre moments of the flux
            Wℓ = angular_klein_nishina(Ei,Ef,L,particle)

            # Compute the Legendre moments of the cross-section
            for ℓ in range(0,L) σℓ[ℓ+1] = σs * Wℓ[ℓ+1] end

        else
            error("Unknown interaction.")
        end

    elseif this.model == "impulse_approximation"

        #----
        # Scattered photon
        #----
        if type == "S"
            σℓ = impulse_approximation(Z,L,δi,Ei,Ef)

        #----
        # Produced electron
        #----
        elseif type == "P"
            _,_,Ui,_,_,_ = electron_subshells(Z)
            Eγ = (Ei-Ui[δi])-Ef
            σℓ = impulse_approximation(Z,L,δi,Ei,Eγ)
        else
            error("Unknown interaction.")
        end
    else
        error("Unknown Compton model.")
    end
    return σℓ
end

"""
    tcs(this::Compton,Ei::Float64,Z::Int64,Eout::Vector{Float64},iz::Int64)

Gives the total cross-section for Compton. 

# Input Argument(s)
- `this::Compton` : Compton structure.
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.
- `Eout::Vector{Float64}` : energy boundaries associated with the outgoing particle.
- `iz::Int64` : element index.

# Output Argument(s)
- `σt::Float64` : total cross-section.

"""
function tcs(this::Compton,Ei::Float64,Z::Int64,Eout::Vector{Float64})

    if this.model ∈ ["klein-nishina","waller-hartree"]
        σt = 0.0
        Ngf = length(Eout)-1
        is_dirac, Np, q_type = out_distribution(this)
        if is_dirac Np = 1; u = [0]; w = [2] else u,w = quadrature(Np,q_type) end
        for gf in range(1,Ngf+1)
            Ef⁻ = Eout[gf]
            if (gf != Ngf+1) Ef⁺ = Eout[gf+1] else Ef⁺ = 0.0 end
            Ef⁻,Ef⁺,isSkip = bounds(this,Ef⁻,Ef⁺,Ei,"S")
            if isSkip continue end
            ΔEf = Ef⁻ - Ef⁺
            for n in range(1,Np)
                Ef = (u[n]*ΔEf + (Ef⁻+Ef⁺))/2

                # Compute the differential scattering cross section
                if this.model == "klein-nishina"
                    σs = Z * klein_nishina(Ei,Ef)
                elseif this.model == "waller-hartree"
                    σs = waller_hartree(Z,Ei,Ef)
                end
                σt += ΔEf/2 * w[n] * σs
            end
        end
        return σt

    elseif this.model == "impulse_approximation"

        mₑc² = 0.510999
        Ei = Ei * 1e6 / 27.211386245988 * mₑc²
        Eout2 = Eout .* (1e6 / 27.211386245988 * mₑc²)
        Nshells,_,_,_,_,_ = electron_subshells(Z)
        σt = 0.0
        Ngf = length(Eout2)-1
        is_dirac, Np, q_type = out_distribution(this)
        if is_dirac Np = 1; u = [0]; w = [2] else u,w = quadrature(Np,q_type) end
        for gf in range(1,Ngf+1)
            Ef⁻ = Eout2[gf]
            if (gf != Ngf+1) Ef⁺ = Eout2[gf+1] else Ef⁺ = 0.0 end
            for δi in range(1,Nshells)
                Ef⁻,Ef⁺,isSkip = bounds(this,Ef⁻,Ef⁺,Ei,"S")
                if isSkip continue end
                ΔEf = Ef⁻ - Ef⁺
                for n in range(1,Np)
                    Ef = (u[n]*ΔEf + (Ef⁻+Ef⁺))/2
                    σs = impulse_approximation(Z,0,δi,Ei,Ef)[1]
                    σt += ΔEf/2 * w[n] * σs
                end
            end
        end
        return σt

    else
        error("Unknown Compton model.")
    end
end

"""
    acs(this::Compton,Ei::Float64,Z::Int64,Ecutoff::Float64)

Gives the absorption cross-section for Compton. 

# Input Argument(s)
- `this::Compton` : Compton structure.
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.
- `Ecutoff::Float64` : cutoff energy.

# Output Argument(s)
- `σa::Float64` : absorption cross-section.

"""
function acs(this::Compton,Ei::Float64,Z::Int64,Ecutoff::Float64)
    return  tcs(this,Ei,Z,[Ecutoff])
end