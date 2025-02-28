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
    Cℓk::Array{Float64}
    Cℓki::Array{Float64}
    is_CSD::Bool
    is_AFP::Bool
    is_AFP_decomposition::Bool
    is_elastic::Bool
    is_preload_data::Bool
    is_subshells_dependant::Bool
    is_waller_hartree_factor::Bool
    incoherent_scattering_factor::Function
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
        this.is_preload_data = true
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

# Output Argument(s)
- `σℓ::Vector{Float64}` : Legendre moments of the scattering cross-sections.

"""
function dcs(this::Compton,L::Int64,Ei::Float64,Ef::Float64,type::String,Z::Int64,iz::Int64,δi::Int64)

    if this.model ∈ ["klein-nishina","waller-hartree"]

        # Initialization
        rₑ = 2.81794092e-13 # (in cm)
        σs = 0
        σℓ = zeros(L+1)

        # Scattered photon
        if type == "S"
            # Compute the differential scattering cross section
            if Ei/(1+2*Ei) ≤ Ef ≤ Ei
                σs = π * rₑ^2 / Ei^2 * (Ei/Ef + Ef/Ei - 2*(1/Ef-1/Ei) + (1/Ef-1/Ei)^2)
            end
            # Compute the Legendre moments of the flux
            if σs != 0
                μ = max(min(1 + 1/Ei - 1/Ef,1),-1)
                if this.model == "waller-hartree"
                    S = this.incoherent_scattering_factor(iz,Ei,μ)
                else
                    S = Z
                end
                Pℓμ = legendre_polynomials(L,μ)
                for ℓ in range(0,L) σℓ[ℓ+1] += Pℓμ[ℓ+1] * σs * S end
            end
        # Produced electron
        elseif type == "P"
            # Change of variable
            Eₑ = copy(Ef)
            Ef = Ei - Eₑ
            # Compute the differential scattering cross section
            if 0 ≤ Eₑ ≤ 2*Ei^2/(1+2*Ei)
                σs = π * rₑ^2 / Ei^2 * (Ei/Ef + Ef/Ei - 2*(1/Ef-1/Ei) + (1/Ef-1/Ei)^2)
            end
            # Compute the Legendre moments of the flux
            if σs != 0
                μ = max(min((1 + Ei)/Ei * 1/sqrt(2/Eₑ+1),1),-1)
                μγ = max(min(1 + 1/Ei - 1/Ef,1),-1)
                if this.model == "waller-hartree"
                    S = this.incoherent_scattering_factor(iz,Ei,μγ)
                else
                    S = Z
                end
                Pℓμ = legendre_polynomials(L,μ)
                for ℓ in range(0,L) σℓ[ℓ+1] += Pℓμ[ℓ+1] * σs * S end
            end
        else
            error("Unknown interaction.")
        end
        return σℓ

    elseif this.model == "impulse_approximation"

        mₑc² = 0.510999
        σs = 0
        σℓ = zeros(L+1)
        J₀i = orbital_compton_profiles(Z)[δi]
        _,Zi,Ui,_,_,_ = electron_subshells(Z)
        mₑ = 1                        # (a₀)
        c = 137.03599908388762        # (a₀×Eₕ/ħ)
        rₑ = 5.325135459237564e-5     # (mₑ)
        a₀ = 5.29177210903e-11        # (a₀)

        # Conversion to atomic units (MeV) -> (Eₕ)
        Ei = Ei * 1e6 / 27.211386245988 * mₑc²
        Ef = Ef * 1e6 / 27.211386245988 * mₑc²
        Ui .= Ui * 1e6 ./ 27.211386245988 * mₑc²

        # Scattered photon
        if type == "S"
            Nμ = 80
            μ,w = quadrature(Nμ,"gauss-lobatto")
            if Ei - Ef - Ui[δi] ≥ 0
                for n in range(1,Nμ)
                    Pℓμ = legendre_polynomials(L,μ[n])
                    Ec = Ei*mₑ*c^2/(mₑ*c^2+Ei*(1-μ[n]))
                    pz = (Ei*Ef*(1-μ[n])-mₑ*c^2*(Ei-Ef))/(c*sqrt(Ei^2+Ef^2-2*Ei*Ef*μ[n]))
                    Ji = J₀i*(1+2*J₀i*abs(pz))*exp(1/2-1/2*(1+2*J₀i*abs(pz))^2)
                    σs = w[n] * π * rₑ^2 * Ef/Ei * (Ec/Ei+Ei/Ec+μ[n]^2-1) * 1/sqrt(Ei^2+Ef^2-2*Ei*Ef*μ[n]+Ei^2*(Ef-Ec)^2/Ec^2) * Zi[δi]*Ji *mₑ*c
                    for ℓ in range(0,L) σℓ[ℓ+1] += Pℓμ[ℓ+1] * σs * (a₀ ^ 2) / 27.211386245988 * 100^2 * 1e6 * mₑc² end
                end
            end
        # Produced electron
        elseif type == "P"
            Eₑ = copy(Ef)
            Ef = (Ei - Ui[δi]) - Eₑ
            Nμ = 80
            μ,w = quadrature(Nμ,"gauss-lobatto")
            if Ei - Ef - Ui[δi] ≥ 0
                for n in range(1,Nμ)
                    Pℓμ = legendre_polynomials(L,μ[n])
                    Ec = Ei*mₑ*c^2/(mₑ*c^2+Ei*(1-μ[n]))
                    pz = (Ei*Ef*(1-μ[n])-mₑ*c^2*(Ei-Ef))/(c*sqrt(Ei^2+Ef^2-2*Ei*Ef*μ[n]))
                    Ji = J₀i*(1+2*J₀i*abs(pz))*exp(1/2-1/2*(1+2*J₀i*abs(pz))^2)
                    σs = w[n] * π * rₑ^2 * Ef/Ei * (Ec/Ei+Ei/Ec+μ[n]^2-1) * 1/sqrt(Ei^2+Ef^2-2*Ei*Ef*μ[n]+Ei^2*(Ef-Ec)^2/Ec^2) * Zi[δi]*Ji *mₑ*c
                    for ℓ in range(0,L) σℓ[ℓ+1] += Pℓμ[ℓ+1] * σs * (a₀ ^ 2) / 27.211386245988 * 100^2 * 1e6 * mₑc² end
                end
            end
        else
            error("Unknown interaction.")
        end
        return σℓ

    else
        error("Unknown Compton model.")
    end
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
function tcs(this::Compton,Ei::Float64,Z::Int64,Eout::Vector{Float64},iz::Int64)

    if this.model ∈ ["klein-nishina","waller-hartree"]
        σt = 0.0
        rₑ = 2.81794092E-13 # (in cm)
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
                σs = 0.0
                if Ei/(1+2*Ei) ≤ Ef ≤ Ei
                    μ = max(min(1 + 1/Ei - 1/Ef,1),-1)
                    if this.model == "waller-hartree"
                        S = this.incoherent_scattering_factor(iz,Ei,μ)
                    else
                        S = Z
                    end
                    σs = S * π * rₑ^2 / Ei^2 * (Ei/Ef + Ef/Ei - 2*(1/Ef-1/Ei) + (1/Ef-1/Ei)^2)
                end
                
                # Cross-sections
                σt += ΔEf/2 * w[n] * σs
            end
        end
        return σt

    elseif this.model == "impulse_approximation"

        mₑc² = 0.510999
        mₑ = 1                        # (a₀)
        c = 137.03599908388762        # (a₀×Eₕ/ħ)
        rₑ = 5.325135459237564e-5     # (mₑ)
        a₀ = 5.29177210903e-11        # (a₀)
        Ei = Ei * 1e6 / 27.211386245988 * mₑc²
        Eout2 = Eout .* (1e6 / 27.211386245988 * mₑc²)
        Nμ = 80
        μ,w2 = quadrature(Nμ,"gauss-lobatto")
        Nshells,Zi,Ui,Ti,ri,subshells = electron_subshells(Z)
        Ui .= Ui * 1e6 ./ 27.211386245988 * mₑc²
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
                    σs = 0.0
                    J₀i = orbital_compton_profiles(Z)[δi]
                    if Ei - Ef - Ui[δi] ≥ 0
                        for n2 in range(1,Nμ)
                            Ec = Ei*mₑ*c^2/(mₑ*c^2+Ei*(1-μ[n2]))
                            pz = (Ei*Ef*(1-μ[n2])-mₑ*c^2*(Ei-Ef))/(c*sqrt(Ei^2+Ef^2-2*Ei*Ef*μ[n2]))
                            Ji = J₀i*(1+2*J₀i*abs(pz))*exp(1/2-1/2*(1+2*J₀i*abs(pz))^2)
                            σs += w2[n2] * π * rₑ^2 * Ef/Ei * (Ec/Ei+Ei/Ec+μ[n2]^2-1) * 1/sqrt(Ei^2+Ef^2-2*Ei*Ef*μ[n2]+Ei^2*(Ef-Ec)^2/Ec^2) * Zi[δi]*Ji *mₑ*c
                        end
                    end
                    # Cross-sections
                    σt += ΔEf/2 * w[n] * σs * (a₀ ^ 2) * 100^2
                end
            end
        end
        return σt

    else
        error("Unknown Compton model.")
    end
end

"""
    preload_data(this::Annihilation,Z::Vector{Int64},Emax::Float64,Emin::Float64,L::Int64,
    type::String,Eout::Vector{Float64},interactions::Vector{Interaction})

Preload data for multigroup annihilation calculations. 

# Input Argument(s)
- `this::Annihilation` : annihilation structure.  
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.

# Output Argument(s)
N/A

"""
function preload_data(this::Compton,Z::Vector{Int64})
    
    # Incoherent scattering factor
    if this.model == "waller-hartree"
        Nz = length(Z)
        x = Vector{Vector{Float64}}(undef,Nz)
        S = Vector{Vector{Float64}}(undef,Nz)
        S_spline = Vector{Function}(undef,Nz)
        for iz in range(1,Nz)
            path = joinpath(find_package_root(), "data", "compton_factors_JENDL5.jld2")
            data = load(path)
            x[iz] = data["x"][Z[iz]]
            S[iz] = data["F"][Z[iz]]
            S_spline[iz] = cubic_hermite_spline(x[iz],S[iz])
        end
        this.incoherent_scattering_factor = function incoherent_scattering_factor(iz::Int64,Ei::Float64,μ::Float64)
            hc = 1/20.60744 # (hc in mₑc² × Å)
            xi = 2*Ei/(hc)*sqrt((1-μ)/2) * sqrt((1+(Ei^2+2*Ei)*(1-μ)/2))/(1+Ei*(1-μ))
            Si = S_spline[iz](xi)
            return Si
        end
    end
    
end
