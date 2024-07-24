"""
    Elastic_Leptons

Structure used to define parameters for production of multigroup elastic cross-sections for leptons.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{String,String},Vector{String}} = Dict(("electrons","electrons") => ["S"],("positrons","positrons") => ["S"])`: Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `("electrons","electrons") => ["S"]`: elastic interaction of electrons.
    - `("positrons","positrons") => ["S"]`: elastic interaction of positrons.

"""
mutable struct Elastic_Leptons <: Interaction

    # Variable(s)
    name::String
    interaction_types::Dict{Tuple{String,String},Vector{String}}
    incoming_particle::Vector{String}
    interaction_particles::Vector{String}
    is_CSD::Bool
    is_elastic::Bool
    is_ETC::Bool
    is_AFP::Bool
    is_preload_data::Bool
    is_subshells_dependant::Bool
    is_kawrakow_correction::Bool
    subshell_dependant_inelastic::Bool
    model::String
    plasma_energy::Float64
    effective_mean_excitation_energy::Float64
    bjk_boschini::Vector{Array{Float64}}
    Cℓk::Array{Float64}
    Cℓki::Array{Float64}
    scattering_model::String

    # Constructor(s)
    function Elastic_Leptons(;
        ### Initial values ###
        model="mott",
        is_kawrakow_correction=true,
        subshell_dependant_inelastic=true,
        is_ETC=true,
        is_AFP=true,
        interaction_types = Dict(("electrons","electrons") => ["S"],("positrons","positrons") => ["S"])
        ######################
        )
        this = new()
        this.name = "elastic_leptons"
        this.set_interaction_types(interaction_types)
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.is_CSD = false
        this.is_elastic = true
        this.set_transport_correction(is_ETC)
        this.set_angular_fokker_planck(is_AFP)
        this.is_preload_data = true
        this.is_subshells_dependant = false
        this.set_kawrakow_correction(is_kawrakow_correction,subshell_dependant_inelastic)
        this.set_model(model)
        this.scattering_model = "BFP"
        return this
    end

end

# Method(s)
"""
    set_interaction_types(this::Elastic_Leptons,interaction_types::Dict{Tuple{String,String},Vector{String}})

To define the interaction types for Elastic_Leptons processes.

# Input Argument(s)
- `this::Elastic_Leptons`: elastic leptons structure.
- `interaction_types::Dict{Tuple{String,String},Vector{String}}`: Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `("electrons","electrons") => ["S"]`: elastic interaction of electrons.
    - `("positrons","positrons") => ["S"]`: elastic interaction of positrons.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_leptons = Elastic_Leptons()
julia> elastic_leptons.set_interaction_types( Dict(("positrons","positrons") => ["S"]) ) # Elastic only for positrons
```
"""
function set_interaction_types(this::Elastic_Leptons,interaction_types::Dict{Tuple{String,String},Vector{String}})
    this.interaction_types = interaction_types
end

"""
    set_model(this::Elastic_Leptons,model::String)

To define the elastic scattering model.

# Input Argument(s)
- `this::Elastic_Leptons`: elastic leptons structure.
- `model::String`: model of elastic scattering:
    - `rutherford`: screened Rutherford cross-sections.
    - `mott`: screened Mott cross-sections.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_leptons = Elastic_Leptons()
julia> elastic_leptons.set_model("rutherford")
```
"""
function set_model(this::Elastic_Leptons,model::String)
    if lowercase(model) ∉ ["rutherford","mott"] error("Unkown elastic model: '$model'.") end
    this.model = lowercase(model)
end

"""
    set_kawrakow_correction(this::Elastic_Leptons,is_kawrakow_correction::Bool)

Apply Kawrakow's correction to elastic cross-sections.

# Input Argument(s)
- `this::Elastic_Leptons`: elastic leptons structure.
- `is_kawrakow_correction::Bool`: Apply Karakow's correction (true) or not (false).

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_leptons = Elastic_Leptons()
julia> elastic_leptons.set_kawrakow_correction(false)
```
"""
function set_kawrakow_correction(this::Elastic_Leptons,is_kawrakow_correction::Bool,subshell_dependant_inelastic::Bool=true)
    this.is_kawrakow_correction = is_kawrakow_correction
    this.subshell_dependant_inelastic = subshell_dependant_inelastic
end

"""
    set_transport_correction(this::Elastic_Leptons,is_ETC::Bool)

Enable or not extended transport correcton.

# Input Argument(s)
- `this::Elastic_Leptons`: elastic leptons structure.
- `is_ETC::Bool`: Enable (true) or not (false) extended transport correcton.

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
    set_angular_fokker_planck(this::Elastic_Leptons,is_AFP::Bool)

Enable or not the extraction of the angular Fokker-Planck.

# Input Argument(s)
- `this::Elastic_Leptons`: elastic leptons structure.
- `is_AFP::Bool`: Enable (true) or not (false) t the extraction of the angular Fokker-Planck.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_leptons = Elastic_Leptons()
julia> elastic_leptons.is_AFP(false)
```
"""
function set_angular_fokker_planck(this::Elastic_Leptons,is_AFP::Bool)
    this.is_AFP = is_AFP
end

function in_distribution(this::Elastic_Leptons)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

function out_distribution(this::Elastic_Leptons)
    is_dirac = true
    N = 1
    quadrature = "dirac"
    return is_dirac, N, quadrature
end

function bounds(this::Elastic_Leptons,Ef⁻::Float64,Ef⁺::Float64,gi::Int64,gf::Int64)
    if (gf != gi) isSkip = true else isSkip = false end # Elastic scattering only
    return Ef⁻,Ef⁺,isSkip
end

function dcs(this::Elastic_Leptons,L::Int64,Ei::Float64,Z::Int64,particle::String,Ecutoff::Float64,iz::Int64)

    # Initialization
    rₑ = 2.81794092E-13 # (in cm)
    β² = Ei*(Ei+2)/(Ei+1)^2
    β = sqrt(β²)
    β₀ = 0.7181287
    α = 1/137

    σℓ = zeros(L+1)
    η = Z^(2/3) * α^2 * (1.13 + 3.76 * (Z*α)^2/β² ) / (4 * (9*π^2/128)^(2/3) * Ei * (Ei+2) )
    ai = zeros(5)
    if this.model == "mott"
        b = this.bjk_boschini[iz]
        for j in range(0,4), k in range(1,6)
            ai[j+1] += b[j+1,k] * (β-β₀)^(k-1)
        end
    elseif this.model == "rutherford"
        ai[1] = 1
    else
        error("Unknown elastic formula.")
    end
    αi = [ai[1]+ai[3]+ai[5],-(ai[3]+2*ai[5]),ai[5],ai[2]+ai[4],-ai[4]]

    # Integral pre-calculations
    μmax = 1; μmin = -1
    𝒢₁⁺ = 𝒢₁(μmax,2+L,1+2*η,-1)
    𝒢₁⁻ = 𝒢₁(μmin,2+L,1+2*η,-1)
    𝒢₂⁺ = 𝒢₂(sqrt(1-μmin),2*(2+L),2*η,1)
    𝒢₂⁻ = 𝒢₂(sqrt(1-μmax),2*(2+L),2*η,1)

    @inbounds for ℓ in range(0,L)
        for k in range(0,div(ℓ,2))
            σℓk = 0.0
            # Compute I₁ -----
            I₁ = zeros(3)
            for i in range(0,2)
                I₁[i+1] += 𝒢₁⁺[ℓ-2*k+i+1] - 𝒢₁⁻[ℓ-2*k+i+1]
                I₁[i+1] *= αi[i+1]
            end 
            σℓk += sum(I₁)
            # Compute I₂ -----
            I₂ = zeros(2)
            for i in range(0,1)
                for g in range(0,ℓ-2*k+i)
                    I₂[i+1] += this.Cℓki[ℓ+1,k+1,i+1,g+1] * ( 𝒢₂⁺[2*(1+g)] - 𝒢₂⁻[2*(1+g)])
                end
                I₂[i+1] *= αi[i+4]
            end
            σℓk += sum(I₂)
            σℓ[ℓ+1] += this.Cℓk[ℓ+1,k+1] * σℓk
        end
        σℓ[ℓ+1] *= 1/(2^ℓ)
    end

    # Kawrakow correction
    if this.is_kawrakow_correction
        if this.subshell_dependant_inelastic
            Nshells,Zi,Ui,Ti,_,_ = electron_subshells(Z)
        else
            Ui = [0.0]; Nshells = 1; Zi = Z[i]; Ti = [0.0]; # Free atomic electron
        end
        gM = 0
        for δi in range(1,Nshells)
            Wc = Ecutoff # Knock-on production cutoff with Møller or Bhabha
            if particle == "electrons"
                Wmax = (Ei-Ui[δi])/2
                if Wc < Wmax
                    Jm₁(x) = ((Ei+2)*((Ei+2)*(x+Ui[δi])*log(x+Ui[δi])-(Ei+2)*(x+Ui[δi])*log(Ei-x+2)+Ui[δi]^2+(Ei+2)*Ui[δi]))/((Ui[δi]+Ei+2)^2*(x+Ui[δi]))
                    Jm₂(x) = ((Ei+2)^2*log(Ei-x)-(Ei+2)^2*log(Ei-x+2))/4-(Ei*(Ei+2))/(2*(x-Ei))
                    Jm₃(x) = -((Ei+2)*((Ei+2)*log(Ei-x+2)+x))/(Ei+1)^2
                    Jm₄(x) = ((Ei+2)*(2*Ei+1)*(2*Ui[δi]*log(x+Ui[δi])+Ei*(Ui[δi]+Ei+2)*log(Ei-x)-(Ei+2)*(Ui[δi]+Ei)*log(Ei-x+2)))/(2*(Ei+1)^2*(Ui[δi]+Ei)*(Ui[δi]+Ei+2))
                    Jm(x) = Jm₁(x) + Jm₂(x) + Jm₃(x) + Jm₄(x)
                    gM += Zi[δi]/Z * (Jm(Wmax) - Jm(Wc))
                end
            elseif particle == "positrons"
                γ = Ei+1
                b = ((γ-1)/γ)^2
                b1 = b * (2*(γ+1)^2-1)/(γ^2-1)
                b2 = b * (3*(γ+1)^2+1)/(γ+1)^2
                b3 = b * (2*(γ-1)*γ)/(γ+1)^2
                b4 = b * (γ-1)^2/(γ+1)^2
                Wmax = Ei-Ui[δi]
                if Wc < Wmax
                    Jb₁(x) = log(x/(Ei-x+2))
                    Jb₂(x) = (Ei+2)/Ei * (b1-(Ei+2)/Ei*b2+(Ei+2)^2/Ei^2*b3-(Ei+2)^3/Ei^3*b4)*log(Ei-x+2)
                    Jb₃(x) = -(Ei+2)/Ei^2 * (b2-(Ei+2)/Ei*b3+(Ei+2)^2/Ei^2*b4)*x
                    Jb₄(x) = (Ei+2)/(2*Ei^3) * (b3-(Ei+2)/Ei*b4)*x^2
                    Jb₅(x) = -(Ei+2)/(3*Ei^4) * b4*x^3
                    Jb(x) = Jb₁(x) + Jb₂(x) + Jb₃(x) + Jb₄(x) + Jb₅(x)
                    gM += Zi[δi]/Z * (Jb(Wmax) - Jb(Wc))
                end
            else
                error("Unkown particle.")
            end

        end
        if this.model == "mott"
            gR = 1/3 * ((30*(16*ai[5]*η^3+12*ai[5]*η^2-6*ai[3]*η^2-4*ai[3]*η+2*ai[1]*η+ai[1])*log(2*(η+1))-30*(16*ai[5]*η^3+12*ai[5]*η^2-6*ai[3]*η^2-4*ai[3]*η+2*ai[1]*η+ai[1])*log(2*η)+15*2^(3/2)*atan(1/sqrt(η))*sqrt(η)*(14*ai[4]*η^2+10*ai[4]*η-5*ai[2]*η-3*ai[2])-15*(32*ai[5]+7*2^(5/2)*ai[4])*η^2-5*(24*ai[5]+2^(11/2)*ai[4]-36*ai[3]-15*2^(3/2)*ai[2])*η+20*ai[5]+2^(9/2)*ai[4]+30*ai[3]+5*2^(7/2)*ai[2]-60*ai[1])/10)
        elseif this.model == "rutherford"
            gR = (1+2*η)*log(1+1/η)-2
        end
        ξ = 1 - gM/gR
    else
        ξ = 0
    end

    Γ = 2*π*rₑ^2*Z*(Z+ξ)/(β²*Ei*(Ei+2))
    σℓ .*= Γ

    # Correction to deal with high-order Legendre moments
    for ℓ in range(1,L)
        if abs(σℓ[1]) < abs(σℓ[ℓ+1])
            σℓ[ℓ+1:end] .= 0.0
            break
        end
    end

    return σℓ
end

function tcs(this::Elastic_Leptons,Ei::Float64,Z::Int64,particle::String,Ecutoff::Float64,iz::Int64)
    σt = dcs(this,0,Ei,Z,particle,Ecutoff,iz)[1]
    return σt
end

function preload_data(this::Elastic_Leptons,Z::Vector{Int64},L::Int64,particle::String)

    # Load Boschini data for Mott cross-sections
    if this.model == "mott"
        path = joinpath(find_package_root(), "data", "mott_data_boschini_2013.jld2")
        data = load(path)
        if ~haskey(data,particle) error("Mott cross-sections are only available for electrons and positrons.") end
        Nz = length(Z)
        this.bjk_boschini = Vector{Array{Float64}}(undef,Nz)
        for iz in range(1,Nz)
            this.bjk_boschini[iz] = data[particle][Z[iz]]
        end
    end

    # Precompute angular integration factors
    this.Cℓk = zeros(L+1,div(L,2)+1)
    this.Cℓki = zeros(L+1,div(L,2)+1,L+1,L+2)
    for ℓ in range(0,L), k in range(0,div(L,2))
        this.Cℓk[ℓ+1,k+1] = (-1)^k * exp( sum(log.(1:2*ℓ-2*k)) - sum(log.(1:k)) - sum(log.(1:ℓ-k)) - sum(log.(1:ℓ-2*k)) )
        for i in range(0,1), g in range(0,ℓ-2*k+i)
            this.Cℓki[ℓ+1,k+1,i+1,g+1] = 2 * exp( sum(log.(1:ℓ-2*k+i)) - sum(log.(1:g)) - sum(log.(1:ℓ-2*k+i-g)) ) * (-1)^g
        end
    end
    
end