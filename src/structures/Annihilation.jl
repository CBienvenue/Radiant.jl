"""
    Annihilation

Structure used to define parameters for production of multigroup annihilation cross-sections.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{String,String},Vector{String}} = Dict(("positrons","positrons") => ["A"],("positrons","photons") => ["P₋","P₊","P_inel","P_brems"],("photons","photons") => ["P_pp"])`: Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `("positrons","positrons") => ["A"]`: absorption of the incoming positron.
    - `("positrons","photons") => ["P₋"]`: production of the lowest energy photon following annihilation.
    - `("positrons","photons") => ["P₊"]`: production of the highest energy photon following annihilation.
    - `("positrons","photons") => ["P_inel"]`: production of annihilation photons from inelastic collisional positrons absorption following scattering under the cutoff energy.
    - `("positrons","photons") => ["P_brems"]`: production of annihilation photons from Bremsstrahlung positrons absorption following scattering under the cutoff energy.
    - `("photons","photons") => ["P_pp"]`: production of annihilation photons from absorption of positrons following their production under the cutoff energy.

"""
mutable struct Annihilation <: Interaction

    # Variable(s)
    name::String
    incoming_particle::Vector{String}
    interaction_particles::Vector{String}
    interaction_types::Dict{Tuple{String,String},Vector{String}}
    is_CSD::Bool
    is_AFP::Bool
    is_elastic::Bool
    is_preload_data::Bool
    is_subshells_dependant::Bool
    prior_interaction::Interaction
    scattering_model::String

    # Constructor(s)
    function Annihilation(;
        ### Initial values ###
        interaction_types = Dict(("positrons","positrons") => ["A"],("positrons","photons") => ["P₋","P₊","P_inel","P_brems"],("photons","photons") => ["P_pp"])
        ######################
        )
        this = new()
        this.name = "annihilation"
        this.set_interaction_types(interaction_types)
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.is_CSD = false
        this.is_AFP = false
        this.is_elastic = false
        this.is_preload_data = true
        this.is_subshells_dependant = false
        this.scattering_model = "BTE"
        return this
    end
end

# Method(s)
"""
    set_interaction_types(this::Annihilation,interaction_types::Dict{Tuple{String,String},Vector{String}})

To define the interaction types for annihilation processes.

# Input Argument(s)
- `this::Annihilation`: annihilation structure.
- `interaction_types::Dict{Tuple{String,String},Vector{String}}`: Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `("positrons","positrons") => ["A"]`: absorption of the incoming positron.
    - `("positrons","photons") => ["P₋"]`: production of the lowest energy photon following annihilation.
    - `("positrons","photons") => ["P₊"]`: production of the highest energy photon following annihilation.
    - `("positrons","photons") => ["P_inel"]`: production of annihilation photons from inelastic collisional positrons absorption following scattering under the cutoff energy.
    - `("positrons","photons") => ["P_brems"]`: production of annihilation photons from Bremsstrahlung positrons absorption following scattering under the cutoff energy.
    - `("photons","photons") => ["P_pp"]`: production of annihilation photons from absorption of positrons following their production under the cutoff energy.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> annihilation = Annihilation()
julia> annihilation.set_interaction_types( Dict(("positrons","positrons") => ["A"]) ) # Annihilation is set to be only absorption of positrons without any production of photons
```
"""
function set_interaction_types(this::Annihilation,interaction_types::Dict{Tuple{String,String},Vector{String}})
    this.interaction_types = interaction_types
end

function in_distribution(this::Annihilation)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

function out_distribution(this::Annihilation,type::String)
    if type ∈ ["A","P₋","P₊"]
        is_dirac = false
        N = 8
        quadrature = "gauss-legendre"
    else
        is_dirac = true
        N = 1
        quadrature = "dirac"
    end
    return is_dirac, N, quadrature
end

function bounds(this::Annihilation,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String)
    γ = Ei+1
    if type == "P₋"
        Ef⁻ = min(Ef⁻,(γ+1)/2)
        Ef⁺ = max(Ef⁺,(γ+1)/(γ+1+sqrt(γ^2-1)))
        if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    elseif type == "P₊"
        Ef⁻ = min(Ef⁻,(γ+1)*(1-1/(γ+1+sqrt(γ^2-1))))
        Ef⁺ = max(Ef⁺,(γ+1)/2)
        if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    elseif type ∈ ["P_inel","P_brems","P_pp"]
        if (Ef⁻-Ef⁺ < 0 || ~(Ef⁺ < 1 < Ef⁻) ) isSkip = true else isSkip = false end
    else
        error("Unknown type of method for annihilation.")
    end
    return Ef⁻,Ef⁺,isSkip
end

function dcs(this::Annihilation,L::Int64,Ei::Float64,Ef::Float64,type::String,gi::Int64,Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,iz::Int64,Ein::Vector{Float64},Ec::Float64)

    # Initialization
    rₑ = 2.81794092e-13 # (in cm)
    γ = Ei+1
    σs = 0.0
    σℓ = zeros(L+1)
    if type ∈ ["P₋","P₊"] S_Heitler(x) = -(γ+1)^2 + (γ^2+4γ+1)/x - 1/x^2 end

    # Lowest energy photon
    if type == "P₋"
        ζ = Ef/(Ei+1)
        ζmin = 1/(γ+1+sqrt(γ^2-1))
        ζmax = 1/2
        if ζmin ≤ ζ ≤ ζmax

            # Scattering cross-sections
            σs = Z[iz]*π*rₑ^2/((γ+1)^2*(γ^2-1)) * (S_Heitler(ζ)+S_Heitler(1-ζ))

            # Legendre moments
            μ = (γ+1-1/ζ)/sqrt(γ^2-1)
            Pℓμ = legendre_polynomials(L,μ)
            for ℓ in range(0,L) σℓ[ℓ+1] += Pℓμ[ℓ+1] * σs end
        end

    # Highest energy photon
    elseif type == "P₊"

        ζ = (Ei+2-Ef)/(Ei+1)
        ζmin = 1/(γ+1+sqrt(γ^2-1))
        ζmax = 1/2
        if ζmin ≤ ζ ≤ ζmax
            
            # Scattering cross-sections
            σs = Z[iz]*π*rₑ^2/((γ+1)^2*(γ^2-1)) * (S_Heitler(ζ)+S_Heitler(1-ζ))

            # Legendre moments
            μ = (γ+1-1/(1-ζ))/sqrt(γ^2-1)
            Pℓμ = legendre_polynomials(L,μ)
            for ℓ in range(0,L) σℓ[ℓ+1] += Pℓμ[ℓ+1] * σs end
        end

    # Annihilation of positrons scattered under the cutoff from inelastic collisionnal interaction
    elseif type == "P_inel"
        σa = tcs(this.prior_interaction,Ei,min(Ein[end],Ec),"positrons",Z[iz])
        σℓ[1] += 2 * σa

    # Annihilation of positrons scattered under the cutoff from Bremsstrahlung interaction
    elseif type == "P_brems"
        σa = tcs(this.prior_interaction,Ei,Z[iz],min(Ein[end],Ec),iz,"positrons","S",[Ein[end]])
        σℓ[1] += 2 * σa

    # Annihilation of positrons produced under the cutoff following pair production event
    elseif type == "P_pp"
        σa = tcs(this.prior_interaction,Ei,Z[iz],iz,[Ein[end]],"A")
        σℓ[1] += 2 * σa
    else
        error("Unknown interaction.")
    end
    return σℓ
end

function dcs_cutoff(this::Annihilation,L::Int64,Ngf::Int64,E_in::Vector{Float64},E_out::Vector{Float64},incoming_particle::String,Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64)
    
    # Initialization
    Σs = zeros(Ngf,L+1)
    Σsₑ = zeros(Ngf)

    # Annihilation for positrons slowing-down under the cutoff
    if incoming_particle == "positrons"
        S_cutoff = sp_dispatch(this.prior_interaction,Z,ωz,ρ,"solid",E_in[end],0.0,incoming_particle,E_in[end],E_out)

        σa = 0.0
        σa2 = 0.0
        σa3 = 0.0
        if typeof(this.prior_interaction) == Inelastic_Leptons
            for iz in range(1,length(Z))
                Ei = E_in[end]
                γ = Ei+1
                rₑ = 2.81794092E-13       # (in cm)
                β² = Ei*(Ei+2)/(Ei+1)^2
                Nshells,Zi,Ui,Ti,ri,subshells = electron_subshells(Z[iz])
                for δi in range(1,Nshells)
                    if Ei ≥ Ui[δi]
                        βᵢ² = Ui[δi]*(Ui[δi]+2)/(Ui[δi]+1)^2
                        fv = βᵢ²/β² * (β²/(β²+βᵢ²-β²*βᵢ²))^(3/2)
                        σa += π*rₑ^2/Ui[δi]^2 * Zi[δi] * fv * (1 + 2/3*(1-Ui[δi]/(2*Ei))*log(2.7+sqrt(Ei/Ui[δi]-1))) * (1-Ui[δi]/Ei)^(3/2) * ωz[iz] * nuclei_density(Z[iz],ρ)

                        b = ((γ-1)/γ)^2
                        b1 = b * (2*(γ+1)^2-1)/(γ^2-1)
                        b2 = b * (3*(γ+1)^2+1)/(γ+1)^2
                        b3 = b * (2*(γ-1)*γ)/(γ+1)^2
                        b4 = b * (γ-1)^2/(γ+1)^2
                        J₀⁺(x) = -1/(x+Ui[δi]) - b1*log(x+Ui[δi])/Ei + b2*x/Ei^2 - b3*(x^2/2+Ui[δi]*x)/Ei^3 + b4*(x^3/3+Ui[δi]*x^2+Ui[δi]^2*x)/Ei^4
                        σa2 += 2*π*rₑ^2/β² * Zi[δi] * (J₀⁺(Ei)-J₀⁺(Ui[δi])) * ωz[iz] * nuclei_density(Z[iz],ρ)


                        Eq = Ei - Ti[δi] - Ui[δi]
                        if Eq > Ti[δi] + Ui[δi]
                            σa3 += π*rₑ^2/Ui[δi]^2 * Zi[δi] * Ui[δi]/Eq * (1+2/3*Ti[δi]/Ui[δi]-Ui[δi]/(Eq-Ti[δi])) * ωz[iz] * nuclei_density(Z[iz],ρ)
                        elseif Ui[δi] ≤ Eq ≤ Ti[δi] + Ui[δi]
                            σa3 += π*rₑ^2/Ui[δi]^2 * Zi[δi] * Ui[δi]/Eq * 2/3*sqrt(Ui[δi]/Ti[δi])*(Eq/Ui[δi]-1)^(3/2) * ωz[iz] * nuclei_density(Z[iz],ρ)
                        end

                    end
                end
            end
        end
        println([σa,σa2,σa3,S_cutoff/(E_in[end-1]-E_in[end])])

        for ig in range(1,Ngf)
            if E_out[ig] ≥ 1 ≥ E_out[ig+1]
                Σs[ig,1] = 2 * σa3 # (S_cutoff+Sa) / E_in[end] #(E_in[end-1]-E_in[end])
                Σsₑ[ig] = Σs[ig,1] * ((E_out[ig]+E_out[ig+1])/2-1)
                break
            end
        end
    end
    return Σs,Σsₑ
end

function tcs(this::Annihilation,Ei::Float64,Z::Int64)
    γ = Ei+1
    rₑ = 2.81794092e-13 # (in cm)
    ζmin = 1/(γ+1+sqrt(γ^2-1))
    σt = Z*π*rₑ^2/((γ+1)*(γ^2-1)) * ( (γ^2+4*γ+1)*log(1/(4*ζmin*(1-ζmin))) -2*(γ+1)^2*(1/2-ζmin) + 4 - 1/ζmin - 1/(1-ζmin) )
    return σt
end

function preload_data(this::Annihilation,Z::Vector{Int64},Emax::Float64,Emin::Float64,L::Int64,type::String,Eout::Vector{Float64},Ein::Vector{Float64},interactions::Vector{Interaction})

    # Preload interaction prior to positron scattering under the cutoff
    interaction = missing
    if type ∈ ["P_inel","P_brems","P_pp"]
        if type == "P_inel"
            # Get impact ionization information for inelastic scattering
            for i in interactions
                if typeof(i) == Inelastic_Leptons
                    interaction = i
                    break
                end
            end
            if ismissing(interaction) interaction = Inelastic_Leptons() end
        elseif type == "P_brems"
            # Get impact ionization information for inelastic scattering
            for i in interactions
                if typeof(i) == Bremsstrahlung
                    interaction = i
                    break
                end
            end
            if ismissing(interaction) interaction = Bremsstrahlung() end
            interaction.preload_data(Z,Emax,Emin,L)
        elseif type == "P_pp"
            # Get impact ionization information for inelastic scattering
            for i in interactions
                if typeof(i) == Pair_Production
                    interaction = i
                    break
                end
            end
            if ismissing(interaction) interaction = Pair_Production() end
            interaction.preload_data(Z,Emax,Emin,Eout,L)
        end
        this.prior_interaction = interaction
    end
end