"""
Inelastic_Leptons

Structure used to define parameters for production of multigroup inelastic collisional cross-sections for leptons.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{String,String},Vector{String}} = Dict(("positrons","positrons") => ["S"],("positrons","electrons") => ["P"],("electrons","electrons") => ["S","P"])`: Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `("positrons","positrons") => ["S"]`: scattering of incident positrons
    - `("electrons","electrons") => ["S"]`: scattering of incident electrons
    - `("positrons","electrons") => ["P"]`: production of electrons by incident positrons
    - `("electrons","electrons") => ["P"]`: production of electrons by incident electrons

"""
mutable struct Inelastic_Leptons <: Interaction

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
    is_shell_correction::Bool
    is_density_correction::Bool
    density_correction_type::String
    plasma_energy::Float64
    effective_mean_excitation_energy::Float64
    shell_correction::Function

    # Constructor(s)
    function Inelastic_Leptons()
        this = new()
        this.name = "inelastic_leptons"
        this.interaction_types = Dict(("positrons","positrons") => ["S"],("positrons","electrons") => ["P"],("electrons","electrons") => ["S","P"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.density_correction_type = "fano"
        this.is_CSD = true
        this.is_AFP = true
        this.is_elastic = false
        this.is_preload_data = true
        this.is_subshells_dependant = true
        this.is_shell_correction = true
        this.is_density_correction = true
        return this
    end

end

# Method(s)
"""
    set_interaction_types(this::Inelastic_Leptons,interaction_types::Dict{Tuple{String,String},Vector{String}})

To define the interaction types for inelastic leptons processes.

# Input Argument(s)
- `this::Inelastic_Leptons`: inelastic leptons structure.
- `interaction_types::Dict{Tuple{String,String},Vector{String}}`: Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `("positrons","positrons") => ["S"]`: scattering of incident positrons
    - `("electrons","electrons") => ["S"]`: scattering of incident electrons
    - `("positrons","electrons") => ["P"]`: production of electrons by incident positrons
    - `("electrons","electrons") => ["P"]`: production of electrons by incident electrons

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_leptons = Inelastic_Leptons()
julia> elastic_leptons.set_interaction_types( Dict(("electrons","electrons") => ["S"]) ) # No knock-on electrons.
```
"""
function set_interaction_types(this::Inelastic_Leptons,interaction_types::Dict{Tuple{String,String},Vector{String}})
    this.interaction_types = interaction_types
end

function in_distribution(this::Inelastic_Leptons)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

function out_distribution(this::Inelastic_Leptons)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

function bounds(this::Inelastic_Leptons,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String,Ec::Float64,Ui::Float64,particle::String)
    if particle == "electrons"
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
    elseif particle == "positrons"
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

function dcs(this::Inelastic_Leptons,L::Int64,Ei::Float64,Ef::Float64,type::String,particle::String,Ui::Float64,Zi::Float64,Ti::Float64)

    # Initialization
    rₑ = 2.81794092e-13 # (in cm)
    γ = Ei+1
    β² = Ei*(Ei+2)/(Ei+1)^2
    σs = 0.0
    σℓ = zeros(L+1)
    if type == "S"
        Ep = Ef
        W = Ei-Ep-Ui
    elseif type == "P"
        W = Ef
        Ep = Ei-W-Ui
    else
        error("Unknown type")
    end

    # Close collisions
    if W ≥ 0 #Ui
        if particle == "electrons"
           Pi = 1 #Ei/(Ei+Ui+Ti)
           Gi = 0 #4*Ti/3 * (1/(W+Ui)^3 + 1/(Ei-W)^3)
           F = Pi * ( 1/(W+Ui)^2 + 1/(Ei-W)^2 + 1/(Ei+1)^2 - (2*Ei+1)/(Ei+1)^2 * 1/((Ei-W)*(W+Ui)) + Gi )
        elseif particle == "positrons"
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


function tcs(this::Inelastic_Leptons,Ei::Float64,Ec::Float64,particle::String,Z::Int64)

    # Inititalisation
    rₑ = 2.81794092E-13       # (in cm)
    γ = Ei+1
    β² = Ei*(Ei+2)/(Ei+1)^2
    Nshells,Zi,Ui,Ti,ri,subshells = electron_subshells(Z)
    σt = 0.0

    # Close collisions
    if particle == "electrons"
        for δi in range(1,Nshells)
            Wmax = (Ei-Ui[δi])/2
            Wmin = Ei-Ec
            if Wmax > Wmin
                Pi = 1 #Ei/(Ei+Ui[δi]+Ti[δi])
                J₀⁻(x) = -1/(x+Ui[δi]) + 1/(Ei-x) + x/(Ei+1)^2 - (2*Ei+1)/(Ei+1)^2 * (log(x+Ui[δi])-log(Ei-x))/(Ei+Ui[δi])
                #+ 4*Ti[δi]/3 * (1/(2*(Ei-x)^2)-1/(2*(x+Ui[δi])^2))
                σt += 2*π*rₑ^2/β² * Zi[δi] * Pi * (J₀⁻(Wmax)-J₀⁻(Wmin))
            end
        end
    elseif particle == "positrons"
        b = ((γ-1)/γ)^2
        b1 = b * (2*(γ+1)^2-1)/(γ^2-1)
        b2 = b * (3*(γ+1)^2+1)/(γ+1)^2
        b3 = b * (2*(γ-1)*γ)/(γ+1)^2
        b4 = b * (γ-1)^2/(γ+1)^2
        for δi in range(1,Nshells)
            Wmax = Ei-Ui[δi]
            Wmin = Ei-Ec
            if Wmax > Wmin
                J₀⁺(x) = -1/(x+Ui[δi]) - b1*log(x+Ui[δi])/Ei + b2*x/Ei^2 - b3*(x^2/2+Ui[δi]*x)/Ei^3 + b4*(x^3/3+Ui[δi]*x^2+Ui[δi]^2*x)/Ei^4   #-1/x - b1*log(x)/Ei + b2*x/Ei^2 - b3*x^2/(2*Ei^3) + b4*x^3/(3*Ei^4)
                σt += 2*π*rₑ^2/β² * Zi[δi] * (J₀⁺(Wmax)-J₀⁺(Wmin))
            end
        end
    else
        error("Unknown particle")
    end
    return σt
end

function sp(this::Inelastic_Leptons,Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,state_of_matter::String,Ei::Float64,Ec::Float64,particle::String)

    # Initialization
    rₑ = 2.81794092e-13 # (in cm)
    γ = Ei+1
    β² = Ei*(Ei+2)/(Ei+1)^2
    𝒩ₑ = Z.*nuclei_density.(Z,ρ)        # (in cm⁻³)
    𝒩ₑ_eff = sum(ωz.*𝒩ₑ)               # (in cm⁻³)
    I = this.effective_mean_excitation_energy
    if (this.is_density_correction) δF = fermi_density_effect(Z,ωz,ρ,Ei,state_of_matter,this.density_correction_type) else δF = 0 end
    if (this.is_shell_correction) Cz = this.shell_correction(Z,ωz,Ei) else Cz = 0 end

    # Compute the total stopping power
    if particle == "electrons"
        f = (1-β²) - (2*γ-1)/γ^2*log(2) +((γ-1)/γ)^2/8
    elseif particle == "positrons"
        f = 2*log(2) - β²/12 * (23 + 14/(γ+1) + 10/(γ+1)^2 + 4/(γ+1)^3)
    else
        error("Unknown particle")
    end
    Stot = 2*π*rₑ^2/β² * 𝒩ₑ_eff * ( log((Ei/I)^2*(γ+1)/2) + f - δF - 2*Cz)
    
    # Compute the catastrophic Møller- or Bhabha- derived stopping power
    Sc = 0
    Nz = length(Z)
    for i in range(1,Nz)
        Nshells,Zi,Ui,Ti,ri,subshells = electron_subshells(Z[i])
        if particle == "electrons"
            for δi in range(1,Nshells)
                Wmax = (Ei-Ui[δi])/2
                Wmin = Ei-Ec
                if Wmax > Wmin
                    Pi = 1 #Ei/(Ei+Ui[δi]+Ti[δi])
                    J₁⁻(x) = log(x+Ui[δi]) + log(Ei-x) + (Ei+Ui[δi])/(Ei-x) + x*(x+2*Ui[δi])/(2*(Ei+1)^2) + (2*Ei+1)/(Ei+1)^2 * log(Ei-x) #+ 4*Ti[δi]/3 * ( (Ui[δi]+Ei)*(3*x+Ui[δi]-2*Ei)/(2*(x^3+(Ui[δi]-2*Ei)*x^2+(Ei^2-2*Ei*Ui[δi])*x+Ei^2*Ui[δi])) )
                    Sc += 2*π*rₑ^2/β² * ωz[i] * nuclei_density(Z[i],ρ) * Zi[δi] * Pi * (J₁⁻(Wmax)-J₁⁻(Wmin))
                end
            end
        elseif particle == "positrons"
            b = ((γ-1)/γ)^2
            b1 = b * (2*(γ+1)^2-1)/(γ^2-1)
            b2 = b * (3*(γ+1)^2+1)/(γ+1)^2
            b3 = b * (2*(γ-1)*γ)/(γ+1)^2
            b4 = b * (γ-1)^2/(γ+1)^2
            for δi in range(1,Nshells)
                Wmax = Ei-Ui[δi]
                Wmin = Ei-Ec
                if Wmax > Wmin
                    J₁⁺(x) = log(x+Ui[δi]) - b1*x/Ei + b2*(x^2/2+Ui[δi]*x)/Ei^2 - b3*(x^3/3+Ui[δi]*x^2+Ui[δi]^2*x)/Ei^3 + b4*(x^4/4+Ui[δi]*x^3+3*Ui[δi]^2*x^2/2+Ui[δi]^3*x)/Ei^4 #log(x)-Ui[δi]/x - b1*(Ui[δi]*log(x)+x)/Ei + b2*(x^2+2*Ui[δi]*x)/(2*Ei^2) - b3*x^2*(2*x+3*Ui[δi])/(6*Ei^3) + b4*x^3*(3*x+4*Ui[δi])/(12*Ei^4)
                    Sc += 2*π*rₑ^2/β² * ωz[i] * nuclei_density(Z[i],ρ) * Zi[δi] * (J₁⁺(Wmax)-J₁⁺(Wmin))
                end
            end
        end
    end

    # Compute the soft stopping power
    Sr = Stot - Sc

    return Sr
end

function mt(this::Inelastic_Leptons)
    return 0
end

function preload_data(this::Inelastic_Leptons,Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,particle::String)
    this.plasma_energy = plasma_energy(Z,ωz,ρ)
    this.effective_mean_excitation_energy = effective_mean_excitation_energy(Z,ωz)
    this.preload_shell_corrections(Z,ωz,particle)
end

"""
    preload_shell_corrections(this::Inelastic_Leptons,Z::Vector{Int64},ωz::Vector{Float64},particle::String)

Compute the shell corrections term for the collisional stopping powers.

# Input Argument(s)
- 'Z::Vector{Int64}': atomic number of the element(s) composing the material.
- 'ωz::Vector{Float64}': weight fraction of the element(s) composing the material.
- 'particle::String': particle type.

# Output Argument(s)
- 'shell_correction::Function': shell correction function.

# Author(s)
Charles Bienvenue

# Reference(s)
- Salvat (2023), SBETHE: Stopping powers of materials for swift charged particles from
  the corrected Bethe formula.

"""
function preload_shell_corrections(this::Inelastic_Leptons,Z::Vector{Int64},ωz::Vector{Float64},particle::String)
    
    # Initialization
    Nz = length(Z)
    S₀ = zeros(Nz)
    Ec = zeros(Nz)
    pn = zeros(Nz,6)
    E = zeros(Nz,153)
    Cz_prime = zeros(Nz,153)
    spline_Cz = Vector{Function}(undef,Nz)

    # Read shell correction data
    path = joinpath(find_package_root(), "data", "shell_corrections_salvat_2023.jld2")
    data = load(path)
    for iz in range(1,Nz)
        datai = data[particle][Z[iz]]
        for n in range(1,153)
            Cz_prime[iz,n] = datai["modified_shell_corrections"][n]
            E[iz,n] = datai["energies"][n]
        end
        Ec[iz] = datai["cutoff_energy"]
        for n in range(1,6)
            pn[iz,n] = datai["fit_parameters"][n]
        end
        S₀[iz] = datai["S₀"]

        spline_Cz[iz] = cubic_hermite_spline(E[iz,:],Cz_prime[iz,:])
    end
    
    this.shell_correction = function shell_correction(Z::Vector{Int64},ωz::Vector{Float64},Ei::Float64)
        Cz = 0.0
        β² = Ei*(Ei+2)/(Ei+1)^2
        Nz = length(Z)
        Zeff = sum(ωz.*Z) 
        for iz in range(1,Nz)
            if Ei < E[iz,end]
                Czprime = spline_Cz[iz](Ei)
            elseif γ < 2
                Czprime = sum(pn[iz,:] .* Ei.^(1:6))
            else
                Czprime = sum(pn[iz,:] .* 2 .^(1:6))
            end
            Cz += Z[iz]/Zeff * (Czprime + (S₀[iz]-Z[iz])/(2*Z[iz])*(log(β²*(Ei+1)^2)-β²))
        end
        return Cz
    end

end


