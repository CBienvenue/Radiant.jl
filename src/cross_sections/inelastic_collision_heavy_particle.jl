"""
    inelastic_collision_heavy_particle(Ei::Float64,W::Float64,particle::Particle)

Gives the inelastic differential cross-sections with atomic electrons for heavy particles.

# Input Argument(s)
- `Zi::Real` : number of electron(s).
- `Ei::Float64` : incoming particle energy.
- `W::Float64` : energy lost by the positron.

# Output Argument(s)
- `σs::Float64` : inelastic differential cross-sections with atomic electrons for heavy
  particles.

# Reference(s)
- Salvat and Heredia (2024), Electromagnetic interaction models for Monte-Carlo simulation
  of protons and alpha particles.

"""
function inelastic_collision_heavy_particle(Zi::Real,Ei::Float64,W::Float64,particle::Particle)

    #----
    # Initialization
    #----
    mₑ = 0.51099895069
    rₑ = 2.81794092e-13 # (in cm)
    Zc = get_charge(particle)
    ratio_mass = get_mass(particle)/mₑ
    γ = (Ei + ratio_mass)/ratio_mass
    β² = (γ^2-1)/γ^2
    R = 1/(1+(1/ratio_mass)^2+2*γ/ratio_mass)
    W_ridge = 2*β²*γ^2*R

    #----
    # Cross-section
    #----
    σs = 0.0
    if W_ridge > W
        σs = 2*π*rₑ^2/β² * Zc^2 * Zi * 1/W^2 * (1 - β²*W/W_ridge + (1-β²)/(2*ratio_mass^2)*W^2)
    end
    return σs
end

function integrate_inelastic_collision_heavy_particle(Z::Int64,Ei::Float64,n::Int64,particle::Particle,E⁻min::Float64=0.0)
    Nshells,Zi,Ui,_,_,_ = electron_subshells(Z)
    σn = 0
    for δi in range(1,1)
        σn += Zi[δi] * integrate_inelastic_collision_heavy_particle_per_subshell(Ei,n,Ui[δi],particle,E⁻min)
    end
    return σn
  end

function integrate_inelastic_collision_heavy_particle_per_subshell(Ei::Float64,n::Int64,Ui::Float64,particle::Particle,E⁻min::Float64=0.0)
    
    #----
    # Initialization
    #----
    mₑ = 0.51099895069
    rₑ = 2.81794092e-13 # (in cm)
    Zc = get_charge(particle)
    ratio_mass = get_mass(particle)/mₑ
    γ = (Ei + ratio_mass)/ratio_mass
    β² = (γ^2-1)/γ^2
    R = 1/(1+(1/ratio_mass)^2+2*γ/ratio_mass)
    W_ridge = 2*β²*γ^2*R
    
    #----
    # Cross-section
    #----
    σni = 0
    W⁺ = min(Ei,W_ridge)
    W⁻ = E⁻min+Ui
    if W⁺ > W⁻
        if n == 0
            J₀(W) = -1/W + β²/W_ridge*log(W) + (1-β²)/(2*ratio_mass^2)*W
            σni = J₀(W⁺) - J₀(W⁻)
        elseif n == 1
            J₁(W) = log(W) + β²/W_ridge*W + (1-β²)/(2*ratio_mass^2)*W^2/2
            σni = J₁(W⁺) - J₁(W⁻)
        else
            error("Integral is given only for n=0 or n=1.")
        end
    end
    σni *= 2*π*rₑ^2/β² * Zc^2
    return σni
end