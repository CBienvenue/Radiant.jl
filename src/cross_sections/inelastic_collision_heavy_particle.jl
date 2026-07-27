"""
    HeavyInelasticCache

Stores reusable constants for heavy-particle inelastic collisions with one subshell
electron population at one incident energy.

# Field(s)
- `Zi::Float64` : number of electrons in the subshell.
- `Ei::Float64` : incoming particle energy [in m_e*c^2].
- `mₑ::Float64` : electron rest energy [MeV].
- `rₑ::Float64` : classical electron radius [cm].
- `Zc::Float64` : incoming particle charge.
- `ratio_mass::Float64` : incoming particle mass divided by the electron rest mass.
- `γ::Float64` : Lorentz factor of the incoming particle.
- `β²::Float64` : squared reduced speed of the incoming particle.
- `R::Float64` : recoil factor used in the maximum energy transfer.
- `W_ridge::Float64` : maximum reduced energy transfer.
- `A::Float64` : leading cross-section coefficient.
- `mu_lin::Float64` : coefficient used in the angular small-energy-transfer expansion.
"""
struct HeavyInelasticCache
    Zi::Float64
    Ei::Float64
    mₑ::Float64
    rₑ::Float64
    Zc::Float64
    ratio_mass::Float64
    γ::Float64
    β²::Float64
    R::Float64
    W_ridge::Float64
    A::Float64  # 2*π*rₑ^2/β² * Zc^2 * Zi
    mu_lin::Float64  # For angular derivatives: 1 / (Ei * (Ei + 2*ratio_mass))
end

"""
    HeavyInelasticCache(Zi::Real, Ei::Float64, particle::Particle)

Builds a cache of heavy-particle inelastic-collision constants.

# Input Argument(s)
- `Zi::Real` : number of electrons in the subshell.
- `Ei::Float64` : incoming particle energy [in m_e*c^2].
- `particle::Particle` : incoming heavy charged particle.

# Output Argument(s)
- `cache::HeavyInelasticCache` : reusable constants for cross-section and feed integrations.
"""
function HeavyInelasticCache(Zi::Real, Ei::Float64, particle::Particle)
    mₑ = 0.51099895069
    rₑ = 2.81794092e-13
    Zc = get_charge(particle)
    ratio_mass = get_mass(particle)/mₑ
    γ = (Ei + ratio_mass)/ratio_mass
    β² = (γ^2-1)/γ^2
    R = 1/(1+(1/ratio_mass)^2+2*γ/ratio_mass)
    W_ridge = 2*β²*γ^2*R
    A = 2*π*rₑ^2/β² * Zc^2 * Zi
    mu_lin = 1.0 / (Ei * (Ei + 2.0 * ratio_mass))
    return HeavyInelasticCache(Float64(Zi), Ei, mₑ, rₑ, Zc, ratio_mass, γ, β², R, W_ridge, A, mu_lin)
end

"""
    inelastic_collision_heavy_particle(Zi::Real, Ei::Float64, W::Float64,
    particle::Particle)

Gives the inelastic differential cross-sections with atomic electrons for heavy particles.

# Input Argument(s)
- `Zi::Real` : number of electron(s).
- `Ei::Float64` : incoming particle energy [in m_e*c^2].
- `W::Float64` : energy transferred to the atomic electron [in m_e*c^2].
- `particle::Particle` : incoming heavy charged particle.

# Output Argument(s)
- `σs::Float64` : inelastic differential cross-section with atomic electrons for heavy particles.

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

"""
    heavy_sigma_components(cache::HeavyInelasticCache, W::Float64)

Splits the heavy-particle inelastic differential cross-section into its singular
`A / W^2` part and the remaining regular part.

# Input Argument(s)
- `cache::HeavyInelasticCache` : reusable constants for this subshell and particle.
- `W::Float64` : energy transferred to the atomic electron [in m_e*c^2].

# Output Argument(s)
- `sigma_reg::Float64` : regular contribution, equal to `sigma_full(W) - A / W^2`.
- `A_over_W2::Float64` : leading singular contribution.
"""
function heavy_sigma_components(cache::HeavyInelasticCache, W::Float64)
    if W <= 0
        return 0.0, 0.0
    end
    if W > cache.W_ridge
        return 0.0, 0.0
    end

    # Full sigma
    σ_full = cache.A * (1/W^2) * (1 - cache.β²*W/cache.W_ridge + (1 - cache.β²)/(2*cache.ratio_mass^2)*W^2)

    A_over_W2 = cache.A / W^2
    sigma_reg = σ_full - A_over_W2
    return sigma_reg, A_over_W2
end

"""
    integrate_A_over_W2_per_subshell(cache::HeavyInelasticCache,
    Ef_minus::Float64, Ef_plus::Float64)

Integrates the leading `A / W^2` part analytically over an outgoing-energy group.
The integration bounds are converted with `W = Ei - Ef` and clipped to
`0 < W <= W_ridge`.

# Input Argument(s)
- `cache::HeavyInelasticCache` : reusable constants for this subshell and particle.
- `Ef_minus::Float64` : upper bound of the outgoing particle energy group [in m_e*c^2].
- `Ef_plus::Float64` : lower bound of the outgoing particle energy group [in m_e*c^2].

# Output Argument(s)
- `σ::Float64` : integral of `A / W^2` over the clipped group.
"""
function integrate_A_over_W2_per_subshell(cache::HeavyInelasticCache, Ef_minus::Float64, Ef_plus::Float64)
    # Energy bounds in W = Ei - Ef
    W_plus = cache.Ei - Ef_plus
    W_minus = cache.Ei - Ef_minus

    # Clip to (0, W_ridge]
    W_lo = max(W_minus, 1e-300)
    W_hi = min(W_plus, cache.W_ridge)
    if W_hi <= W_lo
        return 0.0
    end

    # Numerically stable form: (1/W_lo - 1/W_hi) = (W_hi - W_lo) / (W_lo * W_hi)
    ΔW_clip = W_hi - W_lo
    if ΔW_clip <= 0.0
        return 0.0
    end
    denom = W_lo * W_hi
    # Guard against underflow/overflow by ensuring denom is not subnormal
    if denom == 0.0
        # If one bound is effectively zero, return the dominant 1/W_lo term
        return cache.A * (1.0 / W_lo)
    end
    return cache.A * (ΔW_clip / denom)
end

"""
    heavy_leading_1overW_coeff(cache::HeavyInelasticCache, l::Int64)

Computes the coefficient of the leading `1 / W` angular correction for Legendre order `l`.

# Input Argument(s)
- `cache::HeavyInelasticCache` : reusable constants for this subshell and particle.
- `l::Int64` : Legendre order.

# Output Argument(s)
- `coef::Float64` : coefficient such that the leading correction is `coef / W`.
"""
function heavy_leading_1overW_coeff(cache::HeavyInelasticCache, l::Int64)
    if l == 0
        return 0.0
    end
    # P_l'(1) = l*(l+1)/2
    Plp1 = l * (l + 1) / 2.0
    # coef for 1/W term: -A * P_l'(1) * mu_lin
    coef = -cache.A * Plp1 * cache.mu_lin
    return coef
end

"""
    integrate_leading_1overW_per_subshell(cache::HeavyInelasticCache,
    Ef_minus::Float64, Ef_plus::Float64, l::Int64)

Integrates the leading `1 / W` angular correction analytically over an outgoing-energy
group for one subshell and one Legendre order.

# Input Argument(s)
- `cache::HeavyInelasticCache` : reusable constants for this subshell and particle.
- `Ef_minus::Float64` : upper bound of the outgoing particle energy group [in m_e*c^2].
- `Ef_plus::Float64` : lower bound of the outgoing particle energy group [in m_e*c^2].
- `l::Int64` : Legendre order.

# Output Argument(s)
- `σl::Float64` : integral of the leading angular correction over the clipped group.
"""
function integrate_leading_1overW_per_subshell(cache::HeavyInelasticCache, Ef_minus::Float64, Ef_plus::Float64, l::Int64)
    if l == 0
        return 0.0
    end

    W_plus = cache.Ei - Ef_plus
    W_minus = cache.Ei - Ef_minus
    W_lo = max(W_minus, 1e-300)
    W_hi = min(W_plus, cache.W_ridge)
    if W_hi <= W_lo
        return 0.0
    end

    coef = heavy_leading_1overW_coeff(cache, l)
    ΔW = W_hi - W_lo
    # stable log: log(W_hi/W_lo) = log1p(ΔW/W_lo) when ΔW small
    if ΔW / W_lo < 1e-6
        log_ratio = log1p(ΔW / W_lo)
    else
        log_ratio = log(W_hi) - log(W_lo)
    end
    return coef * log_ratio
end

"""
    angular_inelastic_collision_heavy_particle(Ei::Float64, Ef::Float64, L::Int64,
    particle::Particle)

Computes Legendre polynomial values for the inelastic angular deflection of a heavy
particle after an energy loss.

# Input Argument(s)
- `Ei::Float64` : incoming particle energy [in m_e*c^2].
- `Ef::Float64` : outgoing particle energy [in m_e*c^2].
- `L::Int64` : maximum Legendre order.
- `particle::Particle` : incoming heavy charged particle.

# Output Argument(s)
- `Pl::Vector{Float64}` : Legendre polynomial values from order 0 to `L`.
"""
function angular_inelastic_collision_heavy_particle(Ei::Float64, Ef::Float64, L::Int64, particle::Particle)
    mₑ = 0.51099895069
    M  = get_mass(particle) / mₑ  # ratio_mass

    W = Ei - Ef
    if W <= 0.0 || Ef <= 0.0 || Ei <= 0.0
        return legendre_polynomials_up_to_L(L, 1.0)  # forward
    end

    pi2 = Ei*(Ei + 2.0*M)
    pf2 = Ef*(Ef + 2.0*M)
    pe2 = W*(W + 2.0)

    denom = 2.0*sqrt(pi2*pf2)
    if denom <= 0.0
        μ = 1.0
    else
        μ = (pi2 + pf2 - pe2)/denom
        μ = clamp(μ, -1.0, 1.0)
    end

    return legendre_polynomials_up_to_L(L, μ)
end

"""
    integrate_inelastic_collision_heavy_particle(Z::Int64, Ei::Float64, n::Int64,
    particle::Particle, W⁻min::Float64=0.0, is_subshells::Bool=true)

Computes the subshell-summed analytical energy-loss moment of the heavy-particle
inelastic collision cross-section for one element.

# Input Argument(s)
- `Z::Int64` : atomic number.
- `Ei::Float64` : incoming particle energy [in m_e*c^2].
- `n::Int64` : energy-loss moment order, either 0 or 1.
- `particle::Particle` : incoming heavy charged particle.
- `W⁻min::Float64=0.0` : minimum energy loss W of the incoming particle [in m_e*c^2].
- `is_subshells::Bool=true` : boolean indicating if bounded subshells are resolved.

# Output Argument(s)
- `σn::Float64` : analytical moment summed over all subshells.
"""
function integrate_inelastic_collision_heavy_particle(Z::Int64,Ei::Float64,n::Int64,particle::Particle,W⁻min::Float64=0.0,is_subshells::Bool=true)
    Nshells,Zi,Ui,_,_,_ = electron_subshells(Z,~is_subshells)
    σn = 0
    for δi in range(1,Nshells)
        σn += Zi[δi] * integrate_inelastic_collision_heavy_particle_per_subshell(Ei,n,Ui[δi],particle,W⁻min)
    end
    return σn
end

"""
    integrate_inelastic_collision_heavy_particle_per_subshell(Ei::Float64, n::Int64,
    Ui::Float64, particle::Particle, W⁻min::Float64=0.0)

Computes the analytical energy-loss moment of the heavy-particle inelastic collision
cross-section for one subshell.

# Input Argument(s)
- `Ei::Float64` : incoming particle energy [in m_e*c^2].
- `n::Int64` : energy-loss moment order, either 0 or 1.
- `Ui::Float64` : subshell binding energy [in m_e*c^2].
- `particle::Particle` : incoming heavy charged particle.
- `W⁻min::Float64=0.0` : minimum energy loss W of the incoming particle [in m_e*c^2].

# Output Argument(s)
- `σni::Float64` : analytical moment for the subshell.
"""
function integrate_inelastic_collision_heavy_particle_per_subshell(Ei::Float64,n::Int64,Ui::Float64,particle::Particle,W⁻min::Float64=0.0)

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
    W⁻ = max(W⁻min,Ui)
    if W⁺ > W⁻
        ΔW = W⁺ - W⁻
        if n == 0
            if ΔW / W⁺ < 1e-6
                σni = ΔW/W⁺^2 - (β²/W_ridge)*(ΔW/W⁺ - 0.5*(ΔW/W⁺)^2) + ((1-β²)/(2*ratio_mass^2))*ΔW
            else
                σni = ΔW / (W⁺ * W⁻) - (β² / W_ridge) * log1p(ΔW / W⁻) + ((1 - β²)/(2*ratio_mass^2)) * ΔW
            end
        elseif n == 1
            if ΔW / W⁻ < 1e-6
                log_term = ΔW / W⁻ - 0.5*(ΔW/W⁻)^2
            else
                log_term = log1p(ΔW / W⁻)
            end
            σni = log_term - (β² / W_ridge) * ΔW + ((1 - β²)/(4 * ratio_mass^2)) * ΔW * (W⁺ + W⁻)
        else
            error("Integral is given only for n=0 or n=1.")
        end
    end
    σni *= 2*π*rₑ^2/β² * Zc^2
    return σni
end

"""
    feed_analytical_heavy_particle(cache::HeavyInelasticCache, Ef⁻::Float64,
    Ef⁺::Float64)

Computes the scalar analytical feed integral for heavy-particle inelastic collisions over
an outgoing-energy group.

# Input Argument(s)
- `cache::HeavyInelasticCache` : reusable constants for this subshell and particle.
- `Ef⁻::Float64` : upper bound of the outgoing particle energy group [in m_e*c^2].
- `Ef⁺::Float64` : lower bound of the outgoing particle energy group [in m_e*c^2].

# Output Argument(s)
- `σs::Float64` : scalar feed integral over the group.

"""
function feed_analytical_heavy_particle(cache::HeavyInelasticCache, Ef⁻::Float64, Ef⁺::Float64)
    # Compute the scalar analytical integral of σ over the group (l=0 moment)
    σs = 0.0
    ΔEf = Ef⁻ - Ef⁺

    # Only integrate if energy transfer is within physical limits
    if cache.W_ridge <= 0 || ΔEf <= 0
        return σs
    end

    # Energy bounds for integration (W = Ei - Ef)
    W⁺ = cache.Ei - Ef⁺
    W⁻ = cache.Ei - Ef⁻

    if W⁺ <= W⁻
        return σs
    end

    #----
    # Analytical integration over the energy group for the l=0 moment
    #----
    ΔW = W⁺ - W⁻

    if ΔW / W⁺ < 1e-6
        σs = ΔW/W⁺^2 - (cache.β²/cache.W_ridge)*(ΔW/W⁺ - 0.5*(ΔW/W⁺)^2) + ((1-cache.β²)/(2*cache.ratio_mass^2))*ΔW
    else
        σs = ΔW / (W⁺ * W⁻) - (cache.β² / cache.W_ridge) * log1p(ΔW / W⁻) + ((1 - cache.β²)/(2*cache.ratio_mass^2)) * ΔW
    end
    σs *= cache.A

    return σs
end

"""
    feed_first_moment_heavy_particle(cache::HeavyInelasticCache, Ef⁻::Float64,
    Ef⁺::Float64)

Computes the first energy-loss moment of the heavy-particle inelastic feed integral over
an outgoing-energy group.

# Input Argument(s)
- `cache::HeavyInelasticCache` : reusable constants for this subshell and particle.
- `Ef⁻::Float64` : upper bound of the outgoing particle energy group [in m_e*c^2].
- `Ef⁺::Float64` : lower bound of the outgoing particle energy group [in m_e*c^2].

# Output Argument(s)
- `M₁::Float64` : first energy-loss moment over the group.

"""
function feed_first_moment_heavy_particle(cache::HeavyInelasticCache, Ef⁻::Float64, Ef⁺::Float64)
    # Energy bounds for integration (W = Ei - Ef)
    W⁺ = cache.Ei - Ef⁺
    W⁻ = cache.Ei - Ef⁻

    if W⁺ <= W⁻
        return 0.0
    end

    #----
    # First moment (n=1): ∫ W * σ(W) dW
    #----
    ΔW = W⁺ - W⁻

    if ΔW / W⁻ < 1e-6
        log_term = ΔW / W⁻ - 0.5*(ΔW/W⁻)^2
    else
        log_term = log1p(ΔW / W⁻)
    end
    M₁ = log_term - (cache.β² / cache.W_ridge) * ΔW + ((1 - cache.β²)/(4 * cache.ratio_mass^2)) * ΔW * (W⁺ + W⁻)
    M₁ *= cache.A

    return M₁
end