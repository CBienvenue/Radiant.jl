"""
    momentum_cm(E::Float64, m::Float64, M::Float64;
    kinematics::String="standard")

Gives the projectile center-of-mass momentum using either the ENDF
nonrelativistic convention (`"standard"`) or relativistic two-body kinematics
(`"relativistic"`). `E`, `m`, and `M` must use compatible energy units; the
returned momentum is in the same energy/c units.
"""
function momentum_cm(E::Float64, m::Float64, M::Float64;
                     kinematics::String="standard")::Float64
    kinematics_lc = lowercase(kinematics)
    if kinematics_lc == "standard"
        return (M / (m + M)) * sqrt(2.0 * m * E)
    elseif kinematics_lc == "relativistic"
        E1_total = m + E
        s = m^2 + M^2 + 2.0 * M * E1_total
        p2 = (s - (m + M)^2) * (s - (m - M)^2) / (4.0 * s)
        return sqrt(max(p2, 0.0))
    end
    error("Unknown momentum_cm kinematics $(kinematics) (should be standard or relativistic).")
end

"""
    momentum_lab(E::Float64, m::Float64; kinematics::String="standard")

Gives the initial projectile momentum in the laboratory frame using either
nonrelativistic (`"standard"`) or relativistic kinematics. `E` and `m` must use
compatible energy units; the returned momentum is in the same energy/c units.
"""
function momentum_lab(E::Float64, m::Float64; kinematics::String="standard")::Float64
    kinematics_lc = lowercase(kinematics)
    if kinematics_lc == "standard"
        return sqrt(2.0 * m * E)
    elseif kinematics_lc == "relativistic"
        return sqrt(max((E + m)^2 - m^2, 0.0))
    end
    error("Unknown momentum_lab kinematics $(kinematics) (should be standard or relativistic).")
end

"""
    cm_frame(E_in_mev::Float64, m_mev::Float64, M_mev::Float64)

Gives the center-of-mass frame of a two-body collision of a projectile with kinetic energy
`E_in_mev` on a target at rest (relativistic kinematics).

# Input Argument(s)
- `E_in_mev::Float64` : incident kinetic energy [MeV].
- `m_mev::Float64` : projectile rest energy [MeV].
- `M_mev::Float64` : target rest energy [MeV].

# Output Argument(s)
- `frame::NamedTuple` : `gamma` and `beta` of the center of mass, and the momentum
  `p_star` [MeV/c] and total energy `E1_star` [MeV] of the projectile in it.
"""
function cm_frame(E_in_mev::Float64, m_mev::Float64, M_mev::Float64)
    E1 = E_in_mev + m_mev
    p1 = sqrt(max(E1^2 - m_mev^2, 0.0))
    s = m_mev^2 + M_mev^2 + 2.0 * M_mev * E1
    E1_star = (s + m_mev^2 - M_mev^2) / (2.0 * sqrt(s))
    p_star = sqrt(max(E1_star^2 - m_mev^2, 0.0))
    beta = p1 / (E1 + M_mev)
    return (gamma=1.0 / sqrt(1.0 - beta^2), beta=beta, p_star=p_star, E1_star=E1_star)
end

"""
    cm_frame(cache::AbstractElasticDCSCache)

Gives the center-of-mass frame of an elastic collision from the incident energy and the
masses of an elastic DCS cache (`cm_frame(E_in_mev, m_mev, M_mev)`).

# Input Argument(s)
- `cache::AbstractElasticDCSCache` : elastic DCS cache.

# Output Argument(s)
- `frame::NamedTuple` : center-of-mass frame.
"""
cm_frame(cache::AbstractElasticDCSCache) = cm_frame(cache.E_in_mev, cache.m_mev, cache.M_mev)

"""
    scattered_energy_affine(cache::AbstractElasticDCSCache)

Gives the coefficients of the scattered-projectile energy, which is affine in the
center-of-mass cosine for elastic scattering: E_out = e0 + e1 μ_cm. Since E_out = E_in at
μ_cm = 1, the energy lost is exactly e1 (1 - μ_cm).

# Input Argument(s)
- `cache::AbstractElasticDCSCache` : elastic DCS cache providing the incident energy and
  masses.

# Output Argument(s)
- `e0::Float64` : scattered energy at μ_cm = 0 [MeV].
- `e1::Float64` : slope of the scattered energy in μ_cm [MeV].
"""
function scattered_energy_affine(cache::AbstractElasticDCSCache)
    f = cm_frame(cache)
    return f.gamma * f.E1_star - cache.m_mev, f.gamma * f.beta * f.p_star
end

"""
    one_minus_mu_lab(frame::NamedTuple, x::Float64)

Gives 1 - μ_lab of the scattered projectile for a center-of-mass cosine 1 - x, written
without cancellation for small angles.

# Input Argument(s)
- `frame::NamedTuple` : center-of-mass frame from `cm_frame`.
- `x::Float64` : 1 - μ_cm.

# Output Argument(s)
- `w::Float64` : 1 - μ_lab.
"""
function one_minus_mu_lab(frame::NamedTuple, x::Float64)
    p_perp2 = frame.p_star^2 * x * (2.0 - x)
    p_par = frame.gamma * (frame.p_star * (1.0 - x) + frame.beta * frame.E1_star)
    p_lab = sqrt(p_par^2 + p_perp2)
    p_lab > 0.0 || return 0.0
    return p_par > 0.0 ? p_perp2 / (p_lab * (p_lab + p_par)) : 1.0 - p_par / p_lab
end

"""
Relativistic scattered projectile energy for elastic scattering.
Assumes mu = cos(theta_cm) and target at rest.

Used as the forward map for mu inversion in
`mu_cm_from_Eout_relativistic` when computing mu bounds.
"""
function scattered_energy_relativistic(E_in_mev::Float64, mu::Float64,
                                       m_mev::Float64, M_mev::Float64)::Float64
    f = cm_frame(E_in_mev, m_mev, M_mev)
    E1_prime = f.gamma * (f.E1_star + f.beta * f.p_star * mu)
    return max(E1_prime - m_mev, 0.0)
end

"""
Invert mu_cm from scattered energy using bisection on [mu_lo, 1].
Assumes scattered energy is monotonic in mu_cm.

Used to convert Ef bounds into mu_cm bounds for elastic scattering
integration.
"""
function mu_cm_from_Eout_relativistic(E_in_mev::Float64, E_out_mev::Float64,
                                      m_mev::Float64, M_mev::Float64)

    energy_at_mu(mu_val::Float64) = scattered_energy_relativistic(E_in_mev, mu_val, m_mev, M_mev)
    E_lo = energy_at_mu(-1.0)
    E_hi = energy_at_mu(1.0)

    target = E_out_mev
    if target <= E_lo
        return -1.0
    elseif target >= E_hi
        return 1.0
    end

    lo = -1.0
    hi = 1.0
    for _ in 1:80
        mid = 0.5 * (lo + hi)
        E_mid = energy_at_mu(mid)
        if abs(E_mid - target) <= 1.0e-12 || (hi - lo) < 1.0e-12
            return mid
        end
        if E_mid < target
            lo = mid
        else
            hi = mid
        end
    end
    return 0.5 * (lo + hi)
end

"""
Lab-frame mu (cos(theta_lab)) from CM mu for elastic scattering of projectile.
"""
function mu_lab_from_mu_cm(mu_cm::Float64, E_in_mev::Float64,
                           m_mev::Float64, M_mev::Float64)
    f = cm_frame(E_in_mev, m_mev, M_mev)
    p_par = f.gamma * (f.p_star * mu_cm + f.beta * f.E1_star)
    p_perp = f.p_star * sqrt(max(1.0 - mu_cm^2, 0.0))
    p_lab = sqrt(p_par^2 + p_perp^2)
    return p_lab > 0.0 ? p_par / p_lab : 1.0
end

"""
    CMKinematics

Center-of-mass quantities of a collision with relativistic kinematics, as used by the
optical-model description of elastic scattering.

# Field(s)
- `p0c::Float64` : center-of-mass momentum × c [MeV].
- `mur::Float64` : relativistic reduced mass × c² [MeV].
- `k::Float64` : center-of-mass wave number [1/fm].
- `eta::Float64` : Sommerfeld parameter Z₁Z₂ e² μ_r / (ħ² k).
"""
struct CMKinematics
    p0c::Float64
    mur::Float64
    k::Float64
    eta::Float64
end

"""
    cm_kinematics(E::Float64, m::Float64, M::Float64, Z1Z2::Real)

Gives the center-of-mass quantities of a projectile of rest energy `m` and laboratory kinetic
energy `E` colliding with a target of rest energy `M` at rest: momentum (`momentum_cm`,
relativistic), relativistic reduced mass μ_r = W₁W₂/(W₁ + W₂), wave number and Sommerfeld
parameter.

# Input Argument(s)
- `E::Float64` : incident kinetic energy [MeV].
- `m::Float64` : projectile rest energy [MeV].
- `M::Float64` : target rest energy [MeV].
- `Z1Z2::Real` : product of the projectile and target charge numbers.

# Output Argument(s)
- `kin::CMKinematics` : center-of-mass quantities.

# Reference(s)
- Salvat and Quesada, Nucl. Instrum. Methods B 475 (2020) 49, Eqs. (6), (7), (27), (51).
"""
function cm_kinematics(E::Float64, m::Float64, M::Float64, Z1Z2::Real)
    p0c = momentum_cm(E, m, M; kinematics="relativistic")
    W1 = sqrt(m^2 + p0c^2)
    W2 = sqrt(M^2 + p0c^2)
    mur = W1 * W2 / (W1 + W2)
    hbarc = HBAR_eVs * C_m_per_s * 1.0e9               # ħc [MeV fm]
    return CMKinematics(p0c, mur, p0c / hbarc, Z1Z2 * ALPHA_FS * mur / p0c)
end
