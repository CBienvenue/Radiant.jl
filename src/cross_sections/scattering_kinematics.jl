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
Relativistic scattered projectile energy for elastic scattering.
Assumes mu = cos(theta_cm) and target at rest.

Used as the forward map for mu inversion in
`mu_cm_from_Eout_relativistic` when computing mu bounds.
"""
function scattered_energy_relativistic(E_in_mev::Float64, mu::Float64,
                                       m_mev::Float64, M_mev::Float64)::Float64

    E1 = E_in_mev + m_mev
    p1 = sqrt(max(E1^2 - m_mev^2, 0.0))
    s = m_mev^2 + M_mev^2 + 2.0 * M_mev * E1
    sqrt_s = sqrt(s)

    E1_star = (s + m_mev^2 - M_mev^2) / (2.0 * sqrt_s)
    p_star = sqrt(max(E1_star^2 - m_mev^2, 0.0))

    beta_cm = p1 / (E1 + M_mev)
    gamma_cm = 1.0 / sqrt(1.0 - beta_cm^2)

    E1_prime = gamma_cm * (E1_star + beta_cm * p_star * mu)
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
    E1 = E_in_mev + m_mev
    p1 = sqrt(max(E1^2 - m_mev^2, 0.0))
    s = m_mev^2 + M_mev^2 + 2.0 * M_mev * E1
    sqrt_s = sqrt(s)
    E1_star = (s + m_mev^2 - M_mev^2) / (2.0 * sqrt_s)
    p_star = sqrt(max(E1_star^2 - m_mev^2, 0.0))
    beta_cm = p1 / (E1 + M_mev)
    gamma_cm = 1.0 / sqrt(1.0 - beta_cm^2)

    p_par = gamma_cm * (p_star * mu_cm + beta_cm * E1_star)
    p_perp = p_star * sqrt(max(1.0 - mu_cm^2, 0.0))
    p_lab = sqrt(p_par^2 + p_perp^2)
    return p_lab > 0.0 ? p_par / p_lab : 1.0
end
