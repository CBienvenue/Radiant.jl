const HBAR_eVs = 6.582119569e-16
const C_m_per_s = 299_792_458.0
const ALPHA_FS = 1.0 / 137.035999177

# LTP marker for isotope entries without nuclear data (pure Coulomb scattering).
# ENDF LAW=5 only defines LTP = 1, 2 and > 2, so 0 never collides with library data.
const LTP_COULOMB_ONLY = 0

# Projectile spins (ENDF SPI) by ZAP, used for the identical-particle Coulomb term.
const COULOMB_PROJECTILE_SPIN = Dict(1001 => 0.5, 1002 => 1.0, 1003 => 0.5, 2004 => 0.0)

"""
    warn_outside_elastic_tables(iso::IsotopeElasticDB, Z::Int, A::Int, E_in_MeV::Float64)

Warns, once per isotope and side, when an incident energy lies outside the energies of
the ENDF elastic tables. The nuclear part of the cross-section is then the one of the
nearest tabulated energy, while the Coulomb part follows the incident energy.

# Input Argument(s)
- `iso::IsotopeElasticDB` : isotope elastic-scattering data.
- `Z::Int` : target atomic number.
- `A::Int` : target mass number.
- `E_in_MeV::Float64` : incident particle kinetic energy [MeV].

# Output Argument(s)
N/A
"""
function warn_outside_elastic_tables(iso::IsotopeElasticDB, Z::Int, A::Int, E_in_MeV::Float64)
    isempty(iso.energies) && return nothing
    E_lo, E_hi = iso.energies[1], iso.energies[end]
    if E_in_MeV > E_hi
        @warn "Elastic scattering on Z=$(Z), A=$(A): the ENDF tables end at $(E_hi) MeV. Above (e.g. $(E_in_MeV) MeV), the nuclear part is held at its $(E_hi) MeV value and only the Coulomb part follows the energy." _id=Symbol("elastic_tables_above_", Z, "_", A) maxlog=1
    elseif E_in_MeV < E_lo
        @warn "Elastic scattering on Z=$(Z), A=$(A): the ENDF tables start at $(E_lo) MeV. Below (e.g. $(E_in_MeV) MeV), the nuclear part is held at its $(E_lo) MeV value and only the Coulomb part follows the energy." _id=Symbol("elastic_tables_below_", Z, "_", A) maxlog=1
    end
    return nothing
end

"""
    prepare_elastic_dsig_cache(db::ElasticScatteringENDFDB, Z::Int, A::Int,
    E_in_MeV::Float64)

Prepares the energy-resolved elastic differential cross-section cache for one isotope and
one incident energy: the Coulomb parameters at that energy and the LAW=5 coefficients
(LTP=1, LTP=2) or the tabulated cross-section Coulomb + σ_NI·P_NI on the ENDF cosines
(LTP>2), interpolated in energy between the ENDF tables. Outside the tables, the nuclear
part of the nearest table is used (`warn_outside_elastic_tables`).

# Input Argument(s)
- `db::ElasticScatteringENDFDB` : elastic-scattering database containing ENDF data by
  isotope, with its Coulomb kinematics convention and screening setting.
- `Z::Int` : target atomic number.
- `A::Int` : target mass number.
- `E_in_MeV::Float64` : incident particle kinetic energy in MeV.

# Output Argument(s)
- `cache::ElasticDCSCache` : cache containing the data needed to evaluate the elastic
  differential cross-section at `E_in_MeV`.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering.
"""
function prepare_elastic_dsig_cache(db::ElasticScatteringENDFDB, Z::Int, A::Int, E_in_MeV::Float64)
    iso = db.isotopes[(Z, A)]
    E_eV = E_in_MeV * 1.0e6
    coulomb_kinematics = db.coulomb_kinematics
    m_MeV = iso.m_mev
    M_MeV = iso.M_mev
    m_eV = m_MeV * 1.0e6
    M_eV = M_MeV * 1.0e6
    k_c = k_coulomb(E_eV, m_eV, M_eV; kinematics=coulomb_kinematics)
    eta = eta_coulomb(E_eV, m_eV, iso.Z1, iso.Z2; kinematics=coulomb_kinematics)
    delta = delta_coulomb(E_eV, m_eV, M_eV, iso.Z1, iso.Z2; kinematics=coulomb_kinematics, screening=db.screening)

    if iso.ltp == LTP_COULOMB_ONLY
        return ElasticDCSCache(iso.ltp, E_in_MeV, k_c, eta, delta,
                               Float64[], Float64[], 0,
                               ComplexF64[], Float64[], Float64[],
                               iso.spin_s, m_MeV, M_MeV, iso.identical)
    end
    warn_outside_elastic_tables(iso, Z, A, E_in_MeV)

    if iso.ltp == 1
        coeffs = iso.ba_tables_ltp1
        coeffs === nothing && error("Missing LTP=1 coefficient tables.")
        energies = sort(collect(keys(coeffs)))

        local NL
        local b
        local a
        if length(energies) == 1 || E_eV <= energies[1]
            NL, b, a = coeffs[energies[1]]
        elseif E_eV >= energies[end]
            NL, b, a = coeffs[energies[end]]
        else
            idx = searchsortedlast(energies, E_eV)
            E1 = energies[idx]
            E2 = energies[idx + 1]
            NL1, b1, a1 = coeffs[E1]
            NL2, b2, a2 = coeffs[E2]
            itp = select_interpolation(iso.mf6_nbt, iso.mf6_int, idx)

            NL = min(NL1, NL2)
            nb = min(length(b1), length(b2))
            na = min(length(a1), length(a2))

            b = Vector{Float64}(undef, nb)
            for i in 1:nb
                b[i] = interp_pair(E_eV, E1, b1[i], E2, b2[i], itp)
            end

            a = Vector{ComplexF64}(undef, na)
            for i in 1:na
                ar = interp_pair(E_eV, E1, real(a1[i]), E2, real(a2[i]), itp)
                ai = interp_pair(E_eV, E1, imag(a1[i]), E2, imag(a2[i]), itp)
                a[i] = complex(ar, ai)
            end
        end

        return ElasticDCSCache(iso.ltp, E_in_MeV, k_c, eta, delta,
                               Float64[], Float64[], NL, a, b, Float64[],
                               iso.spin_s, m_MeV, M_MeV, iso.identical)
    elseif iso.ltp == 2
        coeffs = iso.c_tables_ltp2
        coeffs === nothing && error("Missing LTP=2 coefficient tables.")
        energies = sort(collect(keys(coeffs)))

        local NL
        local c
        if length(energies) == 1 || E_eV <= energies[1]
            NL, c = coeffs[energies[1]]
        elseif E_eV >= energies[end]
            NL, c = coeffs[energies[end]]
        else
            idx = searchsortedlast(energies, E_eV)
            E1 = energies[idx]
            E2 = energies[idx + 1]
            NL1, c1 = coeffs[E1]
            NL2, c2 = coeffs[E2]
            itp = select_interpolation(iso.mf6_nbt, iso.mf6_int, idx)

            NL = min(NL1, NL2)
            nc = min(length(c1), length(c2))
            c = Vector{Float64}(undef, nc)
            for i in 1:nc
                c[i] = interp_pair(E_eV, E1, c1[i], E2, c2[i], itp)
            end
        end

        return ElasticDCSCache(iso.ltp, E_in_MeV, k_c, eta, delta,
                               Float64[], Float64[], NL, ComplexF64[], Float64[], c,
                               iso.spin_s, m_MeV, M_MeV, iso.identical)
    elseif iso.ltp > 2
        tables = iso.mupni_tables
        tables === nothing && error("Missing LTP>2 P_NI tables.")
        if iso.sigmaNI_E === nothing || iso.sigmaNI_S === nothing ||
           iso.sigmaNI_nbt === nothing || iso.sigmaNI_int === nothing
            error("Missing sigmaNI tables.")
        end

        energies = sort(collect(keys(tables)))
        common_grid = tables isa Dict{Float64,Vector{Float64}}
        table_at(Ei) = common_grid ? (iso.mu_grid, tables[Ei]) : tables[Ei]
        local mu_tab
        local pni_tab
        if length(energies) == 1 || E_eV <= energies[1]
            mu_tab, pni_tab = table_at(energies[1])
        elseif E_eV >= energies[end]
            mu_tab, pni_tab = table_at(energies[end])
        else
            idx = searchsortedlast(energies, E_eV)
            E1 = energies[idx]
            E2 = energies[idx + 1]
            itp = select_interpolation(iso.mf6_nbt, iso.mf6_int, idx)
            mu1, pni1 = table_at(E1)
            mu2, pni2 = table_at(E2)
            if length(mu1) != length(mu2) || any(i -> mu1[i] != mu2[i], eachindex(mu1))
                error("Energy interpolation for LTP>2 requires identical mu grids at the bracketing energies.")
            end
            mu_tab = mu1
            pni_tab = [interp_pair(E_eV, E1, pni1[i], E2, pni2[i], itp) for i in eachindex(mu_tab)]
        end

        sigmaNI = interp_TAB1(E_eV, iso.sigmaNI_E, iso.sigmaNI_S,
                              iso.sigmaNI_nbt, iso.sigmaNI_int)
        sigma_grid = [elastic_dsig_ltpgreater2(mu_tab[k]; pni=pni_tab[k], sigmaNI=sigmaNI, identical=iso.identical, spin_s=iso.spin_s, k=k_c, eta=eta, delta=delta)
                      for k in eachindex(mu_tab)]

        return ElasticDCSCache(iso.ltp, E_in_MeV, k_c, eta, delta,
                               copy(mu_tab), sigma_grid, 0,
                               ComplexF64[], Float64[], Float64[],
                               iso.spin_s, m_MeV, M_MeV, iso.identical)
    end
    error("Unsupported LTP=$(iso.ltp) for the elastic DCS cache.")
end

"""
    elastic_dsig_raw(cache::ElasticDCSCache, mu_cm::Float64)

Gives the elastic differential cross-section from an energy-resolved cache without clipping
negative numerical artifacts. For tabulated LTP>2 data evaluated beyond the stored
forward angular grid, the Coulomb contribution is evaluated directly from the cached
`k`, `eta`, `delta`, identical-particle flag, and ENDF projectile spin.

Between the tabulated cosines, the total Coulomb + σ_NI·P_NI is interpolated linearly.
ENDF-102 (LTP=12) defines P_NI as linear in μ with the Coulomb term exact; that
representation dips below zero between the forward-most nodes of the TENDL-2023 tables
(e.g. Al, Si at 100 MeV, C at 50-123 MeV), where the destructive interference is strong
and the nodes are coarse, so the total is interpolated instead. It keeps the
cross-section positive between positive nodes; above the last tabulated cosine,
P_NI = 0 (ENDF-102 eq. 6.18).

# Input Argument(s)
- `cache::ElasticDCSCache` : cache generated by `prepare_elastic_dsig_cache`.
- `mu_cm::Float64` : center-of-mass scattering cosine.

# Output Argument(s)
- `σs::Float64` : elastic differential cross-section in barn/sr.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering.
"""
function elastic_dsig_raw(cache::ElasticDCSCache, mu_cm::Float64)
    if cache.ltp == LTP_COULOMB_ONLY
        return coulomb_dcs(mu_cm, cache.k_coulomb, cache.eta_coulomb, cache.delta_coulomb, cache.identical, cache.spin_s)
    elseif is_tabulated(cache)
        if mu_cm > cache.mu_grid[end]
            return coulomb_dcs(mu_cm, cache.k_coulomb, cache.eta_coulomb, cache.delta_coulomb, cache.identical, cache.spin_s)
        end
        return interp_linear(cache.mu_grid, cache.sigma_grid, mu_cm)
    elseif cache.ltp == 1
        return elastic_dsig_ltp1(mu_cm; a=cache.a, b=cache.b, NL=cache.NL, identical=cache.identical, spin_s=cache.spin_s, k=cache.k_coulomb, eta=cache.eta_coulomb, delta=cache.delta_coulomb)
    elseif cache.ltp == 2
        return elastic_dsig_ltp2(mu_cm; c=cache.c, NL=cache.NL, identical=cache.identical, spin_s=cache.spin_s, k=cache.k_coulomb, eta=cache.eta_coulomb, delta=cache.delta_coulomb)
    else
        error("Elastic DCS cache for LTP=$(cache.ltp) is not evaluable.")
    end
end

"""
    elastic_dsig(cache::ElasticDCSCache, mu_cm::Float64)

Gives the non-negative elastic differential cross-section from an energy-resolved cache.

# Input Argument(s)
- `cache::ElasticDCSCache` : cache generated by `prepare_elastic_dsig_cache`.
- `mu_cm::Float64` : center-of-mass scattering cosine.

# Output Argument(s)
- `σs::Float64` : elastic differential cross-section in barn/sr, clipped to zero when
  roundoff produces a small negative value.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering.
"""
function elastic_dsig(cache::ElasticDCSCache, mu_cm::Float64)
    return max(0.0, elastic_dsig_raw(cache, mu_cm))
end

"""
    integrate_elastic_scattering_numerical_GK(cache::ElasticDCSCache, Ec_MeV::Float64)

Gives the total elastic cross-section above an outgoing-energy cutoff using adaptive
Gauss-Kronrod quadrature over the cached center-of-mass angular distribution. Kept for
debugging: the total cross-section is the sum of the group moments of order zero plus
`integrate_elastic_scattering_below`.

# Input Argument(s)
- `cache::ElasticDCSCache` : cache generated by `prepare_elastic_dsig_cache`.
- `Ec_MeV::Float64` : cutoff energy between soft and catastrophic events.

# Output Argument(s)
- `σt::Float64` : integrated elastic cross-section in barn.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering.
"""
function integrate_elastic_scattering_numerical_GK(cache::ElasticDCSCache, Ec_MeV::Float64)
    E_in_MeV = cache.E_in_mev
    E_out_MeV = max(Ec_MeV, 0.0)
    mu_b = mu_cm_from_Eout_relativistic(E_in_MeV, E_out_MeV, cache.m_mev, cache.M_mev)

    mu_lo = cache.identical ? 0.0 : -1.0
    if mu_b <= mu_lo
        return 0.0
    end

    a = mu_lo
    b = mu_b
    if a >= b
        return 0.0
    end

    quadgk = require_quadgk()
    f = mu -> elastic_dsig(cache, mu)
    val, _ = quadgk.quadgk(f, a, b; rtol=1.0e-10)
    return 2.0 * pi * val
end

"""
    is_tabulated(cache::ElasticDCSCache)

Tells whether the elastic cross-section of a cache is tabulated in angle (ENDF LTP>2),
linear between the tabulated cosines, rather than reconstructed from LTP=1 or LTP=2
coefficients or given by the Coulomb term alone.

# Input Argument(s)
- `cache::ElasticDCSCache` : elastic DCS cache.

# Output Argument(s)
- `tabulated::Bool` : boolean indicating a tabulated cross-section.
"""
is_tabulated(cache::ElasticDCSCache) = cache.ltp > 2

"""
    coulomb_dcs_d(cache::ElasticDCSCache, d::Float64)

Gives the Coulomb differential cross-section of `coulomb_dcs` as a function of
d = 1 - μ + δ instead of the cosine μ. Near μ = 1, d is of the order of the screening
shift δ (about 1e-10 for protons), and recomputing it from μ would leave only about six
significant digits; integrations in ln d (`integrate_in_log_cosine`) know d exactly.

# Input Argument(s)
- `cache::ElasticDCSCache` : elastic DCS cache providing `k`, `eta`, `delta`, the
  identical-particle flag and the projectile spin.
- `d::Float64` : 1 - μ + δ, with μ the center-of-mass scattering cosine.

# Output Argument(s)
- `σc::Float64` : Coulomb differential cross-section in barn/sr.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering, eqs. 6.7, 6.8.
"""
function coulomb_dcs_d(cache::ElasticDCSCache, d::Float64)
    eta2k2 = cache.eta_coulomb^2 / cache.k_coulomb^2
    if !cache.identical
        return eta2k2 / d^2
    end
    delta = cache.delta_coulomb
    dplus = 2.0 * (1.0 + delta) - d
    mu = 1.0 + delta - d
    denom = d * dplus
    spin_factor = ((-1.0)^(2.0 * cache.spin_s)) / (2.0 * cache.spin_s + 1.0)
    bracket = ((1.0 + delta)^2 + mu^2) / denom + spin_factor * cos(cache.eta_coulomb * log(dplus / d))
    return 2.0 * eta2k2 / denom * bracket
end

"""
    elastic_dsig_raw_d(cache::ElasticDCSCache, mu::Float64, d::Float64)

Gives the unclipped elastic differential cross-section of `elastic_dsig_raw`, with the
factors 1 - μ + δ of its Coulomb and interference terms taken as d (`coulomb_dcs_d`,
`elastic_dsig_nuclear`). Between the tabulated
cosines of an LTP>2 cross-section, the tabulation is interpolated in μ as in
`elastic_dsig_raw`.

# Input Argument(s)
- `cache::ElasticDCSCache` : elastic DCS cache.
- `mu::Float64` : center-of-mass scattering cosine.
- `d::Float64` : 1 - μ + δ.

# Output Argument(s)
- `σs::Float64` : elastic differential cross-section in barn/sr.
"""
function elastic_dsig_raw_d(cache::ElasticDCSCache, mu::Float64, d::Float64)
    if cache.ltp == LTP_COULOMB_ONLY
        return coulomb_dcs_d(cache, d)
    elseif is_tabulated(cache)
        return mu > cache.mu_grid[end] ? coulomb_dcs_d(cache, d) : interp_linear(cache.mu_grid, cache.sigma_grid, mu)
    end
    return coulomb_dcs_d(cache, d) + elastic_dsig_nuclear(cache, mu, d)
end

"""
    coulomb_dcs_integral(cache::ElasticDCSCache, mu_lo::Float64, mu_hi::Float64)

Gives the integral over the center-of-mass cosine of the Coulomb differential
cross-section. The direct terms are integrated in closed form; the identical-particle
interference term is integrated by adaptive quadrature in ln(1 - μ + δ).

# Input Argument(s)
- `cache::ElasticDCSCache` : elastic DCS cache providing `k`, `eta`, `delta`, the
  identical-particle flag and the projectile spin.
- `mu_lo::Float64` : lower center-of-mass cosine bound.
- `mu_hi::Float64` : upper center-of-mass cosine bound.

# Output Argument(s)
- `σc::Float64` : ∫ dσ_C/dΩ dμ over [mu_lo, mu_hi], in barn.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering.
"""
function coulomb_dcs_integral(cache::ElasticDCSCache, mu_lo::Float64, mu_hi::Float64)
    mu_hi <= mu_lo && return 0.0
    eta = cache.eta_coulomb
    k = cache.k_coulomb
    delta = cache.delta_coulomb
    if !cache.identical
        return (eta^2 / k^2) * (1.0 / (1.0 - mu_hi + delta) - 1.0 / (1.0 - mu_lo + delta))
    end

    # Antiderivative μ / ((1+δ)² - μ²), with the denominator factorized: at μ = 1,
    # (1+δ)² - 1 would lose about half of the significant digits for δ ~ 1e-10.
    F(mu) = mu / ((1.0 - mu + delta) * (1.0 + mu + delta))
    sigma_c = 2.0 * (eta^2 / k^2) * (F(mu_hi) - F(mu_lo))

    spin_factor = ((-1.0)^(2.0 * cache.spin_s)) / (2.0 * cache.spin_s + 1.0)
    interference(mu, d) = begin
        dplus = 2.0 * (1.0 + delta) - d
        2.0 * eta^2 / (k^2 * dplus * d) * spin_factor * cos(eta * log(dplus / d))
    end
    return sigma_c + integrate_in_log_cosine(interference, cache, mu_lo, mu_hi)
end

"""
    tabulated_segment_values(cache::ElasticDCSCache, x0::Float64, x1::Float64)

Gives the values at both ends of [x0, x1] of the linear interpolation of the tabulated
elastic differential cross-section. The interval must not straddle a tabulated cosine:
the line used is the one of the tabulation interval containing its midpoint. Below the
first tabulated cosine the first value is held, as `interp_linear` does.

# Input Argument(s)
- `cache::ElasticDCSCache` : elastic DCS cache with a tabulated `sigma_grid`.
- `x0::Float64` : lower center-of-mass cosine.
- `x1::Float64` : upper center-of-mass cosine.

# Output Argument(s)
- `y0::Float64` : differential cross-section at `x0` [barn/sr].
- `y1::Float64` : differential cross-section at `x1` [barn/sr].
"""
function tabulated_segment_values(cache::ElasticDCSCache, x0::Float64, x1::Float64)
    x = cache.mu_grid
    y = cache.sigma_grid
    xm = 0.5 * (x0 + x1)
    if xm <= x[1]
        return y[1], y[1]
    end
    j = clamp(searchsortedlast(x, xm), 1, length(x) - 1)
    slope = (y[j + 1] - y[j]) / (x[j + 1] - x[j])
    return y[j] + slope * (x0 - x[j]), y[j] + slope * (x1 - x[j])
end

"""
    positive_part_linear(x0::Float64, x1::Float64, y0::Float64, y1::Float64)

Gives the part of [x0, x1] on which the linear function through (x0, y0) and (x1, y1) is
positive. When the function changes sign inside the interval, the interval is cut at
its zero.

# Input Argument(s)
- `x0::Float64`, `x1::Float64` : interval bounds.
- `y0::Float64`, `y1::Float64` : function values at the bounds.

# Output Argument(s)
- `piece::Union{Nothing,NTuple{4,Float64}}` : `(a, b, ya, yb)`, the positive sub-interval
  and the function values at its ends, or `nothing` when the function is nowhere positive.
"""
function positive_part_linear(x0::Float64, x1::Float64, y0::Float64, y1::Float64)
    if y0 >= 0.0 && y1 >= 0.0
        return (y0 == 0.0 && y1 == 0.0) ? nothing : (x0, x1, y0, y1)
    elseif y0 <= 0.0 && y1 <= 0.0
        return nothing
    end
    x_zero = x0 + (x1 - x0) * y0 / (y0 - y1)
    return y0 > 0.0 ? (x0, x_zero, y0, 0.0) : (x_zero, x1, 0.0, y1)
end

"""
    integrate_tabulated_elastic_dcs(cache::ElasticDCSCache, mu_lo::Float64, mu_hi::Float64)

Gives the exact integral over the center-of-mass cosine of the tabulated elastic
differential cross-section, as `elastic_dsig` evaluates it: linear between the tabulated
cosines, with negative values set to zero, and Coulomb above the last tabulated cosine.
Each tabulation interval where the cross-section changes sign is cut at its zero and only
the positive part is kept; the Coulomb part is integrated in closed form.

# Input Argument(s)
- `cache::ElasticDCSCache` : elastic DCS cache with a tabulated `sigma_grid`.
- `mu_lo::Float64` : lower center-of-mass cosine bound.
- `mu_hi::Float64` : upper center-of-mass cosine bound.

# Output Argument(s)
- `σ::Float64` : ∫ dσ/dΩ dμ over [mu_lo, mu_hi], in barn.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering.
"""
function integrate_tabulated_elastic_dcs(cache::ElasticDCSCache, mu_lo::Float64, mu_hi::Float64)
    x = cache.mu_grid
    mu_end = x[end]
    acc = 0.0
    top = min(mu_hi, mu_end)
    if top > mu_lo
        cuts = [mu_lo; filter(xi -> mu_lo < xi < top, x); top]
        for k in 1:(length(cuts) - 1)
            x0, x1 = cuts[k], cuts[k + 1]
            y0, y1 = tabulated_segment_values(cache, x0, x1)
            piece = positive_part_linear(x0, x1, y0, y1)
            isnothing(piece) && continue
            a, b, ya, yb = piece
            acc += 0.5 * (b - a) * (ya + yb)
        end
    end
    if mu_hi > mu_end
        acc += coulomb_dcs_integral(cache, max(mu_lo, mu_end), mu_hi)
    end
    return acc
end

"""
    cm_frame(cache::ElasticDCSCache)

Gives the center-of-mass frame of an elastic collision (relativistic two-body
kinematics, target at rest).

# Input Argument(s)
- `cache::ElasticDCSCache` : elastic DCS cache providing the incident energy and masses.

# Output Argument(s)
- `frame::NamedTuple` : `gamma` and `beta` of the center of mass, and the momentum
  `p_star` [MeV/c] and total energy `E1_star` [MeV] of the projectile in it.
"""
function cm_frame(cache::ElasticDCSCache)
    m_mev = cache.m_mev
    M_mev = cache.M_mev
    E1 = cache.E_in_mev + m_mev
    p1 = sqrt(max(E1^2 - m_mev^2, 0.0))
    s = m_mev^2 + M_mev^2 + 2.0 * M_mev * E1
    E1_star = (s + m_mev^2 - M_mev^2) / (2.0 * sqrt(s))
    p_star = sqrt(max(E1_star^2 - m_mev^2, 0.0))
    beta = p1 / (E1 + M_mev)
    return (gamma=1.0 / sqrt(1.0 - beta^2), beta=beta, p_star=p_star, E1_star=E1_star)
end

"""
    scattered_energy_affine(cache::ElasticDCSCache)

Gives the coefficients of the scattered-projectile energy, which is affine in the
center-of-mass cosine for elastic scattering: E_out = e0 + e1 μ_cm. Since E_out = E_in at
μ_cm = 1, the energy lost is exactly e1 (1 - μ_cm).

# Input Argument(s)
- `cache::ElasticDCSCache` : elastic DCS cache providing the incident energy and masses.

# Output Argument(s)
- `e0::Float64` : scattered energy at μ_cm = 0 [MeV].
- `e1::Float64` : slope of the scattered energy in μ_cm [MeV].
"""
function scattered_energy_affine(cache::ElasticDCSCache)
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

# Gauss-Kronrod (7, 15) rule on [-1, 1] (QUADPACK QK15): Kronrod nodes and weights, and
# the Gauss weights belonging to the even-numbered Kronrod nodes, the last node being 0.
const GK15_NODES = (0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
    0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
    0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
    0.207784955007898467600689403773245, 0.0)
const GK15_WEIGHTS = (0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
    0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
    0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
    0.204432940075298892414161999234649, 0.209482141084727828012999174891714)
const G7_WEIGHTS = (0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
    0.381830050505118944950369775488975, 0.417959183673469387755102040816327)

"""
    gauss_kronrod_15(f::Function, a::Float64, b::Float64)

Applies the Gauss-Kronrod (7, 15) rule to `f` on [a, b].

# Input Argument(s)
- `f::Function` : integrand, returning a scalar or a vector.
- `a::Float64` : lower bound.
- `b::Float64` : upper bound.

# Output Argument(s)
- `K` : 15-point Kronrod estimate.
- `G` : 7-point Gauss estimate.
- `K_abs::Float64` : 15-point Kronrod estimate of the integral of the absolute value of
  the first component of `f`.
"""
function gauss_kronrod_15(f::Function, a::Float64, b::Float64)
    c = 0.5 * (a + b)
    h = 0.5 * (b - a)
    fc = f(c)
    K = GK15_WEIGHTS[8] * fc
    G = G7_WEIGHTS[4] * fc
    K_abs = GK15_WEIGHTS[8] * abs(first(fc))
    for j in 1:7
        x = h * GK15_NODES[j]
        f_minus = f(c - x)
        f_plus = f(c + x)
        fs = f_minus + f_plus
        K = K + GK15_WEIGHTS[j] * fs
        K_abs += GK15_WEIGHTS[j] * (abs(first(f_minus)) + abs(first(f_plus)))
        iseven(j) && (G = G + G7_WEIGHTS[j ÷ 2] * fs)
    end
    return K * h, G * h, K_abs * abs(h)
end

"""
    adaptive_gauss_kronrod(f::Function, a::Float64, b::Float64; rtol::Float64=1.0e-13,
    max_intervals::Int=1000)

Integrates `f` on [a, b] by globally adaptive Gauss-Kronrod (7, 15) quadrature: the
interval with the largest error estimate is bisected until the summed estimate falls
below `rtol` times the integral of the absolute value of the integrand, or until
`max_intervals` intervals are reached, which bounds the work where round-off in the
integrand keeps the estimate above the tolerance. The error is controlled on the first
component when `f` returns a vector.

# Input Argument(s)
- `f::Function` : integrand, returning a scalar or a vector.
- `a::Float64` : lower bound.
- `b::Float64` : upper bound.
- `rtol::Float64` : relative tolerance.
- `max_intervals::Int` : maximum number of intervals.

# Output Argument(s)
- `I` : integral, of the type returned by `f`.
"""
function adaptive_gauss_kronrod(f::Function, a::Float64, b::Float64; rtol::Float64=1.0e-13,
                                max_intervals::Int=1000)
    K, G, A = gauss_kronrod_15(f, a, b)
    lo = [a]; hi = [b]; Ks = [K]; errs = [abs(first(K) - first(G))]; As = [A]
    while sum(errs) > rtol * sum(As) && length(Ks) < max_intervals
        i = argmax(errs)
        mid = 0.5 * (lo[i] + hi[i])
        K_l, G_l, A_l = gauss_kronrod_15(f, lo[i], mid)
        K_r, G_r, A_r = gauss_kronrod_15(f, mid, hi[i])
        push!(lo, mid); push!(hi, hi[i]); push!(Ks, K_r); push!(errs, abs(first(K_r) - first(G_r))); push!(As, A_r)
        hi[i] = mid; Ks[i] = K_l; errs[i] = abs(first(K_l) - first(G_l)); As[i] = A_l
    end
    return sum(Ks)
end

"""
    integrate_in_log_cosine(f::Function, cache::ElasticDCSCache, mu_lo::Float64,
    mu_hi::Float64; rtol::Float64=1.0e-13)

Integrates `f` over the center-of-mass cosine on [mu_lo, mu_hi] in the variable
t = ln d, d = 1 - μ + δ, in which the Coulomb rise towards μ = 1 is smooth, by adaptive
Gauss-Kronrod quadrature. The integrand is called as `f(μ, d)`, with d = e^t exact, so
that it can evaluate the Coulomb term from d (`coulomb_dcs_d`). Without screening
(δ = 0), the range is stopped at μ = 1 - 1e-12.

# Input Argument(s)
- `f::Function` : integrand `f(μ, d)`, returning a scalar or a vector.
- `cache::ElasticDCSCache` : elastic DCS cache providing the screening shift δ.
- `mu_lo::Float64` : lower center-of-mass cosine bound.
- `mu_hi::Float64` : upper center-of-mass cosine bound.
- `rtol::Float64` : relative tolerance.

# Output Argument(s)
- `I` : ∫ f dμ over [mu_lo, mu_hi], of the type returned by `f`.
"""
function integrate_in_log_cosine(f::Function, cache::ElasticDCSCache, mu_lo::Float64,
                                 mu_hi::Float64; rtol::Float64=1.0e-13)
    delta = cache.delta_coulomb
    mu_top = delta > 0.0 ? mu_hi : min(mu_hi, 1.0 - 1.0e-12)
    t_lo = log(1.0 - mu_top + delta)
    t_hi = log(1.0 - mu_lo + delta)
    g(t) = (d = exp(t); f(1.0 + delta - d, d) * d)
    return adaptive_gauss_kronrod(g, t_lo, t_hi; rtol=rtol)
end

"""
    elastic_dsig_nuclear(cache::ElasticDCSCache, mu_cm::Float64, d::Float64=NaN)

Gives the nuclear and Coulomb-nuclear interference part of the elastic differential
cross-section reconstructed from ENDF LTP=1 or LTP=2 coefficients, that is the
cross-section without its Coulomb term.

# Input Argument(s)
- `cache::ElasticDCSCache` : coefficient-based LTP=1 or LTP=2 elastic DCS cache.
- `mu_cm::Float64` : center-of-mass scattering cosine.
- `d::Float64` : 1 - μ + δ when known exactly; computed from `mu_cm` when NaN.

# Output Argument(s)
- `σn::Float64` : nuclear and interference differential cross-section in barn/sr.
"""
function elastic_dsig_nuclear(cache::ElasticDCSCache, mu_cm::Float64, d::Float64=NaN)
    if cache.ltp == 1
        return elastic_dsig_ltp1(mu_cm; a=cache.a, b=cache.b, NL=cache.NL, identical=cache.identical, spin_s=cache.spin_s, k=cache.k_coulomb, eta=cache.eta_coulomb, delta=cache.delta_coulomb, nuclear_only=true, d=d)
    elseif cache.ltp == 2
        return elastic_dsig_ltp2(mu_cm; c=cache.c, NL=cache.NL, identical=cache.identical, spin_s=cache.spin_s, k=cache.k_coulomb, eta=cache.eta_coulomb, delta=cache.delta_coulomb, nuclear_only=true, d=d)
    end
    return 0.0
end

"""
    elastic_dsig_zero_crossings(cache::ElasticDCSCache)

Gives the center-of-mass cosines at which the elastic differential cross-section
reconstructed from ENDF LTP=1 or LTP=2 coefficients changes sign. The cross-section is
sampled over its whole cosine range, uniformly in μ and in ln(1 - μ), always on the same
points, and every sign change is refined by bisection. The result depends on the cache
alone, so that integrals over any sub-range cut at the same cosines.

# Input Argument(s)
- `cache::ElasticDCSCache` : coefficient-based LTP=1 or LTP=2 elastic DCS cache.

# Output Argument(s)
- `zeros::Vector{Float64}` : sign-change cosines, increasing; empty for a Coulomb-only
  cache, whose cross-section is positive.
"""
function elastic_dsig_zero_crossings(cache::ElasticDCSCache)
    roots = Float64[]
    cache.ltp ∈ (1, 2) || return roots
    mu_lo = cache.identical ? 0.0 : -1.0
    mu_top = 1.0 - max(cache.delta_coulomb, 1.0e-12)
    N = 256
    samples = vcat(collect(range(mu_lo, mu_top, length=N)),
                   1.0 .- exp.(range(log(1.0 - mu_lo), log(1.0 - mu_top), length=N)))
    sort!(unique!(samples))
    f(mu) = elastic_dsig_raw(cache, mu)
    x_prev = samples[1]
    f_prev = f(x_prev)
    for x in samples[2:end]
        fx = f(x)
        if (f_prev < 0.0) != (fx < 0.0)
            lo, hi, f_lo = x_prev, x, f_prev
            for _ in 1:200
                mid = 0.5 * (lo + hi)
                (mid <= lo || mid >= hi) && break
                f_mid = f(mid)
                if (f_mid < 0.0) == (f_lo < 0.0)
                    lo, f_lo = mid, f_mid
                else
                    hi = mid
                end
            end
            push!(roots, 0.5 * (lo + hi))
        end
        x_prev, f_prev = x, fx
    end
    return roots
end

"""
    integrate_smooth_elastic_dcs(cache::ElasticDCSCache, mu_lo::Float64, mu_hi::Float64)

Gives the integral over the center-of-mass cosine of the elastic differential
cross-section reconstructed from ENDF LTP=1 or LTP=2 coefficients, or of the Coulomb
cross-section alone, with negative values set to zero, as `elastic_dsig` evaluates it.
The Coulomb term is integrated in closed form and the nuclear and interference terms by
adaptive quadrature in ln(1 - μ + δ). Where the cross-section is negative, between the
cosines from `elastic_dsig_zero_crossings`, its negative part is added back, so that the
result equals the sum of the group moments of `elastic_group_moments!`.

# Input Argument(s)
- `cache::ElasticDCSCache` : coefficient-based LTP=1, LTP=2 or Coulomb-only DCS cache.
- `mu_lo::Float64` : lower center-of-mass cosine bound.
- `mu_hi::Float64` : upper center-of-mass cosine bound.

# Output Argument(s)
- `σ::Float64` : ∫ max(dσ/dΩ, 0) dμ over [mu_lo, mu_hi], in barn.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5, LTP=1 and LTP=2.
"""
function integrate_smooth_elastic_dcs(cache::ElasticDCSCache, mu_lo::Float64, mu_hi::Float64)
    mu_hi <= mu_lo && return 0.0
    acc = coulomb_dcs_integral(cache, mu_lo, mu_hi)
    cache.ltp ∈ (1, 2) || return acc
    cuts = [mu_lo; filter(z -> mu_lo < z < mu_hi, elastic_dsig_zero_crossings(cache)); mu_hi]
    for k in 1:(length(cuts) - 1)
        x0, x1 = cuts[k], cuts[k + 1]
        acc += integrate_in_log_cosine((mu, d) -> elastic_dsig_nuclear(cache, mu, d), cache, x0, x1)
        if elastic_dsig_raw(cache, 0.5 * (x0 + x1)) < 0.0
            acc += integrate_in_log_cosine((mu, d) -> max(-elastic_dsig_raw_d(cache, mu, d), 0.0), cache, x0, x1)
        end
    end
    return acc
end

"""
    elastic_cuts(cache::ElasticDCSCache, mu_lo::Float64, mu_hi::Float64)

Gives the cosines at which an integral of the elastic cross-section over [mu_lo, mu_hi]
is cut so that the cross-section is smooth on every piece: the tabulated cosines (LTP>2),
including the last one above which it is Coulomb, or the cosines where the
reconstructed cross-section changes sign (LTP=1, LTP=2).

# Input Argument(s)
- `cache::ElasticDCSCache` : elastic DCS cache.
- `mu_lo::Float64` : lower center-of-mass cosine bound.
- `mu_hi::Float64` : upper center-of-mass cosine bound.

# Output Argument(s)
- `cuts::Vector{Float64}` : `mu_lo`, the cuts strictly inside, and `mu_hi`, increasing.
"""
function elastic_cuts(cache::ElasticDCSCache, mu_lo::Float64, mu_hi::Float64)
    inner = is_tabulated(cache) ? cache.mu_grid : elastic_dsig_zero_crossings(cache)
    return [mu_lo; filter(x -> mu_lo < x < mu_hi, inner); mu_hi]
end

"""
    elastic_group_moments!(σl::AbstractMatrix{Float64}, σe::AbstractVector{Float64},
    cache::ElasticDCSCache, Eout_mev::Vector{Float64}, mu_lo::Float64, mu_hi::Float64,
    L::Int64; recoil::Bool=false)

Adds, for every outgoing energy group, the Legendre moments and the energy moment of the
elastic cross-section over the center-of-mass cosines [mu_lo, mu_hi] of the scattered
projectile. With `recoil`, they are those of the recoiling target instead, for identical
particles (p + H-1): it leaves at the center-of-mass cosine -μ with the energy the
projectile would have there.

The range is walked once, as for the nuclear-reaction spectra, rather than integrated
group by group. The outgoing energy is affine in the cosine, so the group boundaries map
to cosines in closed form. The range is cut at the group boundaries and at the cuts of
`elastic_cuts`, and every piece is added to the group holding it. Only the positive part
of the cross-section is kept, as in `elastic_dsig`: a tabulated piece where the linear
cross-section changes sign is cut at its zero, and a negative reconstructed piece is
skipped. On a tabulated piece the moment of order zero and the energy moment are exact,
the higher moments by Gauss-Legendre, since the laboratory cosine is not polynomial in
the center-of-mass one. On the Coulomb piece above the tabulation and on reconstructed
pieces, all moments are integrated together by adaptive Gauss-Kronrod in ln(1 - μ + δ)
with the Coulomb term evaluated from 1 - μ + δ (`coulomb_dcs_d`); the moment of order
zero of a pure Coulomb piece is the closed form of `coulomb_dcs_integral`. Pieces outside
the groups, notably below the last group, are skipped.

# Input Argument(s)
- `σl::AbstractMatrix{Float64}` : Legendre moments to add to, size (Ng, L+1), in barn.
- `σe::AbstractVector{Float64}` : energy moments to add to, size Ng, in barn × MeV.
- `cache::ElasticDCSCache` : elastic DCS cache.
- `Eout_mev::Vector{Float64}` : outgoing group boundaries [MeV], decreasing.
- `mu_lo::Float64` : lower center-of-mass cosine bound of the projectile.
- `mu_hi::Float64` : upper center-of-mass cosine bound of the projectile.
- `L::Int64` : Legendre truncation order.
- `recoil::Bool` : boolean to give the moments of the recoiling target.

# Output Argument(s)
N/A

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering.
"""
function elastic_group_moments!(σl::AbstractMatrix{Float64}, σe::AbstractVector{Float64},
        cache::ElasticDCSCache, Eout_mev::Vector{Float64}, mu_lo::Float64, mu_hi::Float64, L::Int64;
        recoil::Bool=false)

    mu_hi <= mu_lo && return nothing
    E_in_mev = cache.E_in_mev
    m_mev = cache.m_mev
    M_mev = cache.M_mev
    Ng = length(Eout_mev) - 1
    e0, e1 = scattered_energy_affine(cache)
    sgn = recoil ? -1.0 : 1.0                 # the particle followed leaves at sgn·μ
    energy(mu) = e0 + sgn * e1 * mu
    lab(mu) = mu_lab_from_mu_cm(sgn * mu, E_in_mev, m_mev, M_mev)
    tabulated = is_tabulated(cache)
    mu_end = tabulated ? cache.mu_grid[end] : Inf
    pure_coulomb = tabulated || cache.ltp == LTP_COULOMB_ONLY

    cuts = elastic_cuts(cache, mu_lo, mu_hi)
    for E_b in Eout_mev
        mu_b = sgn * (E_b - e0) / e1
        mu_lo < mu_b < mu_hi && push!(cuts, mu_b)
    end
    sort!(unique!(cuts))

    u_t, w_t = quadrature(max(L + 2, 4), "gauss-legendre")
    moments(mu, d) = begin
        sigma = max(elastic_dsig_raw_d(cache, mu, d), 0.0)
        vcat(sigma .* legendre_polynomials_up_to_L(L, lab(mu)), sigma * energy(mu))
    end

    for k in 1:(length(cuts) - 1)
        x0, x1 = cuts[k], cuts[k + 1]
        x1 <= x0 && continue
        g = searchsortedfirst(Eout_mev, energy(0.5 * (x0 + x1)); rev=true) - 1
        (1 <= g <= Ng) || continue

        if tabulated && x1 <= mu_end
            # Linear piece of the tabulation, positive part only
            y0, y1 = tabulated_segment_values(cache, x0, x1)
            piece = positive_part_linear(x0, x1, y0, y1)
            isnothing(piece) && continue
            a, b, ya, yb = piece
            h = b - a
            σl[g, 1] += 2π * 0.5 * h * (ya + yb)
            # The outgoing energy is linear too: Simpson's rule is exact for the product
            σe[g] += 2π * h / 6.0 * (ya * energy(a) + 2.0 * (ya + yb) * energy(0.5 * (a + b)) + yb * energy(b))
            if L >= 1
                half = 0.5 * h
                mid = 0.5 * (a + b)
                for q in eachindex(u_t)
                    mu = mid + half * u_t[q]
                    W = 2π * w_t[q] * half * (ya + (yb - ya) * (mu - a) / h)
                    Pl = legendre_polynomials_up_to_L(L, lab(mu))
                    for l in 1:L
                        σl[g, l + 1] += W * Pl[l + 1]
                    end
                end
            end
        else
            # Coulomb above the tabulation, or a reconstructed piece of constant sign
            !pure_coulomb && elastic_dsig_raw(cache, 0.5 * (x0 + x1)) <= 0.0 && continue
            I = integrate_in_log_cosine(moments, cache, x0, x1)
            σl[g, 1] += 2π * (pure_coulomb ? coulomb_dcs_integral(cache, x0, x1) : I[1])
            for l in 1:L
                σl[g, l + 1] += 2π * I[l + 1]
            end
            σe[g] += 2π * I[L + 2]
        end
    end
    return nothing
end

"""
    elastic_soft_integrals(cache::ElasticDCSCache, mu_lo::Float64)

Gives the two integrals of the soft elastic collisions, center-of-mass cosines from
`mu_lo` to 1, on which the stopping power and the momentum transfer rest:
∫ dσ/dΩ (1 - μ) dμ and ∫ dσ/dΩ (1 - μ_lab) dμ. The range is cut as in
`elastic_group_moments!` and only the positive part of the cross-section is kept. On a
tabulated piece the first integral is exact (both factors are linear) and the second is
integrated by adaptive Gauss-Kronrod; on the Coulomb piece above the tabulation and on
reconstructed pieces both are integrated by adaptive Gauss-Kronrod in ln(1 - μ + δ), with
the Coulomb term and 1 - μ evaluated from 1 - μ + δ, and 1 - μ_lab written without
cancellation (`one_minus_mu_lab`).

# Input Argument(s)
- `cache::ElasticDCSCache` : elastic DCS cache.
- `mu_lo::Float64` : lower center-of-mass cosine bound, the cutoff cosine μ_c.

# Output Argument(s)
- `I_cm::Float64` : ∫ dσ/dΩ (1 - μ) dμ over [mu_lo, 1], in barn.
- `I_lab::Float64` : ∫ dσ/dΩ (1 - μ_lab) dμ over [mu_lo, 1], in barn.
"""
function elastic_soft_integrals(cache::ElasticDCSCache, mu_lo::Float64)
    mu_lo >= 1.0 && return 0.0, 0.0
    frame = cm_frame(cache)
    delta = cache.delta_coulomb
    tabulated = is_tabulated(cache)
    mu_end = tabulated ? cache.mu_grid[end] : Inf
    pure_coulomb = tabulated || cache.ltp == LTP_COULOMB_ONLY
    weights(mu, d) = begin
        sigma = max(elastic_dsig_raw_d(cache, mu, d), 0.0)
        x = d - delta
        [sigma * x, sigma * one_minus_mu_lab(frame, x)]
    end

    cuts = elastic_cuts(cache, mu_lo, 1.0)
    I_cm = 0.0
    I_lab = 0.0
    for k in 1:(length(cuts) - 1)
        x0, x1 = cuts[k], cuts[k + 1]
        x1 <= x0 && continue
        if tabulated && x1 <= mu_end
            y0, y1 = tabulated_segment_values(cache, x0, x1)
            piece = positive_part_linear(x0, x1, y0, y1)
            isnothing(piece) && continue
            a, b, ya, yb = piece
            h = b - a
            # dσ/dΩ and 1 - μ both linear: Simpson's rule is exact for the product
            I_cm += h / 6.0 * (ya * (1.0 - a) + 2.0 * (ya + yb) * (1.0 - 0.5 * (a + b)) + yb * (1.0 - b))
            I_lab += adaptive_gauss_kronrod(mu -> (ya + (yb - ya) * (mu - a) / h) * one_minus_mu_lab(frame, 1.0 - mu), a, b)
        else
            !pure_coulomb && elastic_dsig_raw(cache, 0.5 * (x0 + x1)) <= 0.0 && continue
            I = integrate_in_log_cosine(weights, cache, x0, x1)
            I_cm += I[1]
            I_lab += I[2]
        end
    end
    return I_cm, I_lab
end

"""
    integrate_elastic_scattering_below(cache::ElasticDCSCache, Ec_MeV::Float64,
    E_last_mev::Float64)

Gives the part of the total elastic cross-section above an outgoing-energy cutoff whose
scattered projectile ends below the last outgoing group, E_out < `E_last_mev`. It is
integrated from the cosine at which `elastic_group_moments!` stops, by
`integrate_tabulated_elastic_dcs` or `integrate_smooth_elastic_dcs`, so that it completes
the sum of the moments of order zero to the total cross-section (`tcs`).

# Input Argument(s)
- `cache::ElasticDCSCache` : elastic DCS cache.
- `Ec_MeV::Float64` : cutoff energy between soft and catastrophic events.
- `E_last_mev::Float64` : lower boundary of the last outgoing group [MeV].

# Output Argument(s)
- `σb::Float64` : cross-section below the last group, in barn.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering.
"""
function integrate_elastic_scattering_below(cache::ElasticDCSCache, Ec_MeV::Float64,
                                            E_last_mev::Float64)
    mu_lo = cache.identical ? 0.0 : -1.0
    mu_c = mu_cm_from_Eout_relativistic(cache.E_in_mev, max(Ec_MeV, 0.0), cache.m_mev, cache.M_mev)
    mu_c <= mu_lo && return 0.0
    e0, e1 = scattered_energy_affine(cache)
    mu_cut = clamp((E_last_mev - e0) / e1, mu_lo, mu_c)
    below = is_tabulated(cache) ? integrate_tabulated_elastic_dcs(cache, mu_lo, mu_cut) :
                                  integrate_smooth_elastic_dcs(cache, mu_lo, mu_cut)
    return 2.0 * pi * below
end

"""
    integrate_elastic_scattering_analytical(cache::ElasticDCSCache, Ec_MeV::Float64)

Gives the total elastic cross-section above an outgoing-energy cutoff, integrated exactly
or to round-off, with negative values of the reconstructed differential cross-section
set to zero as in `elastic_dsig`. The Coulomb term uses the cached wave number,
Sommerfeld parameter, Moliere screening shift, identical-particle flag and ENDF
projectile spin, and is integrated in closed form.

For a tabulated cross-section (LTP>2) the integral is the exact one of the
piecewise-linear cross-section, from `integrate_tabulated_elastic_dcs`. For LTP=1, LTP=2
and Coulomb-only data it is `integrate_smooth_elastic_dcs`: the nuclear and interference
terms are integrated by adaptive quadrature, and the negative parts are added back.
Either way the result equals the sum of the group moments of `elastic_group_moments!`
over the same range, plus `integrate_elastic_scattering_below`. Kept for debugging.

# Input Argument(s)
- `cache::ElasticDCSCache` : elastic DCS cache.
- `Ec_MeV::Float64` : cutoff energy between soft and catastrophic events.

# Output Argument(s)
- `σt::Float64` : integrated elastic cross-section in barn.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering.
"""
function integrate_elastic_scattering_analytical(cache::ElasticDCSCache, Ec_MeV::Float64)
    mu_b = mu_cm_from_Eout_relativistic(cache.E_in_mev, max(Ec_MeV, 0.0), cache.m_mev, cache.M_mev)
    mu_lo = cache.identical ? 0.0 : -1.0
    if mu_b <= mu_lo
        return 0.0
    end
    if is_tabulated(cache)
        return 2.0 * pi * integrate_tabulated_elastic_dcs(cache, mu_lo, mu_b)
    elseif cache.ltp ∈ (LTP_COULOMB_ONLY, 1, 2)
        return 2.0 * pi * integrate_smooth_elastic_dcs(cache, mu_lo, mu_b)
    end
    error("Elastic DCS cache for LTP=$(cache.ltp) has no analytical integral.")
end

"""
    eta_coulomb(Ein::Float64, m_eV::Float64, Z1::Int, Z2::Int;
    kinematics::String="standard")

Computes the Coulomb Sommerfeld parameter for an incident charged particle.

# Input Argument(s)
- `Ein::Float64` : incident kinetic energy in eV.
- `m_eV::Float64` : projectile rest mass energy in eV.
- `Z1::Int` : projectile charge number.
- `Z2::Int` : target charge number.
- `kinematics::String` : `"standard"` for the ENDF nonrelativistic expression or
  `"relativistic"` for `Z1*Z2*alpha/beta_rel`.

# Output Argument(s)
- `eta::Float64` : Coulomb Sommerfeld parameter.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering.
"""
function eta_coulomb(Ein::Float64, m_eV::Float64, Z1::Int, Z2::Int;
                     kinematics::String="standard")::Float64
    if kinematics == "standard"
        return (Z1 * Z2) * sqrt((ALPHA_FS^2) * m_eV / (2.0 * Ein))
    elseif kinematics == "relativistic"
        gamma = 1.0 + Ein / m_eV
        beta_rel = sqrt(1.0 - 1.0 / gamma^2)
        return Z1 * Z2 * ALPHA_FS / beta_rel
    end
    error("Unknown Coulomb kinematics $(kinematics) (should be standard or relativistic).")
end

"""
    k_coulomb(Ein::Float64, m_eV::Float64, M_eV::Float64;
    kinematics::String="standard")

Computes the Coulomb wave number from the center-of-mass momentum.

# Input Argument(s)
- `Ein::Float64` : incident kinetic energy in eV.
- `m_eV::Float64` : projectile rest mass energy in eV.
- `M_eV::Float64` : target rest mass energy in eV.
- `kinematics::String` : `"standard"` for nonrelativistic center-of-mass momentum or
  `"relativistic"` for relativistic center-of-mass momentum.

# Output Argument(s)
- `k::Float64` : Coulomb wave number in inverse centimeters scaled to the barn/sr
  convention used by the ENDF LAW=5 formulas.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering.
"""
function k_coulomb(Ein::Float64, m_eV::Float64, M_eV::Float64;
                   kinematics::String="standard")::Float64
    hbar_c_eVm = HBAR_eVs * C_m_per_s
    p_cm_eV_c = momentum_cm(Ein, m_eV, M_eV; kinematics=kinematics)
    return (p_cm_eV_c / hbar_c_eVm) * 1.0e-14
end

"""
    delta_coulomb(Ein::Float64, m_eV::Float64, M_eV::Float64, Z1::Int, Z2::Int;
    kinematics::String="standard", screening::Bool=true)

Computes the Coulomb screening shift used in denominators such as
`1 - mu + delta`. The Moliere screening angle is computed in the laboratory frame
and transformed to the center-of-mass frame by `(p_lab/p_cm)^2`.

# Input Argument(s)
- `Ein::Float64` : incident kinetic energy in eV.
- `m_eV::Float64` : projectile rest mass energy in eV.
- `M_eV::Float64` : target rest mass energy in eV.
- `Z1::Int` : projectile charge number.
- `Z2::Int` : target charge number.
- `kinematics::String` : `"standard"` or `"relativistic"` momentum convention.
- `screening::Bool` : when `false`, returns `0.0`.

# Output Argument(s)
- `delta::Float64` : dimensionless screening shift equal to `0.5*chi_a2_cm`.

# Reference(s)
- Bethe (1953), Moliere theory of multiple scattering.
"""
function delta_coulomb(Ein::Float64, m_eV::Float64, M_eV::Float64, Z1::Int, Z2::Int;
                       kinematics::String="standard", screening::Bool=true)::Float64
    if !screening
        return 0.0
    end

    chi_a2_lab = moliere_chi_a2(Ein, m_eV, Z1, Z2)
    p_lab = momentum_lab(Ein, m_eV; kinematics=kinematics)
    p_cm = momentum_cm(Ein, m_eV, M_eV; kinematics=kinematics)
    chi_a2_cm = chi_a2_lab * (p_lab / p_cm)^2
    return 0.5 * chi_a2_cm
end

"""
    coulomb_dcs(mu::Float64, k::Float64, eta::Float64, delta::Float64,
    identical::Bool, spin_s::Float64=0.5)

Gives the Coulomb elastic differential cross-section from precomputed Coulomb
wave number `k`, Sommerfeld parameter `eta`, and screening shift `delta`. For
identical particles it includes the exchange term and the spin-weighted
interference term. The input angle is clamped below `mu = 1` only for the unscreened
case (`delta == 0.0`).

# Input Argument(s)
- `mu::Float64` : center-of-mass scattering cosine.
- `k::Float64` : Coulomb wave number from `k_coulomb`.
- `eta::Float64` : Sommerfeld parameter from `eta_coulomb`.
- `delta::Float64` : screening shift from `delta_coulomb`; `0.0` gives the
  unscreened Coulomb singularity.
- `identical::Bool` : boolean indicating identical-particle scattering.
- `spin_s::Float64` : projectile spin from ENDF `SPI`, used only in the identical
  particle interference term.

# Output Argument(s)
- `σc::Float64` : Coulomb differential cross-section in barn/sr.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering.
"""
function coulomb_dcs(mu::Float64, k::Float64, eta::Float64, delta::Float64,
                     identical::Bool, spin_s::Float64=0.5)::Float64
    eps_mu = 1.0e-12
    mu_hi_safe = delta > 0.0 ? 1.0 : 1.0 - eps_mu
    if identical
        mu_safe = clamp(mu, 0.0, mu_hi_safe)
    else
        mu_safe = clamp(mu, -1.0, mu_hi_safe)
    end

    if !identical
        denom = 1.0 - mu_safe + delta
        return eta^2 / (k^2 * denom^2)
    end

    dminus = 1.0 - mu_safe + delta
    dplus = 1.0 + mu_safe + delta
    denom = dminus*dplus
    phase = eta * log(dplus / dminus)
    spin_factor = ((-1.0)^(2.0 * spin_s)) / (2.0 * spin_s + 1.0)
    bracket = (1.0 + mu_safe^2 + 2 * delta + delta^2) / denom + spin_factor * cos(phase)
    return 2.0 * eta^2 / (k^2 * denom) * bracket
end

"""
    elastic_dsig_ltp1(mu::Float64; a::Vector{ComplexF64}, b::Vector{Float64},
    NL::Int, identical::Bool, spin_s::Float64, k::Float64, eta::Float64, delta::Float64,
    nuclear_only::Bool=false, d::Float64=NaN)

Gives the elastic differential cross-section reconstructed from ENDF LAW=5, LTP=1
coefficients. The Coulomb and Coulomb-nuclear interference terms use precomputed
`k`, `eta`, and `delta`; exact forward scattering is allowed when `delta > 0`, while
the unscreened case is clamped below `mu = 1`.

# Input Argument(s)
- `mu::Float64` : center-of-mass scattering cosine.
- `a::Vector{ComplexF64}` : complex nuclear amplitude coefficients.
- `b::Vector{Float64}` : real nuclear cross-section coefficients.
- `NL::Int` : maximum coefficient order from the ENDF section.
- `identical::Bool` : boolean indicating identical-particle scattering.
- `spin_s::Float64` : projectile spin used in the identical-particle interference term.
- `k::Float64` : Coulomb wave number.
- `eta::Float64` : Sommerfeld parameter.
- `delta::Float64` : screening shift.
- `nuclear_only::Bool` : boolean to leave out the Coulomb term and return only the
  nuclear and Coulomb-nuclear interference terms.
- `d::Float64` : 1 - μ + δ when known exactly (from an integration in ln d); computed
  from `mu` when NaN. Near μ = 1 it is of the order of δ, and recomputing it from μ
  would lose about half of the significant digits.

# Output Argument(s)
- `σs::Float64` : elastic differential cross-section in barn/sr.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5, LTP=1, equations 6.13 and 6.14.
"""
function elastic_dsig_ltp1(mu::Float64;
                           a::Vector{ComplexF64}, b::Vector{Float64}, NL::Int,
                           identical::Bool, spin_s::Float64, k::Float64, eta::Float64, delta::Float64,
                           nuclear_only::Bool=false, d::Float64=NaN)

    eps_mu = 1.0e-12
    mu_hi_safe = delta > 0.0 ? 1.0 : 1.0 - eps_mu
    if identical
        mu_safe = clamp(mu, 0.0, mu_hi_safe)
    else
        mu_safe = clamp(mu, -1.0, mu_hi_safe)
    end


    sigma_c = nuclear_only ? 0.0 : coulomb_dcs(mu_safe, k, eta, delta, identical, spin_s)
    dminus = isnan(d) ? 1.0 - mu_safe + delta : d
    dplus = 2.0 * (1.0 + delta) - dminus

    n_a = min(NL, length(a) - 1)

    if !identical
        phase = exp(im * eta * log(dminus / 2.0))
        sum_a = 0.0 + 0.0im
        for l in 0:n_a
            Pl = legendre_polynomials(l, mu_safe)
            sum_a += ((2 * l + 1) / 2) * a[l + 1] * Pl
        end
        term_a = (2.0 * eta / dminus) * real(phase * sum_a)

        max_b = min(2 * NL, length(b) - 1)
        sum_b = 0.0
        for l in 0:max_b
            sum_b += ((2 * l + 1) / 2) * b[l + 1] * legendre_polynomials(l, mu_safe)
        end

        return sigma_c - term_a + sum_b
    else
        phase_minus = exp(im * eta * log(dminus / 2.0))
        phase_plus = exp(im * eta * log(dplus / 2.0))

        sum_a = 0.0 + 0.0im
        for l in 0:n_a
            Pl = legendre_polynomials(l, mu_safe)
            pref = (dplus * phase_minus + (-1.0)^l * dminus * phase_plus)
            sum_a += pref * ((2 * l + 1) / 2) * a[l + 1] * Pl
        end
        term_a = (2.0 * eta / (dminus * dplus)) * real(sum_a)

        max_b = min(NL, length(b) - 1)
        sum_b = 0.0
        for l in 0:max_b
            sum_b += ((4 * l + 1) / 2) * b[l + 1] * legendre_polynomials(2 * l, mu_safe)
        end

        return sigma_c - term_a + sum_b
    end
end

"""
    elastic_dsig_ltp2(mu::Float64; c::Vector{Float64}, NL::Int,
    identical::Bool, spin_s::Float64, k::Float64, eta::Float64, delta::Float64,
    nuclear_only::Bool=false, d::Float64=NaN)

Gives the elastic differential cross-section reconstructed from ENDF LAW=5, LTP=2
coefficients. The Coulomb term uses precomputed `k`, `eta`, and `delta`; exact
forward scattering is allowed when `delta > 0`, while the unscreened case is clamped
below `mu = 1`.

# Input Argument(s)
- `mu::Float64` : center-of-mass scattering cosine.
- `c::Vector{Float64}` : real nuclear coefficient table.
- `NL::Int` : maximum coefficient order from the ENDF section.
- `identical::Bool` : boolean indicating identical-particle scattering.
- `spin_s::Float64` : projectile spin used in the identical-particle interference term.
- `k::Float64` : Coulomb wave number.
- `eta::Float64` : Sommerfeld parameter.
- `delta::Float64` : screening shift.
- `nuclear_only::Bool` : boolean to leave out the Coulomb term and return only the
  nuclear and Coulomb-nuclear interference terms.
- `d::Float64` : 1 - μ + δ when known exactly (from an integration in ln d); computed
  from `mu` when NaN. Near μ = 1 it is of the order of δ, and recomputing it from μ
  would lose about half of the significant digits.

# Output Argument(s)
- `σs::Float64` : elastic differential cross-section in barn/sr.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5, LTP=2, equations 6.15 to 6.18.
"""
function elastic_dsig_ltp2(mu::Float64;
                           c::Vector{Float64}, NL::Int,
                           identical::Bool, spin_s::Float64, k::Float64, eta::Float64, delta::Float64,
                           nuclear_only::Bool=false, d::Float64=NaN)

    eps_mu = 1.0e-12
    mu_hi_safe = delta > 0.0 ? 1.0 : 1.0 - eps_mu
    if identical
        mu_safe = clamp(mu, 0.0, mu_hi_safe)
    else
        mu_safe = clamp(mu, -1.0, mu_hi_safe)
    end


    sigma_c = nuclear_only ? 0.0 : coulomb_dcs(mu_safe, k, eta, delta, identical, spin_s)
    dminus = isnan(d) ? 1.0 - mu_safe + delta : d
    dplus = 2.0 * (1.0 + delta) - dminus

    max_c = min(NL, length(c) - 1)
    sigma_r = 0.0
    if !identical
        for l in 0:max_c
            sigma_r += ((2 * l + 1) / 2) * c[l + 1] * legendre_polynomials(l, mu_safe)
        end
        return sigma_c + sigma_r / dminus
    else
        for l in 0:max_c
            sigma_r += ((4 * l + 1) / 2) * c[l + 1] * legendre_polynomials(2 * l, mu_safe)
        end
        return sigma_c + sigma_r / (dminus * dplus)
    end
end

"""
    elastic_dsig_ltpgreater2(mu::Float64; pni::Float64, sigmaNI::Float64,
    identical::Bool, spin_s::Float64, k::Float64, eta::Float64, delta::Float64)

Gives the elastic differential cross-section reconstructed from ENDF LAW=5, LTP>2
tabulated nuclear-interference probabilities. The nuclear-interference probability
is combined with the Coulomb term evaluated from precomputed `k`, `eta`, and `delta`.
Exact forward scattering is allowed when `delta > 0`, while the unscreened case is
clamped below `mu = 1`.

# Input Argument(s)
- `mu::Float64` : center-of-mass scattering cosine.
- `pni::Float64` : tabulated nuclear-interference probability density at `mu`.
- `sigmaNI::Float64` : nuclear-interference integrated cross-section in barn.
- `identical::Bool` : boolean indicating identical-particle scattering.
- `spin_s::Float64` : projectile spin used in the identical-particle interference term.
- `k::Float64` : Coulomb wave number.
- `eta::Float64` : Sommerfeld parameter.
- `delta::Float64` : screening shift.

# Output Argument(s)
- `σs::Float64` : elastic differential cross-section in barn/sr.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5, LTP>2, equations 6.19 and 6.20.
"""
function elastic_dsig_ltpgreater2(mu::Float64;
                                  pni::Float64,sigmaNI::Float64,
                                  identical::Bool, spin_s::Float64, k::Float64, eta::Float64, delta::Float64)
    eps_mu = 1.0e-12
    mu_hi_safe = delta > 0.0 ? 1.0 : 1.0 - eps_mu
    if identical
        mu_safe = clamp(mu, 0.0, mu_hi_safe)
    else
        mu_safe = clamp(mu, -1.0, mu_hi_safe)
    end


    sigma_c = coulomb_dcs(mu_safe, k, eta, delta, identical, spin_s)

    return sigma_c + sigmaNI * pni
end

"""
    load_elastic_scattering_endf(endf_path::String, Z::Int, A::Int, particle::Particle)

Loads the ENDF-6 elastic scattering data of one isotope: the LAW=5 coefficients (LTP=1,
LTP=2) or tabulated nuclear-plus-interference distributions P_NI (LTP>2) with the
nuclear-plus-interference cross-section σ_NI of MF=3. The MF=6/MT=2 product subsection is
selected by the projectile `ZAP`, and the projectile spin is read from the ENDF LAW=5
`SPI` field. For LTP>2, a warning is given when the tabulated P_NI do not integrate to 1
over the cosine (ENDF-102 eq. 6.18), which would point to a misread table.

# Input Argument(s)
- `endf_path::String` : path to the ENDF-6 file.
- `Z::Int` : target atomic number.
- `A::Int` : target mass number.
- `particle::Particle` : incident particle type.

# Output Argument(s)
- `data::NamedTuple` : isotope elastic-scattering data: table energies [MeV], common
  cosine grid of the LTP>2 tables, LAW=5 representation, identical-particle flag,
  projectile and target charges, ENDF projectile spin, masses, energy interpolation
  data, and the coefficient or P_NI tables with σ_NI.

# Reference(s)
- ENDF-6 Formats Manual, MF=3 and MF=6, MT=2, LAW=5 charged-particle elastic scattering.
"""
function load_elastic_scattering_endf(endf_path::String, Z::Int, A::Int, particle::Particle)
    lines = readlines(endf_path)
    mf6mt2_subsections = read_mf6_mt(lines, 2; source=endf_path)

    if isempty(mf6mt2_subsections)
        error("MF=6 MT=2 contains no product subsection in $(endf_path).")
    end

    expected_zap = ZAP(particle)
    eladata = nothing
    if isnothing(expected_zap)
        if length(mf6mt2_subsections) != 1
            error("MF=6 MT=2 has $(length(mf6mt2_subsections)) product subsections in $(endf_path), but no ZAP selector is defined for particle $(particle).")
        end
        eladata = mf6mt2_subsections[1]
    else
        for sub in mf6mt2_subsections
            if round(Int, sub.ZAP) == expected_zap
                eladata = sub
                break
            end
        end
        if isnothing(eladata)
            zaps = [round(Int, sub.ZAP) for sub in mf6mt2_subsections]
            error("MF=6 MT=2 does not contain product subsection with ZAP=$(expected_zap) in $(endf_path). Available ZAP values: $(zaps).")
        end
    end

    if eladata.LAW != 5
        error("MF=6 MT=2 product LAW=$(eladata.LAW) (expected 5) in $(endf_path).")
    end

    isnothing(eladata.SPI) && error("MF=6 MT=2 LAW=5 subsection has no SPI projectile spin in $(endf_path).")
    identical = (eladata.LIDP == 1)
    common = (ltp=eladata.LTP, identical=identical, z1=get_charge(particle), z2=Z,
              spin_s=Float64(eladata.SPI), m_MeV=get_mass(particle), M_MeV=get_mass(Z, A),
              mf6_nbt=eladata.mf6_nbt, mf6_int=eladata.mf6_int)
    none = (ba_tables_ltp1=nothing, c_tables_ltp2=nothing, mupni_tables=nothing,
            sigmaNI_E=nothing, sigmaNI_S=nothing, sigmaNI_nbt=nothing, sigmaNI_int=nothing)

    if eladata.LTP == 1 && !isempty(eladata.ba_tables_ltp1)
        energies_eV = sort(collect(keys(eladata.ba_tables_ltp1)))
        return (; common..., none..., energies=energies_eV .* 1.0e-6, mu_grid=Float64[],
                ba_tables_ltp1=eladata.ba_tables_ltp1)
    elseif eladata.LTP == 2 && !isempty(eladata.c_tables_ltp2)
        energies_eV = sort(collect(keys(eladata.c_tables_ltp2)))
        return (; common..., none..., energies=energies_eV .* 1.0e-6, mu_grid=Float64[],
                c_tables_ltp2=eladata.c_tables_ltp2)
    elseif eladata.LTP > 2 && !isempty(eladata.mupni_tables)
        (_ZA3, _AWR3, E3, S3, NBT3, INT3) = read_mf3_mt(lines, 2)
        mf6_energies = sort(collect(keys(eladata.mupni_tables)))

        # P_NI is defined from μ_min = -1 (0 for identical particles) to μ_max, and
        # integrates to 1 there (ENDF-102 eqs. 6.17, 6.18)
        mu_min = identical ? 0.0 : -1.0
        mu_ref, _ = eladata.mupni_tables[mf6_energies[1]]
        same_mu_grid = true
        worst_norm = 0.0
        for Ein in mf6_energies
            mu_tab, pni_tab = eladata.mupni_tables[Ein]
            if isempty(mu_tab) || abs(mu_tab[1] - mu_min) > 1.0e-12
                error("LTP>2 mu grid must start at $(mu_min) in $(endf_path).")
            end
            if length(mu_tab) != length(mu_ref) || any(i -> abs(mu_tab[i] - mu_ref[i]) > 1.0e-8, eachindex(mu_tab))
                same_mu_grid = false
            end
            norm = sum(0.5 * (mu_tab[i + 1] - mu_tab[i]) * (pni_tab[i] + pni_tab[i + 1]) for i in 1:(length(mu_tab) - 1))
            worst_norm = max(worst_norm, abs(norm - 1.0))
        end
        if worst_norm > 0.05
            @warn "LAW=5 LTP=$(eladata.LTP) P_NI tables of $(endf_path) do not integrate to 1 over the cosine (worst deviation $(round(worst_norm, sigdigits=3))); check how the table is read."
        end

        local mu_grid
        local mupni_tables
        if same_mu_grid
            mu_grid = mu_ref
            mupni_tables = Dict{Float64,Vector{Float64}}(Ein => eladata.mupni_tables[Ein][2] for Ein in mf6_energies)
        else
            mu_grid = Float64[]
            mupni_tables = eladata.mupni_tables
        end
        return (; common..., none..., energies=mf6_energies .* 1.0e-6, mu_grid=mu_grid,
                mupni_tables=mupni_tables, sigmaNI_E=E3, sigmaNI_S=S3, sigmaNI_nbt=NBT3,
                sigmaNI_int=INT3)
    end

    error("Unsupported LTP=$(eladata.LTP) or missing elastic scattering data in $(endf_path).")
end

"""
    init_elastic_scattering_endf_db(db_name::String,isotopes::Vector{Tuple{Int,Int}},
    particle::Particle; data_root::Union{Nothing,String}=nothing,
    coulomb_kinematics::String="standard", screening::Bool=true)

Initializes the elastic-scattering ENDF database of a set of target isotopes for one
incident particle.

# Input Argument(s)
- `db_name::String` : ENDF library name (e.g. "TENDL2023").
- `isotopes::Vector{Tuple{Int,Int}}` : list of `(Z, A)` target isotopes.
- `particle::Particle` : incident particle type.
- `data_root::Union{Nothing,String}` : root directory of the ENDF data.
- `coulomb_kinematics::String` : Coulomb DCS kinematics convention, either `"standard"`
  for the ENDF nonrelativistic convention or `"relativistic"` for relativistic beta and
  center-of-mass momentum.
- `screening::Bool` : boolean to include Moliere screening in Coulomb terms.

# Output Argument(s)
- `db::ElasticScatteringENDFDB` : initialized elastic-scattering database.

# Reference(s)
- ENDF-6 Formats Manual, MF=3 and MF=6, MT=2, LAW=5 charged-particle elastic scattering.
"""
function init_elastic_scattering_endf_db(db_name::String,isotopes::Vector{Tuple{Int,Int}},
    particle::Particle; data_root::Union{Nothing,String}=nothing,
    coulomb_kinematics::String="standard", screening::Bool=true)
    root = isnothing(data_root) ? normpath(joinpath(@__DIR__, "..", "..", "data")) : data_root
    db = Dict{Tuple{Int,Int},IsotopeElasticDB}()
    for (Z, A) in isotopes
        endf_path = endf_path_for_isotope(db_name, Z, A, particle; data_root=root)
        if !isfile(endf_path)
            error("ENDF file not found: $(endf_path). Prepare the \"$(db_name)\" ENDF "*
                  "library first by running data/create_radiant_endf_library.jl on your "*
                  "raw ENDF files before starting this simulation.")
        end
        data = load_elastic_scattering_endf(endf_path, Z, A, particle)
        db[(Z, A)] = IsotopeElasticDB(data.energies, data.mu_grid,
              data.ltp, data.identical, data.z1, data.z2, data.spin_s,
              data.m_MeV, data.M_MeV, data.mf6_nbt, data.mf6_int,
              data.ba_tables_ltp1, data.c_tables_ltp2, data.mupni_tables,
              data.sigmaNI_E, data.sigmaNI_S, data.sigmaNI_nbt, data.sigmaNI_int)
    end
    return ElasticScatteringENDFDB(lowercase(coulomb_kinematics), screening, db)
end

"""
    init_elastic_scattering_coulomb_db(isotopes::Vector{Tuple{Int,Int}}, particle::Particle;
    coulomb_kinematics::String="standard", screening::Bool=true)

Initializes an elastic-scattering database containing only the Coulomb (Rutherford)
contribution for a set of target isotopes and one incident particle. No ENDF file is
read: each isotope entry stores the projectile/target charges and masses, the
identical-particle flag and the projectile spin, and is marked with `LTP_COULOMB_ONLY` so
that DCS evaluations and integrations use `coulomb_dcs` alone.

# Input Argument(s)
- `isotopes::Vector{Tuple{Int,Int}}` : list of `(Z, A)` target isotopes.
- `particle::Particle` : incident particle type.
- `coulomb_kinematics::String` : Coulomb DCS kinematics convention, either `"standard"`
  or `"relativistic"`.
- `screening::Bool` : boolean to include Moliere screening in Coulomb terms.

# Output Argument(s)
- `db::ElasticScatteringENDFDB` : initialized Coulomb-only elastic-scattering database.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering (Coulomb term).
"""
function init_elastic_scattering_coulomb_db(isotopes::Vector{Tuple{Int,Int}}, particle::Particle;
    coulomb_kinematics::String="standard", screening::Bool=true)
    Z1 = get_charge(particle)
    Z1 == 0 && error("Coulomb elastic scattering requires a charged incident particle.")
    m_MeV = get_mass(particle)
    zap = ZAP(particle)
    db = Dict{Tuple{Int,Int},IsotopeElasticDB}()
    for (Z, A) in isotopes
        identical = !isnothing(zap) && zap == 1000 * Z + A
        spin_s = identical ? COULOMB_PROJECTILE_SPIN[zap] : 0.0
        db[(Z, A)] = IsotopeElasticDB(Float64[], Float64[],
              LTP_COULOMB_ONLY, identical, Z1, Z, spin_s, m_MeV, get_mass(Z, A),
              Int[], Int[], nothing, nothing, nothing, nothing, nothing, nothing, nothing)
    end
    return ElasticScatteringENDFDB(lowercase(coulomb_kinematics), screening, db)
end
