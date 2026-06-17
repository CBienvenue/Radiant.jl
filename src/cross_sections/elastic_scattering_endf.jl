using Printf

const HBAR_eVs = 6.582119569e-16
const C_m_per_s = 299_792_458.0
const ALPHA_FS = 1.0 / 137.035999177

"""
    prepare_elastic_dsig_cache(db::ElasticScatteringENDFDB, Z::Int, A::Int,
    E_in_MeV::Float64)

Prepares the energy-resolved elastic differential cross-section cache for one isotope and
one incident energy.

# Input Argument(s)
- `db::ElasticScatteringENDFDB` : elastic-scattering database containing ENDF data by
  isotope.
- `Z::Int` : target atomic number.
- `A::Int` : target mass number.
- `E_in_MeV::Float64` : incident particle kinetic energy in MeV.
- `db.dcs_source::String` : source of the differential cross-section, either
  `"database"` for tabulated values or `"coefficients"` for ENDF coefficient
  reconstruction.
- `db.coulomb_kinematics::String` : Coulomb DCS kinematics convention, either
  `"standard"` for the ENDF nonrelativistic convention or `"relativistic"` for
  relativistic beta and center-of-mass momentum.
- `db.screening::Bool` : common screening setting stored in the database; when true,
  Moliere screening is included in Coulomb terms.

# Output Argument(s)
- `cache::ElasticDCSCache` : cache containing the data needed to evaluate the elastic
  differential cross-section at `E_in_MeV`.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering.
"""
function prepare_elastic_dsig_cache(db::ElasticScatteringENDFDB, Z::Int, A::Int, E_in_MeV::Float64)
    iso = db.isotopes[(Z, A)]
    E_eV = E_in_MeV * 1.0e6
    source = db.dcs_source
    coulomb_kinematics = db.coulomb_kinematics
    screening = db.screening
    m_MeV = iso.m_mev
    M_MeV = iso.M_mev
    m_eV = m_MeV * 1.0e6
    M_eV = M_MeV * 1.0e6
    k_c = k_coulomb(E_eV, m_eV, M_eV; kinematics=coulomb_kinematics)
    eta = eta_coulomb(E_eV, m_eV, iso.Z1, iso.Z2; kinematics=coulomb_kinematics)
    delta = delta_coulomb(E_eV, m_eV, M_eV, iso.Z1, iso.Z2; kinematics=coulomb_kinematics, screening=screening)

    if source == "database"
        isempty(iso.energies) && error("No tabulated elastic DCS was stored for isotope (Z=$(Z), A=$(A)). Reinitialize with dcs_source=\"database\".")
        isempty(iso.dcs) && error("No tabulated elastic DCS was stored for isotope (Z=$(Z), A=$(A)). Reinitialize with dcs_source=\"database\".")

        function mu_grid_for_energy(Ei::Float64)
            if iso.ltp > 2 && (iso.mupni_tables isa Dict{Float64,Tuple{Vector{Float64},Vector{Float64}}})
                return iso.mupni_tables[Ei][1]
            end
            return iso.mu_grid
        end

        function sigma_on_grid(Ei_MeV::Float64, mu_grid::Vector{Float64}, dcs_row::Vector{Float64}, mu_cm::Float64)
            if iso.ltp > 2 && mu_cm > mu_grid[end]
                Ei_eV = Ei_MeV * 1.0e6
                k_i = k_coulomb(Ei_eV, m_eV, M_eV; kinematics=coulomb_kinematics)
                eta_i = eta_coulomb(Ei_eV, m_eV, iso.Z1, iso.Z2; kinematics=coulomb_kinematics)
                delta_i = delta_coulomb(Ei_eV, m_eV, M_eV, iso.Z1, iso.Z2; kinematics=coulomb_kinematics, screening=screening)
                return coulomb_dcs(mu_cm, k_i, eta_i, delta_i, iso.identical, iso.spin_s)
            end
            return interp_linear(mu_grid, dcs_row, mu_cm)
        end

        energies = iso.energies
        local mu_grid
        local sigma_grid
        if E_in_MeV <= energies[1]
            mu_grid = copy(mu_grid_for_energy(energies[1]))
            sigma_grid = copy(iso.dcs[1])
        elseif E_in_MeV >= energies[end]
            mu_grid = copy(mu_grid_for_energy(energies[end]))
            sigma_grid = copy(iso.dcs[end])
        else
            i = searchsortedlast(energies, E_in_MeV)
            E1 = energies[i]
            E2 = energies[i + 1]
            mu1 = mu_grid_for_energy(E1)
            mu2 = mu_grid_for_energy(E2)
            row1 = iso.dcs[i]
            row2 = iso.dcs[i + 1]
            itp = select_interpolation(iso.mf6_nbt, iso.mf6_int, i)

            mu_grid = unique(sort(vcat(mu1, mu2)))
            sigma_grid = Vector{Float64}(undef, length(mu_grid))
            for k in eachindex(mu_grid)
                mu = mu_grid[k]
                s1 = sigma_on_grid(E1, mu1, row1, mu)
                s2 = sigma_on_grid(E2, mu2, row2, mu)
                sigma_grid[k] = interp_pair(E_in_MeV, E1, s1, E2, s2, itp)
            end
        end

        return ElasticDCSCache(iso.ltp, E_in_MeV, k_c, eta, delta,
                               mu_grid, sigma_grid, Float64[], 0.0, 0,
                               ComplexF64[], Float64[], Float64[],
                               iso.spin_s, m_MeV, M_MeV, iso.identical)
    elseif source == "coefficients"
        if iso.ltp == 1
            coeffs = iso.ba_tables_ltp1
            coeffs === nothing && error("Missing LTP=1 coefficient tables for coefficient-based DCS.")
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
                                   Float64[], Float64[], Float64[], 0.0, NL, a, b, Float64[],
                                   iso.spin_s, m_MeV, M_MeV, iso.identical)
        elseif iso.ltp == 2
            coeffs = iso.c_tables_ltp2
            coeffs === nothing && error("Missing LTP=2 coefficient tables for coefficient-based DCS.")
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
                                   Float64[], Float64[], Float64[], 0.0, NL, ComplexF64[], Float64[], c,
                                   iso.spin_s, m_MeV, M_MeV, iso.identical)
        elseif iso.ltp > 2
            tables = iso.mupni_tables
            tables === nothing && error("Missing LTP>2 pni tables for coefficient-based DCS.")
            if iso.sigmaNI_E === nothing || iso.sigmaNI_S === nothing ||
               iso.sigmaNI_nbt === nothing || iso.sigmaNI_int === nothing
                error("Missing sigmaNI tables for coefficient-based DCS.")
            end

            energies = sort(collect(keys(tables)))
            local mu_tab
            local pni_tab
            if length(energies) == 1 || E_eV <= energies[1]
                Eref = energies[1]
                if tables isa Dict{Float64,Vector{Float64}}
                    isempty(iso.mu_grid) && error("Missing common mu grid for LTP>2 coefficient-based DCS.")
                    mu_tab = iso.mu_grid
                    pni_tab = tables[Eref]
                else
                    mu_tab, pni_tab = tables[Eref]
                end
            elseif E_eV >= energies[end]
                Eref = energies[end]
                if tables isa Dict{Float64,Vector{Float64}}
                    isempty(iso.mu_grid) && error("Missing common mu grid for LTP>2 coefficient-based DCS.")
                    mu_tab = iso.mu_grid
                    pni_tab = tables[Eref]
                else
                    mu_tab, pni_tab = tables[Eref]
                end
            else
                idx = searchsortedlast(energies, E_eV)
                E1 = energies[idx]
                E2 = energies[idx + 1]
                itp = select_interpolation(iso.mf6_nbt, iso.mf6_int, idx)

                local mu1
                local mu2
                local pni1
                local pni2
                if tables isa Dict{Float64,Vector{Float64}}
                    isempty(iso.mu_grid) && error("Missing common mu grid for LTP>2 coefficient-based DCS.")
                    mu1 = iso.mu_grid
                    mu2 = iso.mu_grid
                    pni1 = tables[E1]
                    pni2 = tables[E2]
                else
                    mu1, pni1 = tables[E1]
                    mu2, pni2 = tables[E2]
                end

                if length(mu1) != length(mu2) || any(i -> mu1[i] != mu2[i], eachindex(mu1))
                    error("Energy interpolation for LTP>2 requires identical mu grids at the bracketing energies.")
                end

                mu_tab = mu1
                pni_tab = Vector{Float64}(undef, length(mu_tab))
                for i in eachindex(mu_tab)
                    pni_tab[i] = interp_pair(E_eV, E1, pni1[i], E2, pni2[i], itp)
                end
            end

            sigmaNI = interp_TAB1(E_eV, iso.sigmaNI_E, iso.sigmaNI_S,
                                  iso.sigmaNI_nbt, iso.sigmaNI_int)
            sigma_grid = Vector{Float64}(undef, length(mu_tab))
            for k in eachindex(mu_tab)
                sigma_grid[k] = elastic_dsig_ltpgreater2(mu_tab[k]; pni=pni_tab[k], sigmaNI=sigmaNI, identical=iso.identical, spin_s=iso.spin_s, k=k_c, eta=eta, delta=delta)
            end

            return ElasticDCSCache(iso.ltp, E_in_MeV, k_c, eta, delta,
                                   copy(mu_tab), sigma_grid, copy(pni_tab), sigmaNI, 0,
                                   ComplexF64[], Float64[], Float64[],
                                   iso.spin_s, m_MeV, M_MeV, iso.identical)
        else
            error("Unsupported LTP=$(iso.ltp) for coefficient-based DCS cache.")
        end
    else
        error("Unknown elastic DCS source $(db.dcs_source) (should be database or coefficients).")
    end
end

"""
    elastic_dsig_raw(cache::ElasticDCSCache, mu_cm::Float64)

Gives the elastic differential cross-section from an energy-resolved cache without clipping
negative numerical artifacts. For tabulated LTP>2 data evaluated beyond the stored
forward angular grid, the Coulomb contribution is evaluated directly from the cached
`k`, `eta`, `delta`, identical-particle flag, and ENDF projectile spin.

# Input Argument(s)
- `cache::ElasticDCSCache` : cache generated by `prepare_elastic_dsig_cache`.
- `mu_cm::Float64` : center-of-mass scattering cosine.

# Output Argument(s)
- `σs::Float64` : elastic differential cross-section in barn/sr.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering.
"""
function elastic_dsig_raw(cache::ElasticDCSCache, mu_cm::Float64)
    if !isempty(cache.sigma_grid)
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
    integrate_elastic_scattering_numerical_GL(cache::ElasticDCSCache, Ec_MeV::Float64)

Gives the total elastic cross-section above an outgoing-energy cutoff using
Gauss-Legendre quadrature over the cached center-of-mass angular distribution.

# Input Argument(s)
- `cache::ElasticDCSCache` : cache generated by `prepare_elastic_dsig_cache`.
- `Ec_MeV::Float64` : cutoff energy between soft and catastrophic events.

# Output Argument(s)
- `σt::Float64` : integrated elastic cross-section in barn.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering.
"""
function integrate_elastic_scattering_numerical_GL(cache::ElasticDCSCache, Ec_MeV::Float64)
    E_in_MeV = cache.E_in_mev
    E_out_MeV = max(Ec_MeV, 0.0)
    mu_b = mu_cm_from_Eout_relativistic(E_in_MeV, E_out_MeV, cache.m_mev, cache.M_mev)

    mu_lo = cache.identical ? 0.0 : -1.0
    if mu_b <= mu_lo
        return 0.0
    end

    Nq = 64
    μ_nodes, w = quadrature(Nq, "gauss-legendre")
    Δ = 0.5 * (mu_b - mu_lo)
    μ_mid = 0.5 * (mu_b + mu_lo)

    acc = 0.0
    for k in eachindex(μ_nodes)
        mu_cm = Δ * μ_nodes[k] + μ_mid
        sigma = elastic_dsig(cache, mu_cm)
        acc += w[k] * sigma
    end
    return 2.0 * pi * acc * Δ
end

"""
    integrate_elastic_scattering_numerical_GK(cache::ElasticDCSCache, Ec_MeV::Float64)

Gives the total elastic cross-section above an outgoing-energy cutoff using adaptive
Gauss-Kronrod quadrature over the cached center-of-mass angular distribution.

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
    integrate_elastic_scattering_numerical_TS(cache::ElasticDCSCache, Ec_MeV::Float64)

Gives the total elastic cross-section above an outgoing-energy cutoff using tanh-sinh
quadrature over the cached center-of-mass angular distribution.

# Input Argument(s)
- `cache::ElasticDCSCache` : cache generated by `prepare_elastic_dsig_cache`.
- `Ec_MeV::Float64` : cutoff energy between soft and catastrophic events.

# Output Argument(s)
- `σt::Float64` : integrated elastic cross-section in barn.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering.
"""
function integrate_elastic_scattering_numerical_TS(cache::ElasticDCSCache, Ec_MeV::Float64)
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

    f = mu -> elastic_dsig(cache, mu)
    val = tanh_sinh_integral(f, a, b)
    return 2.0 * pi * val
end

"""
    negative_ltpgt2_correction(cache::ElasticDCSCache, mu_lo::Float64, mu_hi::Float64)

Computes the integral correction needed when an LTP>2 reconstructed differential
cross-section contains negative tabulation artifacts.

# Input Argument(s)
- `cache::ElasticDCSCache` : elastic DCS cache for an LTP>2 ENDF section.
- `mu_lo::Float64` : lower center-of-mass cosine bound.
- `mu_hi::Float64` : upper center-of-mass cosine bound.

# Output Argument(s)
- `correction::Float64` : positive correction that removes negative DCS contributions
  from the integrated cross-section.
"""
function negative_ltpgt2_correction(cache::ElasticDCSCache, mu_lo::Float64, mu_hi::Float64)
    mu_end = min(mu_hi, cache.mu_grid[end])
    if mu_end <= mu_lo || isempty(cache.pni_grid)
        return 0.0
    end

    correction = 0.0
    for i in 1:(length(cache.mu_grid) - 1)
        seg_lo = max(mu_lo, cache.mu_grid[i])
        seg_hi = min(mu_end, cache.mu_grid[i + 1])
        if seg_hi <= seg_lo
            continue
        end
        negative_part = mu -> min(elastic_dsig_raw(cache, mu), 0.0)
        correction -= tanh_sinh_integral(negative_part, seg_lo, seg_hi; tol=1.0e-10)
    end
    return correction
end

"""
    integrate_elastic_scattering_analytical(cache::ElasticDCSCache, Ec_MeV::Float64)

Gives the total elastic cross-section above an outgoing-energy cutoff using the analytical
ENDF coefficient representation. Coulomb terms use the cached wave number, Sommerfeld
parameter, Moliere screening shift, identical-particle flag, and ENDF projectile spin.

# Input Argument(s)
- `cache::ElasticDCSCache` : coefficient-based cache generated by
  `prepare_elastic_dsig_cache` from an ENDF database initialized with
  `dcs_source="coefficients"`.
- `Ec_MeV::Float64` : cutoff energy between soft and catastrophic events.

# Output Argument(s)
- `σt::Float64` : integrated elastic cross-section in barn.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering.
"""
function integrate_elastic_scattering_analytical(cache::ElasticDCSCache, Ec_MeV::Float64)
    m_MeV = cache.m_mev
    M_MeV = cache.M_mev
    E_in_MeV = cache.E_in_mev
    E_out_MeV = max(Ec_MeV, 0.0)
    mu_b = mu_cm_from_Eout_relativistic(E_in_MeV, E_out_MeV, m_MeV, M_MeV)
    identical = cache.identical

    mu_lo = identical ? 0.0 : -1.0
    if mu_b <= mu_lo
        return 0.0
    end
    
    eta = cache.eta_coulomb
    k = cache.k_coulomb
    delta = cache.delta_coulomb
    spin_s = cache.spin_s
    eps_mu = 1.0e-12
    mu_hi_safe = delta > 0.0 ? 1.0 : 1.0 - eps_mu
    mu_hi = mu_b

    sigma_c = 0.0
    sigma_int = 0.0
    if !identical
        denom_lo = 1.0 - mu_lo + delta
        denom_hi = 1.0 - mu_hi + delta
        sigma_c = (eta^2 / k^2) * (1.0 / denom_hi - 1.0 / denom_lo)
    else
        a2 = 1.0 + delta^2 + 2.0 * delta
        I0 = mu_hi / (a2 - mu_hi^2) - mu_lo / (a2 - mu_lo^2)
        sigma_c = 2.0 * (eta^2 / k^2) * I0

        Nq = 64
        μ_nodes, w = quadrature(Nq, "gauss-legendre")
        Δ = 0.5 * (mu_hi - mu_lo)
        μ_mid = 0.5 * (mu_hi + mu_lo)
        spin_factor = ((-1.0)^(2.0 * spin_s)) / (2.0 * spin_s + 1.0)
        acc = 0.0
        for kq in eachindex(μ_nodes)
            mu = Δ * μ_nodes[kq] + μ_mid
            dplus = 1.0 + mu + delta
            dminus = 1.0 - mu + delta
            denom = dplus * dminus
            phase = eta * log(dplus / dminus)
            sigma_i = 2.0 * (eta^2) / (k^2 * denom) * spin_factor * cos(phase)
            acc += w[kq] * sigma_i
        end
        sigma_int = acc * Δ
    end

    sigma_terms = 0.0
    if cache.ltp == 1
        NL = cache.NL
        b = cache.b
        a = cache.a
        nb = length(b)
        na = length(a)

        Nq = 64
        μ_nodes, w = quadrature(Nq, "gauss-legendre")
        Δ = 0.5 * (mu_hi - mu_lo)
        μ_mid = 0.5 * (mu_hi + mu_lo)
        acc_ta = 0.0
        for kq in eachindex(μ_nodes)
            mu = Δ * μ_nodes[kq] + μ_mid
            mu_safe = clamp(mu, mu_lo, mu_hi_safe)
            if !identical
                phase = exp(im * eta * log((1.0 - mu_safe + delta) / 2.0))
                sum_a = 0.0 + 0.0im
                for l in 0:min(NL, na - 1)
                    sum_a += ((2 * l + 1) / 2) * a[l + 1] * legendre_polynomials(l, mu_safe)
                end
                term_a = (2.0 * eta / (1.0 - mu_safe + delta)) * real(phase * sum_a)
                acc_ta += w[kq] * (-term_a)
            else
                phase_minus = exp(im * eta * log((1.0 - mu_safe + delta) / 2.0))
                phase_plus = exp(im * eta * log((1.0 + mu_safe + delta) / 2.0))
                sum_a = 0.0 + 0.0im
                for l in 0:min(NL, na - 1)
                    Pl = legendre_polynomials(l, mu_safe)
                    pref_ta = ((1.0 + mu_safe + delta) * phase_minus + (-1.0)^l * (1.0 - mu_safe + delta) * phase_plus)
                    sum_a += pref_ta * ((2 * l + 1) / 2) * a[l + 1] * Pl
                end
                denom = (1.0 - mu_safe + delta) * (1.0 + mu_safe + delta)
                term_a = (2.0 * eta / denom) * real(sum_a)
                acc_ta += w[kq] * (-term_a)
            end
        end
        sigma_terms += acc_ta * Δ

        if !identical
            max_b = min(2 * NL, nb - 1)
            for l in 0:max_b
                coeff = ((2 * l + 1) / 2) * b[l + 1]
                sigma_terms += coeff * int_P_l(l, mu_lo, mu_hi)
            end
        else
            max_b = min(NL, nb - 1)
            for l in 0:max_b
                coeff = ((4 * l + 1) / 2) * b[l + 1]
                sigma_terms += coeff * int_P_l(2 * l, mu_lo, mu_hi)
            end
        end
    elseif cache.ltp == 2
        NL = cache.NL
        c = cache.c
        nc = length(c)

        Nq = 64
        μ_nodes, w = quadrature(Nq, "gauss-legendre")
        Δ = 0.5 * (mu_hi - mu_lo)
        μ_mid = 0.5 * (mu_hi + mu_lo)
        acc = 0.0
        for kq in eachindex(μ_nodes)
            mu = Δ * μ_nodes[kq] + μ_mid
            mu_safe = clamp(mu, mu_lo, mu_hi_safe)
            sigma_r = 0.0
            if !identical
                max_c = min(NL, nc - 1)
                for l in 0:max_c
                    sigma_r += ((2 * l + 1) / 2) * c[l + 1] * legendre_polynomials(l, mu_safe)
                end
                acc += w[kq] * (sigma_r / (1.0 - mu_safe + delta))
            else
                max_c = min(NL, nc - 1)
                for l in 0:max_c
                    sigma_r += ((4 * l + 1) / 2) * c[l + 1] * legendre_polynomials(2 * l, mu_safe)
                end
                acc += w[kq] * (sigma_r / ((1.0 - mu_safe + delta) * (1.0 + mu_safe + delta)))
            end
        end
        sigma_terms += acc * Δ
    else
        mu_tab = cache.mu_grid
        pni_tab = cache.pni_grid
        sigmaNI = cache.sigmaNI

        mu_end = min(mu_hi, mu_tab[end])
        if mu_end > mu_lo
            acc = 0.0
            for i in 1:(length(mu_tab) - 1)
                mu1s = mu_tab[i]
                mu2s = mu_tab[i + 1]
                if mu2s < mu_lo || mu1s > mu_end
                    continue
                end
                a = clamp(mu1s, mu_lo, mu_end)
                b = clamp(mu2s, mu_lo, mu_end)
                if b == a
                    continue
                end
                p1 = pni_tab[i]
                p2 = pni_tab[i + 1]
                slope = (p2 - p1) / (mu2s - mu1s)
                intercept = p1 - slope * mu1s
                acc += 0.5 * slope * (b^2 - a^2) + intercept * (b - a)
            end
            sigma_terms += sigmaNI * acc
            sigma_terms += negative_ltpgt2_correction(cache, mu_lo, mu_end)
        end
    end

    return 2.0 * pi * (sigma_c + sigma_int + sigma_terms)
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
    NL::Int, identical::Bool, spin_s::Float64, k::Float64, eta::Float64, delta::Float64)

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

# Output Argument(s)
- `σs::Float64` : elastic differential cross-section in barn/sr.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5, LTP=1, equations 6.13 and 6.14.
"""
function elastic_dsig_ltp1(mu::Float64;
                           a::Vector{ComplexF64}, b::Vector{Float64}, NL::Int,
                           identical::Bool, spin_s::Float64, k::Float64, eta::Float64, delta::Float64)

    eps_mu = 1.0e-12
    mu_hi_safe = delta > 0.0 ? 1.0 : 1.0 - eps_mu
    if identical
        mu_safe = clamp(mu, 0.0, mu_hi_safe)
    else
        mu_safe = clamp(mu, -1.0, mu_hi_safe)
    end

    
    sigma_c = coulomb_dcs(mu_safe, k, eta, delta, identical, spin_s)

    n_a = min(NL, length(a) - 1)

    if !identical
        phase = exp(im * eta * log((1.0 - mu_safe + delta) / 2.0))
        sum_a = 0.0 + 0.0im
        for l in 0:n_a
            Pl = legendre_polynomials(l, mu_safe)
            sum_a += ((2 * l + 1) / 2) * a[l + 1] * Pl
        end
        term_a = (2.0 * eta / (1.0 - mu_safe + delta)) * real(phase * sum_a)

        max_b = min(2 * NL, length(b) - 1)
        sum_b = 0.0
        for l in 0:max_b
            sum_b += ((2 * l + 1) / 2) * b[l + 1] * legendre_polynomials(l, mu_safe)
        end

        return sigma_c - term_a + sum_b
    else
        phase_minus = exp(im * eta * log((1.0 - mu_safe + delta) / 2.0))
        phase_plus = exp(im * eta * log((1.0 + mu_safe + delta) / 2.0))

        sum_a = 0.0 + 0.0im
        for l in 0:n_a
            Pl = legendre_polynomials(l, mu_safe)
            pref = ((1.0 + mu_safe + delta) * phase_minus + (-1.0)^l * (1.0 - mu_safe + delta) * phase_plus)
            sum_a += pref * ((2 * l + 1) / 2) * a[l + 1] * Pl
        end
        denom = (1.0 - mu_safe + delta) * (1.0 + mu_safe + delta)
        term_a = (2.0 * eta / denom) * real(sum_a)

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
    identical::Bool, spin_s::Float64, k::Float64, eta::Float64, delta::Float64)

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

# Output Argument(s)
- `σs::Float64` : elastic differential cross-section in barn/sr.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5, LTP=2, equations 6.15 to 6.18.
"""
function elastic_dsig_ltp2(mu::Float64;
                           c::Vector{Float64}, NL::Int,
                           identical::Bool, spin_s::Float64, k::Float64, eta::Float64, delta::Float64)

    eps_mu = 1.0e-12
    mu_hi_safe = delta > 0.0 ? 1.0 : 1.0 - eps_mu
    if identical
        mu_safe = clamp(mu, 0.0, mu_hi_safe)
    else
        mu_safe = clamp(mu, -1.0, mu_hi_safe)
    end

    
    sigma_c = coulomb_dcs(mu_safe, k, eta, delta, identical, spin_s)

    max_c = min(NL, length(c) - 1)
    sigma_r = 0.0
    if !identical
        for l in 0:max_c
            sigma_r += ((2 * l + 1) / 2) * c[l + 1] * legendre_polynomials(l, mu_safe)
        end
        return sigma_c + sigma_r / (1.0 - mu_safe + delta)
    else
        for l in 0:max_c
            sigma_r += ((4 * l + 1) / 2) * c[l + 1] * legendre_polynomials(2 * l, mu_safe)
        end
        return sigma_c + sigma_r / ((1.0 - mu_safe + delta)*(1.0 + mu_safe + delta))
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
    load_elastic_scattering_endf(endf_path::String, Z::Int, A::Int,
    particle::Particle; analytical::Bool=false, dcs_source::String="database",
    coulomb_kinematics::String="standard", screening::Bool=true)

Loads ENDF-6 elastic scattering data for a single isotope and builds the tabulated or
coefficient-side data used by the elastic-scattering database. The MF=6/MT=2 product
subsection is selected by the projectile `ZAP`, and the projectile spin is read from
the ENDF LAW=5 `SPI` field.

# Input Argument(s)
- `endf_path::String` : path to the ENDF-6 file.
- `Z::Int` : target atomic number.
- `A::Int` : target mass number.
- `particle::Particle` : incident particle type.
- `analytical::Bool` : boolean to keep coefficient data for analytical integration.
- `dcs_source::String` : source of the differential cross-section, either `"database"`
  or `"coefficients"`.
- `coulomb_kinematics::String` : Coulomb DCS kinematics convention, either `"standard"`
  for the ENDF nonrelativistic convention or `"relativistic"` for relativistic beta and
  center-of-mass momentum.
- `screening::Bool` : boolean to include Moliere screening in Coulomb terms. When
  enabled, the laboratory Moliere angle is transformed to the center-of-mass frame.

# Output Argument(s)
- `data::NamedTuple` : isotope elastic-scattering data containing energy grids,
  angular grids, tabulated cross-sections, interpolation data, projectile/target charge,
  ENDF projectile spin, projectile/target masses, and optional coefficient tables.

# Reference(s)
- ENDF-6 Formats Manual, MF=3 and MF=6, MT=2, LAW=5 charged-particle elastic scattering.
"""
function load_elastic_scattering_endf(endf_path::String, Z::Int, A::Int,
                                      particle::Particle;
                                      analytical::Bool=false,
                                      dcs_source::String="database",
                                      coulomb_kinematics::String="standard",
                                      screening::Bool=true)
    lines = readlines(endf_path)
    mf6mt2_subsections = read_mf6_mt(lines, 2; source=endf_path)
    store_database = dcs_source == "database"
    store_coefficients = analytical || dcs_source == "coefficients"

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

    Z1 = get_charge(particle)
    Z2 = Z
    isnothing(eladata.SPI) && error("MF=6 MT=2 LAW=5 subsection has no SPI projectile spin in $(endf_path).")
    spin_s = Float64(eladata.SPI)
    m_MeV = get_mass(particle)
    m_eV = m_MeV * 1.0e6
    M_MeV = get_mass(Z, A)
    M_eV = M_MeV * 1.0e6
    identical = (eladata.LIDP == 1)

    """
        simple_forward_mu_grid(mu_min, mu_max; n_coarse=50, n_dense=100, x_min=1e-8)

    Builds a non-uniform center-of-mass angular grid with a coarse linear part and a
    log-spaced forward-angle part near `mu = 1`.

    # Input Argument(s)
    - `mu_min` : lower center-of-mass scattering cosine.
    - `mu_max` : upper center-of-mass scattering cosine.
    - `n_coarse` : number of points in the coarse linear grid.
    - `n_dense` : number of points in the dense forward-angle grid.
    - `x_min` : smallest value of `1 - mu` in the dense forward-angle grid.

    # Output Argument(s)
    - `mu::Vector{Float64}` : sorted, unique center-of-mass angular grid.
    """
    function simple_forward_mu_grid(mu_min, mu_max;
                                    n_coarse=50,
                                    n_dense=100,
                                    x_min=1e-8)

        # --- coarse linear grid ---
        mu_coarse = range(mu_min, mu_max, length=n_coarse)

        # --- dense grid near mu = 1 ---
        # x = 1 - mu → log spacing in x
        x_max = 1.0 - mu_min
        x_dense = exp.(range(log(x_max), log(x_min), length=n_dense))
        mu_dense = 1.0 .- x_dense

        # keep only values inside range
        mu_dense = filter(m -> mu_min <= m <= mu_max, mu_dense)

        # merge + sort
        mu = unique(sort(vcat(mu_coarse, mu_dense)))
        return mu
    end

    if eladata.LTP == 1 && !isempty(eladata.ba_tables_ltp1)
        energies_eV = sort(collect(keys(eladata.ba_tables_ltp1)))
        mu_grid = Float64[]
        dcs = Vector{Vector{Float64}}()
        if store_database
            if identical
                mu_grid = simple_forward_mu_grid(0.0, 1.0, n_coarse=200, n_dense=200, x_min=1e-8)
            else
                mu_grid = simple_forward_mu_grid(-1.0, 1.0, n_coarse=300, n_dense=300, x_min=1e-8)
            end
            dcs = Vector{Vector{Float64}}(undef, length(energies_eV))
            for (idx, Ein) in enumerate(energies_eV)
                k_c = k_coulomb(Ein, m_eV, M_eV; kinematics=coulomb_kinematics)
                eta = eta_coulomb(Ein, m_eV, Z1, Z2; kinematics=coulomb_kinematics)
                delta = delta_coulomb(Ein, m_eV, M_eV, Z1, Z2; kinematics=coulomb_kinematics, screening=screening)
                NL, b, a = eladata.ba_tables_ltp1[Ein]
                dcs[idx] = [max(0.0, elastic_dsig_ltp1(mu; a=a, b=b, NL=NL, identical=identical, spin_s=spin_s, k=k_c, eta=eta, delta=delta))
                            for mu in mu_grid]
            end
        end
        energies_MeV = [E * 1.0e-6 for E in energies_eV]
        return (energies=energies_MeV, mu_grid=mu_grid, dcs=dcs,
            ltp=eladata.LTP, identical=identical, z1=Z1, z2=Z2,
            spin_s=spin_s, m_MeV=m_MeV, M_MeV=M_MeV,
            mf6_nbt=eladata.mf6_nbt, mf6_int=eladata.mf6_int,
            ba_tables_ltp1=store_coefficients ? eladata.ba_tables_ltp1 : nothing,
            c_tables_ltp2=nothing,
            mupni_tables=nothing,
            sigmaNI_E=nothing,
            sigmaNI_S=nothing,
            sigmaNI_nbt=nothing,
            sigmaNI_int=nothing)
    elseif eladata.LTP == 2 && !isempty(eladata.c_tables_ltp2)
        energies_eV = sort(collect(keys(eladata.c_tables_ltp2)))
        mu_grid = Float64[]
        dcs = Vector{Vector{Float64}}()
        if store_database
            if identical
                mu_grid = simple_forward_mu_grid(0.0, 1.0, n_coarse=200, n_dense=200, x_min=1e-8)
            else
                mu_grid = simple_forward_mu_grid(-1.0, 1.0, n_coarse=300, n_dense=300, x_min=1e-8)
            end
            dcs = Vector{Vector{Float64}}(undef, length(energies_eV))
            for (idx, Ein) in enumerate(energies_eV)
                k_c = k_coulomb(Ein, m_eV, M_eV; kinematics=coulomb_kinematics)
                eta = eta_coulomb(Ein, m_eV, Z1, Z2; kinematics=coulomb_kinematics)
                delta = delta_coulomb(Ein, m_eV, M_eV, Z1, Z2; kinematics=coulomb_kinematics, screening=screening)
                NL, c = eladata.c_tables_ltp2[Ein]
                dcs[idx] = [max(0.0, elastic_dsig_ltp2(mu; c=c, NL=NL, identical=identical, spin_s=spin_s, k=k_c, eta=eta, delta=delta))
                            for mu in mu_grid]
            end
        end
        energies_MeV = [E * 1.0e-6 for E in energies_eV]
        return (energies=energies_MeV, mu_grid=mu_grid, dcs=dcs,
            ltp=eladata.LTP, identical=identical, z1=Z1, z2=Z2,
            spin_s=spin_s, m_MeV=m_MeV, M_MeV=M_MeV,
            mf6_nbt=eladata.mf6_nbt, mf6_int=eladata.mf6_int,
            ba_tables_ltp1=nothing,
            c_tables_ltp2=store_coefficients ? eladata.c_tables_ltp2 : nothing,
            mupni_tables=nothing,
            sigmaNI_E=nothing,
            sigmaNI_S=nothing,
            sigmaNI_nbt=nothing,
            sigmaNI_int=nothing)
    elseif eladata.LTP > 2 && !isempty(eladata.mupni_tables)
        (_ZA3, _AWR3, E3, S3, NBT3, INT3) = read_mf3_mt(lines, 2)
        mf6_energies = sort(collect(keys(eladata.mupni_tables)))
        isempty(mf6_energies) &&
            error("MF=6 MT=2 contains no LTP>2 mu/PNI tables in $(endf_path).")
        mu_ref, _ = eladata.mupni_tables[mf6_energies[1]]
        if isempty(mu_ref) || abs(mu_ref[1] + 1.0) > 1.0e-12
            error("LTP>2 mu grid must start at -1.0 in $(endf_path).")
        end

        same_mu_grid = true
        for Ein in mf6_energies[2:end]
            mu_check, _ = eladata.mupni_tables[Ein]
            if isempty(mu_check) || abs(mu_check[1] + 1.0) > 1.0e-12
                error("LTP>2 mu grid must start at -1.0 in $(endf_path).")
            end
            if length(mu_check) != length(mu_ref)
                same_mu_grid = false
                break
            end
            for i in eachindex(mu_check)
                if abs(mu_check[i] - mu_ref[i]) > 1.0e-8
                    same_mu_grid = false
                    break
                end
            end
        end

        if same_mu_grid
            mu_grid = mu_ref
            pni_tables = Dict{Float64,Vector{Float64}}()
            for Ein in mf6_energies
                _, pni_tab = eladata.mupni_tables[Ein]
                pni_tables[Ein] = pni_tab
            end
            mupni_tables = pni_tables
        else
            mu_grid = Float64[]
            mupni_tables = eladata.mupni_tables
        end

        energies_eV = sort(union(mf6_energies, E3))

        function pni_mu_at_energy(Ein::Float64)
            if same_mu_grid
                return mu_grid, mupni_tables[Ein]
            end
            return mupni_tables[Ein]
        end

        dcs = Vector{Vector{Float64}}()
        if store_database
            dcs = Vector{Vector{Float64}}(undef, length(energies_eV))
            for (idx, Ein) in enumerate(energies_eV)
                k_c = k_coulomb(Ein, m_eV, M_eV; kinematics=coulomb_kinematics)
                eta = eta_coulomb(Ein, m_eV, Z1, Z2; kinematics=coulomb_kinematics)
                delta = delta_coulomb(Ein, m_eV, M_eV, Z1, Z2; kinematics=coulomb_kinematics, screening=screening)
                if haskey(mupni_tables, Ein)
                    mu_tab, pni_tab = pni_mu_at_energy(Ein)
                    sigmaNI = interp_TAB1(Ein, E3, S3, NBT3, INT3)
                    dcs[idx] = [max(0.0, elastic_dsig_ltpgreater2(mu; pni=pni_tab[imu], sigmaNI=sigmaNI, identical=identical, spin_s=spin_s, k=k_c, eta=eta, delta=delta))
                                for (imu, mu) in enumerate(mu_tab)]
                    continue
                end

                idx_hi = searchsortedfirst(mf6_energies, Ein)
                if idx_hi == 1 || idx_hi > length(mf6_energies)
                    error("No bracketing LTP>2 tables for E=$(Ein) in $(endf_path).")
                end
                E1 = mf6_energies[idx_hi - 1]
                E2 = mf6_energies[idx_hi]
                mu1, pni1 = pni_mu_at_energy(E1)
                mu2, pni2 = pni_mu_at_energy(E2)
                if !same_mu_grid
                    if length(mu1) != length(mu2)
                        error("LTP>2 bracketing mu grids differ for interpolation in $(endf_path).")
                    end
                    for i in eachindex(mu1)
                        if abs(mu1[i] - mu2[i]) > 1.0e-12
                            error("LTP>2 bracketing mu grids differ for interpolation in $(endf_path).")
                        end
                    end
                end

                itp = select_interpolation(eladata.mf6_nbt, eladata.mf6_int, idx_hi - 1)
                pni_tab = Vector{Float64}(undef, length(mu1))
                for k in eachindex(mu1)
                    pni_tab[k] = interp_pair(Ein, E1, pni1[k], E2, pni2[k], itp)
                end

                sigmaNI = interp_TAB1(Ein, E3, S3, NBT3, INT3)
                dcs[idx] = [max(0.0, elastic_dsig_ltpgreater2(mu; pni=pni_tab[imu], sigmaNI=sigmaNI, identical=identical, spin_s=spin_s, k=k_c, eta=eta, delta=delta))
                            for (imu, mu) in enumerate(mu1)]

                if mupni_tables isa Dict{Float64,Vector{Float64}}
                    mupni_tables[Ein] = pni_tab
                else
                    mupni_tables[Ein] = (mu1, pni_tab)
                end
            end
        end

        energies_MeV = [E * 1.0e-6 for E in energies_eV]
        return (energies=energies_MeV, mu_grid=mu_grid, dcs=dcs,
            ltp=eladata.LTP, identical=identical, z1=Z1, z2=Z2,
            spin_s=spin_s, m_MeV=m_MeV, M_MeV=M_MeV,
            mf6_nbt=eladata.mf6_nbt, mf6_int=eladata.mf6_int,
            ba_tables_ltp1=nothing,
            c_tables_ltp2=nothing,
            mupni_tables=(store_coefficients || (store_database && !same_mu_grid)) ? mupni_tables : nothing,
            sigmaNI_E=store_coefficients ? E3 : nothing,
            sigmaNI_S=store_coefficients ? S3 : nothing,
            sigmaNI_nbt=store_coefficients ? NBT3 : nothing,
            sigmaNI_int=store_coefficients ? INT3 : nothing)
    end

    error("Unsupported LTP=$(eladata.LTP) or missing elastic scattering data in $(endf_path).")
end


"""
    endf_path_for_isotope(db_name::String,Z::Int,A::Int; data_root::Union{Nothing,String}=nothing)

Builds the ENDF file path for a target isotope in a Radiant ENDF database directory.

# Input Argument(s)
- `db_name::String` : ENDF database directory name.
- `Z::Int` : target atomic number.
- `A::Int` : target mass number.
- `data_root::Union{Nothing,String}` : optional root directory containing ENDF databases.

# Output Argument(s)
- `path::String` : expected ENDF file path for isotope `(Z, A)`.

# Reference(s)
- ENDF database file naming convention used by Radiant.
"""
function endf_path_for_isotope(db_name::String,Z::Int,A::Int; data_root::Union{Nothing,String}=nothing)
    root = isnothing(data_root) ? normpath(joinpath(@__DIR__, "..", "..", "data")) : data_root
    symbol = atomic_symbol(Z)
    filename = @sprintf("p_%03d-%s-%03d.endf", Z, symbol, A)
    return joinpath(root, db_name, filename)
end
"""
    init_elastic_scattering_endf_db(db_name::String, isotopes::Vector{Tuple{Int,Int}},
    particle::Particle; data_root::Union{Nothing,String}=nothing, analytical::Bool=false,
    dcs_source::String="database", coulomb_kinematics::String="standard",
    screening::Bool=true)

Initializes the ENDF elastic-scattering database for a set of target isotopes and one
incident particle. The database stores the Coulomb settings and common screening flag;
each isotope entry stores the selected LAW=5 data, projectile/target masses, and
projectile spin read from ENDF `SPI`.

# Input Argument(s)
- `db_name::String` : ENDF database directory name.
- `isotopes::Vector{Tuple{Int,Int}}` : list of `(Z, A)` target isotopes.
- `particle::Particle` : incident particle type.
- `data_root::Union{Nothing,String}` : optional root directory containing ENDF databases.
- `analytical::Bool` : boolean to keep coefficient data for analytical integration.
- `dcs_source::String` : source of the differential cross-section, either `"database"`
  or `"coefficients"`.
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
    particle::Particle; data_root::Union{Nothing,String}=nothing, analytical::Bool=false,
    dcs_source::String="database", coulomb_kinematics::String="standard",
    screening::Bool=true)
    root = isnothing(data_root) ? normpath(joinpath(@__DIR__, "..", "..", "data")) : data_root
    dcs_source_lc = lowercase(dcs_source)
    coulomb_kinematics_lc = lowercase(coulomb_kinematics)
    db = Dict{Tuple{Int,Int},IsotopeElasticDB}()
    for (Z, A) in isotopes
        endf_path = endf_path_for_isotope(db_name, Z, A; data_root=root)
        if !isfile(endf_path)
            error("ENDF file not found: $(endf_path)")
        end
        data = load_elastic_scattering_endf(endf_path, Z, A, particle;
                    analytical=analytical, dcs_source=dcs_source_lc,
                    coulomb_kinematics=coulomb_kinematics_lc, screening=screening)
        db[(Z, A)] = IsotopeElasticDB(data.energies, data.mu_grid, data.dcs,
              data.ltp, data.identical, data.z1, data.z2, data.spin_s,
              data.m_MeV, data.M_MeV, data.mf6_nbt, data.mf6_int,
              data.ba_tables_ltp1, data.c_tables_ltp2, data.mupni_tables,
              data.sigmaNI_E, data.sigmaNI_S, data.sigmaNI_nbt, data.sigmaNI_int)
    end
    return ElasticScatteringENDFDB(dcs_source_lc, coulomb_kinematics_lc, screening, db)
end

