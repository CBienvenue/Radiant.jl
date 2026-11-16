# Elastic scattering of protons by neutral atoms as Geant4 11.4 simulates it (model "G4" of
# Elastic_Scattering), with the reference physics lists (FTFP_BERT, QGSP_BIC, QBBC, ...) in
# their `_SS` variants, where Coulomb scattering is sampled collision by collision. Geant4 has
# two independent processes, so the whole cross-section is their incoherent sum, without
# Coulomb-nuclear interference:
#     dσ/dΩ_cm = dσ_C/dΩ_lab · dΩ_lab/dΩ_cm + σ_BGG · f_CHIPS(Q²) · p_cm²/π.
#   * Coulomb: G4CoulombScattering → G4eCoulombScatteringModel → G4WentzelOKandVIxSection.
#     Wentzel cross-section in z = 1 - cos θ_lab with the Molière screening parameter A_s,
#         dσ/dz = Z kinFactor / (z + A_s)² · g(z),
#     times the rejection function of SampleSingleScattering, clamped to [0, 1]:
#         g(z) = [1 - z β²/2 + (πα/2) Z β √(z/2) (2 - z)] · (1 + F p² z)⁻² / (1 + z p/M),
#     a spin (McKinley-Feshbach-type) factor, the exponential nuclear form factor of
#     Butkevich et al. and a recoil factor. Hydrogen is restricted to θ_lab ≤ 90°.
#     Scattering off the atomic electrons is part of ionisation and is left out.
#   * Hadronic: G4HadronElasticProcess with the cross-section of G4BGGNucleonElasticXS
#     (below 91 GeV: A σ_el(pp) for Z = 1; for Z > 1 the Barashenkov element cross-section,
#     scaled by the Coulomb-barrier factor at or below 14 MeV) and the Q² = -t distribution of
#     G4ChipsElasticModel (G4ChipsProtonElasticXS::GetExchangeT, written here as a density).
# References: Geant4 Physics Reference Manual 11.4, Sec. 9.2.1 and Ch. 20-21, and the
# classes above in the Geant4 11.4 source. The tables of G4BarashenkovData.hh,
# G4IsotopeList.hh (aeff) and G4NuclearRadii (r0) are in
# `data/elastic_scattering_g4_11_4.jld2`.
#
# The cross-section of one isotope at one energy is tabulated in x = 1 - μ_cm as for the OP
# description (`dcs_table_start_grid`, `tabulate_power_law_dcs`), in a
# `TabulatedElasticDCSCache`, so that the integration kernels shared with OP (end of
# structures/Elastic_Scattering.jl) apply. For p + ¹H the cross-section is folded onto μ_cm ≥ 0: the faster proton is the
# scattered one.

#----
# Constants
#----

const G4_DATA_FILE = "elastic_scattering_g4_11_4.jld2"
const G4_BGG_LOW_ENERGY = 14.0                     # [MeV]
const G4_BGG_GLAUBER_ENERGY = 91.0e3               # [MeV]
const G4_PROTON_RADIUS = 0.895                     # G4NuclearRadii [fm]
const G4_GEV2 = 1.0e6                              # GeV² in MeV²

# CHIPS tabulation of the parameters in ln p (p in GeV/c)
const G4_CHIPS_NPOINTS = 128
const G4_CHIPS_LPMIN = -8.0
const G4_CHIPS_LPMAX = 8.0
const G4_CHIPS_DLNP = (G4_CHIPS_LPMAX - G4_CHIPS_LPMIN) / (G4_CHIPS_NPOINTS - 1)

#----
# Tables
#----

"Mean atomic mass of element Z [u] used by Geant4 (G4IsotopeList aeff, Z ≤ 100)."
g4_aeff(Z::Int) = fast_load(G4_DATA_FILE)["aeff"][min(Z, 100) + 1]

"Nuclear radius constant r0 of G4NuclearRadii [fm] (Z ≤ 92)."
g4_r0(Z::Int) = fast_load(G4_DATA_FILE)["r0"][min(Z, 92) + 1]

#----
# Coulomb scattering (G4WentzelOKandVIxSection)
#----

"""
    g4_screen_rsquare(Z::Int, me::Float64)

Gives the Thomas-Fermi screening factor ScreenRSquare of G4WentzelOKandVIxSection::InitialiseA,
½ α² (mₑc²/0.88534)² (1 + e^(-Z²/1000)) Z^(2/3), with ScreeningFactor = 1. For hydrogen the
factor (1 + e^(-Z²/1000)) is left out. The screening parameter is A_s = ScreenRSquare/p²
times the Molière factor.

# Input Argument(s)
- `Z::Int` : target atomic number.
- `me::Float64` : electron rest energy [MeV].

# Output Argument(s)
- `R2::Float64` : screening factor [MeV²].
"""
function g4_screen_rsquare(Z::Int, me::Float64)
    a0 = me / 0.88534
    afact = 0.5 * ALPHA_FS^2 * a0^2
    Z == 1 && return afact
    x = Z^(1.0 / 3.0)
    return afact * (1.0 + exp(-Z * Z * 0.001)) * x * x
end

"Exponential nuclear form-factor coefficient F of G4WentzelOKandVIxSection [MeV⁻²]."
function g4_formfactor(Z::Int)
    Z == 1 && return 3.097e-6
    x = g4_aeff(Z)^0.27                            # G4NistManager::GetA27
    return 6.937e-6 * x * x
end

"""
    g4_wentzel(iso::IsotopeG4Model, E::Float64, screening::Bool)

Gives the parameters of the Geant4 Coulomb single-scattering cross-section of a proton on one
isotope (G4WentzelOKandVIxSection::SetupKinematic, SetupTarget and SetTargetMass).

# Input Argument(s)
- `iso::IsotopeG4Model` : isotope model.
- `E::Float64` : laboratory kinetic energy [MeV].
- `screening::Bool` : whether the Coulomb field is screened (A_s = 0 otherwise).

# Output Argument(s)
- `w::NamedTuple` : `Z`, p² `mom2` [MeV²], 1/β² `invbeta2`, screening parameter `screenZ`,
  form-factor term F p² `formfactA`, spin term β²/2 `factB`, recoil term p/M `factD`,
  `kinFactor` [fm²] and the largest z, `zmax`.
"""
function g4_wentzel(iso::IsotopeG4Model, E::Float64, screening::Bool)
    Z = iso.Z
    mom2 = momentum_lab(E, iso.m_mev; kinematics="relativistic")^2
    invbeta2 = 1.0 + iso.m_mev^2 / mom2
    screenZ = if !screening
        0.0
    elseif Z == 1
        iso.screen_rsquare / mom2
    else
        min(Z * 1.13, 1.13 + 3.76 * Z * Z * invbeta2 * ALPHA_FS^2) * iso.screen_rsquare / mom2
    end
    return (Z=Z, mom2=mom2, invbeta2=invbeta2, screenZ=screenZ, formfactA=iso.formfactor * mom2,
            factB=0.5 / invbeta2, factD=sqrt(mom2) / iso.M_mev,
            kinFactor=2.0 * π * OP_E2^2 * Z * invbeta2 / mom2, zmax=(Z == 1 ? 1.0 : 2.0))
end

"Rejection function g(z) of G4WentzelOKandVIxSection::SampleSingleScattering."
function g4_wentzel_rejection(w::NamedTuple, z::Float64)
    fm = 1.0 / (1.0 + w.formfactA * z)^2
    return (1.0 - z * w.factB + 0.5 * π * ALPHA_FS * w.Z * sqrt(z * w.factB) * (2.0 - z)) * fm /
           (1.0 + z * w.factD)
end

"Coulomb cross-section dσ/dΩ_lab [b/sr] at z = 1 - cos θ_lab: Wentzel times g clamped to [0, 1]."
function g4_coulomb_dcs_lab(w::NamedTuple, z::Float64)
    (0.0 <= z <= w.zmax) || return 0.0
    dσdz = w.Z * w.kinFactor / (z + w.screenZ)^2 * clamp(g4_wentzel_rejection(w, z), 0.0, 1.0)
    return dσdz / (2.0 * π) * OP_FM2_TO_BARN
end

#----
# Hadronic elastic cross-section (G4BGGNucleonElasticXS)
#----

"G4NuclearRadii::ExplicitRadius [fm] (0 when not explicit)."
function g4_explicit_radius(Z::Int, A::Int)
    Z <= 4 || return 0.0
    A == 1 && return G4_PROTON_RADIUS
    A == 2 && return 2.13
    Z == 1 && A == 3 && return 1.80
    Z == 2 && A == 3 && return 1.96
    Z == 2 && A == 4 && return 1.68
    Z == 3 && return 2.40
    Z == 4 && return 2.51
    return 0.0
end

"G4NuclearRadii::RadiusCB [fm]."
function g4_radius_cb(Z::Int, A::Int)
    R = g4_explicit_radius(Z, A)
    R == 0.0 && (R = g4_r0(Z) * A^(1.0 / 3.0))
    return R
end

"""
    g4_coulomb_factor(Z::Int, A::Int, E::Float64, m::Float64, M::Float64; pp::Bool=false)

Gives the Coulomb-barrier factor of G4NuclearRadii::CoulombFactor, 1 - B_C/T_cm with
B_C = ½ e² Z/(R_p + R_A), or 0 below the barrier.

# Input Argument(s)
- `Z::Int`, `A::Int` : target.
- `E::Float64` : laboratory kinetic energy [MeV].
- `m::Float64`, `M::Float64` : projectile and target rest energies [MeV].
- `pp::Bool` : proton target (radius of the proton).

# Output Argument(s)
- `CF::Float64` : Coulomb factor.
"""
function g4_coulomb_factor(Z::Int, A::Int, E::Float64, m::Float64, M::Float64; pp::Bool=false)
    tR = pp ? G4_PROTON_RADIUS : g4_radius_cb(Z, A)
    totTcm = sqrt(m^2 + M^2 + 2.0 * (E + m) * M) - m - M
    bC = 0.5 * OP_E2 * Z / (G4_PROTON_RADIUS + tR)
    return totTcm > bC ? 1.0 - bC / totTcm : 0.0
end

"Elastic p + p cross-section [b] of G4HadronNucleonXsc::HadronNucleonXscNS (p_lab < 10 GeV/c)."
function g4_pp_elastic(E::Float64, m::Float64)
    ekin = max(E, 0.1)
    pLab = momentum_lab(ekin, m; kinematics="relativistic") / 1000.0
    pLab < 10.0 || error("The Geant4 p + p elastic cross-section is implemented below 10 GeV/c only.")
    if pLab < 0.73
        tot = 23.0 + 50.0 * sqrt(log(0.73 / pLab)^7)
        el = tot
    elseif pLab < 1.05
        x = log(pLab / 0.73)
        tot = 23.0 + 40.0 * x^2
        el = 23.0 + 20.0 * x^2
    else
        tot = 39.0 + 75.0 * (pLab - 1.2) / (pLab^3 + 0.15)
        el = 6.0 + 20.0 / ((log(pLab) - 0.182)^2 + 1.0)
    end
    el = min(el, tot) * 1.0e-3                     # mb -> b
    ekin < 100.0 && (el *= g4_coulomb_factor(1, 1, ekin, m, m; pp=true))
    return el
end

"G4PiData::TotalXSection / ReactionXSection: linear in E, constant outside the table."
function g4_pidata_interp(E::AbstractVector{Float64}, X::AbstractVector{Float64}, e::Float64)
    e <= E[1] && return X[1]
    e >= E[end] && return X[end]
    i = max(searchsortedfirst(E, e), 2)
    return max(0.0, X[i - 1] + (e - E[i - 1]) * (X[i] - X[i - 1]) / (E[i] - E[i - 1]))
end

"""
    g4_barashenkov_elastic(Z::Int, E::Float64)

Gives the proton elastic cross-section σ_tot - σ_inel of G4ComponentBarNucleonNucleusXsc:
the Barashenkov tables of 17 elements, interpolated in Z with an A^(2/3) scaling.

# Input Argument(s)
- `Z::Int` : target atomic number (92 above).
- `E::Float64` : laboratory kinetic energy [MeV].

# Output Argument(s)
- `σ::Float64` : elastic cross-section [b].
"""
function g4_barashenkov_elastic(Z::Int, E::Float64)
    tables = fast_load(G4_DATA_FILE)["barashenkov"]
    Z = min(Z, 92)
    e = E / 1000.0                                 # [GeV]
    zs = sort(collect(keys(tables)))
    it = something(findfirst(>=(Z), zs), length(zs))
    tot(z) = g4_pidata_interp(view(tables[z], :, 1), view(tables[z], :, 2), e)
    inel(z) = g4_pidata_interp(view(tables[z], :, 1), view(tables[z], :, 3), e)
    if zs[it] == Z
        σt, σi = tot(Z), inel(Z)
    else
        it == 1 && (it = 2)
        Z1, Z2 = zs[it - 1], zs[it]
        A75(z) = g4_aeff(z)^(2.0 / 3.0)
        alp1 = g4_aeff(Z) - g4_aeff(Z1)
        alp2 = g4_aeff(Z2) - g4_aeff(Z)
        interp(x1, x2) = (x1 * A75(Z) / A75(Z1) * alp2 + x2 * A75(Z) / A75(Z2) * alp1) / (alp1 + alp2)
        σt = interp(tot(Z1), tot(Z2))
        σi = interp(inel(Z1), inel(Z2))
    end
    return max(σt - σi, 0.0) * 1.0e-3              # mb -> b
end

"""
    g4_bgg_elastic(iso::IsotopeG4Model, E::Float64)

Gives the hadronic elastic cross-section of a proton on one isotope as G4BGGNucleonElasticXS
gives it to the elastic process below 91 GeV: A σ_el(pp) for Z = 1 (isotope cross-section);
otherwise the element cross-section, the same for every isotope of Z, from the Barashenkov
tables above 14 MeV and their value at 14 MeV scaled by the Coulomb factor below.

# Input Argument(s)
- `iso::IsotopeG4Model` : isotope model.
- `E::Float64` : laboratory kinetic energy [MeV].

# Output Argument(s)
- `σ::Float64` : hadronic elastic cross-section [b].
"""
function g4_bgg_elastic(iso::IsotopeG4Model, E::Float64)
    E > G4_BGG_GLAUBER_ENERGY && error("The Geant4 hadronic elastic cross-section is implemented below 91 GeV only.")
    iso.Z == 1 && return iso.A * g4_pp_elastic(E, iso.m_mev)
    Z = min(iso.Z, 92)
    E > G4_BGG_LOW_ENERGY && return g4_barashenkov_elastic(Z, E)
    # BGG target: theA[Z] = lrint(mean atomic mass), which may be unstable (⁵⁹Ni, ⁷⁹Se, ...)
    Ab = round(Int, g4_aeff(Z))
    Mb = get_mass(Z, Ab)
    fac = g4_barashenkov_elastic(Z, G4_BGG_LOW_ENERGY) /
          g4_coulomb_factor(Z, Ab, G4_BGG_LOW_ENERGY, iso.m_mev, Mb)
    return fac * g4_coulomb_factor(Z, Ab, E, iso.m_mev, Mb)
end

#----
# Hadronic angular distribution (G4ChipsElasticModel, G4ChipsProtonElasticXS)
#----

"Target (Z, N) of the CHIPS parameterisation for (Z, A): ³H → ²H and ³He → ⁴He (G4ChipsElasticModel)."
function g4_chips_target(Z::Int, A::Int)
    N = A - Z
    if Z == 1 && N == 2
        N = 1
    elseif Z == 2 && N == 1
        N = 2
    end
    return Z, N
end

"""
    g4_chips_parameters(tgZ::Int, tgN::Int)

Gives the parameters of the CHIPS proton elastic parameterisation for one target, lastPAR of
G4ChipsProtonElasticXS::GetPTables (index i + 1 for the source index i).

# Input Argument(s)
- `tgZ::Int`, `tgN::Int` : CHIPS target.

# Output Argument(s)
- `P::Vector{Float64}` : parameters.
"""
function g4_chips_parameters(tgZ::Int, tgN::Int)
    P = zeros(G4_CHIPS_NPOINTS)
    if tgZ == 1 && tgN == 0
        P[1:32] .= [2.865, 18.9, .6461, 3., 9., .425, .4276, .0022, 5., 74., 3., 3.4, .2, .17,
                    .001, 8., .055, 3.64, 5.e-5, 4000., 1500., .46, 1.2e6, 3.5e6, 5.e-5, 1.e10,
                    8.5e8, 1.e10, 1.1, 3.4e6, 6.8e6, 0.]
        return P
    end
    tgZ == 0 && error("The CHIPS proton elastic parameterisation of neutron targets is not implemented.")
    a = Float64(tgZ + tgN)
    sa = sqrt(a); ssa = sqrt(sa); asa = a * sa
    a2 = a * a; a3 = a2 * a; a4 = a3 * a; a5 = a4 * a; a6 = a4 * a2; a7 = a6 * a; a8 = a7 * a
    a9 = a8 * a; a10 = a5 * a5; a12 = a6 * a6; a14 = a7 * a7; a16 = a8 * a8; a17 = a16 * a
    a20 = a16 * a4; a32 = a16 * a16
    p(i, v) = (P[i + 1] = v)
    p(0, 5 / (1 + 22 / asa))
    p(1, 4.8 * a^1.14 / (1 + 3.6 / a3))
    p(2, 1 / (1 + 4.E-3 * a4) + 2.E-6 * a3 / (1 + 1.3E-6 * a3))
    p(3, 1.3 * a)
    p(4, 3.E-8 * a3 / (1 + 4.E-7 * a4))
    p(5, .07 * asa / (1 + .009 * a2))
    p(6, (3 + 3.E-16 * a20) / (1 + a20 * (2.E-16 / a + 3.E-19 * a)))
    p(7, (5.E-9 * a4 * sa + .27 / a) / (1 + 5.E16 / a20) / (1 + 6.E-9 * a4) + .015 / a2)
    p(8, (.001 * a + .07 / a) / (1 + 5.E13 / a16 + 5.E-7 * a3) + .0003 / sa)
    if a < 6.5
        a28 = a16 * a12
        p(9, 4000 * a)
        p(10, 1.2e7 * a8 + 380 * a17)
        p(11, .7 / (1 + 4.e-12 * a16))
        p(12, 2.5 / a8 / (a4 + 1.e-16 * a32))
        p(13, .28 * a)
        p(14, 1.2 * a2 + 2.3)
        p(15, 3.8 / a)
        p(16, .01 / (1 + .0024 * a5))
        p(17, .2 * a)
        p(18, 9.e-7 / (1 + .035 * a5))
        p(19, (42 + 2.7e-11 * a16) / (1 + .14 * a))
        p(20, 2.25 * a3)
        p(21, 18.)
        p(22, 2.4e-3 * a8 / (1 + 2.6e-4 * a7))
        p(23, 3.5e-36 * a32 * a8 / (1 + 5.e-15 * a32 / a))
        p(24, 1.e5 / (a8 + 2.5e12 / a16))
        p(25, 8.e7 / (a12 + 1.e-27 * a28 * a28))
        p(26, .0006 * a3)
        p(27, 10 + 4.e-8 * a12 * a)
        p(28, .114)
        p(29, .003)
        p(30, 2.e-23)
        p(31, 1 / (1 + .0001 * a8))
        p(32, 1.5e-4 / (1 + 5.e-6 * a12))
        p(33, .03)
        p(34, a / 2)
        p(35, 2.e-7 * a4)
        p(36, 4.)
        p(37, 64 / a3)
        p(38, 1.e8 * exp(.32 * asa))
        p(39, 20 * exp(.45 * asa))
        p(40, 7.e3 + 2.4e6 / a5)
        p(41, 2.5e5 * exp(.085 * a3))
        p(42, 2.5 * a)
        p(43, 920 + .03 * a8 * a3)
        p(44, 93 + .0023 * a12)
    else
        p1a10 = 2.2e-28 * a10
        r4a16 = 6.e14 / a16
        s4a16 = r4a16 * r4a16
        p(9, 4.5 * a^1.15)
        p(10, .06 * a^.6)
        p(11, .6 * a / (1 + 2.e15 / a16))
        p(12, .17 / (a + 9.e5 / a3 + 1.5e33 / a32))
        p(13, (.001 + 7.e-11 * a5) / (1 + 4.4e-11 * a5))
        p(14, (p1a10 * p1a10 + 2.e-29) / (1 + 2.e-22 * a12))
        p(15, 400 / a12 + 2.e-22 * a9)
        p(16, 1.e-32 * a12 / (1 + 5.e22 / a14))
        p(17, 1000 / a2 + 9.5 * sa * ssa)
        p(18, 4.e-6 * a * asa + 1.e11 / a16)
        p(19, (120 / a + .002 * a2) / (1 + 2.e14 / a16))
        p(20, 9 + 100 / a)
        p(21, .002 * a3 + 3.e7 / a6)
        p(22, 7.e-15 * a4 * asa)
        p(23, 9000 / a4)
        p(24, .0011 * asa / (1 + 3.e34 / a32 / a4))
        p(25, 1.e-5 * a2 + 2.e14 / a16)
        p(26, 1.2e-11 * a2 / (1 + 1.5e19 / a12))
        p(27, .016 * asa / (1 + 5.e16 / a16))
        p(28, .002 * a4 / (1 + 7.e7 / (a - 6.83)^14))
        p(29, 2.e6 / a6 + 7.2 / a^.11)
        p(30, 11 * a3 / (1 + 7.e23 / a16 / a8))
        p(31, 100 / asa)
        p(32, (.1 + 4.4e-5 * a2) / (1 + 5.e5 / a4))
        p(33, 3.5e-4 * a2 / (1 + 1.e8 / a8))
        p(34, 1.3 + 3.e5 / a4)
        p(35, 500 / (a2 + 50) + 3)
        p(36, 1.e-9 / a + s4a16 * s4a16)
        p(37, .4 * asa + 3.e-9 * a6)
        p(38, .0005 * a5)
        p(39, .002 * a5)
        p(40, 10.)
        p(41, .05 + .005 * a)
        p(42, 7.e-8 / sa)
        p(43, .8 * sa)
        p(44, .02 * sa)
        p(45, 1.e8 / a3)
        p(46, 3.e32 / (a32 + 1.e32))
        p(47, 24.)
        p(48, 20 / sa)
        p(49, 7.e3 * a / (sa + 1))
        p(50, 900 * sa / (1 + 500 / a3))
    end
    p(51, 1.e15 + 2.e27 / a4 / (1 + 2.e-18 * a16))
    return P
end

"(σ_CHIPS [mb], SS, S1, B1, S2, B2, S3, B3, S4, B4) of G4ChipsProtonElasticXS::GetTabValues at ln p (p in GeV/c)."
function g4_chips_tab_values(P::Vector{Float64}, lp::Float64, tgZ::Int, tgN::Int)
    L(i) = P[i + 1]
    p = exp(lp); sp = sqrt(p); p2 = p * p; p3 = p2 * p; p4 = p3 * p
    if tgZ == 1 && tgN == 0
        p2s = p2 * sp
        dl2 = lp - L(8)
        SS = L(31)
        S1 = (L(9) + L(10) * dl2 * dl2) / (1 + L(11) / p4 / p) + (L(12) / p2 + L(13) * p) / (p4 + L(14) * sp)
        B1 = L(15) * p^L(16) / (1 + L(17) / p3)
        S2 = L(18) + L(19) / (p4 + L(20) * p)
        B2 = L(21) + L(22) / (p4 + L(23) / sp)
        S3 = L(24) + L(25) / (p4 * p4 + L(26) * p2 + L(27))
        B3 = L(28) + L(29) / (p4 + L(30))
        dl1 = lp - L(3)
        sig = L(0) / p2s / (1 + L(7) / p2s) + (L(1) + L(2) * dl1 * dl1 + L(4) / p) / (1 + L(5) * lp) / (1 + L(6) / p4)
        return (sig, SS, S1, B1, S2, B2, S3, B3, 0.0, 0.0)
    end
    p5 = p4 * p; p6 = p5 * p; p8 = p6 * p2; p10 = p8 * p2; p12 = p10 * p2; p16 = p8 * p8
    dl = lp - 5.0
    a = Float64(tgZ + tgN)
    if a < 6.5
        pah = p^(a / 2); pa = pah * pah; pa2 = pa * pa
        S1 = L(9) / (1 + L(10) * p4 * pa) + L(11) / (p4 + L(12) * p4 / pa2) + (L(13) * dl * dl + L(14)) / (1 + L(15) / p2)
        B1 = (L(16) + L(17) * p2) / (p4 + L(18) / pah) + L(19)
        SS = L(20) / (1 + L(21) / p2) + L(22) / (p6 / pa + L(23) / p16)
        S2 = L(24) / (pa / p2 + L(25) / p4) + L(26)
        B2 = L(27) * p^L(28) + L(29) / (p8 + L(30) / p16)
        S3 = L(31) / (pa * p + L(32) / pa) + L(33)
        B3 = L(34) / (p3 + L(35) / p6) + L(36) / (1 + L(37) / p2)
        S4 = p2 * (pah * L(38) * exp(-pah * L(39)) + L(40) / (1 + L(41) * p^L(42)))
        B4 = L(43) * pa / p2 / (1 + pa * L(44))
    else
        S1 = L(9) / (1 + L(10) / p4) + L(11) / (p4 + L(12) / p2) + L(13) / (p5 + L(14) / p16)
        B1 = (L(15) / p8 + L(19)) / (p + L(16) / p^L(20)) + L(17) / (1 + L(18) / p4)
        SS = L(21) / (p4 / p^L(23) + L(22) / p4)
        S2 = L(24) / p4 / (p^L(25) + L(26) / p12) + L(27)
        B2 = L(28) / p^L(29) + L(30) / p^L(31)
        S3 = L(32) / p^L(35) / (1 + L(36) / p12) + L(33) / (1 + L(34) / p6)
        B3 = L(37) / p8 + L(38) / p2 + L(39) / (1 + L(40) / p8)
        S4 = (L(41) / p4 + L(46) / p) / (1 + L(42) / p10) + (L(43) + L(44) * dl * dl) / (1 + L(45) / p12)
        B4 = L(47) / (1 + L(48) / p) + L(49) * p4 / (1 + L(50) * p5)
    end
    sig = (L(0) * dl * dl + L(1)) / (1 + L(2) / p + L(5) / p6) + L(3) / (p3 + L(4) / p3) + L(7) / (p4 + (L(8) / p)^L(6))
    return (sig, SS, S1, B1, S2, B2, S3, B3, S4, B4)
end

"Q²_max [GeV²] of G4ChipsProtonElasticXS::GetQ2max; `pP` [GeV/c], `mt` and `mp` the target and proton rest energies [GeV]."
function g4_chips_q2max(tgZ::Int, tgN::Int, pP::Float64, mt::Float64, mp::Float64)
    pP2 = pP * pP
    if tgZ == 1 && tgN == 0
        tMid = sqrt(pP2 + mp^2) * mp - mp^2
        return tMid + tMid
    end
    dmt = mt + mt
    mds = dmt * sqrt(pP2 + mp^2) + mp^2 + mt * mt
    return dmt * dmt * pP2 / mds
end

"""
    g4_chips_state(iso::IsotopeG4Model, E::Float64, frame::NamedTuple)

Gives the Q² distribution of G4ChipsElasticModel for a proton on one isotope: the CHIPS
parameters, interpolated linearly on the 128 nodes in ln p as GetChipsCrossSection does, the
Q²_max of the CHIPS target (`TM`) and 4 p_cm² of the true target (`tmax`).

# Input Argument(s)
- `iso::IsotopeG4Model` : isotope model.
- `E::Float64` : laboratory kinetic energy [MeV].
- `frame::NamedTuple` : center-of-mass frame from `cm_frame`.

# Output Argument(s)
- `c::NamedTuple` : `tgZ`, `tgN`, true `A`, `lp` = ln p [GeV/c], CHIPS cross-section `sig`
  [mb] (only its sign is used), `SS`, `S1`...`B4`, `TM` and `tmax` [GeV²].
"""
function g4_chips_state(iso::IsotopeG4Model, E::Float64, frame::NamedTuple)
    tgZ, tgN = iso.chips_Z, iso.chips_N
    plab = momentum_lab(E, iso.m_mev; kinematics="relativistic") / 1000.0
    lp = log(plab)
    vals = if G4_CHIPS_LPMIN < lp < G4_CHIPS_LPMAX
        sh = (lp - G4_CHIPS_LPMIN) / G4_CHIPS_DLNP
        b = clamp(floor(Int, sh), 0, G4_CHIPS_NPOINTS - 2)
        f = sh - b
        lo = g4_chips_tab_values(iso.chips_par, G4_CHIPS_LPMIN + b * G4_CHIPS_DLNP, tgZ, tgN)
        hi = g4_chips_tab_values(iso.chips_par, G4_CHIPS_LPMIN + (b + 1) * G4_CHIPS_DLNP, tgZ, tgN)
        lo .+ f .* (hi .- lo)
    else
        g4_chips_tab_values(iso.chips_par, lp, tgZ, tgN)
    end
    TM = g4_chips_q2max(tgZ, tgN, plab, iso.chips_M_mev / 1000.0, iso.m_mev / 1000.0)
    return (tgZ=tgZ, tgN=tgN, A=iso.A, lp=lp, sig=max(vals[1], 0.0), SS=vals[2], S1=vals[3],
            B1=vals[4], S2=vals[5], B2=vals[6], S3=vals[7], B3=vals[8], S4=vals[9], B4=vals[10],
            TM=TM, tmax=4.0 * frame.p_star^2 / G4_GEV2)
end

"Weights I_k of the components of G4ChipsProtonElasticXS::GetExchangeT."
function g4_chips_components(c::NamedTuple)
    T = c.TM
    if c.tgZ == 1 && c.tgN == 0
        return (c.S1 / c.B1 * (1 - exp(-T * c.B1)), c.S2 * (1 - exp(-(T * c.B2)^3)), c.S3 * (1 - exp(-T * c.B3)))
    end
    heavy = c.tgZ + c.tgN > 6.5
    E2 = T^3 * c.B2 * (heavy ? T^2 : 1.0)
    E3 = T * c.B3 * (heavy ? T^6 : 1.0)
    return (c.S1 * (1 - exp(-T * (c.B1 + T * c.SS))), c.S2 * (1 - exp(-E2)), c.S3 * (1 - exp(-E3)),
            c.S4 * (1 - exp(-T * c.B4)))
end

"""
    g4_chips_density_raw(c::NamedTuple, q::Float64)

Gives the unnormalised density of Q² that G4ChipsProtonElasticXS::GetExchangeT samples, the
sum of its component densities on [0, Q²_max] (T = Q²_max):
- p + ¹H: S1 e^(-B1 Q²) + 3 S2 B2³ Q⁴ e^(-(B2 Q²)³) + S3 B3 e^(-B3 Q²);
- A < 6.5: S1 (B1 + 2 SS Q²) e^(-B1 Q² - SS Q⁴) + 3 S2 B2 Q⁴ e^(-B2 Q⁶) + S3 B3 e^(-B3 Q²)
  + S4 B4 e^(-B4 (T - Q²));
- A > 6.5: S1 (B1 + 2 SS Q²) e^(-B1 Q² - SS Q⁴) + 5 S2 B2 Q⁸ e^(-B2 Q¹⁰)
  + 7 S3 B3 Q¹² e^(-B3 Q¹⁴) + S4 B4 e^(-B4 Q²).

# Input Argument(s)
- `c::NamedTuple` : CHIPS state from `g4_chips_state`.
- `q::Float64` : Q² [GeV²].

# Output Argument(s)
- `f::Float64` : unnormalised density [1/GeV²].
"""
function g4_chips_density_raw(c::NamedTuple, q::Float64)
    (0.0 <= q <= c.TM) || return 0.0
    if c.tgZ == 1 && c.tgN == 0
        return c.S1 * exp(-c.B1 * q) + 3 * c.S2 * c.B2^3 * q^2 * exp(-(c.B2 * q)^3) + c.S3 * c.B3 * exp(-c.B3 * q)
    end
    tss = 2 * c.SS
    f1 = abs(tss) > 1.0e-7 ? c.S1 * (c.B1 + tss * q) * exp(-q * (c.B1 + c.SS * q)) : c.S1 * c.B1 * exp(-c.B1 * q)
    if c.tgZ + c.tgN > 6.5
        f2 = 5 * c.S2 * c.B2 * q^4 * exp(-c.B2 * q^5)
        f3 = 7 * c.S3 * c.B3 * q^6 * exp(-c.B3 * q^7)
        f4 = c.S4 * c.B4 * exp(-c.B4 * q)
    else
        f2 = 3 * c.S2 * c.B2 * q^2 * exp(-c.B2 * q^3)
        f3 = c.S3 * c.B3 * exp(-c.B3 * q)
        f4 = c.S4 * c.B4 * exp(-c.B4 * (c.TM - q))
    end
    return f1 + f2 + f3 + f4
end

"Coefficients (a₁, b₁, a₂, b₂) of the two-exponential form of G4HadronElastic::SampleInvariantT."
function g4_gheisha_coefficients(A::Int)
    if A <= 62
        bb = 14.5 * A^(2.0 / 3.0); dd = 20.0
        aa = A * A / bb; cc = 1.4 * A^(1.0 / 3.0) / dd
    else
        bb = 60.0 * A^(1.0 / 3.0); dd = 25.0
        aa = A^1.33 / bb; cc = 0.2 * A^0.4 / dd
    end
    return aa, bb, cc, dd
end

"Normalised density of Q² [1/GeV²] of G4HadronElastic::SampleInvariantT on [0, tmax]."
function g4_gheisha_density(A::Int, tmax::Float64, t::Float64)
    (0.0 <= t <= tmax) || return 0.0
    aa, bb, cc, dd = g4_gheisha_coefficients(A)
    s1 = aa * (1 - exp(-min(bb * tmax, 18.0)))
    s2 = cc * (1 - exp(-min(dd * tmax, 18.0)))
    f = (t <= 18.0 / bb ? aa * bb * exp(-bb * t) : 0.0) + (t <= 18.0 / dd ? cc * dd * exp(-dd * t) : 0.0)
    return f / (s1 + s2)
end

"Probability that CHIPS samples Q² above 4 p_cm² of the true target (³H, ³He only), trapezoidal rule."
function g4_chips_tail(c::NamedTuple, norm::Float64)
    n = 400
    h = (c.TM - c.tmax) / n
    acc = 0.0
    for k in 0:n
        acc += ((k == 0 || k == n) ? 0.5 : 1.0) * g4_chips_density_raw(c, c.tmax + k * h)
    end
    return acc * h / norm
end

"""
    g4_chips_t_pdf(c::NamedTuple, t::Float64)

Gives the normalised density of Q² = -t on [0, 4 p_cm²] with which G4HadronElastic and
G4ChipsElasticModel scatter a proton: the CHIPS density, uniform below ln p = -4.3 (S-wave,
T < 0.1 MeV), with the values above 4 p_cm² (³H, ³He) resampled from the two-exponential form
of G4HadronElastic, which is used throughout when the CHIPS cross-section is zero.

# Input Argument(s)
- `c::NamedTuple` : CHIPS state from `g4_chips_state`.
- `t::Float64` : Q² [GeV²].

# Output Argument(s)
- `f::Float64` : density [1/GeV²].
"""
function g4_chips_t_pdf(c::NamedTuple, t::Float64)
    (0.0 <= t <= c.tmax) || return 0.0
    c.sig > 0 || return g4_gheisha_density(c.A, c.tmax, t)
    if c.lp < -4.3
        f = 1.0 / c.TM
        pout = max(c.TM - c.tmax, 0.0) / c.TM
    else
        norm = sum(g4_chips_components(c))
        f = g4_chips_density_raw(c, t) / norm
        pout = c.TM > c.tmax ? g4_chips_tail(c, norm) : 0.0
    end
    return (t <= c.TM ? f : 0.0) + (pout > 0 ? pout * g4_gheisha_density(c.A, c.tmax, t) : 0.0)
end

#----
# Whole cross-section in the center-of-mass frame
#----

"Geant4 Coulomb cross-section dσ/dΩ_cm [b/sr] at x = 1 - μ_cm (projectile, not folded)."
g4_coulomb_dcs_cm(w::NamedTuple, frame::NamedTuple, x::Float64) =
    g4_coulomb_dcs_lab(w, one_minus_mu_lab(frame, x)) * lab_solid_angle_ratio(frame, x)

"Geant4 hadronic cross-section dσ/dΩ_cm [b/sr] at x = 1 - μ_cm, σ f(Q²) dQ²/dΩ, Q² = tmax x/2."
g4_nuclear_dcs_cm(c::NamedTuple, σ_nuc::Float64, x::Float64) =
    σ_nuc * g4_chips_t_pdf(c, 0.5 * c.tmax * x) * c.tmax / (4.0 * π)

"""
    g4_dcs_function(iso::IsotopeG4Model, E::Float64, screening::Bool)

Gives the Geant4 elastic cross-section of one isotope at one energy as a function of
x = 1 - μ_cm: the Coulomb and hadronic cross-sections summed, folded onto μ_cm ≥ 0 for p + ¹H.

# Input Argument(s)
- `iso::IsotopeG4Model` : isotope model.
- `E::Float64` : laboratory kinetic energy [MeV].
- `screening::Bool` : whether the Coulomb field is screened.

# Output Argument(s)
- `f::Function` : x -> dσ/dΩ [b/sr].
- `frame::NamedTuple` : center-of-mass frame from `cm_frame`.
- `w::NamedTuple` : Coulomb parameters from `g4_wentzel`.
"""
function g4_dcs_function(iso::IsotopeG4Model, E::Float64, screening::Bool)
    frame = cm_frame(E, iso.m_mev, iso.M_mev)
    w = g4_wentzel(iso, E, screening)
    c = g4_chips_state(iso, E, frame)
    σ_nuc = g4_bgg_elastic(iso, E)
    g(x) = g4_coulomb_dcs_cm(w, frame, x) + g4_nuclear_dcs_cm(c, σ_nuc, x)
    iso.identical && return (x -> x <= 1.0 ? g(x) + g(2.0 - x) : 0.0), frame, w
    return g, frame, w
end

#----
# Database and tabulated cross-section
#----

"""
    init_elastic_scattering_g4_db(isotopes::Vector{Tuple{Int,Int}}, particle::Particle;
    screening::Bool=true)

Builds the Geant4 description of elastic scattering of protons on the given isotopes.

# Input Argument(s)
- `isotopes::Vector{Tuple{Int,Int}}` : target isotope pairs `(Z, A)`.
- `particle::Particle` : incident particle (proton).
- `screening::Bool` : whether the Coulomb field is screened by the atomic electrons.

# Output Argument(s)
- `db::ElasticScatteringG4DB` : Geant4 description.
"""
function init_elastic_scattering_g4_db(isotopes::Vector{Tuple{Int,Int}}, particle::Particle;
    screening::Bool=true)
    is_proton(particle) || error("The Geant4 elastic scattering is defined for protons only.")
    m = get_mass(particle)
    me = get_mass(Electron())
    db = Dict{Tuple{Int,Int},IsotopeG4Model}()
    for (Z, A) in isotopes
        1 <= Z <= 99 || error("The Geant4 elastic scattering is defined for Z = 1-99 (Z = $(Z)).")
        tgZ, tgN = g4_chips_target(Z, A)
        M_chips = (tgZ == 1 && tgN == 0) ? m : get_mass(tgZ, tgZ + tgN)
        db[(Z, A)] = IsotopeG4Model(Z, A, m, get_mass(Z, A), Z == 1 && A == 1, g4_screen_rsquare(Z, me),
                                    g4_formfactor(Z), tgZ, tgN, M_chips, g4_chips_parameters(tgZ, tgN))
    end
    return ElasticScatteringG4DB(screening, db)
end

"""
    prepare_g4_dsig_cache(db::ElasticScatteringG4DB, Z::Int, A::Int, E_in_MeV::Float64)

Tabulates the Geant4 elastic cross-section of one isotope at one incident energy in
x = 1 - μ_cm, as `prepare_op_dsig_cache` does: start grid (`dcs_table_start_grid`) from
x₂ = 10⁻⁴ A_s (p_lab/p_cm)², the Coulomb screening parameter in the center-of-mass frame
(10⁻¹² without screening), with the diffraction nodes of R = 1.2 (A^(1/3) + 1) fm, then
power-law refinement (`tabulate_power_law_dcs`).

# Input Argument(s)
- `db::ElasticScatteringG4DB` : Geant4 description.
- `Z::Int` : target atomic number.
- `A::Int` : target mass number.
- `E_in_MeV::Float64` : incident kinetic energy [MeV].

# Output Argument(s)
- `cache::TabulatedElasticDCSCache` : tabulated cross-section.
"""
function prepare_g4_dsig_cache(db::ElasticScatteringG4DB, Z::Int, A::Int, E_in_MeV::Float64)
    iso = db.isotopes[(Z, A)]
    f, frame, w = g4_dcs_function(iso, E_in_MeV, db.screening)
    x_max = iso.identical ? 1.0 : 2.0
    x_first = db.screening ? 1.0e-4 * w.screenZ * w.mom2 / frame.p_star^2 : 1.0e-12
    R = 1.2 * (A^(1.0 / 3.0) + (iso.m_mev / AMU_MEV)^(1.0 / 3.0))
    nodes = dcs_table_start_grid(x_first, x_max, frame.p_star / OP_HBARC * R)
    X, Σ, p = tabulate_power_law_dcs(f, nodes, db.screening;
        label="Geant4 elastic cross-section of Z=$(Z), A=$(A) at $(E_in_MeV) MeV")
    return TabulatedElasticDCSCache(E_in_MeV, iso.m_mev, iso.M_mev, iso.identical, X, Σ, p)
end
