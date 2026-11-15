# Elastic scattering of protons and alphas by neutral atoms described by optical-model
# potentials (model "OP" of Elastic_Scattering), after
#   F. Salvat and J. M. Quesada, Nucl. Instrum. Methods B 475 (2020) 49, and
#   F. Salvat and L. Heredia, Nucl. Instrum. Methods B 546 (2024) 165157:
#     dσ/dΩ = dσ_scr/dΩ · F_nuc(θ),
# with dσ_scr the eikonal cross-section for the screened Coulomb field of a point nucleus
# (DHFS screening) and F_nuc = dσ_nuc/dσ_R the nuclear factor from partial waves in the
# optical-model potential (Koning-Delaroche, Watson et al., Su-Han). Protons on A ≤ 6 use the
# empirical parameterisation of Galyuzov and Kosov instead. Everything is in the center-of-mass
# frame with relativistic kinematics.
#
# For each isotope and incident energy, the cross-section is tabulated in x = 1 - μ_cm on
# nodes chosen adaptively (`prepare_op_dsig_cache`). It is interpolated as a power law between
# nodes, so that its integrals weighted by polynomials in x are closed forms; the moments that
# involve the laboratory cosine are integrated in ln x. The integration kernels of the ENDF
# description have methods here for `OpticalDCSCache`.

#----
# Constants
#----

const OP_HBARC = HBAR_eVs * C_m_per_s * 1.0e9      # ħc [MeV fm]
const OP_E2 = ALPHA_FS * OP_HBARC                  # e² [MeV fm]
const OP_A0_FM = A0_M * 1.0e15                     # Bohr radius [fm]
const OP_LAMBDA_PI2 = 1.429502^2                   # (ħ/m_π c)² [fm²], Koning-Delaroche convention
const OP_FM2_TO_BARN = 0.01
const OP_DATA_FILE = "elastic_scattering_op_salvat_2024.jld2"

# Upper energy of the Watson et al. potential [MeV]. Salvat & Heredia (2024) switch to
# Koning-Delaroche at 35 MeV; the Watson potential is used here up to the upper end of its
# fitted range, 50 MeV, where it agrees better with the experimental data of their Fig. 2.
const OP_WATSON_EMAX = 50.0
const OP_KD_EFREEZE = 300.0                        # Koning-Delaroche frozen above [MeV]
const OP_SU_HAN_EFREEZE = 386.0                    # Su-Han frozen above [MeV]

# Tabulation of the cross-section: relative tolerance of the power-law interpolation at the
# log-midpoint of every piece, and maximum number of nodes.
const OP_TABLE_RTOL = 1.0e-3
const OP_TABLE_MAX_NODES = 4000

# Quadrature of the moments involving the laboratory cosine on each piece, in ln x:
# :adaptive (Gauss-Kronrod) or :gauss (fixed Gauss-Legendre).
const OP_QUADRATURE = Ref(:adaptive)

#----
# Point-Coulomb cross-section (kinematics: `cm_kinematics`, scattering_kinematics.jl)
#----

"Point-Coulomb (Rutherford) cross-section [fm²/sr] at x = 1 - cos θ."
op_rutherford(kin::CMKinematics, x::Float64) = (kin.eta / (kin.k * x))^2

#----
# Electronic screening: eikonal cross-section for the screened field of a point nucleus
#----

"""
    op_dhfs_screening(Z::Int)

Gives the coefficients Aᵢ and the screening constants aᵢ [1/fm] of the analytical DHFS
screening function Φ(r) = Σ Aᵢ exp(-aᵢ r), read from `data/elastic_scattering_op_salvat_2024.jld2`
(key "dhfs": Z => [A1, A2, α1, α2, α3] with α in units of 1/a₀, A3 = 1 - A1 - A2, and α3 = NaN
for the elements with two terms).

# Input Argument(s)
- `Z::Int` : atomic number (1 to 92).

# Output Argument(s)
- `A::Vector{Float64}` : coefficients.
- `a::Vector{Float64}` : screening constants [1/fm].

# Reference(s)
- Salvat, Martínez, Mayol and Parellada, Phys. Rev. A 36 (1987) 467, Table I.
"""
function op_dhfs_screening(Z::Int)
    data = fast_load(OP_DATA_FILE)["dhfs"]
    haskey(data, Z) || error("DHFS screening parameters are tabulated for Z = 1-92 only (Z = $(Z)).")
    A1, A2, α1, α2, α3 = data[Z]
    A, a = isnan(α3) ? ([A1, A2], [α1, α2]) : ([A1, A2, 1.0 - A1 - A2], [α1, α2, α3])
    return A, a ./ OP_A0_FM
end

"""
    op_charge_radius(Z::Int, A::Int)

Gives the rms nuclear charge radius [fm] of a light isotope, read from
`data/elastic_scattering_op_salvat_2024.jld2` (key "charge_radius": Z => A => r).

# Input Argument(s)
- `Z::Int` : atomic number.
- `A::Int` : mass number.

# Output Argument(s)
- `r::Float64` : rms charge radius [fm].

# Reference(s)
- Angeli and Marinova, At. Data Nucl. Data Tables 99 (2013) 69.
"""
function op_charge_radius(Z::Int, A::Int)
    data = fast_load(OP_DATA_FILE)["charge_radius"]
    (haskey(data, Z) && haskey(data[Z], A)) || error("No charge radius for Z = $(Z), A = $(A).")
    return data[Z][A]
end

"""
    OPAtomGrid

Impact-parameter quadrature of the eikonal amplitude for one atom, independent of the
energy: nodes and weights, and the screening sums at the nodes, with which the eikonal phase
of the Wallace approximation is χ(b) = -2η [Φ₀(b) - (η/k) Ψ₀(b)] and its derivative
χ'(b) = 2η [Φ₁(b) - (η/k) Ψ₁(b)]:
Φ₀ = Σᵢ Aᵢ K₀(aᵢb), Ψ₀ = Σᵢⱼ AᵢAⱼaⱼ K₀((aᵢ+aⱼ)b), Φ₁ = Σᵢ Aᵢaᵢ K₁(aᵢb),
Ψ₁ = Σᵢⱼ AᵢAⱼaⱼ(aᵢ+aⱼ) K₁((aᵢ+aⱼ)b).

# Field(s)
- `q::Float64` : largest momentum transfer the quadrature is built for [1/fm].
- `b::Vector{Float64}` : nodes [fm].
- `w::Vector{Float64}` : weights [fm].
- `Φ0`, `Ψ0`, `Φ1`, `Ψ1` : screening sums at the nodes.
"""
struct OPAtomGrid
    q::Float64
    b::Vector{Float64}
    w::Vector{Float64}
    Φ0::Vector{Float64}
    Ψ0::Vector{Float64}
    Φ1::Vector{Float64}
    Ψ1::Vector{Float64}
end

# Grids by (Z, η bucket, level), built once per session
const OP_ATOM_GRIDS = Dict{Tuple{Int,Int,Int},OPAtomGrid}()

"""
    op_atom_grid(Z::Int, η_bucket::Int, level::Int)

Gives the impact-parameter quadrature of the eikonal amplitude for the atom of atomic number
`Z`, built once and cached. The b range extends to 30 decay lengths of the slowest screening
exponential, in panels of 20 Gauss-Legendre points that are geometric near the origin and
limited to a change of the eikonal phase of 4π (for Sommerfeld parameters up to
η = 6·2^η_bucket) and to 4 oscillations of J₁(q b) up to q = q_ref/4^level, where
q_ref = 200/R_at bounds the momentum transfer at θ_eik for every energy.

# Input Argument(s)
- `Z::Int` : atomic number.
- `η_bucket::Int` : bucket of the Sommerfeld parameter.
- `level::Int` : level of the largest momentum transfer.

# Output Argument(s)
- `grid::OPAtomGrid` : quadrature.
"""
function op_atom_grid(Z::Int, η_bucket::Int, level::Int)
    key = (Z, η_bucket, level)
    haskey(OP_ATOM_GRIDS, key) && return OP_ATOM_GRIDS[key]
    A, a = op_dhfs_screening(Z)
    n = length(A)
    Φ1f(b) = sum(A[i] * a[i] * besselk(1, a[i] * b) for i in 1:n)
    q_ref = 200.0 / (Z^(-1.0 / 3.0) * OP_A0_FM)
    q = q_ref / 4.0^level
    η_max = 6.0 * 2.0^η_bucket
    a_min = minimum(a[abs.(A) .> 0.0])
    b_max = 30.0 / a_min
    edges = [0.0, 1.0e-4 / q_ref]
    b = edges[end]
    while b < b_max
        dχ = 2.4 * η_max * abs(Φ1f(b))
        b += min(b, 4.0π / max(dχ, 1.0e-300), 8.0π / q, 0.5 / a_min)
        push!(edges, b)
    end
    u, wq = quadrature(20, "gauss-legendre")
    N = 20 * (length(edges) - 1)
    nodes = zeros(N); weights = zeros(N)
    Φ0 = zeros(N); Ψ0 = zeros(N); Φ1 = zeros(N); Ψ1 = zeros(N)
    K0 = zeros(n, n); K1 = zeros(n, n)
    m = 0
    for p in 1:(length(edges) - 1)
        c, h = 0.5 * (edges[p + 1] + edges[p]), 0.5 * (edges[p + 1] - edges[p])
        for j in eachindex(u)
            m += 1
            bb = c + h * u[j]
            nodes[m] = bb; weights[m] = h * wq[j]
            for i in 1:n, l in i:n
                K0[i, l] = K0[l, i] = besselk(0, (a[i] + a[l]) * bb)
                K1[i, l] = K1[l, i] = besselk(1, (a[i] + a[l]) * bb)
            end
            for i in 1:n
                Φ0[m] += A[i] * besselk(0, a[i] * bb)
                Φ1[m] += A[i] * a[i] * besselk(1, a[i] * bb)
                for l in 1:n
                    Ψ0[m] += A[i] * A[l] * a[l] * K0[i, l]
                    Ψ1[m] += A[i] * A[l] * a[l] * (a[i] + a[l]) * K1[i, l]
                end
            end
        end
    end
    grid = OPAtomGrid(q, nodes, weights, Φ0, Ψ0, Φ1, Ψ1)
    OP_ATOM_GRIDS[key] = grid
    return grid
end

"""
    OPScreening

Eikonal description of the scattering by the screened Coulomb field of a point nucleus at
one energy. For each level of the impact-parameter quadrature (`op_atom_grid`), the weights
times χ'(b) e^{iχ(b)} b, which do not depend on the angle, are computed once; an amplitude at
momentum transfer q then costs one sum of J₁(q b) on the coarsest level built for q. Beyond
θ_eik, the cross-section is continued with the form of Salvat and Quesada (2020), Eq. (45),
matched to the eikonal one at θ_eik. Without screening, the cross-section is point Coulomb.

# Field(s)
- `screened::Bool` : whether the field is screened.
- `kin::CMKinematics` : center-of-mass quantities.
- `x_eik::Float64` : 1 - cos θ_eik.
- `q::Vector{Float64}` : largest momentum transfer of each level [1/fm], decreasing.
- `b::Vector{Vector{Float64}}` : impact-parameter nodes of each level [fm].
- `g::Vector{Vector{ComplexF64}}` : weights w χ'(b) e^{iχ(b)} b of each level.
- `f0::ComplexF64` : forward amplitude [fm].
- `ABC::NTuple{3,Float64}` : coefficients of the large-angle form.
"""
struct OPScreening
    screened::Bool
    kin::CMKinematics
    x_eik::Float64
    q::Vector{Float64}
    b::Vector{Vector{Float64}}
    g::Vector{Vector{ComplexF64}}
    f0::ComplexF64
    ABC::NTuple{3,Float64}
end

"Eikonal amplitude f(q) [fm] (Salvat and Quesada 2020, Eqs. (39), (42))."
function op_eikonal_amplitude(S::OPScreening, q::Float64)
    q == 0.0 && return S.f0
    lev = findlast(>=(q), S.q)
    isnothing(lev) && (lev = 1)
    b, g = S.b[lev], S.g[lev]
    acc = zero(ComplexF64)
    @inbounds for n in eachindex(b)
        acc += besselj1(q * b[n]) * g[n]
    end
    return -(S.kin.k / q) * acc
end

"""
    op_screening(Z::Int, kin::CMKinematics, screening::Bool)

Builds the eikonal description of the scattering by the screened Coulomb field of a point
nucleus of atomic number `Z`, from the cached impact-parameter quadratures of the atom.

# Input Argument(s)
- `Z::Int` : atomic number of the target.
- `kin::CMKinematics` : center-of-mass quantities.
- `screening::Bool` : whether the field is screened; otherwise point Coulomb.

# Output Argument(s)
- `S::OPScreening` : eikonal description.

# Reference(s)
- Salvat and Quesada (2020), Eqs. (36)-(46); Salvat and Heredia (2024), Eqs. (27)-(36).
"""
function op_screening(Z::Int, kin::CMKinematics, screening::Bool)
    screening || return OPScreening(false, kin, 0.0, Float64[], Vector{Float64}[], Vector{ComplexF64}[],
                                    0.0im, (0.0, 0.0, 0.0))
    A, a = op_dhfs_screening(Z)
    η, k = kin.eta, kin.k
    θ_eik = min(200.0 / (k * Z^(-1.0 / 3.0) * OP_A0_FM), 0.1 * π)
    q_max = 2.0 * k * sin(θ_eik / 2.0)
    bucket = max(0, ceil(Int, log2(abs(η) / 6.0)))
    a_min = minimum(a[abs.(A) .> 0.0])

    qs = Float64[]; bs = Vector{Float64}[]; gs = Vector{ComplexF64}[]
    f0 = 0.0im
    level = 0
    while true
        grid = op_atom_grid(Z, bucket, level)
        χ = -2.0 .* η .* (grid.Φ0 .- (η / k) .* grid.Ψ0)
        dχ = 2.0 .* η .* (grid.Φ1 .- (η / k) .* grid.Ψ1)
        push!(qs, grid.q); push!(bs, grid.b)
        push!(gs, grid.w .* dχ .* cis.(χ) .* grid.b)
        # Coarser levels until the oscillation limit no longer sets the panels
        if 8.0π / grid.q >= 0.5 / a_min || grid.q < 1.0e-3 * q_max
            f0 = -im * k * sum(grid.w .* (cis.(χ) .- 1.0) .* grid.b)
            break
        end
        level += 1
    end
    S = OPScreening(true, kin, 2.0 * sin(θ_eik / 2.0)^2, qs, bs, gs, f0, (0.0, 0.0, 0.0))

    # Large-angle form matched to the eikonal cross-section at θ_eik
    qm = [q_max * 0.90, q_max * 0.95, q_max]
    ys = [2.0 * η * k / abs(op_eikonal_amplitude(S, q)) for q in qm]
    Mq = vcat([[1.0 q^(2.0 / 3.0) q^(4.0 / 3.0)] for q in qm]...)
    ABC = Mq \ (ys .- qm .^ 2)
    return OPScreening(true, kin, S.x_eik, qs, bs, gs, f0, (ABC[1], ABC[2], ABC[3]))
end

"""
    op_screened_dcs(S::OPScreening, x::Float64)

Gives the cross-section [fm²/sr] for the screened Coulomb field of a point nucleus (spinless)
at x = 1 - cos θ.

# Input Argument(s)
- `S::OPScreening` : eikonal description.
- `x::Float64` : 1 - cos θ (center of mass).

# Output Argument(s)
- `σ::Float64` : cross-section [fm²/sr].
"""
function op_screened_dcs(S::OPScreening, x::Float64)
    S.screened || return op_rutherford(S.kin, x)
    q = S.kin.k * sqrt(2.0 * x)
    x <= S.x_eik && return abs2(op_eikonal_amplitude(S, q))
    A, B, C = S.ABC
    return (2.0 * S.kin.eta * S.kin.k)^2 / (A + B * q^(2.0 / 3.0) + C * q^(4.0 / 3.0) + q^2)^2
end

#----
# Optical-model potentials
#----

"""
    OpticalPotential

Optical-model potential for one target isotope at one energy (Koning-Delaroche sign
convention, depths positive [MeV], absolute radii [fm]):
U(r) = -Vv f(r;Rv,av) - i Wv f(r;Rw,aw) + 4 avd Vd f'(r;Rvd,avd) + 4i ad Wd f'(r;Rd,ad)
       + (ħ/m_π c)² [Vso f'(r;Rso,aso) + i Wso f'(r;Rwso,awso)] (l·σ)/r_so + V_C(r),
with f the Woods-Saxon form factor, r_so = r (Thomas form) or R_so (form of Watson et al.),
and V_C the potential of a uniformly charged sphere of radius Rc.
"""
Base.@kwdef struct OpticalPotential
    Vv::Float64 = 0.0;  Rv::Float64 = 1.0;  av::Float64 = 0.5
    Wv::Float64 = 0.0;  Rw::Float64 = 1.0;  aw::Float64 = 0.5
    Vd::Float64 = 0.0;  Rvd::Float64 = 1.0; avd::Float64 = 0.5
    Wd::Float64 = 0.0;  Rd::Float64 = 1.0;  ad::Float64 = 0.5
    Vso::Float64 = 0.0; Rso::Float64 = 1.0; aso::Float64 = 0.5
    Wso::Float64 = 0.0; Rwso::Float64 = 1.0; awso::Float64 = 0.5
    Rc::Float64 = 1.0
    so_watson::Bool = false
end

"Range of the short-range part of the potential [fm]: largest radius plus 16 diffusenesses."
function op_nuclear_range(U::OpticalPotential)
    r = U.Rc
    for (R, a, V) in ((U.Rv, U.av, U.Vv), (U.Rw, U.aw, U.Wv), (U.Rvd, U.avd, U.Vd),
                      (U.Rd, U.ad, U.Wd), (U.Rso, U.aso, U.Vso), (U.Rwso, U.awso, U.Wso))
        V != 0.0 && (r = max(r, R + 16.0 * a))
    end
    return r
end

"""
    op_koning_delaroche_proton(Z::Int, A::Int, E::Float64)

Gives the global proton optical-model potential of Koning and Delaroche, with its energy
dependence frozen above 300 MeV as in Salvat and Heredia (2024).

# Input Argument(s)
- `Z::Int` : atomic number of the target.
- `A::Int` : mass number of the target.
- `E::Float64` : laboratory kinetic energy [MeV].

# Output Argument(s)
- `U::OpticalPotential` : potential.

# Reference(s)
- Koning and Delaroche, Nucl. Phys. A 713 (2003) 231, Tables 10-11.
"""
function op_koning_delaroche_proton(Z::Int, A::Int, E::Float64)
    N = A - Z
    A13 = A^(1.0 / 3.0)
    f = min(E, OP_KD_EFREEZE) - (-8.4075 + 0.01378 * A)
    v1 = 59.30 + 21.0 * (N - Z) / A - 0.024 * A
    v2 = 0.007067 + 4.23e-6 * A
    v3 = 1.729e-5 + 1.136e-8 * A
    v4 = 7.0e-9
    w1 = 14.667 + 0.009629 * A
    w2 = 73.55 + 0.0795 * A
    d1 = 16.0 + 16.0 * (N - Z) / A
    d2 = 0.0180 + 0.003802 / (1.0 + exp((A - 156.0) / 8.0))
    d3 = 11.5
    rv = 1.3039 - 0.4054 / A13;  av = 0.6778 - 1.487e-4 * A
    rd = 1.3424 - 0.01585 * A13; ad = 0.5187 + 5.205e-4 * A
    rso = 1.1854 - 0.647 / A13;  aso = 0.59
    rc = 1.198 + 0.697 * A^(-2.0 / 3.0) + 12.994 * A^(-5.0 / 3.0)
    Vc = 1.73 / rc * Z / A13
    Vv = v1 * (1.0 - v2 * f + v3 * f^2 - v4 * f^3) + Vc * v1 * (v2 - 2.0 * v3 * f + 3.0 * v4 * f^2)
    return OpticalPotential(
        Vv=Vv, Rv=rv * A13, av=av, Wv=w1 * f^2 / (f^2 + w2^2), Rw=rv * A13, aw=av,
        Wd=d1 * f^2 * exp(-d2 * f) / (f^2 + d3^2), Rd=rd * A13, ad=ad,
        Vso=(5.922 + 0.0030 * A) * exp(-0.0040 * f), Rso=rso * A13, aso=aso,
        Wso=-3.1 * f^2 / (f^2 + 160.0^2), Rwso=rso * A13, awso=aso,
        Rc=rc * A13)
end

"""
    op_watson_proton(Z::Int, A::Int, E::Float64, m::Float64, M::Float64)

Gives the proton optical-model potential of Watson, Singh and Segel for 1p-shell nuclei. The
energies of the prescription are center-of-mass energies; the spin-orbit term is their
modification of the Thomas form, (ħ/m_π c)² R_so⁻¹ df/dr.

# Input Argument(s)
- `Z::Int` : atomic number of the target.
- `A::Int` : mass number of the target.
- `E::Float64` : laboratory kinetic energy [MeV].
- `m::Float64` : proton rest energy [MeV].
- `M::Float64` : target rest energy [MeV].

# Output Argument(s)
- `U::OpticalPotential` : potential.

# Reference(s)
- Watson, Singh and Segel, Phys. Rev. 182 (1969) 977, Sec. IV and Table II.
"""
function op_watson_proton(Z::Int, A::Int, E::Float64, m::Float64, M::Float64)
    A13 = A^(1.0 / 3.0)
    asym = (A - 2Z) / A
    Ecm = E * M / (m + M)
    R = (1.15 - 0.001 * Ecm) * A13
    Ws = (Ecm < 13.8 ? 0.64 * Ecm : 9.60 - 0.06 * Ecm) + 10.0 * asym
    Wv = Ecm < 32.7 ? 0.0 : (Ecm <= 39.3 ? 1.15 * (Ecm - 32.7) : 7.5)
    return OpticalPotential(
        Vv=60.0 + 0.4 * Z / A13 + 27.0 * asym - 0.3 * Ecm, Rv=R, av=0.57, Wv=Wv, Rw=R, aw=0.50,
        Wd=max(Ws, 0.0), Rd=R, ad=0.50, Vso=5.5, Rso=R, aso=0.57, Rc=R, so_watson=true)
end

"""
    op_su_han_alpha(Z::Int, A::Int, E::Float64)

Gives the global alpha optical-model potential of Su and Han, with its energy dependence
frozen above 386 MeV as in Salvat and Heredia (2024).

# Input Argument(s)
- `Z::Int` : atomic number of the target.
- `A::Int` : mass number of the target.
- `E::Float64` : laboratory kinetic energy [MeV].

# Output Argument(s)
- `U::OpticalPotential` : potential.

# Reference(s)
- Su and Han, Int. J. Mod. Phys. E 24 (2015) 1550092, Table 1.
"""
function op_su_han_alpha(Z::Int, A::Int, E::Float64)
    A13 = A^(1.0 / 3.0)
    asym = (A - 2Z) / A
    Ee = min(E, OP_SU_HAN_EFREEZE)
    return OpticalPotential(
        Vv=175.0881 - 0.6236 * Ee + 0.0006 * Ee^2 + 30.0 * asym - 0.2360 * Z / A13,
        Rv=1.3421 * A13, av=0.6578,
        Wv=max(0.0, -4.0174 + 0.1409 * Ee), Rw=1.4259 * A13, aw=0.5587,
        Wd=max(0.0, 27.5816 - 0.0797 * Ee + 48.0 * asym), Rd=1.2928 * A13, ad=0.6359,
        Rc=1.35 * A13)
end

"""
    op_potential(iso::IsotopeOpticalModel, E::Float64)

Gives the optical-model potential selected by Salvat and Heredia (2024): Watson et al. for
protons on 6 < A < 24 up to `OP_WATSON_EMAX`, Koning-Delaroche for the other protons, Su-Han
for alphas.

# Input Argument(s)
- `iso::IsotopeOpticalModel` : isotope model.
- `E::Float64` : laboratory kinetic energy [MeV].

# Output Argument(s)
- `U::OpticalPotential` : potential.
"""
function op_potential(iso::IsotopeOpticalModel, E::Float64)
    iso.Z1 == 2 && return op_su_han_alpha(iso.Z2, iso.A, E)
    (iso.A < 24 && E <= OP_WATSON_EMAX) && return op_watson_proton(iso.Z2, iso.A, E, iso.m_mev, iso.M_mev)
    return op_koning_delaroche_proton(iso.Z2, iso.A, E)
end

#----
# Partial waves (nuclear S-matrix) and nuclear factor
#----

"""
    OPPartialWaves

Nuclear S-matrix elements relative to point-Coulomb scattering for one isotope and energy.
`Splus[l+1]` belongs to j = l + 1/2 (or to l for spin 0), `Sminus[l+1]` to j = l - 1/2.
`coulomb_only` means that the classical turning point lies so far outside the nucleus that
the S-matrix is one.
"""
struct OPPartialWaves
    kin::CMKinematics
    spin_half::Bool
    σ::Vector{Float64}
    Splus::Vector{ComplexF64}
    Sminus::Vector{ComplexF64}
    coulomb_only::Bool
end

op_ws(r, R, a) = 1.0 / (1.0 + exp((r - R) / a))
op_dws(r, R, a) = (f = op_ws(r, R, a); -f * (1.0 - f) / a)

function op_central(U::OpticalPotential, Zc::Float64, r::Float64)
    V = -U.Vv * op_ws(r, U.Rv, U.av) + 4.0 * U.avd * U.Vd * op_dws(r, U.Rvd, U.avd)
    W = -U.Wv * op_ws(r, U.Rw, U.aw) + 4.0 * U.ad * U.Wd * op_dws(r, U.Rd, U.ad)
    Vc = r < U.Rc ? Zc / (2.0 * U.Rc) * (3.0 - (r / U.Rc)^2) : Zc / r
    return complex(V + Vc, W)
end

function op_spin_orbit(U::OpticalPotential, r::Float64)
    rso = U.so_watson ? U.Rso : r
    return OP_LAMBDA_PI2 * complex(U.Vso * op_dws(r, U.Rso, U.aso), U.Wso * op_dws(r, U.Rwso, U.awso)) / rso
end

"Coulomb phase shifts σ_l = arg Γ(l + 1 + iη), l = 0..Lmax."
function op_coulomb_phases(η::Float64, Lmax::Int)
    σ = zeros(Lmax + 1)
    σ[1] = imag(loggamma(complex(1.0, η)))
    for l in 1:Lmax
        σ[l + 1] = σ[l] + atan(η / l)
    end
    return σ
end

op_cf_S(L, η, x) = L / x + η / L
op_cf_R2(L, η) = 1.0 + (η / L)^2

# Steed's continued fraction CF1: F'_L / F_L (modified Lentz)
function op_coulomb_cf1(η::Float64, x::Float64, L::Int)
    tiny = 1.0e-300
    f = op_cf_S(L + 1, η, x)
    f == 0.0 && (f = tiny)
    C = f; D = 0.0
    for n in 1:1_000_000
        an = -op_cf_R2(L + n, η)
        bn = op_cf_S(L + n, η, x) + op_cf_S(L + n + 1, η, x)
        D = bn + an * D; D == 0.0 && (D = tiny)
        C = bn + an / C; C == 0.0 && (C = tiny)
        D = 1.0 / D
        Δ = C * D
        f *= Δ
        abs(Δ - 1.0) < 1.0e-16 && return f
    end
    error("Coulomb CF1 did not converge (η = $(η), x = $(x), L = $(L)).")
end

# Steed's continued fraction CF2: p + iq = H'_0 / H_0
function op_coulomb_cf2(η::Float64, x::Float64)
    xi = 1.0 / x
    p = 0.0; q = 1.0 - η * xi
    ar = -η^2; ai = η
    br = 2.0 * (x - η); bi = 2.0
    wi = 2.0 * η
    den = br^2 + bi^2
    dr = br / den; di = -bi / den
    dp = -xi * (ar * di + ai * dr)
    dq = xi * (ar * dr - ai * di)
    pk = 0.0
    for _ in 1:1_000_000
        p += dp; q += dq
        pk += 2.0
        ar += pk; ai += wi
        bi += 2.0
        d = ar * dr - ai * di + br
        di = ai * dr + ar * di + bi
        c = 1.0 / (d^2 + di^2)
        dr = c * d; di = -c * di
        a = br * dr - bi * di - 1.0
        b = bi * dr + br * di
        c = dp * a - dq * b
        dq = dp * b + dq * a
        dp = c
        abs(dp) + abs(dq) < (abs(p) + abs(q)) * 1.0e-16 && return p, q
    end
    error("Coulomb CF2 did not converge (η = $(η), x = $(x)).")
end

"""
    op_coulomb_FG(η::Float64, x::Float64, Lmax::Int)

Gives the regular and irregular Coulomb functions F_L(η, x) and G_L(η, x), L = 0..Lmax, by
Steed's method. x must lie outside the classical turning point of L = 0.

# Input Argument(s)
- `η::Float64` : Sommerfeld parameter.
- `x::Float64` : k r.
- `Lmax::Int` : largest angular momentum.

# Output Argument(s)
- `F::Vector{Float64}` : F_L, index L + 1.
- `G::Vector{Float64}` : G_L, index L + 1.

# Reference(s)
- Barnett et al., Comput. Phys. Commun. 8 (1974) 377; Thompson and Barnett, J. Comput.
  Phys. 64 (1986) 490.
"""
function op_coulomb_FG(η::Float64, x::Float64, Lmax::Int)
    F = zeros(Lmax + 1); Fp = zeros(Lmax + 1)
    G = zeros(Lmax + 1); Gp = zeros(Lmax + 1)
    F[end] = 1.0e-200; Fp[end] = op_coulomb_cf1(η, x, Lmax) * F[end]
    for L in Lmax:-1:1
        S, R = op_cf_S(L, η, x), sqrt(op_cf_R2(L, η))
        F[L] = (S * F[L + 1] + Fp[L + 1]) / R
        Fp[L] = S * F[L] - R * F[L + 1]
        if abs(F[L]) > 1.0e200
            F[L:end] ./= 1.0e200; Fp[L:end] ./= 1.0e200
        end
    end
    p, q = op_coulomb_cf2(η, x)
    f0 = Fp[1] / F[1]
    scale = sign(F[1]) * sqrt(q / ((f0 - p)^2 + q^2)) / F[1]
    F .*= scale; Fp .*= scale
    G[1] = (f0 - p) * F[1] / q
    Gp[1] = p * G[1] - q * F[1]
    for L in 0:(Lmax - 1)
        S, R = op_cf_S(L + 1, η, x), sqrt(op_cf_R2(L + 1, η))
        G[L + 2] = (S * G[L + 1] - Gp[L + 1]) / R
        Gp[L + 2] = R * G[L + 1] - S * G[L + 2]
    end
    return F, G
end

# u(r_a)/u(r_b) of the regular radial solution of one (l, l·σ) channel, Numerov's method
function op_radial_ratio(l::Int, ls::Float64, Uc::Vector{ComplexF64}, Uso::Vector{ComplexF64},
                         κ::Float64, k::Float64, h::Float64, na::Int, nb::Int)
    ll = l * (l + 1)
    Q(n) = ll / (n * h)^2 + κ * (Uc[n] + ls * Uso[n]) - k^2
    n0 = max(1, ceil(Int, 0.5 * (l + 1)))
    r0, r1 = n0 * h, (n0 + 1) * h
    q0 = κ * (Uc[1] + ls * Uso[1]) - k^2
    um = complex(1.0)
    u = complex((r1 / r0)^(l + 1)) * (1.0 + q0 * r1^2 / (2 * (2l + 3))) / (1.0 + q0 * r0^2 / (2 * (2l + 3)))
    logscale = 0.0
    h12 = h^2 / 12.0
    Qm = Q(n0); Qc = Q(n0 + 1)
    ua = zero(ComplexF64); la = 0.0
    n0 + 1 == na && (ua = u; la = logscale)
    for n in (n0 + 1):(nb - 1)
        Qp = Q(n + 1)
        up = (2.0 * (1.0 + 5.0 * h12 * Qc) * u - (1.0 - h12 * Qm) * um) / (1.0 - h12 * Qp)
        um, u = u, up
        Qm, Qc = Qc, Qp
        if abs(u) > 1.0e100
            u /= 1.0e100; um /= 1.0e100; logscale += 100.0
        end
        n + 1 == na && (ua = u; la = logscale)
    end
    return ua / u * 10.0^(la - logscale)
end

"""
    op_partial_waves(iso::IsotopeOpticalModel, E::Float64, U::OpticalPotential)

Solves the radial equations of the projectile in the optical-model potential and gives the
nuclear S-matrix relative to point-Coulomb scattering. The radial equation, with the
relativistic reduced mass, is integrated outwards by Numerov's method and matched to the
Coulomb functions outside the potential and the classical turning point.

# Input Argument(s)
- `iso::IsotopeOpticalModel` : isotope model.
- `E::Float64` : laboratory kinetic energy [MeV].
- `U::OpticalPotential` : optical-model potential.

# Output Argument(s)
- `res::OPPartialWaves` : S-matrix and Coulomb phase shifts.

# Reference(s)
- Salvat and Quesada (2020), Eqs. (47)-(56); Salvat and Heredia (2024), Eqs. (15)-(26).
"""
function op_partial_waves(iso::IsotopeOpticalModel, E::Float64, U::OpticalPotential)
    kin = cm_kinematics(E, iso.m_mev, iso.M_mev, iso.Z1 * iso.Z2)
    k, η = kin.k, kin.eta
    spin_half = iso.spin == 0.5
    κ = 2.0 * kin.mur / OP_HBARC^2
    Zc = iso.Z1 * iso.Z2 * OP_E2
    r_nuc = op_nuclear_range(U)
    r_t = 2.0 * η / k
    if r_t > 4.0 * r_nuc
        return OPPartialWaves(kin, spin_half, op_coulomb_phases(η, 0), [1.0 + 0im], [1.0 + 0im], true)
    end

    depth = abs(U.Vv) + abs(U.Wv) + abs(U.Vd) + abs(U.Wd) + 1.5 * Zc / U.Rc
    h = min(0.02, 0.08 / sqrt(k^2 + κ * depth))
    na = ceil(Int, max(r_nuc, 1.25 * r_t) / h)
    nb = na + max(2, round(Int, (π / 2.0) / (k * h)))
    r = (1:nb) .* h
    Uc = [op_central(U, Zc, ri) for ri in r]
    Uso = spin_half ? [op_spin_orbit(U, ri) for ri in r] : zeros(ComplexF64, nb)

    Lcap = ceil(Int, k * r[nb] + 4.0 * (k * r[nb])^(1.0 / 3.0) + 40.0)
    Fa, Ga = op_coulomb_FG(η, k * na * h, Lcap)
    Fb, Gb = op_coulomb_FG(η, k * nb * h, Lcap)
    σ = op_coulomb_phases(η, Lcap)

    Splus = ComplexF64[]; Sminus = ComplexF64[]
    nconv = 0
    for l in 0:Lcap
        Hpa, Hma = complex(Ga[l + 1], Fa[l + 1]), complex(Ga[l + 1], -Fa[l + 1])
        Hpb, Hmb = complex(Gb[l + 1], Fb[l + 1]), complex(Gb[l + 1], -Fb[l + 1])
        Smat(ρ) = (Hma - Hmb * ρ) / (Hpa - Hpb * ρ)
        if spin_half
            Sp = Smat(op_radial_ratio(l, Float64(l), Uc, Uso, κ, k, h, na, nb))
            Sm = l == 0 ? Sp : Smat(op_radial_ratio(l, -(l + 1.0), Uc, Uso, κ, k, h, na, nb))
        else
            Sp = Smat(op_radial_ratio(l, 0.0, Uc, Uso, κ, k, h, na, nb)); Sm = Sp
        end
        push!(Splus, Sp); push!(Sminus, Sm)
        if max(abs(Sp - 1.0), abs(Sm - 1.0)) < 1.0e-11 && l > k * r_nuc
            nconv += 1
            nconv >= 4 && break
        else
            nconv = 0
        end
    end
    return OPPartialWaves(kin, spin_half, σ[1:length(Splus)], Splus, Sminus, false)
end

"""
    op_amplitudes(res::OPPartialWaves, x::Float64)

Gives the Coulomb amplitude, the nuclear part of the direct amplitude and the spin-flip
amplitude [fm] at x = 1 - cos θ (center of mass).

# Input Argument(s)
- `res::OPPartialWaves` : S-matrix and Coulomb phase shifts.
- `x::Float64` : 1 - cos θ.

# Output Argument(s)
- `fC::ComplexF64`, `fN::ComplexF64`, `g::ComplexF64` : amplitudes [fm].

# Reference(s)
- Salvat and Quesada (2020), Eqs. (53)-(56).
"""
function op_amplitudes(res::OPPartialWaves, x::Float64)
    η, k = res.kin.eta, res.kin.k
    s2 = 0.5 * x
    fC = -η * cis(2.0 * res.σ[1] - η * log(s2)) / (2.0 * k * s2)
    res.coulomb_only && return fC, zero(ComplexF64), zero(ComplexF64)
    c = 1.0 - x
    s = sqrt(max(x * (2.0 - x), 0.0))
    Pm, P = 0.0, 1.0
    P1m, P1 = 0.0, 0.0
    fN = zero(ComplexF64); g = zero(ComplexF64)
    for l in 0:(length(res.Splus) - 1)
        e2σ = cis(2.0 * res.σ[l + 1])
        if res.spin_half
            fN += e2σ * ((l + 1) * (res.Splus[l + 1] - 1.0) + l * (res.Sminus[l + 1] - 1.0)) * P
            g += e2σ * (res.Sminus[l + 1] - res.Splus[l + 1]) * P1
        else
            fN += e2σ * (2l + 1) * (res.Splus[l + 1] - 1.0) * P
        end
        Pn = ((2l + 1) * c * P - l * Pm) / (l + 1)
        P1n = l == 0 ? -s : ((2l + 1) * c * P1 - (l + 1) * P1m) / l
        Pm, P = P, Pn
        P1m, P1 = P1, P1n
    end
    return fC, fN / (2im * k), g / (2im * k)
end

"""
    op_nuclear_factor(res::OPPartialWaves, x::Float64, identical::Bool)

Gives the nuclear factor F_nuc = dσ_nuc/dσ_R at x = 1 - cos θ: the cross-section for the bare
nucleus relative to point-Coulomb scattering. For identical spin-0 particles (α + ⁴He), the
amplitude is symmetrized, f(θ) + f(π - θ), Coulomb part included.

# Input Argument(s)
- `res::OPPartialWaves` : S-matrix and Coulomb phase shifts.
- `x::Float64` : 1 - cos θ.
- `identical::Bool` : whether projectile and target are identical spin-0 particles.

# Output Argument(s)
- `F::Float64` : nuclear factor.
"""
function op_nuclear_factor(res::OPPartialWaves, x::Float64, identical::Bool)
    fC, fN, g = op_amplitudes(res, x)
    if identical
        res.spin_half && error("Identical spin-1/2 particles are described by the Galyuzov-Kosov parameterisation.")
        fCr, fNr, _ = op_amplitudes(res, 2.0 - x)
        return abs2(fC + fN + fCr + fNr) / abs2(fC)
    end
    return 1.0 + (2.0 * real(conj(fC) * fN) + abs2(fN) + abs2(g)) / abs2(fC)
end

#----
# Galyuzov-Kosov parameterisation (protons on A ≤ 6)
#----

"""
    OPGalyuzovKosov

Galyuzov-Kosov parameterisation for one light isotope at one energy: lab momentum `p`
[GeV/c], center-of-mass momentum squared `pcm2` [GeV²] and |t|max = 4 p_cm² [GeV²]. The
tables of the paper state MeV/c, but the coefficients reproduce the published curves with p
in GeV/c, t in GeV² and dσ/dt in mb/GeV².
"""
struct OPGalyuzovKosov
    kin::CMKinematics
    Z::Int
    A::Int
    identical::Bool
    p::Float64
    pcm2::Float64
    tmax::Float64
end

function op_gk(iso::IsotopeOpticalModel, E::Float64)
    kin = cm_kinematics(E, iso.m_mev, iso.M_mev, iso.Z2)
    pcm2 = (kin.p0c / 1000.0)^2
    return OPGalyuzovKosov(kin, iso.Z2, iso.A, iso.identical, sqrt(E * (E + 2.0 * iso.m_mev)) / 1000.0,
                           pcm2, 4.0 * pcm2)
end

# p + p, strong-interaction part of Eq. (16) for one channel, y = -t (or -u) [GeV²] -> mb/GeV²
function op_gk_pp_strong(p::Float64, y::Float64)
    lp = log(p)
    Ad = (74.0 + 3.0 * (lp - 5.0)^2) / (1.0 + 3.4 / p^5) + (0.2 / p^2 + 17.0 * p) / (p^4 + 0.001 * sqrt(p))
    Bd = 8.0 * p^0.055 / (1.0 + 3.64 / p^3)
    Am = 5.0e-5 + 4000.0 / (p^4 + 1500.0 * p)
    Bm = 0.46 + 1.2e6 / (p^4 + 3.5e6 / sqrt(p))
    Ah = 5.0e-5 + 1.0e10 / (p^8 + 8.5e8 * p^2 + 1.0e10)
    Bh = 1.1 + 3.4e6 / (p^4 + 6.8e6)
    return Ad * exp(-Bd * y) + 2.0 * Am * Bm^2 * y * exp(-(Bm * y)^2) + Ah * Bh * exp(-Bh * y)
end

# p + p, empirical electromagnetic term of Eq. (16), A_e corrected by the erratum
function op_gk_pp_em(p::Float64, y::Float64)
    Ae = 1.71 + 5.24 / p + 2.78 / p^3
    Be = 10.0 * log(p) + 105.0 / sqrt(p)
    return Ae * Be * exp(-Be * sqrt(y)) / (2.0 * sqrt(y))
end

# p + A, 2 ≤ A ≤ 6: coefficients of Eq. (17), Table 2 with the erratum
function op_gk_light_coefficients(p::Float64, Ai::Int)
    A = Float64(Ai)
    a1 = 1.2e7 * A^8 + 380.0 * A^17
    a2 = 0.7 / (1.0 + 4.0e-12 * A^16)
    a3 = 2.5 / (A^12 + 1.0e-16 * A^40)
    a4 = 1.2 * A^2 + 2.3
    Ad = 4000.0 * A / (1.0 + a1 * p^(A + 4.0)) + a2 / (p^4 + a3 * p^(4.0 - 2.0 * A)) +
         (0.28 * A * (log(p) - 5.0)^2 + a4) / (1.0 + 3.8 * p^(-2.0 / A))
    b1 = 0.01 / (1.0 + 0.0024 * A^5)
    b2 = 9.0e-7 / (1.0 + 0.035 * A^5)
    b3 = (42.0 + 2.7e-11 * A^16) / (1.0 + 0.14 * A)
    Bd = (b1 + 0.2 * A * p^2) / (p^4 + b2 / p^(A / 2.0)) + b3
    c1 = 0.0024 * A^8 / (1.0 + 2.6e-4 * A^7)
    c2 = 3.5e-36 * A^40 / (1.0 + 5.0e-15 * A^31)
    Cd = 2.25 * A^3 / (1.0 + 18.0 / p^2) + c1 / (p^(6.0 - A) + c2 / p^16)
    d = 8.0e7 * A^-12 / (1.0 + 1.0e-28 * A^44)
    Am = 1.0e5 / ((A^8 + 2.5e12 / A^16) * (p^(A - 2.0) + d / p^4)) + 0.0006 * A^3
    Bm = (10.0 + 4.0e-8 * A^13) * p^0.114 + 0.003 / (p^8 + 2.0e-23 / p^16)
    Ah = 0.03 + 1.0 / ((p^(A + 1.0) + 1.5e-4 * p^(-A) / (1.0 + 5.0e-6 * A^12)) * (1.0 + 1.0e-4 * A^8))
    Bh = (A / 2.0) / (p^3 + 2.0e-7 * A^4 / p^6) + 4.0 / (1.0 + 64.0 * A^-3 / p^2)
    Au = p^2 * (1.0e8 * p^(A / 2.0) * exp(0.32 * A * sqrt(A) - 20.0 * p^(A / 2.0) * exp(0.45 * A * sqrt(A))) +
                (7000.0 + 2.4e6 / A^5) / (1.0 + 2.5e5 * exp(0.085 * A^3) * p^(2.5 * A)))
    Bu = p^(A - 2.0) * (920.0 + 0.03 * A^11) / (1.0 + p^A * (93.0 + 0.0023 * A^12))
    Ae = 1.5 + 0.2 * (A - 1.0) / p^4
    Be = 10.0 * log(p) + 210.0 / (A * sqrt(p))
    return (; Ad, Bd, Cd, Am, Bm, Ah, Bh, Au, Bu, Ae, Be)
end

"""
    op_gk_dsdt(g::OPGalyuzovKosov, x::Float64; em::Bool)

Gives dσ/dt [mb/GeV²] of the Galyuzov-Kosov parameterisation at x = 1 - cos θ: the
strong-interaction terms, and with `em` also their empirical electromagnetic term (Eq. 16 for
p + ¹H, symmetrized in t and u; Eq. 17 for 2 ≤ A ≤ 6).

# Input Argument(s)
- `g::OPGalyuzovKosov` : parameterisation at one energy.
- `x::Float64` : 1 - cos θ (center of mass).
- `em::Bool` : whether to add the electromagnetic term.

# Output Argument(s)
- `dsdt::Float64` : dσ/dt [mb/GeV²].

# Reference(s)
- Galyuzov and Kosov, Phys. At. Nucl. 84 (2021) 615, Eqs. (16)-(17), Tables 1-2; erratum,
  Phys. At. Nucl. 86 (2023) 593.
"""
function op_gk_dsdt(g::OPGalyuzovKosov, x::Float64; em::Bool)
    y = g.tmax * 0.5 * x
    yu = g.tmax * (1.0 - 0.5 * x)
    if g.identical
        v = op_gk_pp_strong(g.p, y) + op_gk_pp_strong(g.p, yu)
        em && (v += op_gk_pp_em(g.p, y) + op_gk_pp_em(g.p, yu))
    else
        c = op_gk_light_coefficients(g.p, g.A)
        v = c.Ad * (c.Bd + 2.0 * c.Cd * y) * exp(-c.Bd * y - c.Cd * y^2) +
            3.0 * c.Am * c.Bm * y^2 * exp(-c.Bm * y^3) + c.Ah * c.Bh * exp(-c.Bh * y) +
            c.Au * c.Bu * exp(-c.Bu * yu)
        em && (v += c.Ae * c.Be * exp(-c.Be * sqrt(y)) / (2.0 * sqrt(y)))
    end
    return v
end

"Mott factor of identical particles relative to the projectile Rutherford cross-section (GK Eq. 14)."
function op_mott_factor(kin::CMKinematics, x::Float64, spin::Float64)
    t2 = x / (2.0 - x)
    return 1.0 + t2^2 + (-1.0)^round(Int, 2spin) / (spin + 0.5) * t2 * cos(kin.eta * log(t2))
end

"""
    op_gk_bounded_dcs(S::OPScreening, g::OPGalyuzovKosov, spin::Float64, x::Float64)

Gives the cross-section [fm²/sr] of protons on A ≤ 6 at x = 1 - cos θ. The complete
Galyuzov-Kosov formula (strong terms plus their electromagnetic term) is kept within the
bounds (√C - √N)² ≤ dσ ≤ (√C + √N)² of Coulomb-nuclear interference, with N the strong part
and C the Coulomb cross-section of the extended charges: the screened point-Coulomb
cross-section (times the Mott factor for p + ¹H) times F(q)², F(q) = exp(-q²(r_p² + r_T²)/6).
The electromagnetic term of the formula falls off much more slowly than Rutherford's at small
|t|, so that towards small angles the lower bound takes over; at large |t| both fall off like
the charge form factor. For p + ¹H, this reproduces Fig. 1 of Salvat and Heredia (2024).

# Input Argument(s)
- `S::OPScreening` : eikonal description of the screened Coulomb field.
- `g::OPGalyuzovKosov` : parameterisation at one energy.
- `spin::Float64` : projectile spin.
- `x::Float64` : 1 - cos θ (center of mass).

# Output Argument(s)
- `σ::Float64` : cross-section [fm²/sr].
"""
function op_gk_bounded_dcs(S::OPScreening, g::OPGalyuzovKosov, spin::Float64, x::Float64)
    q2 = 2.0 * g.kin.k^2 * x
    F2 = exp(-q2 * (op_charge_radius(1, 1)^2 + op_charge_radius(g.Z, g.A)^2) / 3.0)
    C = op_screened_dcs(S, x) * (g.identical ? op_mott_factor(g.kin, x, spin) : 1.0) * F2
    N = g.pcm2 / π * op_gk_dsdt(g, x; em=false) * 0.1
    full = g.pcm2 / π * op_gk_dsdt(g, x; em=true) * 0.1
    return clamp(full, (sqrt(C) - sqrt(N))^2, (sqrt(C) + sqrt(N))^2)
end

#----
# Database and tabulated cross-section
#----

"""
    init_elastic_scattering_op_db(isotopes::Vector{Tuple{Int,Int}}, particle::Particle;
    screening::Bool=true)

Builds the optical-model description of elastic scattering for one incident particle (proton
or alpha) on the given isotopes. The model of each isotope is the Galyuzov-Kosov
parameterisation for protons on A ≤ 6 and partial waves otherwise.

# Input Argument(s)
- `isotopes::Vector{Tuple{Int,Int}}` : target isotope pairs `(Z, A)`.
- `particle::Particle` : incident particle.
- `screening::Bool` : whether the Coulomb field is screened by the atomic electrons.

# Output Argument(s)
- `db::ElasticScatteringOPDB` : optical-model description.
"""
function init_elastic_scattering_op_db(isotopes::Vector{Tuple{Int,Int}}, particle::Particle;
    screening::Bool=true)
    (is_proton(particle) || is_alpha(particle)) || error("The optical-model elastic scattering is defined for protons and alphas.")
    Z1 = get_charge(particle)
    spin = is_proton(particle) ? 0.5 : 0.0
    m = get_mass(particle)
    zap = ZAP(particle)
    db = Dict{Tuple{Int,Int},IsotopeOpticalModel}()
    for (Z, A) in isotopes
        screening && Z > 92 && error("DHFS screening parameters are tabulated for Z = 1-92 only (Z = $(Z)).")
        method = is_proton(particle) && A <= 6 ? :galyuzov_kosov : :partial_waves
        db[(Z, A)] = IsotopeOpticalModel(Z1, Z, A, spin, m, get_mass(Z, A), zap == 1000 * Z + A, method)
    end
    return ElasticScatteringOPDB(screening, db)
end

"""
    op_dcs_function(iso::IsotopeOpticalModel, E::Float64, screening::Bool)

Gives the cross-section of one isotope at one energy as a function of x = 1 - cos θ, in
barn/sr: the screened Coulomb cross-section times the nuclear factor from partial waves, or
the bounded Galyuzov-Kosov cross-section for protons on A ≤ 6. The S-matrix and the eikonal
quadrature are computed once here.

# Input Argument(s)
- `iso::IsotopeOpticalModel` : isotope model.
- `E::Float64` : laboratory kinetic energy [MeV].
- `screening::Bool` : whether the Coulomb field is screened.

# Output Argument(s)
- `f::Function` : x -> dσ/dΩ [b/sr].
- `kin::CMKinematics` : center-of-mass quantities.
"""
function op_dcs_function(iso::IsotopeOpticalModel, E::Float64, screening::Bool)
    if iso.method == :galyuzov_kosov
        g = op_gk(iso, E)
        S = op_screening(iso.Z2, g.kin, screening)
        return x -> OP_FM2_TO_BARN * op_gk_bounded_dcs(S, g, iso.spin, x), g.kin
    end
    res = op_partial_waves(iso, E, op_potential(iso, E))
    S = op_screening(iso.Z2, res.kin, screening)
    return x -> OP_FM2_TO_BARN * op_screened_dcs(S, x) * op_nuclear_factor(res, x, iso.identical), res.kin
end

"""
    prepare_op_dsig_cache(db::ElasticScatteringOPDB, Z::Int, A::Int, E_in_MeV::Float64)

Tabulates the optical-model cross-section of one isotope at one incident energy in
x = 1 - μ_cm. The start grid has 10 nodes per decade from x₂ = 10⁻⁴ x_scr (x_scr = θ_scr²/2,
θ_scr = 1/(k R_at) the screening angle; 10⁻¹² without screening) to x_max (2, or 1 for
identical particles), plus 8 nodes per diffraction period, Δθ = π/(8 k R), beyond
θ = 1/(k R). Each piece is then bisected at its log-midpoint until the power-law interpolation
there agrees with the cross-section within `OP_TABLE_RTOL`. The first node is x₁ = 0: on
[0, x₂] the screened cross-section is flat and is taken constant, equal to its value at x₂;
without screening it is taken as zero there, as the ENDF description stops the Coulomb term
at μ = 1 - 10⁻¹².

# Input Argument(s)
- `db::ElasticScatteringOPDB` : optical-model description.
- `Z::Int` : target atomic number.
- `A::Int` : target mass number.
- `E_in_MeV::Float64` : incident kinetic energy [MeV].

# Output Argument(s)
- `cache::OpticalDCSCache` : tabulated cross-section.
"""
function prepare_op_dsig_cache(db::ElasticScatteringOPDB, Z::Int, A::Int, E_in_MeV::Float64)
    iso = db.isotopes[(Z, A)]
    f, kin = op_dcs_function(iso, E_in_MeV, db.screening)
    x_max = iso.identical ? 1.0 : 2.0
    x_first = db.screening ? 1.0e-4 * 0.5 / (kin.k * Z^(-1.0 / 3.0) * OP_A0_FM)^2 : 1.0e-12

    # Start grid
    nodes = collect(exp10.(range(log10(x_first), log10(x_max), length=ceil(Int, 10 * log10(x_max / x_first)) + 1)))
    R = 1.2 * (A^(1.0 / 3.0) + (iso.m_mev / AMU_MEV)^(1.0 / 3.0))
    θ_max = iso.identical ? π / 2.0 : π
    for θ in (1.0 / (kin.k * R)):(π / (8.0 * kin.k * R)):θ_max
        push!(nodes, 2.0 * sin(θ / 2.0)^2)
    end
    sort!(nodes)
    nodes = nodes[[true; diff(nodes) .> 1.0e-9 .* nodes[2:end]]]
    nodes[end] = x_max

    # Refinement at the log-midpoints, relative to the larger of the true and interpolated
    # values (the cross-section spans some 25 decades, so no absolute floor; an exact zero is
    # isolated by the bisection limit and the node cap)
    X = Float64[]; Σ = Float64[]
    stack = [(nodes[i], f(nodes[i])) for i in length(nodes):-1:1]
    xl, σl = pop!(stack)
    push!(X, xl); push!(Σ, σl)
    while !isempty(stack)
        xr, σr = stack[end]
        xm = sqrt(xl * xr)
        σm = f(xm)
        σi = sqrt(max(σl, 0.0) * max(σr, 0.0))
        refine = abs(σm - σi) > OP_TABLE_RTOL * max(σm, σi) &&
                 length(X) + length(stack) < OP_TABLE_MAX_NODES && xr / xl > 1.0 + 1.0e-9
        if refine
            push!(stack, (xm, σm))
        else
            pop!(stack)
            push!(X, xr); push!(Σ, σr)
            xl, σl = xr, σr
        end
    end
    length(X) >= OP_TABLE_MAX_NODES &&
        @warn "Optical-model elastic cross-section of Z=$(Z), A=$(A) at $(E_in_MeV) MeV: the tabulation stopped at $(OP_TABLE_MAX_NODES) nodes." maxlog=5

    # Power laws need positive values
    Σ .= max.(Σ, floatmin(Float64))
    pushfirst!(X, 0.0)
    pushfirst!(Σ, db.screening ? Σ[1] : 0.0)
    p = zeros(length(X) - 1)
    for i in 2:(length(X) - 1)
        p[i] = log(Σ[i + 1] / Σ[i]) / log(X[i + 1] / X[i])
    end
    return OpticalDCSCache(E_in_MeV, iso.m_mev, iso.M_mev, iso.identical, X, Σ, p)
end

#----
# Interpolation and integration of the tabulated cross-section
#----

"Piece of the table holding x: 1 for [x₁, x₂], i for [xᵢ, xᵢ₊₁]."
op_piece(cache::OpticalDCSCache, x::Float64) = clamp(searchsortedlast(cache.x_grid, x), 1, length(cache.x_grid) - 1)

"Interpolated cross-section [b/sr] at x = 1 - μ_cm on piece i."
function op_sigma(cache::OpticalDCSCache, i::Int, x::Float64)
    i == 1 && return cache.sigma_grid[1]
    return cache.sigma_grid[i] * (x / cache.x_grid[i])^cache.p[i]
end

"""
    elastic_dsig(cache::OpticalDCSCache, mu_cm::Float64)

Gives the optical-model elastic differential cross-section, interpolated from its
tabulation.

# Input Argument(s)
- `cache::OpticalDCSCache` : tabulated cross-section.
- `mu_cm::Float64` : center-of-mass scattering cosine.

# Output Argument(s)
- `σs::Float64` : elastic differential cross-section in barn/sr.
"""
function elastic_dsig(cache::OpticalDCSCache, mu_cm::Float64)
    x = 1.0 - mu_cm
    return op_sigma(cache, op_piece(cache, x), x)
end

# ∫ σ xⁿ dx over [a, b] inside piece i, in closed form (n = 0, 1)
function op_piece_moment(cache::OpticalDCSCache, i::Int, a::Float64, b::Float64, n::Int)
    b <= a && return 0.0
    i == 1 && return cache.sigma_grid[1] * (b^(n + 1) - a^(n + 1)) / (n + 1)
    p = cache.p[i]
    Fa = op_sigma(cache, i, a) * a^(n + 1)
    Fb = op_sigma(cache, i, b) * b^(n + 1)
    q = p + n + 1
    abs(q) < 1.0e-10 && return Fa * log(b / a)
    return (Fb - Fa) / q
end

"""
    op_moment(cache::OpticalDCSCache, xa::Float64, xb::Float64, n::Int)

Gives ∫ σ(x) xⁿ dx over [xa, xb], x = 1 - μ_cm, in closed form on every piece of the
tabulation (power law, or constant on the first piece).

# Input Argument(s)
- `cache::OpticalDCSCache` : tabulated cross-section.
- `xa::Float64` : lower bound.
- `xb::Float64` : upper bound.
- `n::Int` : power of x (0 or 1).

# Output Argument(s)
- `I::Float64` : integral [b].
"""
function op_moment(cache::OpticalDCSCache, xa::Float64, xb::Float64, n::Int)
    xb <= xa && return 0.0
    X = cache.x_grid
    acc = 0.0
    for i in op_piece(cache, xa):op_piece(cache, xb)
        acc += op_piece_moment(cache, i, max(xa, X[i]), min(xb, X[i + 1]), n)
    end
    return acc
end

"""
    op_quadrature(cache::OpticalDCSCache, w::Function, xa::Float64, xb::Float64, n_gauss::Int)

Gives ∫ σ(x) w(x) dx over [xa, xb], x = 1 - μ_cm, for a weight w (scalar or vector) that is
not polynomial in x, piece by piece: in ln x on the power-law pieces, in x on the first one.
`OP_QUADRATURE[]` selects adaptive Gauss-Kronrod (`:adaptive`) or fixed Gauss-Legendre with
`n_gauss` points (`:gauss`).

# Input Argument(s)
- `cache::OpticalDCSCache` : tabulated cross-section.
- `w::Function` : weight w(x).
- `xa::Float64` : lower bound.
- `xb::Float64` : upper bound.
- `n_gauss::Int` : number of Gauss-Legendre points per piece.

# Output Argument(s)
- `I` : integral, of the type returned by `w`.
"""
function op_quadrature(cache::OpticalDCSCache, w::Function, xa::Float64, xb::Float64, n_gauss::Int)
    X = cache.x_grid
    adaptive = OP_QUADRATURE[] == :adaptive
    u_g, w_g = quadrature(n_gauss, "gauss-legendre")
    acc = nothing
    for i in op_piece(cache, xa):op_piece(cache, xb)
        a, b = max(xa, X[i]), min(xb, X[i + 1])
        b <= a && continue
        if i == 1 || a == 0.0
            g = x -> op_sigma(cache, i, x) * w(x)
            lo, hi = a, b
        else
            g = t -> (x = exp(t); op_sigma(cache, i, x) * w(x) * x)
            lo, hi = log(a), log(b)
        end
        if adaptive
            I = quadgk(g, lo, hi; rtol=1.0e-12, maxevals=30_000)[1]
        else
            c, h = 0.5 * (hi + lo), 0.5 * (hi - lo)
            I = h * w_g[1] * g(c + h * u_g[1])
            for q in 2:n_gauss
                I = I + h * w_g[q] * g(c + h * u_g[q])
            end
        end
        acc = isnothing(acc) ? I : acc + I
    end
    return isnothing(acc) ? zero(w(max(xa, 1.0e-300))) : acc
end

"""
    elastic_group_moments!(σl::AbstractMatrix{Float64}, σe::AbstractVector{Float64},
    cache::OpticalDCSCache, Eout_mev::Vector{Float64}, mu_lo::Float64, mu_hi::Float64,
    L::Int64; recoil::Bool=false)

Adds, for every outgoing energy group, the Legendre moments and the energy moment of the
optical-model cross-section over the center-of-mass cosines [mu_lo, mu_hi] of the scattered
projectile, or with `recoil` of the recoiling target of identical particles. The range is cut
at the nodes of the tabulation and at the cosines of the group boundaries (the outgoing
energy is affine in the cosine). On every piece the moment of order zero and the energy
moment are closed forms; the moments of order l ≥ 1 are integrated by `op_quadrature`, since
the laboratory cosine is not polynomial in x. Pieces outside the groups are skipped.

# Input Argument(s)
- `σl::AbstractMatrix{Float64}` : Legendre moments to add to, size (Ng, L+1), in barn.
- `σe::AbstractVector{Float64}` : energy moments to add to, size Ng, in barn × MeV.
- `cache::OpticalDCSCache` : tabulated cross-section.
- `Eout_mev::Vector{Float64}` : outgoing group boundaries [MeV], decreasing.
- `mu_lo::Float64` : lower center-of-mass cosine bound of the projectile.
- `mu_hi::Float64` : upper center-of-mass cosine bound of the projectile.
- `L::Int64` : Legendre truncation order.
- `recoil::Bool` : boolean to give the moments of the recoiling target.

# Output Argument(s)
N/A
"""
function elastic_group_moments!(σl::AbstractMatrix{Float64}, σe::AbstractVector{Float64},
        cache::OpticalDCSCache, Eout_mev::Vector{Float64}, mu_lo::Float64, mu_hi::Float64, L::Int64;
        recoil::Bool=false)
    mu_hi <= mu_lo && return nothing
    Ng = length(Eout_mev) - 1
    e0, e1 = scattered_energy_affine(cache)
    sgn = recoil ? -1.0 : 1.0                 # the particle followed leaves at sgn·μ
    energy(x) = e0 + sgn * e1 * (1.0 - x)
    lab(x) = mu_lab_from_mu_cm(sgn * (1.0 - x), cache.E_in_mev, cache.m_mev, cache.M_mev)

    xa, xb = 1.0 - mu_hi, 1.0 - mu_lo
    cuts = [xa; filter(x -> xa < x < xb, cache.x_grid); xb]
    for E_b in Eout_mev
        # x of the group boundary, written without cancellation (e0 + e1 = E_in)
        x_b = recoil ? 1.0 + (E_b - e0) / e1 : (cache.E_in_mev - E_b) / e1
        xa < x_b < xb && push!(cuts, x_b)
    end
    sort!(unique!(cuts))

    for k in 1:(length(cuts) - 1)
        x0, x1 = cuts[k], cuts[k + 1]
        x1 <= x0 && continue
        g = searchsortedfirst(Eout_mev, energy(0.5 * (x0 + x1)); rev=true) - 1
        (1 <= g <= Ng) || continue
        I0 = op_moment(cache, x0, x1, 0)
        I1 = op_moment(cache, x0, x1, 1)
        σl[g, 1] += 2π * I0
        σe[g] += 2π * ((e0 + sgn * e1) * I0 - sgn * e1 * I1)
        if L >= 1
            I = op_quadrature(cache, x -> legendre_polynomials_up_to_L(L, lab(x)), x0, x1, max(L + 2, 8))
            for l in 1:L
                σl[g, l + 1] += 2π * I[l + 1]
            end
        end
    end
    return nothing
end

"""
    elastic_soft_integrals(cache::OpticalDCSCache, mu_lo::Float64)

Gives the two integrals of the soft elastic collisions, center-of-mass cosines from `mu_lo`
to 1, on which the stopping power and the momentum transfer rest: ∫ dσ/dΩ (1 - μ) dμ, in
closed form, and ∫ dσ/dΩ (1 - μ_lab) dμ, by `op_quadrature` with 1 - μ_lab written without
cancellation (`one_minus_mu_lab`).

# Input Argument(s)
- `cache::OpticalDCSCache` : tabulated cross-section.
- `mu_lo::Float64` : lower center-of-mass cosine bound, the cutoff cosine μ_c.

# Output Argument(s)
- `I_cm::Float64` : ∫ dσ/dΩ (1 - μ) dμ over [mu_lo, 1], in barn.
- `I_lab::Float64` : ∫ dσ/dΩ (1 - μ_lab) dμ over [mu_lo, 1], in barn.
"""
function elastic_soft_integrals(cache::OpticalDCSCache, mu_lo::Float64)
    mu_lo >= 1.0 && return 0.0, 0.0
    xb = 1.0 - mu_lo
    frame = cm_frame(cache)
    I_cm = op_moment(cache, 0.0, xb, 1)
    I_lab = op_quadrature(cache, x -> one_minus_mu_lab(frame, x), 0.0, xb, 8)
    return I_cm, I_lab
end

"""
    integrate_elastic_scattering_below(cache::OpticalDCSCache, Ec_MeV::Float64,
    E_last_mev::Float64)

Gives the part of the total optical-model elastic cross-section above an outgoing-energy
cutoff whose scattered projectile ends below the last outgoing group, E_out < `E_last_mev`,
in closed form. It completes the sum of the moments of order zero of
`elastic_group_moments!` to the total cross-section (`tcs`).

# Input Argument(s)
- `cache::OpticalDCSCache` : tabulated cross-section.
- `Ec_MeV::Float64` : cutoff energy between soft and catastrophic events.
- `E_last_mev::Float64` : lower boundary of the last outgoing group [MeV].

# Output Argument(s)
- `σb::Float64` : cross-section below the last group, in barn.
"""
function integrate_elastic_scattering_below(cache::OpticalDCSCache, Ec_MeV::Float64,
                                            E_last_mev::Float64)
    mu_lo = cache.identical ? 0.0 : -1.0
    mu_c = mu_cm_from_Eout_relativistic(cache.E_in_mev, max(Ec_MeV, 0.0), cache.m_mev, cache.M_mev)
    mu_c <= mu_lo && return 0.0
    e0, e1 = scattered_energy_affine(cache)
    mu_cut = clamp((E_last_mev - e0) / e1, mu_lo, mu_c)
    return 2.0 * pi * op_moment(cache, 1.0 - mu_cut, 1.0 - mu_lo, 0)
end
