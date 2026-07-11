function gn_weights_spherical_harmonics(L::Int64,Nv::Int64,Ndims::Int64;tiling::String="polar-anchored")
    if L < 0 error("Legendre order is greater or equal to zero.") end
    if Nv <= 0 error("Number of direction cosine patches should be greater than zero.") end
    if ~(1 ≤ Ndims ≤ 3) error("Number of dimensions should be between 1 and 3.") end
    if tiling == "symmetric"
        return gn_weights_spherical_harmonics_symmetric(L,Nv,Ndims)
    end
    Np = spherical_harmonics_number_basis(L)
    pl,pm = spherical_harmonics_indices(L)
    𝒩 = zeros(Np,Np,Ndims,8,Nv,Nv)
    N = 32
    x,weight = gauss_legendre(N)
    t = 0.5 .* (x .+ 1.0)

    sx = [1,1,1,1,-1,-1,-1,-1]
    sy = [1,1,-1,-1,1,1,-1,-1]
    sz = [1,-1,1,-1,1,-1,1,-1]

    # Offsets azimutaux par octant (ne dépend que de u)
    ϕ_offset_u = Vector{Float64}(undef, 8)
    for u in 1:8
        ϕ_offset_u[u] = (π / 2) * (2 + (sy[u] + 1) / 2 - (sz[u] + 1) / 2 - (sy[u] + 1) * (sz[u] + 1) / 2)
    end

    # Pré-calcul des valeurs 𝒯m(m,4*φp) avec 4*φp = π*(x+1)
    # et tables azimutales sur les m uniques
    m_unique = sort!(collect(Set(pm)))
    m_to_idx = Dict{Int64,Int}(m => i for (i, m) in enumerate(m_unique))
    m_idx_p = [m_to_idx[pm[p]] for p in 1:Np]

    Tm = Matrix{Float64}(undef, length(m_unique), N)
    for (im, m) in enumerate(m_unique)
        if m ≥ 0
            for n in 1:N
                Tm[im, n] = cos(m * (π * (x[n] + 1.0)))
            end
        else
            for n in 1:N
                Tm[im, n] = sin((-m) * (π * (x[n] + 1.0)))
            end
        end
    end
    TmW = Matrix{Float64}(undef, length(m_unique), N)
    for n in 1:N
        wn = weight[n]
        for im in 1:length(m_unique)
            TmW[im, n] = Tm[im, n] * wn
        end
    end
    az0_table = Matrix{Float64}(undef, length(m_unique), length(m_unique))
    LinearAlgebra.mul!(az0_table, TmW, LinearAlgebra.transpose(Tm))
    az0_pq = Matrix{Float64}(undef, Np, Np)
    for p in 1:Np
        imp = m_idx_p[p]
        for q in 1:Np
            imq = m_idx_p[q]
            az0_pq[p, q] = az0_table[imp, imq]
        end
    end

    # Pré-calcul des polynômes associés Pnm(l,|m|,±x)
    P_pos = Matrix{Float64}(undef, Np, N)
    P_neg = Matrix{Float64}(undef, Np, N)
    for p in 1:Np
        lp = pl[p]
        mp = pm[p]
        amp = abs(mp)
        for n in 1:N
            P_pos[p, n] = Pnm(lp, amp, x[n])
            P_neg[p, n] = Pnm(lp, amp, -x[n])
        end
    end

    # Factorisation des constantes Cp*Cq (en conservant la dépendance originale en mq==0)
    cp_base = Vector{Float64}(undef, Np)
    for p in 1:Np
        lp = pl[p]
        mp = pm[p]
        cp_base[p] = sqrt((2 * lp + 1) * factorial_factor([lp - abs(mp)], [lp + abs(mp)]))
    end
    cq_base = Vector{Float64}(undef, Np)
    mq_factor = Vector{Float64}(undef, Np)
    for q in 1:Np
        lq = pl[q]
        mq = pm[q]
        cq_base[q] = sqrt((2 * lq + 1) * factorial_factor([lq - abs(mq)], [lq + abs(mq)]))
        mq_factor[q] = 2 * (2 - (mq == 0)) / π
    end
    sq = Vector{Float64}(undef, Np)
    for q in 1:Np
        sq[q] = mq_factor[q] * cq_base[q]
    end

    # Buffers réutilisés
    Pw = Matrix{Float64}(undef, Np, N)
    cosine_x = Matrix{Float64}(undef, Np, Np)
    cosine_yz = Matrix{Float64}(undef, Np, Np)

    cw = Vector{Float64}(undef, N)
    sw = Vector{Float64}(undef, N)
    cos_dt = Vector{Float64}(undef, N)
    sin_dt = Vector{Float64}(undef, N)

    TmC = Matrix{Float64}(undef, length(m_unique), N)
    az_tmp = Matrix{Float64}(undef, length(m_unique), length(m_unique))

    # Helper μ(v) : évite les closures μ_uv/ϕ_uw
    denom = Float64(Nv * (Nv + 1))
    mu_of_v(s::Int, v::Int) = (s == 1) ? (1.0 - ((Nv + 1 - v) * (Nv + 2 - v)) / denom) : (-1.0 + ((v - 1) * v) / denom)

    # Simplification: (Δμ*Δϕ)/4 * (π/(2ΔμΔϕ)) = π/8
    pref = π / 8

    for u in 1:8
        su = sx[u]
        Psign = (su == 1) ? P_pos : P_neg
        ϕ_offset = ϕ_offset_u[u]

        for v in 1:Nv
            Nw = (su == 1) ? (Nv + 1 - v) : v
            Δϕ = (π / 2) / Nw

            # Cos(Δϕ * t[n]) et Sin(Δϕ * t[n]) (utiles pour Ndims>=2)
            if Ndims ≥ 2
                for n in 1:N
                    θ = Δϕ * t[n]
                    cos_dt[n] = cos(θ)
                    sin_dt[n] = sin(θ)
                end
            end

            # Poids pour les intégrales cosinus (dépend de v, pas de w)
            μ0 = mu_of_v(su, v)
            μ1 = mu_of_v(su, v + 1)
            Δμ = μ1 - μ0

            # cosine_integral_x : poids = weight[n] * μi
            for p in 1:Np
                for n in 1:N
                    μi = μ0 + Δμ * t[n]
                    Pw[p, n] = Psign[p, n] * (weight[n] * μi)
                end
            end
            LinearAlgebra.mul!(cosine_x, Pw, LinearAlgebra.transpose(Psign))

            if Ndims ≥ 2
                # cosine_integral_yz : poids = weight[n] * sqrt(1-μi^2)
                for p in 1:Np
                    for n in 1:N
                        μi = μ0 + Δμ * t[n]
                        Pw[p, n] = Psign[p, n] * (weight[n] * sqrt(1 - μi^2))
                    end
                end
                LinearAlgebra.mul!(cosine_yz, Pw, LinearAlgebra.transpose(Psign))
            end

            # Ndims == 1 : tout est indépendant de w (par (u,v)) sauf la sortie elle-même
            if Ndims == 1
                for w in 1:Nw
                    for p in 1:Np
                        sp = pref * cp_base[p]
                        for q in 1:Np
                            𝒩[p, q, 1, u, v, w] += (sp * sq[q]) * az0_pq[p, q] * cosine_x[p, q]
                        end
                    end
                end
                continue
            end

            # Ndims >= 2 : partie azimutale dépend de w via cos(φi)/sin(φi)
            cosΔ = cos(Δϕ)
            sinΔ = sin(Δϕ)
            cosϕ0 = cos(ϕ_offset)
            sinϕ0 = sin(ϕ_offset)

            for w in 1:Nw
                # cw/sw = weight[n] * cos(φi_n) / sin(φi_n) avec φi_n = ϕ0 + Δϕ*t[n]
                for n in 1:N
                    c = cosϕ0 * cos_dt[n] - sinϕ0 * sin_dt[n]
                    s = sinϕ0 * cos_dt[n] + cosϕ0 * sin_dt[n]
                    cw[n] = weight[n] * c
                    sw[n] = weight[n] * s
                end

                # azimutal_integral_x est constant (az0_pq), azimutal_integral_y/z via cw/sw
                for n in 1:N
                    cwn = cw[n]
                    swn = sw[n]
                    for im in 1:length(m_unique)
                        TmC[im, n] = Tm[im, n] * cwn
                    end
                end
                LinearAlgebra.mul!(az_tmp, TmC, LinearAlgebra.transpose(Tm))

                # Remplissage dim 1 & 2
                for p in 1:Np
                    sp = pref * cp_base[p]
                    imp = m_idx_p[p]
                    for q in 1:Np
                        imq = m_idx_p[q]
                        common = (sp * sq[q])
                        𝒩[p, q, 1, u, v, w] += common * az0_pq[p, q] * cosine_x[p, q]
                        𝒩[p, q, 2, u, v, w] += common * az_tmp[imp, imq] * cosine_yz[p, q]
                    end
                end

                if Ndims == 3
                    for n in 1:N
                        swn = sw[n]
                        for im in 1:length(m_unique)
                            TmC[im, n] = Tm[im, n] * swn
                        end
                    end
                    LinearAlgebra.mul!(az_tmp, TmC, LinearAlgebra.transpose(Tm))

                    for p in 1:Np
                        sp = pref * cp_base[p]
                        imp = m_idx_p[p]
                        for q in 1:Np
                            imq = m_idx_p[q]
                            𝒩[p, q, 3, u, v, w] += (sp * sq[q]) * az_tmp[imp, imq] * cosine_yz[p, q]
                        end
                    end
                end

                # ϕ0 += Δϕ (update trig via récurrence)
                cosϕ1 = cosϕ0 * cosΔ - sinϕ0 * sinΔ
                sinϕ1 = sinϕ0 * cosΔ + cosϕ0 * sinΔ
                cosϕ0 = cosϕ1
                sinϕ0 = sinϕ1
            end
        end
    end
    return 𝒩
end

"""
    gn_weights_spherical_harmonics_2D_quarter(L_elem::Int64)

Compute the per-quadrant streaming-weight matrices `𝒩[p, q, d, u, 1, 1]` (`d = 1, 2`
for x- and y-streaming) for the 2D spherical-harmonics GN basis over four quadrants,
the four-quadrant counterpart of the octant streaming-weight builder. The 2D angular flux is
symmetric under `μ_z → -μ_z`, so each domain is a full-`μ_z` quadrant of the sphere
selected by `(sign μ_x, sign μ_y)`. The four quadrants are stored in octant slots
`u ∈ {1, 3, 5, 7}` (one representative per `(sx, sy)` sign pair); the remaining slots
stay zero. Each quadrant spans the azimuthal range `Δϕ = π`. With `L_elem = L` the
`u = 1` slice reproduces the DPN four-quadrant (Double-PN) streaming weights.

# Input Argument(s)
- `L_elem::Int64` : per-patch local Legendre order (`Nq = (L_elem+1)(L_elem+2)/2`).

# Output Argument(s)
- `𝒩::Array{Float64,6}` : streaming weights, shape `(Nq, Nq, 2, 8, 1, 1)`.
"""
function gn_weights_spherical_harmonics_2D_quarter(L_elem::Int64)
    if L_elem < 0 error("Local Legendre order must be ≥ 0.") end
    Nq = spherical_harmonics_number_basis(L_elem)
    pl, pm = spherical_harmonics_indices(L_elem)
    𝒩 = zeros(Nq, Nq, 2, 8, 1, 1)
    N = 32
    x, weight = gauss_legendre(N)
    t = 0.5 .* (x .+ 1.0)
    C = Vector{Float64}(undef, Nq)
    for q in 1:Nq
        lq = pl[q]; mq = pm[q]
        C[q] = sqrt(2 * (2 - (mq == 0)) / π * (2 * lq + 1) * factorial_factor([lq - abs(mq)], [lq + abs(mq)]))
    end
    # Patch-local azimuthal basis 𝒯m(m, π(x+1)) (reference over [0,2π]).
    Tref = Matrix{Float64}(undef, Nq, N)
    for q in 1:Nq
        mq = pm[q]
        for n in 1:N
            Tref[q, n] = (mq ≥ 0) ? cos(mq * π * (x[n] + 1.0)) : sin(-mq * π * (x[n] + 1.0))
        end
    end
    P_pos = Matrix{Float64}(undef, Nq, N)
    P_neg = Matrix{Float64}(undef, Nq, N)
    for q in 1:Nq
        lq = pl[q]; amq = abs(pm[q])
        for n in 1:N
            P_pos[q, n] = Pnm(lq, amq, x[n])
            P_neg[q, n] = Pnm(lq, amq, -x[n])
        end
    end
    pref = π / 8
    for u in (1, 3, 5, 7)
        su = _GN_SX[u]; sv = _GN_SY[u]
        Psign = (su == 1) ? P_pos : P_neg
        μ0 = (su == 1) ? 0.0 : -1.0; μ1 = (su == 1) ? 1.0 : 0.0; Δμ = μ1 - μ0
        ϕ0 = (sv == 1) ? (-π / 2) : (π / 2)
        for p in 1:Nq
            for q in 1:Nq
                azx = 0.0; azy = 0.0; cx = 0.0; cyz = 0.0
                for n in 1:N
                    μi = μ0 + Δμ * t[n]
                    PpPq = Psign[p, n] * Psign[q, n]
                    cx += weight[n] * μi * PpPq
                    cyz += weight[n] * sqrt(1 - μi^2) * PpPq
                    TpTq = Tref[p, n] * Tref[q, n]
                    azx += weight[n] * TpTq
                    azy += weight[n] * cos(ϕ0 + π * t[n]) * TpTq
                end
                𝒩[p, q, 1, u, 1, 1] += pref * C[p] * C[q] * azx * cx
                𝒩[p, q, 2, u, 1, 1] += pref * C[p] * C[q] * azy * cyz
            end
        end
    end
    return 𝒩
end

"""
    gn_weights_spherical_harmonics_symmetric(L, Nv, Ndims)

Compute the Galerkin patch-weight matrix `𝒩[p, q, d, u, i, j]` for the symmetric
tiling: each octant is barycentrically subdivided into `Nv²` spherical sub-triangles
treating the three axis-vertices equivalently.

Patches are indexed by `(u, i, j)` with `i ∈ 1:Nv` (barycentric row) and
`j ∈ 1:(2i-1)` (position within the row, alternating up/down). The third axis of
the returned array has length `2*Nv-1`; entries with `j > 2i-1` remain zero.

Each sub-triangle integral is computed via a Duffy-chart Gauss-Legendre quadrature
of order 32 (1024 points per patch). Unlike the polar-anchored variant, the
polar/azimuthal BLAS factorisation does not apply here — the patches are curved
spherical triangles not aligned with `(μ, ϕ)` axes.
"""
function gn_weights_spherical_harmonics_symmetric(L::Int64,Nv::Int64,Ndims::Int64)
    Np = spherical_harmonics_number_basis(L)
    Nw_max = 2*Nv - 1
    𝒩 = zeros(Np,Np,Ndims,8,Nv,Nw_max)
    N = 32
    x,weight = gauss_legendre(N)
    Nquad = N * N

    Ψw = Matrix{Float64}(undef, Np, Nquad)

    for u in 1:8, i in 1:Nv
        for j in 1:(2i - 1)
            Ψ, Px, Py, Pz, JW = symmetric_patch_orthonormal_basis(L, Nv, u, i, j, x, weight)

            # Direction 1: Ω_x = Px
            for k in 1:Nquad
                wkx = Px[k] * JW[k]
                for p in 1:Np
                    Ψw[p, k] = Ψ[p, k] * wkx
                end
            end
            LinearAlgebra.mul!(view(𝒩, :, :, 1, u, i, j), Ψw, LinearAlgebra.transpose(Ψ))

            if Ndims ≥ 2
                for k in 1:Nquad
                    wky = Py[k] * JW[k]
                    for p in 1:Np
                        Ψw[p, k] = Ψ[p, k] * wky
                    end
                end
                LinearAlgebra.mul!(view(𝒩, :, :, 2, u, i, j), Ψw, LinearAlgebra.transpose(Ψ))
            end

            if Ndims == 3
                for k in 1:Nquad
                    wkz = Pz[k] * JW[k]
                    for p in 1:Np
                        Ψw[p, k] = Ψ[p, k] * wkz
                    end
                end
                LinearAlgebra.mul!(view(𝒩, :, :, 3, u, i, j), Ψw, LinearAlgebra.transpose(Ψ))
            end
        end
    end

    return 𝒩
end

"""
    gn_weights_legendre_1D(L_elem::Int64, Nv::Int64)

Compute the per-patch streaming-weight matrix `𝒩[p, q, 1, u, v, 1]` for the 1D
Legendre (azimuthally symmetric) GN basis, the Legendre counterpart of
`gn_weights_spherical_harmonics`. Here
`𝒩[p,q,1,u,v,1] = ∫_band μ · ψ_p^loc(μ) ψ_q^loc(μ) dμ` in the patch-orthonormal
Legendre basis `ψ_q^loc(μ) = sqrt((2 q_l + 1)/Δ) · P_{q_l}(ξ)`,
`ξ = (2μ - μ0 - μ1)/Δ`. Only octants `u ∈ {1, 5}` (sx = ±1) carry μ-band patches.

# Input Argument(s)
- `L_elem::Int64` : per-patch local Legendre order (`Nq = L_elem + 1`).
- `Nv::Int64` : number of μ-bands per hemisphere.

# Output Argument(s)
- `𝒩::Array{Float64,6}` : streaming weights, shape `(Nq, Nq, 1, 8, Nv, 1)`.
"""
function gn_weights_legendre_1D(L_elem::Int64, Nv::Int64)
    if L_elem < 0 error("Local Legendre order must be ≥ 0.") end
    if Nv <= 0 error("Number of μ-bands per hemisphere must be > 0.") end
    Nq = L_elem + 1
    𝒩 = zeros(Nq, Nq, 1, 8, Nv, 1)
    N = 32
    x, weight = gauss_legendre(N)
    t = 0.5 .* (x .+ 1.0)
    # Local orthonormal Legendre at the (band-independent) reference points ξ = x.
    Pξ = [legendre_polynomials_up_to_L(L_elem, x[n]) for n in 1:N]
    for u in (1, 5)
        edges = gn_1D_band_edges(u, Nv, :legendre)
        for v in 1:Nv
            μ0 = edges[v]; μ1 = edges[v + 1]
            Δ = μ1 - μ0
            for n in 1:N
                μ = μ0 + Δ * t[n]
                wn = weight[n] * (Δ / 2)
                for p in 1:Nq
                    ψp = sqrt((2 * (p - 1) + 1) / Δ) * Pξ[n][p]
                    for q in 1:Nq
                        ψq = sqrt((2 * (q - 1) + 1) / Δ) * Pξ[n][q]
                        𝒩[p, q, 1, u, v, 1] += wn * μ * ψp * ψq
                    end
                end
            end
        end
    end
    return 𝒩
end

"""
    gn_weights_spherical_harmonics_1D(L_elem::Int64, Nv::Int64)

Compute the per-patch x-streaming weight matrix `𝒩[p, q, 1, u, v, 1]` for the 1D
spherical-harmonics GN basis over two half-spheres, the SH counterpart of
`gn_weights_legendre_1D`. In 1D the polar axis is the x-axis and the angular
domain is azimuthally symmetric, so each domain is a full-azimuth (`ϕ ∈ [0, 2π]`)
μ-band of a hemisphere. Only octants `u ∈ {1, 5}` (`sx = ±1`) carry patches and
`Nw = 1` per band. Here
`𝒩[p,q,1,u,v,1] = ∫_band ∫_0^{2π} μ · ψ_p^loc(μ,ϕ) ψ_q^loc(μ,ϕ) dμ dϕ` in the
patch-local real-SH basis. By azimuthal orthogonality the matrix is block-diagonal
in `m`; the `m = 0` block (the only one excited by an azimuthally symmetric 1D
problem) reproduces the DPN half-sphere (Double-PN) streaming weights. With `Nv = 1, L_elem = L`
this matches that reference per hemisphere.

# Input Argument(s)
- `L_elem::Int64` : per-patch local Legendre order (`Nq = (L_elem+1)(L_elem+2)/2`).
- `Nv::Int64` : number of μ-bands per hemisphere.

# Output Argument(s)
- `𝒩::Array{Float64,6}` : streaming weights, shape `(Nq, Nq, 1, 8, Nv, 1)`.
"""
function gn_weights_spherical_harmonics_1D(L_elem::Int64, Nv::Int64)
    if L_elem < 0 error("Local Legendre order must be ≥ 0.") end
    if Nv <= 0 error("Number of μ-bands per hemisphere must be > 0.") end
    Nq = spherical_harmonics_number_basis(L_elem)
    pl, pm = spherical_harmonics_indices(L_elem)
    𝒩 = zeros(Nq, Nq, 1, 8, Nv, 1)
    N = 32
    x, weight = gauss_legendre(N)
    t = 0.5 .* (x .+ 1.0)
    # Local-basis azimuthal-times-cosine normalization (same as the octant
    # patch q-basis in patch_to_full_range_matrix_spherical_harmonics).
    C = Vector{Float64}(undef, Nq)
    for q in 1:Nq
        lq = pl[q]; mq = pm[q]
        C[q] = sqrt(2 * (2 - (mq == 0)) / π * (2 * lq + 1) * factorial_factor([lq - abs(mq)], [lq + abs(mq)]))
    end
    # Pre-compute Pnm(l,|m|,±x).
    P_pos = Matrix{Float64}(undef, Nq, N)
    P_neg = Matrix{Float64}(undef, Nq, N)
    for q in 1:Nq
        lq = pl[q]; amq = abs(pm[q])
        for n in 1:N
            P_pos[q, n] = Pnm(lq, amq, x[n])
            P_neg[q, n] = Pnm(lq, amq, -x[n])
        end
    end
    pref = π / 8
    for u in (1, 5)
        su = _GN_SX[u]
        Psign = (su == 1) ? P_pos : P_neg
        edges = gn_1D_band_edges(u, Nv, :spherical)
        for v in 1:Nv
            μ0 = edges[v]; μ1 = edges[v + 1]
            Δμ = μ1 - μ0
            for p in 1:Nq
                mp = pm[p]
                sp = pref * C[p]
                for q in 1:Nq
                    if pm[q] != mp continue end   # azimuthal orthogonality (full 2π)
                    az = 1.0 + (mp == 0)          # ∫_{-1}^1 𝒯m(m,π(x+1))² dx
                    cos_int = 0.0
                    for n in 1:N
                        μi = μ0 + Δμ * t[n]
                        cos_int += weight[n] * μi * Psign[p, n] * Psign[q, n]
                    end
                    𝒩[p, q, 1, u, v, 1] += sp * C[q] * az * cos_int
                end
            end
        end
    end
    return 𝒩
end

# function gn_weights_spherical_harmonics(L::Int64,Nv::Int64,Ndims::Int64)
#     if L < 0 error("Legendre order is greater or equal to zero.") end
#     if Nv <= 0 error("Number of direction cosine patches should be greater than zero.") end
#     if ~(1 ≤ Ndims ≤ 3) error("Number of dimensions should be between 1 and 3.") end
#     Np = spherical_harmonics_number_basis(L)
#     pl,pm = spherical_harmonics_indices(L)
#     𝒩 = zeros(Np,Np,Ndims,8,Nv,Nv)
#     N = 32
#     x,weight = gauss_legendre(N)
#     sx = [1,1,1,1,-1,-1,-1,-1]
#     sy = [1,1,-1,-1,1,1,-1,-1]
#     sz = [1,-1,1,-1,1,-1,1,-1]
#     μ_uv(u,v) = sx[u]*(1-(1-v+(sx[u]+1)/2*Nv)*(-v+(sx[u]+1)/2*(Nv+2))/(Nv*(Nv+1)))
#     ϕ_uw(u,v,w) = (π/2)*(w-1)/(-sx[u]*v + (sx[u]+1)/2*(Nv+1)) + π/2 * (2 + (sy[u]+1)/2 - (sz[u]+1)/2 - (sy[u]+1)*(sz[u]+1)/2)
#     for p in range(1,Np), q in range(1,Np)
#         lp = pl[p]
#         lq = pl[q]
#         mp = pm[p]
#         mq = pm[q]
#         Cp = sqrt(2*(2-(mq==0))/π * (2*lp+1) * factorial_factor([lp-abs(mp)],[lp+abs(mp)]))
#         Cq = sqrt(2*(2-(mq==0))/π * (2*lq+1) * factorial_factor([lq-abs(mq)],[lq+abs(mq)]))
#         for u in range(1,8)
#             for v in range(1,Nv)
#                 Nw = Int(-sx[u]*v + (sx[u]+1)/2*(Nv+1))
#                 for w in range(1,Nw)
#                     Δμ_uv = μ_uv(u,v+1) - μ_uv(u,v)
#                     Δϕ_uw = ϕ_uw(u,v,w+1) - ϕ_uw(u,v,w)

#                     # Azimuthal integral
#                     azimutal_integral_x = 0.0
#                     azimutal_integral_y = 0.0
#                     azimutal_integral_z = 0.0
#                     for n in range(1,N)
#                         φi = ϕ_uw(u,v,w) + 0.5 * Δϕ_uw * (x[n]+1)
#                         φp = π/2 * (φi-ϕ_uw(u,v,w))/Δϕ_uw
#                         azimutal_integral_x += weight[n] * 𝒯m(mp,4*φp) * 𝒯m(mq,4*φp)
#                         if (Ndims ≥ 2) azimutal_integral_y += weight[n] * 𝒯m(mp,4*φp) * 𝒯m(mq,4*φp) * cos(φi) end
#                         if (Ndims == 3) azimutal_integral_z += weight[n] * 𝒯m(mp,4*φp) * 𝒯m(mq,4*φp) * sin(φi) end
#                     end

#                     # Cosine integral
#                     cosine_integral_x = 0.0
#                     cosine_integral_yz = 0.0
#                     for n in range(1,N)
#                         μi = μ_uv(u,v) + 0.5 * Δμ_uv * (x[n]+1)
#                         μp = sx[u] * ((sx[u]-1)/2+(μi-μ_uv(u,v))/Δμ_uv)
#                         cosine_integral_x += weight[n] * Pnm(lp,abs(mp),2*μp-1) * Pnm(lq,abs(mq),2*μp-1) * μi
#                         if (Ndims ≥ 2) cosine_integral_yz += weight[n] * Pnm(lp,abs(mp),2*μp-1) * Pnm(lq,abs(mq),2*μp-1) * sqrt(1-μi^2) end
#                     end

#                     # Factor
#                     Cuvw = π/(2*Δμ_uv*Δϕ_uw)

#                     # Accumulate weights
#                     𝒩[p,q,1,u,v,w] += (Δμ_uv*Δϕ_uw)/4 * Cp * Cq * Cuvw * azimutal_integral_x * cosine_integral_x
#                     if (Ndims ≥ 2) 𝒩[p,q,2,u,v,w] += (Δμ_uv*Δϕ_uw)/4 * Cp * Cq * Cuvw * azimutal_integral_y * cosine_integral_yz end
#                     if (Ndims == 3) 𝒩[p,q,3,u,v,w] += (Δμ_uv*Δϕ_uw)/4 * Cp * Cq * Cuvw * azimutal_integral_z * cosine_integral_yz end

#                 end
#             end
#         end
#     end
#     return 𝒩
# end

