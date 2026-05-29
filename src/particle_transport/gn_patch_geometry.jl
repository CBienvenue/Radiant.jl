"""
    gn_patch_geometry(tiling::String, Nv::Int64)

Compute the centroids, solid angles, and edge adjacency of all angular patches
of the GN finite-element tiling. Used by the finite-difference angular
Fokker-Planck discretization to assemble a Laplace-Beltrami operator on the
patch graph (no Voronoi tessellation needed - the (u,v,w) indexing already
encodes the connectivity).

# Input Argument(s)
- `tiling::String`: either `"polar-anchored"` or `"symmetric"`.
- `Nv::Int64`: subdivision order (polar bands per octant, or barycentric rows).

# Output Argument(s)
- `Nk::Int64`: total number of patches over the unit sphere.
- `lin_index::Array{Int64,3}`: `lin_index[u, v_or_i, w_or_j]` → linear patch index `k`, or `0` if not a valid patch.
- `uvw_of_k::Vector{NTuple{3,Int64}}`: inverse mapping `k → (u, v_or_i, w_or_j)`.
- `centroids::Vector{NTuple{3,Float64}}`: unit-vector spherical centroid of each patch.
- `areas::Vector{Float64}`: solid angle of each patch.
- `edges::Vector{Tuple{Int64,Int64,Float64,Float64}}`: unique edges `(k, kp, edge_len, d_geo)` with `k < kp`.

# Reference(s)
- Bienvenue et al. (2025), A Flexible, Moment-Preserving, and Monotone Discretization of the Multidimensional Angular Fokker-Planck Operator.
"""
function gn_patch_geometry(tiling::String, Nv::Int64)
    if tiling == "polar-anchored"
        return _gn_patch_geometry_polar_anchored(Nv)
    elseif tiling == "symmetric"
        return _gn_patch_geometry_symmetric(Nv)
    else
        error("Unknown tiling \"$tiling\"; expected \"polar-anchored\" or \"symmetric\".")
    end
end

# Octant sign tables (consistent with restricted_to_full_domain_matrix.jl)
const _GN_SX = (1, 1, 1, 1, -1, -1, -1, -1)
const _GN_SY = (1, 1, -1, -1, 1, 1, -1, -1)
const _GN_SZ = (1, -1, 1, -1, 1, -1, 1, -1)

# Octant adjacent in φ at the high-φ end of octant u, within the same hemisphere
# (same sx). Derived from φ_offset[u]: the next octant in increasing φ.
#   1 (φ∈[0,π/2])    → 3 (φ∈[π/2,π])
#   3 (φ∈[π/2,π])    → 4 (φ∈[π,3π/2])
#   4 (φ∈[π,3π/2])   → 2 (φ∈[3π/2,2π])
#   2 (φ∈[3π/2,2π])  → 1 (φ∈[0,π/2])   (wrap)
# Same pattern for the southern hemisphere (octants 5..8).
const _GN_NEXT_PHI = (3, 1, 4, 2, 7, 5, 8, 6)
const _GN_PREV_PHI = (2, 4, 1, 3, 6, 8, 5, 7)

# φ_offset[u] per octant (from restricted_to_full_domain_matrix.jl:184)
function _gn_phi_offset(u::Int64)
    sy = _GN_SY[u]; sz = _GN_SZ[u]
    return (π / 2) * (2 + (sy + 1) / 2 - (sz + 1) / 2 - (sy + 1) * (sz + 1) / 2)
end

# Polar-anchored band edges: μ_uv = sx * (1 - ((Nv+1-v)*(Nv+2-v)) / (Nv*(Nv+1)))
# when sx = +1, and μ_uv = -1 + ((v-1)*v) / (Nv*(Nv+1)) when sx = -1.
function _gn_mu_of_v(sx::Int, v::Int, Nv::Int64)
    denom = Float64(Nv * (Nv + 1))
    return sx == 1 ? (1.0 - ((Nv + 1 - v) * (Nv + 2 - v)) / denom) :
                     (-1.0 + ((v - 1) * v) / denom)
end

# Number of azimuthal slots per polar band (polar-anchored).
_gn_nw(u::Int64, v::Int64, Nv::Int64) = _GN_SX[u] == 1 ? (Nv + 1 - v) : v

# Equal-width μ-band [μ0, μ1] of band `v ∈ 1:Nv` in hemisphere octant `u` for the
# 1D Legendre (azimuthally symmetric) GN tiling. Nv bands per hemisphere:
#   u = 1 (sx = +1): μ ∈ [0, 1]   → band v = [(v-1)/Nv, v/Nv]
#   u = 5 (sx = -1): μ ∈ [-1, 0]  → band v = [-1+(v-1)/Nv, -1+v/Nv]
function _gn_legendre_band(u::Int64, v::Int64, Nv::Int64)
    if _GN_SX[u] == 1
        return ((v - 1) / Nv, v / Nv)
    else
        return ((v - 1) / Nv - 1.0, v / Nv - 1.0)
    end
end

# Antiderivative G(μ) = ∫ √(1-μ²) dμ = (μ √(1-μ²) + asin μ) / 2.
function _gn_G(μ::Float64)
    return 0.5 * (μ * sqrt(max(0.0, 1.0 - μ * μ)) + asin(clamp(μ, -1.0, 1.0)))
end

# Spherical centroid of the polar-anchored quadrilateral
# K = {(μ,φ) : μ ∈ [μ0, μ1], φ ∈ [φa, φb]}.
# Returns the unit-vector spherical centroid.
function _gn_centroid_polar_anchored(μ0::Float64, μ1::Float64, φa::Float64, φb::Float64)
    Ix = 0.5 * (μ1 * μ1 - μ0 * μ0) * (φb - φa)
    Gμ = _gn_G(μ1) - _gn_G(μ0)
    Iy = Gμ * (sin(φb) - sin(φa))
    Iz = Gμ * (cos(φa) - cos(φb))
    r = sqrt(Ix * Ix + Iy * Iy + Iz * Iz)
    if r == 0.0
        # Degenerate; fall back to mid-point.
        μm = 0.5 * (μ0 + μ1); φm = 0.5 * (φa + φb)
        s = sqrt(max(0.0, 1.0 - μm * μm))
        return (μm, s * cos(φm), s * sin(φm))
    end
    return (Ix / r, Iy / r, Iz / r)
end

# Great-circle distance between two unit vectors.
function _gn_d_geo(A::NTuple{3,Float64}, B::NTuple{3,Float64})
    d = A[1] * B[1] + A[2] * B[2] + A[3] * B[3]
    return acos(clamp(d, -1.0, 1.0))
end

# Overlap length of two real intervals [a_lo, a_hi] and [b_lo, b_hi].
_gn_overlap(a_lo::Float64, a_hi::Float64, b_lo::Float64, b_hi::Float64) =
    max(0.0, min(a_hi, b_hi) - max(a_lo, b_lo))

#----
# Polar-anchored tiling
#----
function _gn_patch_geometry_polar_anchored(Nv::Int64)
    # 1. Linear indexing
    lin_index = zeros(Int64, 8, Nv, Nv)
    uvw_of_k = Vector{NTuple{3,Int64}}()
    k = 0
    for u in 1:8, v in 1:Nv
        Nw = _gn_nw(u, v, Nv)
        for w in 1:Nw
            k += 1
            lin_index[u, v, w] = k
            push!(uvw_of_k, (u, v, w))
        end
    end
    Nk = k

    # 2. Centroids and areas
    centroids = Vector{NTuple{3,Float64}}(undef, Nk)
    areas = Vector{Float64}(undef, Nk)
    for k in 1:Nk
        u, v, w = uvw_of_k[k]
        sx = _GN_SX[u]
        μ_a = _gn_mu_of_v(sx, v, Nv)
        μ_b = _gn_mu_of_v(sx, v + 1, Nv)
        μ0, μ1 = min(μ_a, μ_b), max(μ_a, μ_b)
        Nw = _gn_nw(u, v, Nv)
        Δφ = (π / 2) / Nw
        φa = _gn_phi_offset(u) + (w - 1) * Δφ
        φb = φa + Δφ
        centroids[k] = _gn_centroid_polar_anchored(μ0, μ1, φa, φb)
        areas[k] = (μ1 - μ0) * Δφ
    end

    # 3. Adjacency
    edges = Vector{Tuple{Int64,Int64,Float64,Float64}}()

    # 3a. Intra-octant azimuthal neighbors: (u, v, w) ↔ (u, v, w+1).
    # Shared edge: meridian arc at φ = φ_a + w·Δφ, μ ∈ [μ0, μ1].
    for u in 1:8, v in 1:Nv
        Nw = _gn_nw(u, v, Nv)
        if Nw < 2 continue end
        sx = _GN_SX[u]
        μ_a = _gn_mu_of_v(sx, v, Nv)
        μ_b = _gn_mu_of_v(sx, v + 1, Nv)
        μ0, μ1 = min(μ_a, μ_b), max(μ_a, μ_b)
        # Meridian arc length (great-circle along the meridian) between (μ0, φ)
        # and (μ1, φ) is |acos μ0 - acos μ1|.
        edge_len = abs(acos(clamp(μ0, -1.0, 1.0)) - acos(clamp(μ1, -1.0, 1.0)))
        for w in 1:(Nw - 1)
            k = lin_index[u, v, w]
            kp = lin_index[u, v, w + 1]
            d_geo = _gn_d_geo(centroids[k], centroids[kp])
            push!(edges, (k, kp, edge_len, d_geo))
        end
    end

    # 3b. Cross-octant azimuthal neighbors at the φ-boundary between adjacent
    # octants (high-φ end w=Nw connects to w=1 of next_phi[u] octant).
    # To enumerate each edge only once, restrict to the half u → next_phi[u]
    # with u < next_phi[u].
    for u in 1:8
        un = _GN_NEXT_PHI[u]
        if u >= un continue end
        for v in 1:Nv
            Nw = _gn_nw(u, v, Nv)
            # Both octants share the same v-band μ range and the same
            # number of slots in their respective halves.
            sx = _GN_SX[u]
            μ_a = _gn_mu_of_v(sx, v, Nv)
            μ_b = _gn_mu_of_v(sx, v + 1, Nv)
            μ0, μ1 = min(μ_a, μ_b), max(μ_a, μ_b)
            edge_len = abs(acos(clamp(μ0, -1.0, 1.0)) - acos(clamp(μ1, -1.0, 1.0)))
            k = lin_index[u, v, Nw]
            kp = lin_index[un, v, 1]
            d_geo = _gn_d_geo(centroids[k], centroids[kp])
            push!(edges, (k, kp, edge_len, d_geo))
        end
    end

    # 3c. Intra-octant polar neighbors: (u, v, w) ↔ (u, v+1, w') for all w'
    # whose φ-extent in band v+1 overlaps that of w in band v.
    # Shared edge is on the latitude circle μ = μ_uv (between the two bands).
    for u in 1:8, v in 1:(Nv - 1)
        Nw_v = _gn_nw(u, v, Nv)
        Nw_v1 = _gn_nw(u, v + 1, Nv)
        Δφ_v = (π / 2) / Nw_v
        Δφ_v1 = (π / 2) / Nw_v1
        φ_off = _gn_phi_offset(u)
        sx = _GN_SX[u]
        μ_int = _gn_mu_of_v(sx, v + 1, Nv)  # shared band boundary
        sin_int = sqrt(max(0.0, 1.0 - μ_int * μ_int))
        for w in 1:Nw_v
            φa = φ_off + (w - 1) * Δφ_v
            φb = φa + Δφ_v
            for wp in 1:Nw_v1
                φap = φ_off + (wp - 1) * Δφ_v1
                φbp = φap + Δφ_v1
                ov = _gn_overlap(φa, φb, φap, φbp)
                if ov <= 0.0 continue end
                edge_len = sin_int * ov
                if edge_len == 0.0 continue end  # equator case μ_int=0 → handled in 3d
                k = lin_index[u, v, w]
                kp = lin_index[u, v + 1, wp]
                d_geo = _gn_d_geo(centroids[k], centroids[kp])
                push!(edges, (k, kp, edge_len, d_geo))
            end
        end
    end

    # 3d. Cross-octant polar neighbors at the equator μ = 0:
    # (u, Nv, 1) ↔ (u+4, 1, 1) for u in 1..4 (both have a single slot at the
    # equator, Nw=1, covering the same φ-quarter).
    # Note: at μ = 0, sin(arccos(0)) = 1 → edge_len = π/2 (the full φ-quarter).
    for u in 1:4
        u_south = u + 4
        Nw_n = _gn_nw(u, Nv, Nv)
        Nw_s = _gn_nw(u_south, 1, Nv)
        # Both should equal 1 by construction:
        #   sx=+1, v=Nv: Nv+1-Nv = 1; sx=-1, v=1: 1.
        # Be defensive in case of future change.
        for w in 1:Nw_n
            φa_n = _gn_phi_offset(u) + (w - 1) * (π / 2) / Nw_n
            φb_n = φa_n + (π / 2) / Nw_n
            for ws in 1:Nw_s
                φa_s = _gn_phi_offset(u_south) + (ws - 1) * (π / 2) / Nw_s
                φb_s = φa_s + (π / 2) / Nw_s
                ov = _gn_overlap(φa_n, φb_n, φa_s, φb_s)
                if ov <= 0.0 continue end
                edge_len = ov  # sin(π/2) = 1 at the equator
                k = lin_index[u, Nv, w]
                kp = lin_index[u_south, 1, ws]
                d_geo = _gn_d_geo(centroids[k], centroids[kp])
                push!(edges, (k, kp, edge_len, d_geo))
            end
        end
    end

    return Nk, lin_index, uvw_of_k, centroids, areas, edges
end

#----
# Symmetric tiling
#----
function _gn_patch_geometry_symmetric(Nv::Int64)
    # 1. Linear indexing
    Nw_max = 2 * Nv - 1
    lin_index = zeros(Int64, 8, Nv, Nw_max)
    uvw_of_k = Vector{NTuple{3,Int64}}()
    k = 0
    for u in 1:8, i in 1:Nv, j in 1:(2 * i - 1)
        k += 1
        lin_index[u, i, j] = k
        push!(uvw_of_k, (u, i, j))
    end
    Nk = k

    # 2. Centroids, areas, and vertex storage (per patch).
    # Use Vector{Vector{Float64}} for vertices to share with downstream
    # adjacency matching by raw signed-integer signatures.
    centroids = Vector{NTuple{3,Float64}}(undef, Nk)
    areas = Vector{Float64}(undef, Nk)
    # vertex_signatures[k] is a sorted tuple of 3 vertex signatures, each
    # encoded as (sx*a, sy*b, sz*c) signed integers before normalization.
    # This lets us identify shared edges across octants by tuple equality.
    vertex_sigs = Vector{NTuple{3,NTuple{3,Int64}}}(undef, Nk)
    for k in 1:Nk
        u, i, j = uvw_of_k[k]
        A, B, C = octant_subtriangle_vertices(u, i, j, Nv)
        # Spherical area (l'Huilier).
        areas[k] = spherical_triangle_area(A, B, C)
        # Approximate spherical centroid = normalize(A + B + C). Sufficient
        # for TPFA; the small displacement vs. the exact center-of-mass on
        # a curved triangle is O((triangle diameter)²).
        Sx = A[1] + B[1] + C[1]
        Sy = A[2] + B[2] + C[2]
        Sz = A[3] + B[3] + C[3]
        r = sqrt(Sx * Sx + Sy * Sy + Sz * Sz)
        centroids[k] = (Sx / r, Sy / r, Sz / r)
        vertex_sigs[k] = _gn_vertex_signatures_symmetric(u, i, j, Nv)
    end

    # 3. Adjacency: enumerate all unique edges by appariement of shared vertex
    # signatures. Each patch contributes 3 edges. The vertex signature is the
    # signed integer triplet before normalisation; two patches share an edge iff
    # they share two vertex signatures.
    edge_map = Dict{Tuple{NTuple{3,Int64},NTuple{3,Int64}},Vector{Int64}}()
    for k in 1:Nk
        sig = vertex_sigs[k]
        for a in 1:3, b in (a + 1):3
            v1 = sig[a]; v2 = sig[b]
            key = v1 < v2 ? (v1, v2) : (v2, v1)
            push!(get!(edge_map, key, Int64[]), k)
        end
    end

    edges = Vector{Tuple{Int64,Int64,Float64,Float64}}()
    for (key, ks) in edge_map
        if length(ks) != 2 continue end  # boundary; shouldn't happen on closed sphere
        k, kp = ks[1], ks[2]
        if k > kp k, kp = kp, k end
        # Compute the geodesic edge length: the two shared vertices, with
        # the signs of either patch (they coincide on the shared edge).
        # Extract vertices from the patch with k.
        u, i, j = uvw_of_k[k]
        A, B, C = octant_subtriangle_vertices(u, i, j, Nv)
        sig_k = vertex_sigs[k]
        verts = (A, B, C)
        # Pick the two vertices matching the edge key.
        v1_idx = findfirst(==(key[1]), sig_k)
        v2_idx = findfirst(==(key[2]), sig_k)
        if v1_idx === nothing || v2_idx === nothing
            # Edge crosses an axis plane; the shared vertices have at least one
            # zero component, so the sign of the patch doesn't affect them.
            # Recover them via the other patch instead.
            u2, i2, j2 = uvw_of_k[kp]
            A2, B2, C2 = octant_subtriangle_vertices(u2, i2, j2, Nv)
            sig_kp = vertex_sigs[kp]
            verts = (A2, B2, C2)
            v1_idx = findfirst(==(key[1]), sig_kp)
            v2_idx = findfirst(==(key[2]), sig_kp)
        end
        V1 = verts[v1_idx]
        V2 = verts[v2_idx]
        edge_len = _gn_d_geo((V1[1], V1[2], V1[3]), (V2[1], V2[2], V2[3]))
        d_geo = _gn_d_geo(centroids[k], centroids[kp])
        push!(edges, (k, kp, edge_len, d_geo))
    end

    return Nk, lin_index, uvw_of_k, centroids, areas, edges
end

# Build the three signed-integer vertex signatures of sub-triangle (u, i, j)
# of the symmetric tiling. Mirrors the raw-coordinate logic of
# `octant_subtriangle_vertices` but returns the signed integers before
# normalisation so that two triangles in different octants sharing an axis edge
# expose the same two vertex signatures (since the shared edge has at least
# one zero coordinate, the sign on that coordinate is irrelevant).
function _gn_vertex_signatures_symmetric(u::Int64, i::Int64, j::Int64, Nv::Int64)
    sx = _GN_SX[u]; sy = _GN_SY[u]; sz = _GN_SZ[u]
    if isodd(j)
        k = (j + 1) ÷ 2
        a1, b1, c1 = Nv - i + 1, i - k,     k - 1
        a2, b2, c2 = Nv - i,     i - k + 1, k - 1
        a3, b3, c3 = Nv - i,     i - k,     k
    else
        k = j ÷ 2
        a1, b1, c1 = Nv - i + 1, i - k,     k - 1
        a2, b2, c2 = Nv - i + 1, i - k - 1, k
        a3, b3, c3 = Nv - i,     i - k,     k
    end
    return (
        (sx * a1, sy * b1, sz * c1),
        (sx * a2, sy * b2, sz * c2),
        (sx * a3, sy * b3, sz * c3),
    )
end
