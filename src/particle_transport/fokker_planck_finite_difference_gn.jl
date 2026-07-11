"""
    fokker_planck_finite_difference_gn(L, L_elem, Nv, Ndims, tiling, Mll, pl)

Build the angular Fokker-Planck scattering matrix for the GN solver using a
finite-volume (two-point flux approximation) discretization of the
Laplace-Beltrami operator on the patch graph induced by the GN finite-element
tiling. No Voronoi tessellation is needed - the (u, v, w) indexing of the
tiling already provides the connectivity.

The returned matrix `ℳ` has shape `(Np, Np)` with `Np = (L+1)²`, identical to
the Galerkin variant, and is consumed by `fokker_planck_source` exactly like
the SN finite-difference and the Galerkin GN matrices.

# Input Argument(s)
- `L::Int64`: global Legendre order (full-range basis size `Np = (L+1)²`).
- `L_elem::Int64`: per-patch local Legendre order (unused beyond `q=1`).
- `Nv::Int64`: subdivision order of the tiling.
- `Ndims::Int64`: geometry dimension (unused; carried for signature parity).
- `tiling::String`: `"polar-anchored"` or `"symmetric"`.
- `Mll::Array{Float64,5}`: patch-to-full-range moment matrix, shape `(Np, Nq, 8, Nv, Nw_max)`.
- `pl::Vector{Int64}`: Legendre order per full-range basis index (used only for diagnostics).

# Output Argument(s)
- `ℳ::Array{Float64,2}`: Fokker-Planck scattering matrix, shape `(Np, Np)`.
- `λ₀::Float64`: stabilisation factor added to `Σₜ` by the caller.

# Reference(s)
- Bienvenue et al. (2025), A Flexible, Moment-Preserving, and Monotone Discretization of the Multidimensional Angular Fokker-Planck Operator.
"""
function fokker_planck_finite_difference_gn(L::Int64,_L_elem::Int64,Nv::Int64,_Ndims::Int64,tiling::String,Mll::Array{Float64,5},_pl::Vector{Int64})

    # 1. Patch geometry and adjacency.
    Nk, _, uvw_of_k, _, areas, edges = gn_patch_geometry(tiling, Nv)
    sqω = sqrt.(areas)

    # 2. Graph Laplacian (TPFA) on patch averages.
    #    ℒ[k, kp] = γ / ω_k    for kp ∈ N(k)
    #    ℒ[k, k]  = -Σ γ / ω_k
    # γ_{k,kp} = edge_len / d_geo. Anti-symmetric in row-index weighting; the
    # operator is self-adjoint w.r.t. ⟨u,v⟩_ω = Σ_k ω_k u_k v_k.
    ℒ = zeros(Float64, Nk, Nk)
    for (k, kp, edge_len, d_geo) in edges
        if d_geo == 0.0 continue end
        γ = edge_len / d_geo
        ℒ[k,  kp] += γ / areas[k]
        ℒ[kp,  k] += γ / areas[kp]
        ℒ[k,   k] -= γ / areas[k]
        ℒ[kp, kp] -= γ / areas[kp]
    end

    # 3. Projection matrices in the Radiant full-range convention
    #    φ(Ω) = Σ_p (2ℓ_p+1)/(4π) · R_p · φ_p, so the moment-space matrix of an
    #    angular operator A is A_pq = (2ℓ_q+1)/(4π) · ⟨R_p, A · R_q⟩.
    #
    #    Round-trip through the patch-constant subspace:
    #    - patch-mean value of φ on K:   v_K = (Mc[p,K] · (2ℓ_p+1)/(4π)) · φ_p / √|K|
    #    - apply the value-space graph Laplacian ℒ on v
    #    - re-project to full moments:   u_p = √|K| · Mc[p,K] · (ℒv)_K
    #
    #    Hence with Mc[p,K] := Mll[p,1,K]:
    #      𝓓̃[p, k] = √ω_k · Mc[p, k]
    #      𝓜̃[k, q] = ((2ℓ_q+1)/(4π)) · Mc[q, k] / √ω_k
    Np = (L + 1)^2
    inv_4π = 1.0 / (4π)
    𝓜̃ = zeros(Float64, Nk, Np)
    𝓓̃ = zeros(Float64, Np, Nk)
    pl_local, _ = spherical_harmonics_indices(L)
    for k in 1:Nk
        u, v_or_i, w_or_j = uvw_of_k[k]
        s = sqω[k]
        if s == 0.0 continue end
        invs = 1.0 / s
        for p in 1:Np
            mp1 = Mll[p, 1, u, v_or_i, w_or_j]
            𝓜̃[k, p] = (2 * pl_local[p] + 1) * inv_4π * mp1 * invs
            𝓓̃[p, k] = s * mp1
        end
    end

    # 4. Assemble ℳ = 𝓓̃ · ℒ · 𝓜̃ in the full-range moment space.
    # The Laplace–Beltrami operator has eigenvalue -ℓ(ℓ+1) on R_{ℓ,m}; this
    # discretisation reproduces that diagonal asymptotically.
    ℳ = (𝓓̃ * (ℒ * 𝓜̃))

    # 5. Stabilisation: shift diagonal so that all diagonal entries are ≥ 0.
    λ₀ = 0.0
    @inbounds for p in 1:Np
        if ℳ[p, p] < -λ₀
            λ₀ = -ℳ[p, p]
        end
    end
    @inbounds for p in 1:Np
        ℳ[p, p] += λ₀
    end

    return ℳ, λ₀
end

"""
    fokker_planck_finite_difference_gn_legendre_1D(L, L_elem, Nv, Mll)

Build the angular Fokker-Planck scattering matrix `ℳ` for the 1D Legendre
(azimuthally symmetric) GN basis, the Legendre counterpart of
`fokker_planck_finite_difference_gn`. It uses a finite-volume (TPFA)
discretization of the azimuthally-symmetric Laplace-Beltrami operator
`d/dμ[(1-μ²) d/dμ]` on the chain of μ-bands.

The 2·Nv μ-bands form a chain ordered from μ = -1 to μ = +1 (octant `u = 5`
bands first, then octant `u = 1`). Adjacent bands `k, k+1` share the boundary
`μ_int`; the TPFA conductance is `γ = (1-μ_int²) / (μ_c[k+1]-μ_c[k])` and the
value-space measure is the band width `ω_k = Δμ_k`.

# Input Argument(s)
- `L::Int64` : global Legendre order (full-range basis size `Np = L + 1`).
- `_L_elem::Int64` : per-patch local Legendre order (unused beyond `q = 1`).
- `Nv::Int64` : number of μ-bands per hemisphere.
- `Mll::Array{Float64,5}` : patch-to-full-range moment matrix, shape `(Np, Nq, 8, Nv, 1)`.

# Output Argument(s)
- `ℳ::Array{Float64,2}` : Fokker-Planck scattering matrix, shape `(L+1, L+1)`.
- `λ₀::Float64` : stabilisation factor added to `Σₜ` by the caller.
"""
function fokker_planck_finite_difference_gn_legendre_1D(L::Int64, _L_elem::Int64, Nv::Int64, Mll::Array{Float64,5})

    # 1. μ-band chain ordered from μ = -1 to μ = +1.
    edges_of_u = Dict(u => gn_1D_band_edges(u, Nv, :legendre) for u in (1, 5))
    patches = NTuple{2,Int64}[]
    for v in 1:Nv push!(patches, (5, v)) end   # μ ∈ [-1, 0]
    for v in 1:Nv push!(patches, (1, v)) end   # μ ∈ [0, 1]
    Nk = length(patches)
    μ0k = zeros(Nk); μ1k = zeros(Nk); μc = zeros(Nk); ω = zeros(Nk)
    for k in 1:Nk
        u, v = patches[k]
        edges = edges_of_u[u]
        a, b = edges[v], edges[v + 1]
        μ0k[k] = a; μ1k[k] = b; μc[k] = 0.5 * (a + b); ω[k] = b - a
    end

    # 2. Graph Laplacian (TPFA) on patch-mean values, self-adjoint w.r.t. ⟨u,v⟩_ω.
    ℒ = zeros(Nk, Nk)
    for k in 1:(Nk - 1)
        μ_int = μ1k[k]               # shared band boundary (= μ0k[k+1])
        d_geo = μc[k + 1] - μc[k]
        if d_geo == 0.0 continue end
        γ = (1.0 - μ_int^2) / d_geo
        ℒ[k, k + 1]     += γ / ω[k]
        ℒ[k + 1, k]     += γ / ω[k + 1]
        ℒ[k, k]         -= γ / ω[k]
        ℒ[k + 1, k + 1] -= γ / ω[k + 1]
    end

    # 3. Projection matrices in the Radiant 1D Legendre convention:
    #    φ(μ) = Σ_p (2ℓ_p+1)/2 · P_p · φ_p, with Mc[p,k] = Mll[p,1,u,v,1].
    Np = L + 1
    𝓜̃ = zeros(Nk, Np)
    𝓓̃ = zeros(Np, Nk)
    for k in 1:Nk
        u, v = patches[k]
        s = sqrt(ω[k])
        if s == 0.0 continue end
        invs = 1.0 / s
        for p in 1:Np
            mp1 = Mll[p, 1, u, v, 1]
            𝓜̃[k, p] = (2 * (p - 1) + 1) / 2 * mp1 * invs
            𝓓̃[p, k] = s * mp1
        end
    end

    # 4. Assemble ℳ = 𝓓̃ · ℒ · 𝓜̃ in the full-range moment space.
    ℳ = 𝓓̃ * (ℒ * 𝓜̃)

    # 5. Stabilisation: shift the diagonal so all diagonal entries are ≥ 0.
    λ₀ = 0.0
    @inbounds for p in 1:Np
        if ℳ[p, p] < -λ₀ λ₀ = -ℳ[p, p] end
    end
    @inbounds for p in 1:Np
        ℳ[p, p] += λ₀
    end

    return ℳ, λ₀
end
