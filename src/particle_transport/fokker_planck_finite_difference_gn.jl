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
