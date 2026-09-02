"""
    gn_inner_pass!(...)

Apply one in-group source-iteration pass for the generalized-harmonics (GN) solver, the GN
counterpart of `sn_inner_pass!`: build the source from the current flux `𝚽l`, transform to the
restricted-angle patches, sweep every octant/patch, transform back, and update `𝚽l` and the
incoming boundary fluxes (`𝚽x12⁻`, plus `𝚽y12⁻`/`𝚽z12⁻` in 2D/3D). This is the building block
every in-group acceleration scheme is built on (see `gn_one_speed`).

The pass is the affine map `T(z) = A·z + c` over the state `z = (𝚽l, boundary fluxes)`:
- `homogeneous = false` gives `T(z)` (all fixed sources active);
- `homogeneous = true` gives the linear part `A·z`, by dropping the fixed sources (`Qlout = 0`,
  the zeroed surface source `sources_q_zero`, zeroed incoming energy flux `𝚽E12`). The GN surface
  source enters through `sources_q`, not `Np_source`, hence the explicit zeroed clone.

`𝚽E12_temp` receives the outgoing energy flux (when `isCSD`); `Ql`, `*_q` and `*_ws`/`*_buf` are
scratch. Remaining arguments mirror `gn_one_speed`.
"""
function gn_inner_pass!(𝚽l,Qlout,Σt,Σs,mat,Ndims,Ns,Δs,Np,Nq,pl,Np_surf,𝒪,Nm,isFC,C,ω,isCSD,solver,𝚽E12,S⁻,S⁺,S,T,ℳ,𝒲,𝒩,boundary_conditions,Np_source,Nv,Mll,Mll_surf,Rpq,Mll_factored,tiling,is_SPH,fold,Ql,𝚽E12_temp,sources_q,sources_q_zero,𝚽x12⁻,𝚽x12⁺,𝚽y12⁻,𝚽y12⁺,𝚽z12⁻,𝚽z12⁺,Q_q,𝚽_q,𝚽E12_q,𝚽x12_q,𝚽y12_q,𝚽z12_q,𝒮_ws,Q_ws,𝚽_ws,𝚽x12_buf,𝚽y12_buf,𝚽z12_buf;homogeneous::Bool)

    # Octant sign patterns and patch indexing helper (recomputed cheaply; see gn_one_speed).
    # In 1D the azimuth collapses (both Legendre and spherical-harmonics bases): a
    # single slot per octant and only octants u ∈ {1, 5} carry μ-band patches.
    sx = [1,1,1,1,-1,-1,-1,-1]
    if (Ndims > 1) sy = [1,1,-1,-1,1,1,-1,-1] end
    if (Ndims > 2) sz = [1,-1,1,-1,1,-1,1,-1] end
    # In 1D the angular domain collapses to two half-spheres (u ∈ {1,5}, full
    # azimuth) for the Legendre basis and for the folded spherical-harmonics basis;
    # the unfolded spherical-harmonics basis keeps the full octant tiling. In 2D the
    # z-symmetry fold skips the even octants (four quadrants {1,3,5,7}).
    azim_collapsed = (Ndims == 1) && (!is_SPH || fold)
    z_fold_2D = (Ndims == 2) && (Nv == 1) && fold
    Nw_of(u, v) = azim_collapsed ? ((u == 1 || u == 5) ? 1 : 0) :
                  (z_fold_2D && iseven(u)) ? 0 :
                  ((tiling == "symmetric") ? (2*v - 1) : ((sx[u] == 1) ? (Nv + 1 - v) : v))

    # Calculation of the Legendre components of the source (in-scattering). For the homogeneous
    # operator A·z the out-of-group source Qlout is dropped; the scattering/Fokker-Planck builders
    # are linear in 𝚽l and stay active in both cases.
    if homogeneous Ql .= 0.0 else Ql .= Qlout end
    if solver ∉ [4,5,6] Ql = scattering_source(Ql,𝚽l,Σs,mat,Np,pl,Nm[5],Ns) end

    # Finite element treatment of the angular Fokker-Planck term
    if solver ∈ [2,4] Ql = fokker_planck_source(Np,Nm[5],T,𝚽l,Ql,Ns,mat,ℳ) end

    # Fixed-source switches: drop the surface source and incoming energy flux for the homogeneous
    # operator (the surface source enters the sweeps through sources_q, not Np_source).
    sources_q_eff = homogeneous ? sources_q_zero : sources_q
    𝚽E12_eff = (isCSD && homogeneous) ? zero(𝚽E12) : 𝚽E12

    #----
    # Loop over all discrete ordinates
    #----
    𝚽l .= 0
    𝚽E12_temp .= 0
    if Ndims == 1
        fill!(Q_q, 0.0)
        fill!(𝚽_q, 0.0)
        fill!(𝚽E12_q, 0.0)
        fill!(𝚽x12_q, 0.0)

        NS = Nm[5] * Ns[1]
        NSE = Nm[4] * Ns[1]
        NSx = Nm[1]
        Ql_mat = reshape(Ql, Np, NS)

        # Transformation of full-range fluxes to restricted-angle fluxes
        @views for u in 1:8, v in 1:Nv
            Nw = Nw_of(u, v)
            for w in 1:Nw
                Mt = transpose(Mll_factored[:,:,u,v,w])
                mul!(reshape(Q_q[:,:,:,u,v,w], Nq, NS), Mt, Ql_mat)
                if isCSD
                    mul!(reshape(𝚽E12_q[:,:,:,u,v,w], Nq, NSE), Mt, reshape(𝚽E12_eff, Np, NSE))
                end
                for ib in 1:2
                    mul!(reshape(𝚽x12_q[:,:,ib,u,v,w], Nq, NSx),transpose(Mll_surf[:,:,u,v,w,ib,1]),reshape(𝚽x12⁻[:,:,ib], Np_surf, NSx))
                end
            end
        end
        # Computation of the restricted-angle fluxes by sweeping through the spatial grid
        @views for u in 1:8, v in 1:Nv
            Nw = Nw_of(u, v)
            for w in 1:Nw
                gn_sweep_1D!(𝚽_q[:,:,:,u,v,w],
                            𝚽E12_q[:,:,:,u,v,w],
                            𝚽x12_q[:,:,:,u,v,w],
                            sx[u],Σt,mat[:,1,1],Ns[1],Δs[1],
                            Q_q[:,:,:,u,v,w],
                            Nq,Np_source,𝒪,Nm,C,ω,
                            sources_q_eff[:,:,u,v,w],
                            S⁻,S⁺,S,𝒲,isFC,isCSD,
                            𝒩[:,:,1,u,v,w],
                            𝒮_ws,Q_ws,𝚽_ws,
                            𝚽x12_buf)
            end
        end
        # Transformation of restricted-angle fluxes to full-range fluxes (BLAS gemm, accumulate)
        𝚽l_mat = reshape(𝚽l, Np, NS)
        @views for u in 1:8, v in 1:Nv
            Nw = Nw_of(u, v)
            for w in 1:Nw
                M = Mll[:,:,u,v,w]
                mul!(𝚽l_mat, M, reshape(𝚽_q[:,:,:,u,v,w], Nq, NS), 1.0, 1.0)
                if isCSD
                    mul!(reshape(𝚽E12_temp, Np, NSE), M, reshape(𝚽E12_q[:,:,:,u,v,w], Nq, NSE), 1.0, 1.0)
                end
                for ib in 1:2
                    mul!(reshape(𝚽x12⁺[:,:,ib], Np_surf, NSx),Mll_surf[:,:,u,v,w,ib,2],reshape(𝚽x12_q[:,:,ib,u,v,w], Nq, NSx),1.0, 1.0)
                end
            end
        end
        # Boundary conditions treatment
        𝚽x12⁻ .= 0.0
        for ib in range(1,2)
            if boundary_conditions[ib] != 0
                if boundary_conditions[ib] == 1 # Reflective boundary condition
                    # Batched reflection with the stride-1 basis index p innermost (q-summation
                    # order preserved ⇒ bit-identical to the original p-outer form).
                    @inbounds for is in range(1,Nm[1]), q in range(1,Np_surf)
                        s = 𝚽x12⁺[q,is,ib]
                        @simd for p in range(1,Np_surf)
                            𝚽x12⁻[p,is,ib] += Rpq[p,q,ib] * s
                        end
                    end
                elseif boundary_conditions[ib] == 2 # Periodic boundary condition
                    for p in range(1,Np_surf), is in range(1,Nm[1])
                        if ib == 1
                            𝚽x12⁻[p,is,1] += 𝚽x12⁺[p,is,2]
                        else
                            𝚽x12⁻[p,is,2] += 𝚽x12⁺[p,is,1]
                        end
                    end
                else
                    error("Invalid boundary condition type.")
                end
            end
        end
        𝚽x12⁺ .= 0.0

    elseif Ndims == 2

        fill!(Q_q, 0.0)
        fill!(𝚽_q, 0.0)
        fill!(𝚽E12_q, 0.0)
        fill!(𝚽x12_q, 0.0)
        fill!(𝚽y12_q, 0.0)

        NS = Nm[5] * Ns[1] * Ns[2]
        NSE = Nm[4] * Ns[1] * Ns[2]
        NSx = Nm[1] * Ns[2]
        NSy = Nm[2] * Ns[1]
        Ql_mat = reshape(Ql, Np, NS)

        # Transformation of full-range fluxes to restricted-angle fluxes
        @views for u in 1:8, v in 1:Nv
            Nw = Nw_of(u, v)
            for w in 1:Nw
                Mt = transpose(Mll_factored[:,:,u,v,w])
                mul!(reshape(Q_q[:,:,:,:,u,v,w], Nq, NS), Mt, Ql_mat)
                if isCSD
                    mul!(reshape(𝚽E12_q[:,:,:,:,u,v,w], Nq, NSE), Mt, reshape(𝚽E12_eff, Np, NSE))
                end
                for ib in 1:2
                    mul!(reshape(𝚽x12_q[:,:,:,ib,u,v,w], Nq, NSx),transpose(Mll_surf[:,:,u,v,w,ib,1]),reshape(𝚽x12⁻[:,:,:,ib], Np_surf, NSx))
                    mul!(reshape(𝚽y12_q[:,:,:,ib,u,v,w], Nq, NSy),transpose(Mll_surf[:,:,u,v,w,ib+2,1]),reshape(𝚽y12⁻[:,:,:,ib], Np_surf, NSy))
                end
            end
        end
        # Computation of the restricted-angle fluxes by sweeping through the spatial grid
        @views for u in 1:8, v in 1:Nv
            Nw = Nw_of(u, v)
            for w in 1:Nw
                gn_sweep_2D!(𝚽_q[:,:,:,:,u,v,w],
                            𝚽E12_q[:,:,:,:,u,v,w],
                            𝚽x12_q[:,:,:,:,u,v,w],
                            𝚽y12_q[:,:,:,:,u,v,w],
                            sx[u],sy[u],Σt,mat[:,:,1],Ns[1],Ns[2],Δs[1],Δs[2],
                            Q_q[:,:,:,:,u,v,w],
                            Nq,Np_source,𝒪,Nm,C,ω,
                            sources_q_eff[:,:,u,v,w],
                            S⁻,S⁺,S,𝒲,isFC,isCSD,
                            𝒩[:,:,1,u,v,w],𝒩[:,:,2,u,v,w],
                            𝒮_ws,Q_ws,𝚽_ws,
                            𝚽x12_buf,𝚽y12_buf)
            end
        end
        # Transformation of restricted-angle fluxes to full-range fluxes (BLAS gemm, accumulate)
        𝚽l_mat = reshape(𝚽l, Np, NS)
        @views for u in 1:8, v in 1:Nv
            Nw = Nw_of(u, v)
            for w in 1:Nw
                M = Mll[:,:,u,v,w]
                mul!(𝚽l_mat, M, reshape(𝚽_q[:,:,:,:,u,v,w], Nq, NS), 1.0, 1.0)
                if isCSD
                    mul!(reshape(𝚽E12_temp, Np, NSE), M, reshape(𝚽E12_q[:,:,:,:,u,v,w], Nq, NSE), 1.0, 1.0)
                end
                for ib in 1:2
                    mul!(reshape(𝚽x12⁺[:,:,:,ib], Np_surf, NSx),Mll_surf[:,:,u,v,w,ib,2],reshape(𝚽x12_q[:,:,:,ib,u,v,w], Nq, NSx),1.0, 1.0)
                    mul!(reshape(𝚽y12⁺[:,:,:,ib], Np_surf, NSy),Mll_surf[:,:,u,v,w,ib+2,2],reshape(𝚽y12_q[:,:,:,ib,u,v,w], Nq, NSy),1.0, 1.0)
                end
            end
        end
        # Boundary conditions treatment
        𝚽x12⁻ .= 0.0
        𝚽y12⁻ .= 0.0
        for ib in range(1,2)
            # X-axis boundary conditions
            if boundary_conditions[ib] != 0
                if boundary_conditions[ib] == 1 # Reflective boundary condition
                    # Batched reflection with the stride-1 basis index p innermost (q-summation
                    # order preserved ⇒ bit-identical to the original p-outer form).
                    @inbounds for iy in range(1,Ns[2]), is in range(1,Nm[1]), q in range(1,Np_surf)
                        s = 𝚽x12⁺[q,is,iy,ib]
                        @simd for p in range(1,Np_surf)
                            𝚽x12⁻[p,is,iy,ib] += Rpq[p,q,ib] * s
                        end
                    end
                elseif boundary_conditions[ib] == 2 # Periodic boundary condition
                    for p in range(1,Np_surf), is in range(1,Nm[1]), iy in range(1,Ns[2])
                        if ib == 1
                            𝚽x12⁻[p,is,iy,1] += 𝚽x12⁺[p,is,iy,2]
                        else
                            𝚽x12⁻[p,is,iy,2] += 𝚽x12⁺[p,is,iy,1]
                        end
                    end
                else
                    error("Invalid boundary condition type.")
                end
            end
            # Y-axis boundary conditions
            if boundary_conditions[ib+2] != 0
                if boundary_conditions[ib+2] == 1 # Reflective boundary condition
                    @inbounds for ix in range(1,Ns[1]), is in range(1,Nm[2]), q in range(1,Np_surf)
                        s = 𝚽y12⁺[q,is,ix,ib]
                        @simd for p in range(1,Np_surf)
                            𝚽y12⁻[p,is,ix,ib] += Rpq[p,q,ib+2] * s
                        end
                    end
                elseif boundary_conditions[ib+2] == 2 # Periodic boundary condition
                    for p in range(1,Np_surf), is in range(1,Nm[2]), ix in range(1,Ns[1])
                        if ib == 1
                            𝚽y12⁻[p,is,ix,1] += 𝚽y12⁺[p,is,ix,2]
                        else
                            𝚽y12⁻[p,is,ix,2] += 𝚽y12⁺[p,is,ix,1]
                        end
                    end
                else
                    error("Invalid boundary condition type.")
                end
            end
        end
        𝚽x12⁺ .= 0.0
        𝚽y12⁺ .= 0.0

    elseif Ndims == 3

        fill!(Q_q, 0.0)
        fill!(𝚽_q, 0.0)
        fill!(𝚽E12_q, 0.0)
        fill!(𝚽x12_q, 0.0)
        fill!(𝚽y12_q, 0.0)
        fill!(𝚽z12_q, 0.0)

        NS = Nm[5] * Ns[1] * Ns[2] * Ns[3]
        NSE = Nm[4] * Ns[1] * Ns[2] * Ns[3]
        NSx = Nm[1] * Ns[2] * Ns[3]
        NSy = Nm[2] * Ns[1] * Ns[3]
        NSz = Nm[3] * Ns[1] * Ns[2]
        Ql_mat = reshape(Ql, Np, NS)

        # Transformation of full-range fluxes to restricted-angle fluxes (BLAS gemm)
        @views for u in 1:8, v in 1:Nv
            Nw = Nw_of(u, v)
            for w in 1:Nw
                Mt = transpose(Mll_factored[:,:,u,v,w])
                mul!(reshape(Q_q[:,:,:,:,:,u,v,w], Nq, NS), Mt, Ql_mat)
                if isCSD
                    mul!(reshape(𝚽E12_q[:,:,:,:,:,u,v,w], Nq, NSE), Mt, reshape(𝚽E12_eff, Np, NSE))
                end
                for ib in 1:2
                    mul!(reshape(𝚽x12_q[:,:,:,:,ib,u,v,w], Nq, NSx),transpose(Mll_surf[:,:,u,v,w,ib,1]),reshape(𝚽x12⁻[:,:,:,:,ib], Np_surf, NSx))
                    mul!(reshape(𝚽y12_q[:,:,:,:,ib,u,v,w], Nq, NSy),transpose(Mll_surf[:,:,u,v,w,ib+2,1]),reshape(𝚽y12⁻[:,:,:,:,ib], Np_surf, NSy))
                    mul!(reshape(𝚽z12_q[:,:,:,:,ib,u,v,w], Nq, NSz),transpose(Mll_surf[:,:,u,v,w,ib+4,1]),reshape(𝚽z12⁻[:,:,:,:,ib], Np_surf, NSz))
                end
            end
        end
        # Computation of the restricted-angle fluxes by sweeping through the spatial grid
        @views for u in 1:8, v in 1:Nv
            Nw = Nw_of(u, v)
            for w in 1:Nw
                gn_sweep_3D!(𝚽_q[:,:,:,:,:,u,v,w],
                            𝚽E12_q[:,:,:,:,:,u,v,w],
                            𝚽x12_q[:,:,:,:,:,u,v,w],
                            𝚽y12_q[:,:,:,:,:,u,v,w],
                            𝚽z12_q[:,:,:,:,:,u,v,w],
                            sx[u],sy[u],sz[u],Σt,mat,Ns[1],Ns[2],Ns[3],Δs[1],Δs[2],Δs[3],
                            Q_q[:,:,:,:,:,u,v,w],
                            Nq,Np_source,𝒪,Nm,C,ω,
                            sources_q_eff[:,:,u,v,w],
                            S⁻,S⁺,S,𝒲,isFC,isCSD,
                            𝒩[:,:,1,u,v,w],𝒩[:,:,2,u,v,w],𝒩[:,:,3,u,v,w],
                            𝒮_ws,Q_ws,𝚽_ws,
                            𝚽x12_buf,𝚽y12_buf,𝚽z12_buf)
            end
        end
        # Transformation of restricted-angle fluxes to full-range fluxes (BLAS gemm, accumulate)
        𝚽l_mat = reshape(𝚽l, Np, NS)
        @views for u in 1:8, v in 1:Nv
            Nw = Nw_of(u, v)
            for w in 1:Nw
                M = Mll[:,:,u,v,w]
                mul!(𝚽l_mat, M, reshape(𝚽_q[:,:,:,:,:,u,v,w], Nq, NS), 1.0, 1.0)
                if isCSD
                    mul!(reshape(𝚽E12_temp, Np, NSE), M, reshape(𝚽E12_q[:,:,:,:,:,u,v,w], Nq, NSE), 1.0, 1.0)
                end
                for ib in 1:2
                    mul!(reshape(𝚽x12⁺[:,:,:,:,ib], Np_surf, NSx),Mll_surf[:,:,u,v,w,ib,2],reshape(𝚽x12_q[:,:,:,:,ib,u,v,w], Nq, NSx),1.0, 1.0)
                    mul!(reshape(𝚽y12⁺[:,:,:,:,ib], Np_surf, NSy),Mll_surf[:,:,u,v,w,ib+2,2],reshape(𝚽y12_q[:,:,:,:,ib,u,v,w], Nq, NSy),1.0, 1.0)
                    mul!(reshape(𝚽z12⁺[:,:,:,:,ib], Np_surf, NSz),Mll_surf[:,:,u,v,w,ib+4,2],reshape(𝚽z12_q[:,:,:,:,ib,u,v,w], Nq, NSz),1.0, 1.0)
                end
            end
        end
        # Boundary conditions treatment
        𝚽x12⁻ .= 0.0
        𝚽y12⁻ .= 0.0
        𝚽z12⁻ .= 0.0
        for ib in range(1,2)
            # X-axis boundary conditions
            if boundary_conditions[ib] != 0
                if boundary_conditions[ib] == 1 # Reflective boundary condition
                    # Batched reflection 𝚽x12⁻[:,col] += Rpq[:,:,ib]·𝚽x12⁺[:,col]; loop with the
                    # stride-1 basis index p innermost and q just outside (summation order over q
                    # preserved ⇒ bit-identical to the original p-outer form).
                    @inbounds for iz in range(1,Ns[3]), iy in range(1,Ns[2]), is in range(1,Nm[1]), q in range(1,Np_surf)
                        s = 𝚽x12⁺[q,is,iy,iz,ib]
                        @simd for p in range(1,Np_surf)
                            𝚽x12⁻[p,is,iy,iz,ib] += Rpq[p,q,ib] * s
                        end
                    end
                elseif boundary_conditions[ib] == 2 # Periodic boundary condition
                    for p in range(1,Np_surf), is in range(1,Nm[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3])
                        if ib == 1
                            𝚽x12⁻[p,is,iy,iz,1] += 𝚽x12⁺[p,is,iy,iz,2]
                        else
                            𝚽x12⁻[p,is,iy,iz,2] += 𝚽x12⁺[p,is,iy,iz,1]
                        end
                    end
                else
                    error("Invalid boundary condition type.")
                end
            end
            # Y-axis boundary conditions
            if boundary_conditions[ib+2] != 0
                if boundary_conditions[ib+2] == 1 # Reflective boundary condition
                    @inbounds for iz in range(1,Ns[3]), ix in range(1,Ns[1]), is in range(1,Nm[2]), q in range(1,Np_surf)
                        s = 𝚽y12⁺[q,is,ix,iz,ib]
                        @simd for p in range(1,Np_surf)
                            𝚽y12⁻[p,is,ix,iz,ib] += Rpq[p,q,ib+2] * s
                        end
                    end
                elseif boundary_conditions[ib+2] == 2 # Periodic boundary condition
                    for p in range(1,Np_surf), is in range(1,Nm[2]), ix in range(1,Ns[1]), iz in range(1,Ns[3])
                        if ib == 1
                            𝚽y12⁻[p,is,ix,iz,1] += 𝚽y12⁺[p,is,ix,iz,2]
                        else
                            𝚽y12⁻[p,is,ix,iz,2] += 𝚽y12⁺[p,is,ix,iz,1]
                        end
                    end
                else
                    error("Invalid boundary condition type.")
                end
            end
            # Z-axis boundary conditions
            if boundary_conditions[ib+4] != 0
                if boundary_conditions[ib+4] == 1 # Reflective boundary condition
                    @inbounds for iy in range(1,Ns[2]), ix in range(1,Ns[1]), is in range(1,Nm[3]), q in range(1,Np_surf)
                        s = 𝚽z12⁺[q,is,ix,iy,ib]
                        @simd for p in range(1,Np_surf)
                            𝚽z12⁻[p,is,ix,iy,ib] += Rpq[p,q,ib+4] * s
                        end
                    end
                elseif boundary_conditions[ib+4] == 2 # Periodic boundary condition
                    for p in range(1,Np_surf), is in range(1,Nm[3]), ix in range(1,Ns[1]), iy in range(1,Ns[2])
                        if ib == 1
                            𝚽z12⁻[p,is,ix,iy,1] += 𝚽z12⁺[p,is,ix,iy,2]
                        else
                            𝚽z12⁻[p,is,ix,iy,2] += 𝚽z12⁺[p,is,ix,iy,1]
                        end
                    end
                else
                    error("Invalid boundary condition type.")
                end
            end
        end
        𝚽x12⁺ .= 0.0
        𝚽y12⁺ .= 0.0
        𝚽z12⁺ .= 0.0

    else
        error("Geometry dimension is either 1D, 2D or 3D.")
    end

    return 𝚽l
end

"""
    gn_inner_pass_fast!(...)

Optimized counterpart of `gn_inner_pass!`, restricted to 3D — the GN counterpart of one
in-group source-iteration pass. Same three phases as the reference (transform the source
onto the angular patches, sweep every patch, transform back and apply the boundary
conditions) and the same affine map `T(z) = A·z + c`, with `homogeneous = true` selecting
the linear part.

What differs:

- the patches are enumerated on a flat list, so the `(8,Nv,Nw_max)` arrays of the
  reference no longer carry the `8·Nv·Nw_max − 4·Nv·(Nv+1)` slots that never hold a patch
  (24 of 72 at `Nv = 3`);
- each patch sweeps with its own `GNFastWorkspace`, and `gn_3D_*_fast!` stays out of
  LAPACK: `lu!`/`ldiv!` go through OpenBLAS' globally-locked workspace allocator, which
  costs more per call than the arithmetic on systems this small;
- the per-pass `fill!` of the full-grid patch arrays is gone: every element of `Q_q`,
  `𝚽_q`, `𝚽E12_q` and the boundary arrays is overwritten before it is read — by the
  forward `mul!` for the first three, by the sweep for `𝚽_q`;
- the boundary sources reach the sweep as concrete arrays (see `gn_sweep_3D_fast!`).

`wss` holds one workspace per patch, `ctx` the precomputed context (`gn_fast_context.jl`).
The remaining arguments mirror `gn_inner_pass!`.
"""
function gn_inner_pass_fast!(𝚽l,Qlout,Σs,mat,Ndims,Ns,Np,Nq,pl,Np_surf,Nm,isCSD,solver,𝚽E12,T,ℳ,boundary_conditions,Mll,Mll_surf,Rpq,Mll_factored,Mll_all,Mfact_all,ctx,patch_uvw,wss,Ql,𝚽E12_temp,srcx,srcy,srcz,𝚽x12⁻,𝚽x12⁺,𝚽y12⁻,𝚽y12⁺,𝚽z12⁻,𝚽z12⁺,Q_q,𝚽_q,𝚽E12_q,𝚽x12_q,𝚽y12_q,𝚽z12_q,E12_state,need_boundary_flux::Bool,tbufs;homogeneous::Bool,reconstruct::Bool=false)

    Npatch = length(patch_uvw)

    # Calculation of the Legendre components of the source (in-scattering). For the homogeneous
    # operator A·z the out-of-group source Qlout is dropped; the scattering/Fokker-Planck builders
    # are linear in 𝚽l and stay active in both cases.
    if homogeneous Ql .= 0.0 else Ql .= Qlout end
    if solver ∉ [4,5,6] Ql = scattering_source_fast(Ql,𝚽l,Σs,mat,Np,pl,Nm[5],Ns) end

    # Finite element treatment of the angular Fokker-Planck term
    if solver ∈ [2,4] Ql = fokker_planck_source_fast(Np,Nm[5],T,𝚽l,Ql,Ns,mat,ℳ) end

    # The homogeneous operator drops the incoming energy flux; it is switched below by
    # zeroing 𝚽E12_q directly rather than by transforming a zeroed clone of 𝚽E12. (The
    # surface source is switched by the caller, which passes zeroed srcx/srcy/srcz.)
    𝚽E12_eff = 𝚽E12

    NS  = Nm[5] * Ns[1] * Ns[2] * Ns[3]
    # Without the CSD term there is no energy-axis flux at all: the caller allocates
    # 𝚽E12_q with an empty moment axis, and 𝚽E12/𝚽E12_temp are scalar placeholders.
    NSE = isCSD ? Nm[4] * Ns[1] * Ns[2] * Ns[3] : 0
    NSx = Nm[1] * Ns[2] * Ns[3]
    NSy = Nm[2] * Ns[1] * Ns[3]
    NSz = Nm[3] * Ns[1] * Ns[2]
    Nmf = [Nm[1],Nm[2],Nm[3]]
    has_z = Ndims ≥ 3
    has_y = Ndims ≥ 2

    𝚽l .= 0
    𝚽E12_temp .= 0

    # Flat (angular, spatial, patch) views of the patch arrays. The leading axes of each
    # array are contiguous, so these reshapes are free and every `view(·,:,range,k)` below
    # is a strided matrix BLAS can take directly.
    Ql_mat      = reshape(Ql, Np, NS)
    𝚽l_mat      = reshape(𝚽l, Np, NS)
    Q_q_mat     = reshape(Q_q, Nq, NS, Npatch)
    𝚽_q_mat     = reshape(𝚽_q, Nq, NS, Npatch)
    𝚽E12_q_mat  = reshape(𝚽E12_q, Nq, NSE, Npatch)
    𝚽E12_mat    = isCSD ? reshape(𝚽E12_eff , Np, NSE) : Matrix{Float64}(undef,Np,0)
    𝚽E12_t_mat  = isCSD ? reshape(𝚽E12_temp, Np, NSE) : Matrix{Float64}(undef,Np,0)

    # Transformation of full-range fluxes to restricted-angle fluxes, cut along the spatial
    # axis: the ranges are independent (each writes its own columns) and the full-range
    # slice stays in cache across the patch loop.
    #
    # With one basis function per patch (`Nq == 1`, i.e. `legendre_order_local = 0`, the
    # usual GN setting) the whole patch loop collapses into a single `gemm`. The patch axis
    # of `Q_q` is then its slowest axis and the spatial axis its fastest, so `Q_qᵀ = Qlᵀ·M`
    # lands exactly in the array's own layout — one compute-bound `gemm` reading `Ql` once,
    # instead of `Npatch` memory-bound `gemv` passes over it.
    # With `Nq > 1` the same collapse is still possible, but the single `gemm` then produces a
    # patch-major block that has to be scattered back into the sweep's layout. Even paying that
    # copy it beats the per-patch form by 2.59× (measured at Nq = 3, 72 patches): those are
    # `gemm`s of m = 3, where BLAS spends more time packing than multiplying.
    chunks = gn_fast_chunks(NS,1)
    nth = length(tbufs)
    blocks = gn_fast_chunks(NS,max(1,cld(NS,GN_FAST_BLOCK)))
    Q_q_v = vec(Q_q)
    if Nq == 1
        Q_q_2d = reshape(Q_q, NS, Npatch)
        for ci in 1:length(chunks)
            rng = chunks[ci]
            @views mul!(Q_q_2d[rng,:], transpose(Ql_mat[:,rng]), Mfact_all)
        end
    else
        for t in 1:nth
            buf = tbufs[t]
            for bi in t:nth:length(blocks)
                rng = blocks[bi]
                @views mul!(buf[:,1:length(rng)], transpose(Mfact_all), Ql_mat[:,rng])
                gn_fast_scatter!(Q_q_v,buf,rng,Nq,Npatch,NS)
            end
        end
    end
    # Incoming energy flux. 𝚽E12 is inherited from the previous energy group and is not part
    # of the fixed-point state, so its transform onto the patches is the same at every pass —
    # and the sweep now leaves 𝚽E12_q pristine (see the `do_E` guard in gn_3D_BFP_fast!).
    # It therefore only has to be rebuilt when the caller switches between the affine map and
    # its homogeneous part, which happens at most a couple of times per group.
    if isCSD
        want = homogeneous ? 2 : 1
        if E12_state[] != want
            if homogeneous
                fill!(𝚽E12_q, 0.0)
            else
                chunksE = gn_fast_chunks(NSE,1)
                if Nq == 1
                    𝚽E12_q_2d = reshape(𝚽E12_q, NSE, Npatch)
                    for ci in 1:length(chunksE)
                        rng = chunksE[ci]
                        @views mul!(𝚽E12_q_2d[rng,:], transpose(𝚽E12_mat[:,rng]), Mfact_all)
                    end
                else
                    blocksE = gn_fast_chunks(NSE,max(1,cld(NSE,GN_FAST_BLOCK)))
                    𝚽E12_q_v = vec(𝚽E12_q)
                    for t in 1:nth
                        buf = tbufs[t]
                        for bi in t:nth:length(blocksE)
                            rng = blocksE[bi]
                            @views mul!(buf[:,1:length(rng)], transpose(Mfact_all), 𝚽E12_mat[:,rng])
                            gn_fast_scatter!(𝚽E12_q_v,buf,rng,Nq,Npatch,NSE)
                        end
                    end
                end
            end
            E12_state[] = want
        end
    end
    # Incoming boundary fluxes. With every face void these are identically zero for the whole
    # solve — only reflective and periodic conditions feed them — so the transforms would
    # move nothing but zeros, and 𝚽*12_q (allocated zero, and re-zeroed by the sweep) already
    # holds the right values.
    @views for k in (need_boundary_flux ? (1:Npatch) : (1:0))
        (u,v,w) = patch_uvw[k]
        for ib in 1:2
            mul!(reshape(𝚽x12_q[:,:,:,:,ib,k], Nq, NSx),transpose(Mll_surf[:,:,u,v,w,ib,1]),reshape(𝚽x12⁻[:,:,:,:,ib], Np_surf, NSx))
            if has_y
                mul!(reshape(𝚽y12_q[:,:,:,:,ib,k], Nq, NSy),transpose(Mll_surf[:,:,u,v,w,ib+2,1]),reshape(𝚽y12⁻[:,:,:,:,ib], Np_surf, NSy))
            end
            if has_z
                mul!(reshape(𝚽z12_q[:,:,:,:,ib,k], Nq, NSz),transpose(Mll_surf[:,:,u,v,w,ib+4,1]),reshape(𝚽z12⁻[:,:,:,:,ib], Np_surf, NSz))
            end
        end
    end

    # Computation of the restricted-angle fluxes by sweeping through the spatial grid.
    # The patches are independent within a pass — they only meet in the reduction into 𝚽l
    # below — so the sweeps run concurrently, each on its own workspace.
    # Flat views of every patch array: the sweep addresses them linearly from a per-voxel
    # offset, which is what keeps it from building a SubArray per voxel. `vec` on an Array
    # shares storage, so these cost nothing.
    𝚽_q_v = vec(𝚽_q); Q_q_v = vec(Q_q); 𝚽E12_q_v = vec(𝚽E12_q)
    𝚽x12_q_v = vec(𝚽x12_q); 𝚽y12_q_v = vec(𝚽y12_q); 𝚽z12_q_v = vec(𝚽z12_q)
    srcx_v = vec(srcx); srcy_v = vec(srcy); srcz_v = vec(srcz)
    NmEf = isCSD ? Nm[4] : 0

    for k in 1:Npatch
        if has_z
            gn_sweep_3D_fast!(𝚽_q_v,𝚽E12_q_v,𝚽x12_q_v,𝚽y12_q_v,𝚽z12_q_v,k,
                              mat,Ns[1],Ns[2],Ns[3],Q_q_v,
                              srcx_v,srcy_v,srcz_v,
                              ctx.mom,ctx.cells,ctx.patch[k],wss[k],Nmf,NmEf,isCSD,reconstruct,need_boundary_flux)
        elseif has_y
            gn_sweep_2D_fast!(𝚽_q_v,𝚽E12_q_v,𝚽x12_q_v,𝚽y12_q_v,k,
                              mat,Ns[1],Ns[2],Q_q_v,
                              srcx_v,srcy_v,
                              ctx.mom,ctx.cells,ctx.patch[k],wss[k],Nmf,NmEf,isCSD,reconstruct,need_boundary_flux)
        else
            gn_sweep_1D_fast!(𝚽_q_v,𝚽E12_q_v,𝚽x12_q_v,k,
                              mat,Ns[1],Q_q_v,srcx_v,
                              ctx.mom,ctx.cells,ctx.patch[k],wss[k],Nmf,NmEf,isCSD,reconstruct,need_boundary_flux)
        end
    end

    # Transformation of restricted-angle fluxes to full-range fluxes (accumulate). Cut along
    # the spatial axis like the forward transform, so each range keeps its slice of 𝚽l in
    # cache across the whole patch loop.
    if Nq == 1
        𝚽_q_2d = reshape(𝚽_q, NS, Npatch)
        for ci in 1:length(chunks)
            rng = chunks[ci]
            @views mul!(𝚽l_mat[:,rng], Mll_all, transpose(𝚽_q_2d[rng,:]), 1.0, 1.0)
        end
    else
        # Mirror of the forward transform: gather the block into the patch-major staging
        # buffer, then one accumulating gemm. Each range owns its columns of 𝚽l, so the
        # accumulation over patches stays race-free.
        𝚽_q_vv = vec(𝚽_q)
        for t in 1:nth
            buf = tbufs[t]
            for bi in t:nth:length(blocks)
                rng = blocks[bi]
                gn_fast_gather!(buf,𝚽_q_vv,rng,Nq,Npatch,NS)
                @views mul!(𝚽l_mat[:,rng], Mll_all, buf[:,1:length(rng)], 1.0, 1.0)
            end
        end
    end
    # Outgoing energy flux: only the caller's final reconstruction pass needs it — it is not
    # part of the fixed-point state, it is read once gn_one_speed_fast has converged.
    if isCSD && reconstruct
        chunksE = gn_fast_chunks(NSE,1)
        if Nq == 1
            𝚽E12_q_2d = reshape(𝚽E12_q, NSE, Npatch)
            for ci in 1:length(chunksE)
                rng = chunksE[ci]
                @views mul!(𝚽E12_t_mat[:,rng], Mll_all, transpose(𝚽E12_q_2d[rng,:]), 1.0, 1.0)
            end
        else
            blocksE = gn_fast_chunks(NSE,max(1,cld(NSE,GN_FAST_BLOCK)))
            𝚽E12_q_vv = vec(𝚽E12_q)
            for t in 1:nth
                buf = tbufs[t]
                for bi in t:nth:length(blocksE)
                    rng = blocksE[bi]
                    gn_fast_gather!(buf,𝚽E12_q_vv,rng,Nq,Npatch,NSE)
                    @views mul!(𝚽E12_t_mat[:,rng], Mll_all, buf[:,1:length(rng)], 1.0, 1.0)
                end
            end
        end
    end
    # Outgoing boundary fluxes: only the boundary conditions consume them, so they are dead
    # work when every face is void.
    @views for k in (need_boundary_flux ? (1:Npatch) : (1:0))
        (u,v,w) = patch_uvw[k]
        for ib in 1:2
            mul!(reshape(𝚽x12⁺[:,:,:,:,ib], Np_surf, NSx),Mll_surf[:,:,u,v,w,ib,2],reshape(𝚽x12_q[:,:,:,:,ib,k], Nq, NSx),1.0, 1.0)
            if has_y
                mul!(reshape(𝚽y12⁺[:,:,:,:,ib], Np_surf, NSy),Mll_surf[:,:,u,v,w,ib+2,2],reshape(𝚽y12_q[:,:,:,:,ib,k], Nq, NSy),1.0, 1.0)
            end
            if has_z
                mul!(reshape(𝚽z12⁺[:,:,:,:,ib], Np_surf, NSz),Mll_surf[:,:,u,v,w,ib+4,2],reshape(𝚽z12_q[:,:,:,:,ib,k], Nq, NSz),1.0, 1.0)
            end
        end
    end

    # Boundary conditions treatment
    if ~need_boundary_flux return 𝚽l end
    𝚽x12⁻ .= 0.0
    if has_y 𝚽y12⁻ .= 0.0 end
    if has_z 𝚽z12⁻ .= 0.0 end
    for ib in range(1,2)
        # X-axis boundary conditions
        if boundary_conditions[ib] != 0
            if boundary_conditions[ib] == 1 # Reflective boundary condition
                @inbounds for iz in range(1,Ns[3]), iy in range(1,Ns[2]), is in range(1,Nm[1]), q in range(1,Np_surf)
                    s = 𝚽x12⁺[q,is,iy,iz,ib]
                    @simd for p in range(1,Np_surf)
                        𝚽x12⁻[p,is,iy,iz,ib] += Rpq[p,q,ib] * s
                    end
                end
            elseif boundary_conditions[ib] == 2 # Periodic boundary condition
                for p in range(1,Np_surf), is in range(1,Nm[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3])
                    if ib == 1
                        𝚽x12⁻[p,is,iy,iz,1] += 𝚽x12⁺[p,is,iy,iz,2]
                    else
                        𝚽x12⁻[p,is,iy,iz,2] += 𝚽x12⁺[p,is,iy,iz,1]
                    end
                end
            else
                error("Invalid boundary condition type.")
            end
        end
        # Y-axis boundary conditions
        if has_y && boundary_conditions[ib+2] != 0
            if boundary_conditions[ib+2] == 1 # Reflective boundary condition
                @inbounds for iz in range(1,Ns[3]), ix in range(1,Ns[1]), is in range(1,Nm[2]), q in range(1,Np_surf)
                    s = 𝚽y12⁺[q,is,ix,iz,ib]
                    @simd for p in range(1,Np_surf)
                        𝚽y12⁻[p,is,ix,iz,ib] += Rpq[p,q,ib+2] * s
                    end
                end
            elseif boundary_conditions[ib+2] == 2 # Periodic boundary condition
                for p in range(1,Np_surf), is in range(1,Nm[2]), ix in range(1,Ns[1]), iz in range(1,Ns[3])
                    if ib == 1
                        𝚽y12⁻[p,is,ix,iz,1] += 𝚽y12⁺[p,is,ix,iz,2]
                    else
                        𝚽y12⁻[p,is,ix,iz,2] += 𝚽y12⁺[p,is,ix,iz,1]
                    end
                end
            else
                error("Invalid boundary condition type.")
            end
        end
        # Z-axis boundary conditions
        if has_z && boundary_conditions[ib+4] != 0
            if boundary_conditions[ib+4] == 1 # Reflective boundary condition
                @inbounds for iy in range(1,Ns[2]), ix in range(1,Ns[1]), is in range(1,Nm[3]), q in range(1,Np_surf)
                    s = 𝚽z12⁺[q,is,ix,iy,ib]
                    @simd for p in range(1,Np_surf)
                        𝚽z12⁻[p,is,ix,iy,ib] += Rpq[p,q,ib+4] * s
                    end
                end
            elseif boundary_conditions[ib+4] == 2 # Periodic boundary condition
                for p in range(1,Np_surf), is in range(1,Nm[3]), ix in range(1,Ns[1]), iy in range(1,Ns[2])
                    if ib == 1
                        𝚽z12⁻[p,is,ix,iy,1] += 𝚽z12⁺[p,is,ix,iy,2]
                    else
                        𝚽z12⁻[p,is,ix,iy,2] += 𝚽z12⁺[p,is,ix,iy,1]
                    end
                end
            else
                error("Invalid boundary condition type.")
            end
        end
    end
    𝚽x12⁺ .= 0.0
    if has_y 𝚽y12⁺ .= 0.0 end
    if has_z 𝚽z12⁺ .= 0.0 end

    return 𝚽l
end
