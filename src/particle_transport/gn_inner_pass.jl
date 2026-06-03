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
                    for p in range(1,Np_surf), q in range(1,Np_surf), is in range(1,Nm[1])
                        𝚽x12⁻[p,is,ib] += Rpq[p,q,ib] * 𝚽x12⁺[q,is,ib]
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
                    for p in range(1,Np_surf), q in range(1,Np_surf), is in range(1,Nm[1]), iy in range(1,Ns[2])
                        𝚽x12⁻[p,is,iy,ib] += Rpq[p,q,ib] * 𝚽x12⁺[q,is,iy,ib]
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
                    for p in range(1,Np_surf), q in range(1,Np_surf), is in range(1,Nm[2]), ix in range(1,Ns[1])
                        𝚽y12⁻[p,is,ix,ib] += Rpq[p,q,ib+2] * 𝚽y12⁺[q,is,ix,ib]
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
                    for p in range(1,Np_surf), q in range(1,Np_surf), is in range(1,Nm[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3])
                        𝚽x12⁻[p,is,iy,iz,ib] += Rpq[p,q,ib] * 𝚽x12⁺[q,is,iy,iz,ib]
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
                    for p in range(1,Np_surf), q in range(1,Np_surf), is in range(1,Nm[2]), ix in range(1,Ns[1]), iz in range(1,Ns[3])
                        𝚽y12⁻[p,is,ix,iz,ib] += Rpq[p,q,ib+2] * 𝚽y12⁺[q,is,ix,iz,ib]
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
                    for p in range(1,Np_surf), q in range(1,Np_surf), is in range(1,Nm[3]), ix in range(1,Ns[1]), iy in range(1,Ns[2])
                        𝚽z12⁻[p,is,ix,iy,ib] += Rpq[p,q,ib+4] * 𝚽z12⁺[q,is,ix,iy,ib]
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
