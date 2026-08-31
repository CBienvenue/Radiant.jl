
function gn_sweep_3D!(𝚽l::AbstractArray{Float64,5},𝚽E12::AbstractArray{Float64},𝚽x12::AbstractArray{Float64,5},𝚽y12::AbstractArray{Float64,5},𝚽z12::AbstractArray{Float64,5},sx::Int64,sy::Int64,sz::Int64,Σt::Vector{Float64},mat::Array{Int64},Nx::Int64,Ny::Int64,Nz::Int64,Δx::Vector{Float64},Δy::Vector{Float64},Δz::Vector{Float64},Ql::AbstractArray{Float64,5},Np::Int64,Np_source::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},C::Vector{Float64},ω::Vector{Vector{Float64}},sources::AbstractMatrix{Union{Float64,Array{Float64}}},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},𝒲::Array{Float64},isFC::Bool,is_CSD::Bool,𝒩x::AbstractMatrix{Float64},𝒩y::AbstractMatrix{Float64},𝒩z::AbstractMatrix{Float64},𝒮_ws::Matrix{Float64},Q_ws::Vector{Float64},𝚽_ws::Vector{Float64},𝚽x12_buf::Array{Float64,4},𝚽y12_buf::Array{Float64,3},𝚽z12_buf::Matrix{Float64})

    # Sweep ordering and incoming/outgoing boundary slots
    if sx > 0; x_sweep = 1:Nx;     in_x = 1; out_x = 2 else x_sweep = Nx:-1:1; in_x = 2; out_x = 1 end
    if sy > 0; y_sweep = 1:Ny;     in_y = 1; out_y = 2 else y_sweep = Ny:-1:1; in_y = 2; out_y = 1 end
    if sz > 0; z_sweep = 1:Nz;     in_z = 1; out_z = 2 else z_sweep = Nz:-1:1; in_z = 2; out_z = 1 end

    # Reset moving x-boundary workspace
    fill!(𝚽x12_buf, 0.0)

    @inbounds @views for ix in x_sweep
        # Reset moving y-boundary workspace at each ix
        fill!(𝚽y12_buf, 0.0)

        for iy in y_sweep
            # Reset moving z-boundary workspace at each iy
            fill!(𝚽z12_buf, 0.0)

            # Z-boundary initialization (sources + incoming face)
            src_z = sz > 0 ? 5 : 6
            for p in 1:Np
                𝚽z12_buf[p,1] += sources[p,src_z][ix,iy]
                for is in 1:Nm[3]
                    𝚽z12_buf[p,is] += 𝚽z12[p,is,ix,iy,in_z]
                end
            end

            for iz in z_sweep
                # Y-boundary initialization (at iy entrance only)
                if (iy == 1 && sy > 0) || (iy == Ny && sy < 0)
                    src_y = sy > 0 ? 3 : 4
                    for p in 1:Np
                        𝚽y12_buf[p,1,iz] += sources[p,src_y][ix,iz]
                        for is in 1:Nm[2]
                            𝚽y12_buf[p,is,iz] += 𝚽y12[p,is,ix,iz,in_y]
                        end
                    end
                end

                # X-boundary initialization (at ix entrance only)
                if (ix == 1 && sx > 0) || (ix == Nx && sx < 0)
                    src_x = sx > 0 ? 1 : 2
                    for p in 1:Np
                        𝚽x12_buf[p,1,iy,iz] += sources[p,src_x][iy,iz]
                        for is in 1:Nm[1]
                            𝚽x12_buf[p,is,iy,iz] += 𝚽x12[p,is,iy,iz,in_x]
                        end
                    end
                end

                # Flux calculation
                if ~is_CSD
                    gn_3D_BTE!(𝚽l[:,:,ix,iy,iz],
                              𝚽x12_buf[:,:,iy,iz],
                              𝚽y12_buf[:,:,iz],
                              𝚽z12_buf,
                              sx,sy,sz,Σt[mat[ix,iy,iz]],Δx[ix],Δy[iy],Δz[iz],
                              Ql[:,:,ix,iy,iz],
                              𝒮_ws,Q_ws,𝚽_ws,
                              𝒪[1],𝒪[2],𝒪[3],Np,C,ω[1],ω[2],ω[3],
                              𝒩x,𝒩y,𝒩z,isFC)
                else
                    gn_3D_BFP!(𝚽l[:,:,ix,iy,iz],
                              𝚽x12_buf[:,:,iy,iz],
                              𝚽y12_buf[:,:,iz],
                              𝚽z12_buf,
                              𝚽E12[:,:,ix,iy,iz],
                              sx,sy,sz,Σt[mat[ix,iy,iz]],S⁻[mat[ix,iy,iz]],S⁺[mat[ix,iy,iz]],S[mat[ix,iy,iz],:],
                              Δx[ix],Δy[iy],Δz[iz],
                              Ql[:,:,ix,iy,iz],
                              𝒮_ws,Q_ws,𝚽_ws,
                              𝒪[1],𝒪[2],𝒪[3],𝒪[4],Np,C,ω[1],ω[2],ω[3],ω[4],
                              𝒩x,𝒩y,𝒩z,𝒲,isFC)
                end

                # Save x-outgoing boundary at the exit face
                if (ix == Nx && sx > 0) || (ix == 1 && sx < 0)
                    for p in 1:Np, is in 1:Nm[1]
                        𝚽x12[p,is,iy,iz,out_x] = 𝚽x12_buf[p,is,iy,iz]
                    end
                end
            end

            # Save z-outgoing boundary (end of z-sweep for this (ix,iy))
            for p in 1:Np, is in 1:Nm[3]
                𝚽z12[p,is,ix,iy,out_z] = 𝚽z12_buf[p,is]
            end
        end

        # Save y-outgoing boundary (end of y-sweep for this ix)
        for p in 1:Np, is in 1:Nm[2], iz in 1:Nz
            𝚽y12[p,is,ix,iz,out_y] = 𝚽y12_buf[p,is,iz]
        end
    end

    # Zero out the incoming ib slots so the outgoing-transform mul! in
    # gn_one_speed only sees the freshly-computed outgoing contribution.
    fill!(view(𝚽x12,:,:,:,:,in_x), 0.0)
    fill!(view(𝚽y12,:,:,:,:,in_y), 0.0)
    fill!(view(𝚽z12,:,:,:,:,in_z), 0.0)

    return nothing
end

"""
    gn_sweep_3D_fast!(𝚽l,𝚽E12,𝚽x12,𝚽y12,𝚽z12,mat,Nx,Ny,Nz,Ql,src_x,src_y,src_z,
    mom,cells,pat,ws,Nmf,is_CSD)

Optimized counterpart of `gn_sweep_3D!`: sweep the spatial grid for one angular patch.

The sweep ordering and the moving-boundary bookkeeping are identical to the reference; two
things change.

The boundary sources arrive as three concrete `Array{Float64,3}` — the incoming x, y and z
faces already selected for this octant — instead of the
`AbstractMatrix{Union{Float64,Array{Float64}}}` the reference indexes. That container's
element type is abstract, so every `sources[p,face][i,j]` access in the sweep costs a
dynamic dispatch (visible in the profile).

The per-cell work goes to `gn_3D_BTE_fast!`/`gn_3D_BFP_fast!` with the voxel's cell-system
configuration, rather than rebuilding and refactorizing the cell matrix in place.

All scratch comes from `ws`, one `GNFastWorkspace` per thread, so several patches can be
swept concurrently. See `set_fast_path`.
"""
function gn_sweep_3D_fast!(𝚽l::Vector{Float64},𝚽E12::Vector{Float64},𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},𝚽z12::Vector{Float64},k::Int64,mat::Array{Int64,3},Nx::Int64,Ny::Int64,Nz::Int64,Ql::Vector{Float64},src_x::Vector{Float64},src_y::Vector{Float64},src_z::Vector{Float64},mom::GNFastMoments{NMOM},cells::GNFastCells,pat::GNFastPatch{NQ},ws::GNFastWorkspace,Nmf::Vector{Int64},NmEf::Int64,is_CSD::Bool,do_E::Bool,save_out::Bool) where {NMOM,NQ}

    Nq = NQ
    sx = pat.sx; sy = pat.sy; sz = pat.sz
    Nmom = NMOM
    Nm1 = Nmf[1]; Nm2 = Nmf[2]; Nm3 = Nmf[3]

    # Sweep ordering and incoming/outgoing boundary slots
    if sx > 0; x_sweep = 1:Nx;     in_x = 1; out_x = 2 else x_sweep = Nx:-1:1; in_x = 2; out_x = 1 end
    if sy > 0; y_sweep = 1:Ny;     in_y = 1; out_y = 2 else y_sweep = Ny:-1:1; in_y = 2; out_y = 1 end
    if sz > 0; z_sweep = 1:Nz;     in_z = 1; out_z = 2 else z_sweep = Nz:-1:1; in_z = 2; out_z = 1 end

    bufx = ws.bufx; bufy = ws.bufy; bufz = ws.bufz

    # Strides of the flat patch arrays. 𝚽l/Ql are (Nq,Nmom,Nx,Ny,Nz,Npatch), 𝚽E12 is
    # (Nq,NmEf,Nx,Ny,Nz,Npatch), and the three face arrays are (Nq,Nmf[a],·,·,2,Npatch).
    sc  = Nq*Nmom                       # one voxel of 𝚽l / Ql
    scE = Nq*NmEf                       # one voxel of 𝚽E12
    o_k  = sc  * Nx*Ny*Nz * (k-1)       # base of this patch
    o_kE = scE * Nx*Ny*Nz * (k-1)
    sfx = Nq*Nm1; sfy = Nq*Nm2; sfz = Nq*Nm3
    o_x_in  = sfx*Ny*Nz*((in_x -1) + 2*(k-1)); o_x_out = sfx*Ny*Nz*((out_x-1) + 2*(k-1))
    o_y_in  = sfy*Nx*Nz*((in_y -1) + 2*(k-1)); o_y_out = sfy*Nx*Nz*((out_y-1) + 2*(k-1))
    o_z_in  = sfz*Nx*Ny*((in_z -1) + 2*(k-1)); o_z_out = sfz*Nx*Ny*((out_z-1) + 2*(k-1))
    o_sx = Nq*Ny*Nz*(k-1); o_sy = Nq*Nx*Nz*(k-1); o_sz = Nq*Nx*Ny*(k-1)

    # Reset the moving z-boundary workspace (a full x-y plane: z is the outermost sweep axis)
    fill!(bufz, 0.0)

    @inbounds for iz in z_sweep
        # Reset the moving y-boundary workspace at each iz
        fill!(bufy, 0.0)
        ikz = Int64(cells.kz[iz])
        # Entrance/exit-plane tests are invariant inside the sweep: hoist them out.
        is_z_entry = (iz == 1 && sz > 0) || (iz == Nz && sz < 0)
        is_z_exit  = save_out && ((iz == Nz && sz > 0) || (iz == 1 && sz < 0))

        for iy in y_sweep
            # Reset the moving x-boundary workspace at each iy
            fill!(bufx, 0.0)
            iky = Int64(cells.ky[iy])

            # X-boundary initialization (sources + incoming face)
            bsx = o_sx + Nq*((iy-1) + Ny*(iz-1))
            bx0 = o_x_in + sfx*((iy-1) + Ny*(iz-1))
            for p in 1:Nq
                bufx[p] += src_x[bsx+p]
                for is in 1:Nm1
                    bufx[p + Nq*(is-1)] += 𝚽x12[bx0 + p + Nq*(is-1)]
                end
            end

            # Y-boundary initialization (at the iy entrance plane only). Each ix slot is
            # independent and is read by the kernel of that same ix, so the whole loop can
            # run before the x sweep rather than inside it.
            if (iy == 1 && sy > 0) || (iy == Ny && sy < 0)
                for ix in 1:Nx
                    bsy = o_sy + Nq*((ix-1) + Nx*(iz-1))
                    by0 = o_y_in + sfy*((ix-1) + Nx*(iz-1))
                    bb  = sfy*(ix-1)
                    for p in 1:Nq
                        bufy[bb+p] += src_y[bsy+p]
                        for is in 1:Nm2
                            bufy[bb + p + Nq*(is-1)] += 𝚽y12[by0 + p + Nq*(is-1)]
                        end
                    end
                end
            end

            # Z-boundary initialization (at the iz entrance plane only)
            if is_z_entry
                for ix in 1:Nx
                    bsz = o_sz + Nq*((ix-1) + Nx*(iy-1))
                    bz0 = o_z_in + sfz*((ix-1) + Nx*(iy-1))
                    bb  = sfz*((ix-1) + Nx*(iy-1))
                    for p in 1:Nq
                        bufz[bb+p] += src_z[bsz+p]
                        for is in 1:Nm3
                            bufz[bb + p + Nq*(is-1)] += 𝚽z12[bz0 + p + Nq*(is-1)]
                        end
                    end
                end
            end

            for ix in x_sweep
                ikx = Int64(cells.kx[ix])
                icell = (ix-1) + Nx*((iy-1) + Ny*(iz-1))
                ocell = o_k  + sc  * icell
                ocelE = o_kE + scE * icell
                oby = sfy*(ix-1)
                obz = sfz*((ix-1) + Nx*(iy-1))

                # Flux calculation
                m = mat[ix,iy,iz]
                conf = gn_fast_conf(cells,m,ix,iy,iz)
                if ~is_CSD
                    gn_3D_BTE_fast!(𝚽l,ocell,bufx,0,bufy,oby,bufz,obz,Ql,ocell,
                                    mom,pat,ws,conf,ikx,iky,ikz)
                else
                    gn_3D_BFP_fast!(𝚽l,ocell,bufx,0,bufy,oby,bufz,obz,𝚽E12,ocelE,Ql,ocell,
                                    mom,pat,ws,conf,ikx,iky,ikz,m,do_E)
                end

                # Save z-outgoing boundary at the exit plane
                if is_z_exit
                    bz1 = o_z_out + sfz*((ix-1) + Nx*(iy-1))
                    for i in 1:sfz
                        𝚽z12[bz1+i] = bufz[obz+i]
                    end
                end
            end

            # Save x-outgoing boundary (end of the x sweep for this (iy,iz))
            if save_out
                bx1 = o_x_out + sfx*((iy-1) + Ny*(iz-1))
                for i in 1:sfx
                    𝚽x12[bx1+i] = bufx[i]
                end
            end
        end

        # Save y-outgoing boundary (end of the y sweep for this iz)
        if save_out
            for ix in 1:Nx
                by1 = o_y_out + sfy*((ix-1) + Nx*(iz-1))
                bb  = sfy*(ix-1)
                for i in 1:sfy
                    𝚽y12[by1+i] = bufy[bb+i]
                end
            end
        end
    end

    # Zero out the incoming ib slots so the outgoing-transform mul! in
    # gn_one_speed_fast only sees the freshly-computed outgoing contribution.
    if save_out
        @inbounds begin
            for i in 1:sfx*Ny*Nz; 𝚽x12[o_x_in+i] = 0.0 end
            for i in 1:sfy*Nx*Nz; 𝚽y12[o_y_in+i] = 0.0 end
            for i in 1:sfz*Nx*Ny; 𝚽z12[o_z_in+i] = 0.0 end
        end
    end

    return nothing
end
