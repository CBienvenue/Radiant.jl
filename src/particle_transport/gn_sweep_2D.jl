
function gn_sweep_2D!(𝚽l::AbstractArray{Float64,4},𝚽E12::AbstractArray{Float64},𝚽x12::AbstractArray{Float64,4},𝚽y12::AbstractArray{Float64,4},sx::Int64,sy::Int64,Σt::Vector{Float64},mat::AbstractMatrix{Int64},Nx::Int64,Ny::Int64,Δx::Vector{Float64},Δy::Vector{Float64},Ql::AbstractArray{Float64,4},Np::Int64,Np_source::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},C::Vector{Float64},ω::Vector{Vector{Float64}},sources::AbstractMatrix{Union{Float64,Array{Float64}}},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},𝒲::Array{Float64},isFC::Bool,is_CSD::Bool,𝒩x::AbstractMatrix{Float64},𝒩y::AbstractMatrix{Float64},𝒮_ws::Matrix{Float64},Q_ws::Vector{Float64},𝚽_ws::Vector{Float64},𝚽x12_buf::Array{Float64,3},𝚽y12_buf::Matrix{Float64})

    # Sweep ordering and incoming/outgoing boundary slots
    if sx > 0; x_sweep = 1:Nx;     in_x = 1; out_x = 2 else x_sweep = Nx:-1:1; in_x = 2; out_x = 1 end
    if sy > 0; y_sweep = 1:Ny;     in_y = 1; out_y = 2 else y_sweep = Ny:-1:1; in_y = 2; out_y = 1 end

    # Reset moving x-boundary workspace
    fill!(𝚽x12_buf, 0.0)

    @inbounds @views for ix in x_sweep
        # Reset moving y-boundary workspace at each ix
        fill!(𝚽y12_buf, 0.0)

        # Y-boundary initialization (sources + incoming face)
        src_y = sy > 0 ? 3 : 4
        for p in 1:Np
            𝚽y12_buf[p,1] += sources[p,src_y][ix]
            for is in 1:Nm[2]
                𝚽y12_buf[p,is] += 𝚽y12[p,is,ix,in_y]
            end
        end

        for iy in y_sweep
            # X-boundary initialization (at ix entrance only)
            if (ix == 1 && sx > 0) || (ix == Nx && sx < 0)
                src_x = sx > 0 ? 1 : 2
                for p in 1:Np
                    𝚽x12_buf[p,1,iy] += sources[p,src_x][iy]
                    for is in 1:Nm[1]
                        𝚽x12_buf[p,is,iy] += 𝚽x12[p,is,iy,in_x]
                    end
                end
            end

            # Flux calculation
            if ~is_CSD
                gn_2D_BTE!(𝚽l[:,:,ix,iy],
                          𝚽x12_buf[:,:,iy],
                          𝚽y12_buf,
                          sx,sy,Σt[mat[ix,iy]],Δx[ix],Δy[iy],
                          Ql[:,:,ix,iy],
                          𝒮_ws,Q_ws,𝚽_ws,
                          𝒪[1],𝒪[2],Np,C,ω[1],ω[2],
                          𝒩x,𝒩y,isFC)
            else
                gn_2D_BFP!(𝚽l[:,:,ix,iy],
                          𝚽x12_buf[:,:,iy],
                          𝚽y12_buf,
                          𝚽E12[:,:,ix,iy],
                          sx,sy,Σt[mat[ix,iy]],S⁻[mat[ix,iy]],S⁺[mat[ix,iy]],S[mat[ix,iy],:],
                          Δx[ix],Δy[iy],
                          Ql[:,:,ix,iy],
                          𝒮_ws,Q_ws,𝚽_ws,
                          𝒪[1],𝒪[2],𝒪[4],Np,C,ω[1],ω[2],ω[4],
                          𝒩x,𝒩y,𝒲,isFC)
            end

            # Save x-outgoing boundary at the exit face
            if (ix == Nx && sx > 0) || (ix == 1 && sx < 0)
                for p in 1:Np, is in 1:Nm[1]
                    𝚽x12[p,is,iy,out_x] = 𝚽x12_buf[p,is,iy]
                end
            end
        end

        # Save y-outgoing boundary (end of y-sweep for this ix)
        for p in 1:Np, is in 1:Nm[2]
            𝚽y12[p,is,ix,out_y] = 𝚽y12_buf[p,is]
        end
    end

    # Zero out the incoming ib slots so the outgoing-transform mul!
    # in gn_one_speed only sees the freshly-computed outgoing contribution.
    fill!(view(𝚽x12,:,:,:,in_x), 0.0)
    fill!(view(𝚽y12,:,:,:,in_y), 0.0)

    return nothing
end

"""
    gn_sweep_2D_fast!(𝚽l,𝚽E12,𝚽x12,𝚽y12,k,mat,Nx,Ny,Ql,src_x,src_y,mom,cells,pat,ws,
    Nmf,NmEf,is_CSD,do_E,save_out)

Optimized counterpart of `gn_sweep_2D!`: sweep the spatial grid for one angular patch.

Same sweep ordering and moving-boundary bookkeeping as the reference, with the three changes
made in 3D (see `gn_sweep_3D_fast!`): concrete boundary-source arrays instead of the
abstractly-typed container the reference indexes, flat arrays addressed from a per-voxel
offset instead of a `SubArray` per voxel, and the cell system taken from its cached
factorization instead of being rebuilt.

The loop nesting is `iy` outer, `ix` inner, so the innermost index is the fastest one of the
arrays; the moving buffers follow suit — `bufy` spans the `x` line, `bufx` is a single face.
"""
function gn_sweep_2D_fast!(𝚽l::Vector{Float64},𝚽E12::Vector{Float64},𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},k::Int64,mat::Array{Int64,3},Nx::Int64,Ny::Int64,Ql::Vector{Float64},src_x::Vector{Float64},src_y::Vector{Float64},mom::GNFastMoments{NMOM},cells::GNFastCells,pat::GNFastPatch{NQ},ws::GNFastWorkspace,Nmf::Vector{Int64},NmEf::Int64,is_CSD::Bool,do_E::Bool,save_out::Bool) where {NMOM,NQ}

    Nq = NQ
    sx = pat.sx; sy = pat.sy
    Nmom = NMOM
    Nm1 = Nmf[1]; Nm2 = Nmf[2]

    # Sweep ordering and incoming/outgoing boundary slots
    if sx > 0; x_sweep = 1:Nx;     in_x = 1; out_x = 2 else x_sweep = Nx:-1:1; in_x = 2; out_x = 1 end
    if sy > 0; y_sweep = 1:Ny;     in_y = 1; out_y = 2 else y_sweep = Ny:-1:1; in_y = 2; out_y = 1 end

    bufx = ws.bufx; bufy = ws.bufy

    # Strides of the flat patch arrays. 𝚽l/Ql are (Nq,Nmom,Nx,Ny,Npatch), 𝚽E12 is
    # (Nq,NmEf,Nx,Ny,Npatch), and the face arrays are (Nq,Nmf[a],·,2,Npatch).
    sc  = Nq*Nmom
    scE = Nq*NmEf
    o_k  = sc  * Nx*Ny * (k-1)
    o_kE = scE * Nx*Ny * (k-1)
    sfx = Nq*Nm1; sfy = Nq*Nm2
    o_x_in  = sfx*Ny*((in_x -1) + 2*(k-1)); o_x_out = sfx*Ny*((out_x-1) + 2*(k-1))
    o_y_in  = sfy*Nx*((in_y -1) + 2*(k-1)); o_y_out = sfy*Nx*((out_y-1) + 2*(k-1))
    o_sx = Nq*Ny*(k-1); o_sy = Nq*Nx*(k-1)

    # Reset the moving y-boundary workspace (an x line: y is the outer sweep axis)
    fill!(bufy, 0.0)

    @inbounds for iy in y_sweep
        # Reset the moving x-boundary workspace at each iy
        fill!(bufx, 0.0)
        iky = Int64(cells.ky[iy])
        is_y_entry = (iy == 1 && sy > 0) || (iy == Ny && sy < 0)
        is_y_exit  = save_out && ((iy == Ny && sy > 0) || (iy == 1 && sy < 0))

        # X-boundary initialization (sources + incoming face)
        bsx = o_sx + Nq*(iy-1)
        bx0 = o_x_in + sfx*(iy-1)
        for p in 1:Nq
            bufx[p] += src_x[bsx+p]
            for is in 1:Nm1
                bufx[p + Nq*(is-1)] += 𝚽x12[bx0 + p + Nq*(is-1)]
            end
        end

        # Y-boundary initialization (at the iy entrance plane only). Each ix slot is
        # independent and is read by the kernel of that same ix, so the loop can run before
        # the x sweep rather than inside it.
        if is_y_entry
            for ix in 1:Nx
                bsy = o_sy + Nq*(ix-1)
                by0 = o_y_in + sfy*(ix-1)
                bb  = sfy*(ix-1)
                for p in 1:Nq
                    bufy[bb+p] += src_y[bsy+p]
                    for is in 1:Nm2
                        bufy[bb + p + Nq*(is-1)] += 𝚽y12[by0 + p + Nq*(is-1)]
                    end
                end
            end
        end

        for ix in x_sweep
            ikx = Int64(cells.kx[ix])
            icell = (ix-1) + Nx*(iy-1)
            ocell = o_k  + sc  * icell
            ocelE = o_kE + scE * icell
            oby = sfy*(ix-1)

            # Flux calculation
            m = mat[ix,iy,1]
            conf = gn_fast_conf(cells,m,ix,iy,1)
            if ~is_CSD
                gn_2D_BTE_fast!(𝚽l,ocell,bufx,0,bufy,oby,Ql,ocell,
                                mom,pat,ws,conf,ikx,iky)
            else
                gn_2D_BFP_fast!(𝚽l,ocell,bufx,0,bufy,oby,𝚽E12,ocelE,Ql,ocell,
                                mom,pat,ws,conf,ikx,iky,m,do_E)
            end

            # Save y-outgoing boundary at the exit plane
            if is_y_exit
                by1 = o_y_out + sfy*(ix-1)
                for i in 1:sfy
                    𝚽y12[by1+i] = bufy[oby+i]
                end
            end
        end

        # Save x-outgoing boundary (end of the x sweep for this iy)
        if save_out
            bx1 = o_x_out + sfx*(iy-1)
            for i in 1:sfx
                𝚽x12[bx1+i] = bufx[i]
            end
        end
    end

    # Zero out the incoming ib slots so the outgoing-transform mul! in
    # gn_one_speed_fast only sees the freshly-computed outgoing contribution.
    if save_out
        @inbounds begin
            for i in 1:sfx*Ny; 𝚽x12[o_x_in+i] = 0.0 end
            for i in 1:sfy*Nx; 𝚽y12[o_y_in+i] = 0.0 end
        end
    end

    return nothing
end
