function gn_sweep_1D!(𝚽l::AbstractArray{Float64,3},𝚽E12::AbstractArray{Float64},𝚽x12::AbstractArray{Float64,3},sx::Int64,Σt::Vector{Float64},mat::AbstractVector{Int64},Nx::Int64,Δx::Vector{Float64},Ql::AbstractArray{Float64,3},Np::Int64,Np_source::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},C::Vector{Float64},ω::Vector{Vector{Float64}},sources::AbstractMatrix{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},𝒲::Array{Float64},isFC::Bool,is_CSD::Bool,𝒩x::AbstractMatrix{Float64},𝒮_ws::Matrix{Float64},Q_ws::Vector{Float64},𝚽_ws::Vector{Float64},𝚽x12_buf::Matrix{Float64})

    # Sweep ordering and incoming/outgoing boundary slots
    if sx > 0; x_sweep = 1:Nx;     in_x = 1; out_x = 2 else x_sweep = Nx:-1:1; in_x = 2; out_x = 1 end

    # Reset moving x-boundary workspace, then seed with sources + incoming face
    fill!(𝚽x12_buf, 0.0)
    @inbounds @views begin
        src_x = sx > 0 ? 1 : 2
        for p in 1:Np
            𝚽x12_buf[p,1] += sources[p,src_x]
            for is in 1:Nm[1]
                𝚽x12_buf[p,is] += 𝚽x12[p,is,in_x]
            end
        end

        for ix in x_sweep
            if ~is_CSD
                gn_1D_BTE!(𝚽l[:,:,ix],
                          𝚽x12_buf[:,1],
                          sx,Σt[mat[ix]],Δx[ix],
                          Ql[:,:,ix],
                          𝒮_ws,Q_ws,𝚽_ws,
                          𝒪[1],Np,C,ω[1],𝒩x)
            else
                gn_1D_BFP!(𝚽l[:,:,ix],
                          𝚽x12_buf,
                          𝚽E12[:,:,ix],
                          sx,Σt[mat[ix]],S⁻[mat[ix]],S⁺[mat[ix]],S[mat[ix],:],Δx[ix],
                          Ql[:,:,ix],
                          𝒮_ws,Q_ws,𝚽_ws,
                          𝒪[1],𝒪[4],Np,C,ω[1],ω[4],𝒩x,𝒲,isFC)
            end
        end

        # Save outgoing boundary
        for p in 1:Np, is in 1:Nm[1]
            𝚽x12[p,is,out_x] = 𝚽x12_buf[p,is]
        end

        # Zero out the incoming ib slot so the outgoing-transform mul!
        # in gn_one_speed only sees the freshly-computed outgoing contribution.
        fill!(view(𝚽x12,:,:,in_x), 0.0)
    end

    return nothing
end

"""
    gn_sweep_1D_fast!(𝚽l,𝚽E12,𝚽x12,k,mat,Nx,Ql,src_x,mom,cells,pat,ws,Nmf,NmEf,
    is_CSD,do_E,save_out)

Optimized counterpart of `gn_sweep_1D!`: sweep the line for one angular patch.

The 1D sweep has a single entrance face, so the boundary seeding happens once before the
loop rather than at an entrance plane. Otherwise the same three changes as in 2D/3D: a
concrete boundary-source array, flat arrays addressed from a per-voxel offset, and the cell
system taken from its cached factorization.
"""
function gn_sweep_1D_fast!(𝚽l::Vector{Float64},𝚽E12::Vector{Float64},𝚽x12::Vector{Float64},k::Int64,mat::Array{Int64,3},Nx::Int64,Ql::Vector{Float64},src_x::Vector{Float64},mom::GNFastMoments{NMOM},cells::GNFastCells,pat::GNFastPatch{NQ},ws::GNFastWorkspace,Nmf::Vector{Int64},NmEf::Int64,is_CSD::Bool,do_E::Bool,save_out::Bool) where {NMOM,NQ}

    Nq = NQ
    sx = pat.sx
    Nmom = NMOM
    Nm1 = Nmf[1]

    if sx > 0; x_sweep = 1:Nx; in_x = 1; out_x = 2 else x_sweep = Nx:-1:1; in_x = 2; out_x = 1 end

    bufx = ws.bufx

    sc  = Nq*Nmom
    scE = Nq*NmEf
    o_k  = sc  * Nx * (k-1)
    o_kE = scE * Nx * (k-1)
    sfx = Nq*Nm1
    o_x_in  = sfx*((in_x -1) + 2*(k-1))
    o_x_out = sfx*((out_x-1) + 2*(k-1))
    o_sx = Nq*(k-1)

    # Reset the moving boundary workspace, then seed it with the source and incoming face
    fill!(bufx, 0.0)
    @inbounds begin
        for p in 1:Nq
            bufx[p] += src_x[o_sx+p]
            for is in 1:Nm1
                bufx[p + Nq*(is-1)] += 𝚽x12[o_x_in + p + Nq*(is-1)]
            end
        end

        for ix in x_sweep
            ikx = Int64(cells.kx[ix])
            ocell = o_k  + sc  * (ix-1)
            ocelE = o_kE + scE * (ix-1)
            m = mat[ix,1,1]
            conf = gn_fast_conf(cells,m,ix,1,1)
            if ~is_CSD
                gn_1D_BTE_fast!(𝚽l,ocell,bufx,0,Ql,ocell,mom,pat,ws,conf,ikx)
            else
                gn_1D_BFP_fast!(𝚽l,ocell,bufx,0,𝚽E12,ocelE,Ql,ocell,mom,pat,ws,conf,ikx,m,do_E)
            end
        end

        # Save the outgoing boundary, then clear the incoming slot so the outgoing-transform
        # mul! in gn_one_speed_fast only sees the freshly-computed contribution.
        if save_out
            for i in 1:sfx
                𝚽x12[o_x_out+i] = bufx[i]
                𝚽x12[o_x_in+i]  = 0.0
            end
        end
    end

    return nothing
end
