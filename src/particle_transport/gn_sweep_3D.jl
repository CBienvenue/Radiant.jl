
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
