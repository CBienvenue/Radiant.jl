
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
