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
