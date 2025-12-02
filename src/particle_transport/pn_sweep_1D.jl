
function pn_sweep_1D(sx::Int64,𝚽l::Array{Float64,3},Ql::Array{Float64,3},Σt::Vector{Float64},mat::Vector{Int64},Nx::Int64,Δx::Vector{Float64},Np::Int64,Np_surf::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},C::Vector{Float64},ω::Vector{Array{Float64}},sources::Matrix{Union{Float64,Array{Float64}}},is_SPH::Bool,pl::Vector{Int64},pm::Vector{Int64},S⁻,S⁺,S,𝚽E12,𝒲,isFC,is_CSD)

    # Initialization
    𝒪x = 𝒪[1]
    𝒪E = 𝒪[4]
    𝚽x12 = zeros(Np,Nm[1])

    # Monodirectional boundary sources
    if sx > 0
        x_sweep = 1:Nx
        # Surface X-
        for p in range(1,min(Np_surf,Np))
            𝚽x12[p,1] += sources[p,1]
        end
    elseif sx < 0
        x_sweep = Nx:-1:1
        # Surface X+
        for p in range(1,min(Np_surf,Np))
            𝚽x12[p,1] += sources[p,2]
        end
    end

    for ix in x_sweep

        # Flux calculation
        if ~is_CSD
            𝚽l[:,:,ix],𝚽x12[:,1] = pn_1D_BTE(sx,Σt[mat[ix]],Δx[ix],Ql[:,:,ix],𝚽x12[:,1],𝒪x,Np,C,copy(ω[1]),is_SPH,pl,pm)
        else
            𝚽l[:,:,ix],𝚽x12,𝚽E12[:,:,ix] = pn_1D_BFP(sx,Σt[mat[ix]],Δx[ix],Ql[:,:,ix],𝚽x12,S⁻[mat[ix]],S⁺[mat[ix]],S[mat[ix],:],𝚽E12[:,:,ix],𝒪E,𝒪x,Np,C,copy(ω[1]),copy(ω[2]),is_SPH,pl,pm,𝒲,isFC)
        end

    end
    return 𝚽l, 𝚽E12
end