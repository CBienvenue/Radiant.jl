
function pn_sweep_1D(sx::Int64,𝚽l::Array{Float64,3},Ql::Array{Float64,3},Σt::Vector{Float64},mat::Vector{Int64},Nx::Int64,Δx::Vector{Float64},Np::Int64,Np_surf::Int64,𝒪::Vector{Int64},C::Vector{Float64},ω::Vector{Array{Float64}},sources::Matrix{Union{Float64,Array{Float64}}},is_SPH::Bool,pl::Vector{Int64},pm::Vector{Int64})

    # Initialization
    𝒪x = 𝒪[1]
    𝚽x12 = zeros(Np)

    # Monodirectional boundary sources
    if sx > 0
        x_sweep = 1:Nx
        # Surface X-
        for p in range(1,min(Np_surf,Np))
            𝚽x12[p] += sources[p,1]
        end
    elseif sx < 0
        x_sweep = Nx:-1:1
        # Surface X+
        for p in range(1,min(Np_surf,Np))
            𝚽x12[p] += sources[p,2]
        end
    end

    for ix in x_sweep

        # Flux calculation
        𝚽l[:,:,ix],𝚽x12 = pn_1D_BTE(sx,Σt[mat[ix]],Δx[ix],Ql[:,:,ix],𝚽x12,𝒪x,Np,C,copy(ω[1]),is_SPH,pl,pm)

    end
    return 𝚽l
end