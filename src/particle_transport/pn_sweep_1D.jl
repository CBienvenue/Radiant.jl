
function pn_sweep_1D(sx::Int64,𝚽l::Array{Float64,3},Ql::Array{Float64,3},Σt::Vector{Float64},mat::Vector{Int64},Nx::Int64,Δx::Vector{Float64},Np::Int64,Np_source::Int64,Np_surf::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},C::Vector{Float64},ω::Vector{Array{Float64}},sources::Matrix{Union{Float64,Array{Float64}}},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},𝚽E12::Array{Float64},𝒲::Array{Float64},isFC::Bool,is_CSD::Bool,𝒩x::Matrix{Float64},𝚽x12⁺_in::Array{Float64,3},𝚽x12⁻_in::Array{Float64,3},boundary_conditions::Vector{Int64})

    # Initialization
    𝒪x = 𝒪[1]
    𝒪E = 𝒪[4]
    𝚽x12 = zeros(Np,Nm[1])
    𝚽x12_out = zeros(Np,Nm[1],2)

    # Boundary conditions and sources
    if sx > 0
        x_sweep = 1:Nx
        # Surface X-
        for p in range(1,Np_source)
            𝚽x12[p,1] += sources[p,1]
        end
        if boundary_conditions[1] != 0 # Not void
            for p in range(1,Np), is in range(1,Nm[1])
                if boundary_conditions[1] == 1 # Reflective
                    𝚽x12[p,is] += 𝚽x12⁻_in[p,is,1]
                elseif boundary_conditions[1] == 2 # Periodic
                    𝚽x12[p,is] += 𝚽x12⁺_in[p,is,2]
                end
            end
        end
    else
        x_sweep = Nx:-1:1
        # Surface X+
        for p in range(1,Np_source)
            𝚽x12[p,1] += sources[p,2]
        end
        if boundary_conditions[2] != 0 # Not void
            for p in range(1,Np), is in range(1,Nm[1])
                if boundary_conditions[2] == 1 # Reflective
                    𝚽x12[p,is] += 𝚽x12⁺_in[p,is,2]
                elseif boundary_conditions[2] == 2 # Periodic
                    𝚽x12[p,is] += 𝚽x12⁻_in[p,is,1]
                end
            end
        end
    end

    for ix in x_sweep
        # Flux calculation
        if ~is_CSD
            𝚽l[:,:,ix],𝚽x12[:,1] = pn_1D_BTE(sx,Σt[mat[ix]],Δx[ix],Ql[:,:,ix],𝚽x12[:,1],𝒪x,Np,C,ω[1],𝒩x)
        else
            𝚽l[:,:,ix],𝚽x12,𝚽E12[:,:,ix] = pn_1D_BFP(sx,Σt[mat[ix]],Δx[ix],Ql[:,:,ix],𝚽x12,S⁻[mat[ix]],S⁺[mat[ix]],S[mat[ix],:],𝚽E12[:,:,ix],𝒪E,𝒪x,Np,C,ω[1],ω[2],𝒲,isFC,𝒩x)
        end
    end

    # Save boundary fluxes
    for p in range(1,Np), is in range(1,Nm[1])
        if sx > 0 # Surface X+
            𝚽x12_out[p,is,2] += 𝚽x12[p,is]
        else # Surface X-
            𝚽x12_out[p,is,1] += 𝚽x12[p,is]
        end
    end

    return 𝚽l, 𝚽E12, 𝚽x12_out
end