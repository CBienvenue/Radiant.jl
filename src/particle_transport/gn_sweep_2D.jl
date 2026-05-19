
function gn_sweep_2D(sx::Int64,sy::Int64,𝚽l::Array{Float64,4},Ql::Array{Float64,4},Σt::Vector{Float64},mat::Matrix{Int64},Nx::Int64,Ny::Int64,Δx::Vector{Float64},Δy::Vector{Float64},Np::Int64,Np_source::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},C::Vector{Float64},ω::Vector{Vector{Float64}},sources::Matrix{Union{Float64,Array{Float64}}},𝚽x12⁻::Array{Float64,4},𝚽y12⁻::Array{Float64,4},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},𝚽E12::Array{Float64},𝒲::Array{Float64},isFC::Bool,is_CSD::Bool,𝒩x::Matrix{Float64},𝒩y::Matrix{Float64})

    # Initialization
    𝒪x = 𝒪[1]
    𝒪y = 𝒪[2]
    𝒪E = 𝒪[4]
    if (sx > 0) x_sweep = (1:Nx) else x_sweep = (Nx:-1:1) end
    if (sy > 0) y_sweep = (1:Ny) else y_sweep = (Ny:-1:1) end
    𝚽x12⁺ = zeros(Np,Nm[1],Ny,2)
    𝚽y12⁺ = zeros(Np,Nm[2],Nx,2)

    # Sweep over x-axis
    𝚽x12 = zeros(Np,Nm[1],Ny)
    for ix in x_sweep
        𝚽y12 = zeros(Np,Nm[2])

        # Y-boundary initialization (sources + reflective/periodic)
        if sy > 0
            # Surface Y-
            for p in range(1,Np)
                𝚽y12[p,1] += sources[p,3][ix]
                for is in range(1,Nm[2])
                    𝚽y12[p,is] += 𝚽y12⁻[p,is,ix,1]
                end
            end
        else
            # Surface Y+
            for p in range(1,Np)
                𝚽y12[p,1] += sources[p,4][ix]
                for is in range(1,Nm[2])
                    𝚽y12[p,is] += 𝚽y12⁻[p,is,ix,2]
                end
            end
        end

        # Sweep over y-axis
        for iy in y_sweep
            # X-boundary initialization (sources + reflective/periodic)
            if (ix == 1 && sx > 0) || (ix == Nx && sx < 0 )
                if sx > 0
                    # Surface X-
                    for p in range(1,Np)
                        𝚽x12[p,1,iy] += sources[p,1][iy]
                        for is in range(1,Nm[1])
                            𝚽x12[p,is,iy] += 𝚽x12⁻[p,is,iy,1]
                        end
                    end
                else
                    # Surface X+
                    for p in range(1,Np)
                        𝚽x12[p,1,iy] += sources[p,2][iy]
                        for is in range(1,Nm[1])
                            𝚽x12[p,is,iy] += 𝚽x12⁻[p,is,iy,2]
                        end
                    end
                end
            end

            # Flux calculation
            if ~is_CSD
                𝚽l[:,:,ix,iy],𝚽x12[:,:,iy],𝚽y12 = gn_2D_BTE(sx,sy,Σt[mat[ix,iy]],Δx[ix],Δy[iy],Ql[:,:,ix,iy],𝚽x12[:,:,iy],𝚽y12,𝒪x,𝒪y,Np,C,ω[1],ω[2],𝒩x,𝒩y,isFC)
            else
                𝚽l[:,:,ix,iy],𝚽x12[:,:,iy],𝚽y12,𝚽E12[:,:,ix,iy] = gn_2D_BFP(sx,sy,Σt[mat[ix,iy]],S⁻[mat[ix,iy]],S⁺[mat[ix,iy]],S[mat[ix,iy],:],Δx[ix],Δy[iy],Ql[:,:,ix,iy],𝚽x12[:,:,iy],𝚽y12,𝚽E12[:,:,ix,iy],𝒪x,𝒪y,𝒪E,Np,C,ω[1],ω[2],ω[4],𝒩x,𝒩y,𝒲,isFC)
            end

            # Save boundary fluxes along x-axis (far boundary for this sweep)
            if (ix == Nx && sx > 0) || (ix == 1 && sx < 0 )
                for p in range(1,Np), is in range(1,Nm[1])
                    if sx > 0 # Surface X+
                        𝚽x12⁺[p,is,iy,2] = 𝚽x12[p,is,iy]
                    else # Surface X-
                        𝚽x12⁺[p,is,iy,1] = 𝚽x12[p,is,iy]
                    end
                end
            end
        end

        # Save boundary fluxes along y-axis (far boundary for this sweep)
        for p in range(1,Np), is in range(1,Nm[2])
            if sy > 0 # Surface Y+
                𝚽y12⁺[p,is,ix,2] = 𝚽y12[p,is]
            else # Surface Y-
                𝚽y12⁺[p,is,ix,1] = 𝚽y12[p,is]
            end
        end
    end
    return 𝚽l, 𝚽E12, 𝚽x12⁺, 𝚽y12⁺
end