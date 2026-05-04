
function gn_sweep_3D(sx::Int64,sy::Int64,sz::Int64,𝚽l::Array{Float64,5},Ql::Array{Float64,5},Σt::Vector{Float64},mat::Array{Int64},Nx::Int64,Ny::Int64,Nz::Int64,Δx::Vector{Float64},Δy::Vector{Float64},Δz::Vector{Float64},Np::Int64,Np_source::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},C::Vector{Float64},ω::Vector{Vector{Float64}},sources::Matrix{Union{Float64,Array{Float64}}},𝚽x12⁻::Array{Float64,5},𝚽y12⁻::Array{Float64,5},𝚽z12⁻::Array{Float64,5},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},𝚽E12::Array{Float64},𝒲::Array{Float64},isFC::Bool,is_CSD::Bool,𝒩x::Matrix{Float64},𝒩y::Matrix{Float64},𝒩z::Matrix{Float64})

    # Initialization
    𝒪x = 𝒪[1]
    𝒪y = 𝒪[2]
    𝒪z = 𝒪[3]
    𝒪E = 𝒪[4]
    if (sx > 0) x_sweep = (1:Nx) else x_sweep = (Nx:-1:1) end
    if (sy > 0) y_sweep = (1:Ny) else y_sweep = (Ny:-1:1) end
    if (sz > 0) z_sweep = (1:Nz) else z_sweep = (Nz:-1:1) end
    𝚽x12⁺ = zeros(Np,Nm[1],Ny,Nz,2)
    𝚽y12⁺ = zeros(Np,Nm[2],Nx,Nz,2)
    𝚽z12⁺ = zeros(Np,Nm[3],Nx,Ny,2)

    # Sweep over x-axis
    𝚽x12 = zeros(Np,Nm[1],Ny,Nz)
    for ix in x_sweep
        𝚽y12 = zeros(Np,Nm[2],Nz)

        # Sweep over y-axis
        for iy in y_sweep
            𝚽z12 = zeros(Np,Nm[3])

            # Z-boundary initialization (sources + reflective/periodic)
            if sz > 0
                # Surface Z-
                for p in range(1,Np)
                    𝚽z12[p,1] += sources[p,5][ix,iy]
                    for is in range(1,Nm[3])
                        𝚽z12[p,is] += 𝚽z12⁻[p,is,ix,iy,1]
                    end
                end
            else
                # Surface Z+
                for p in range(1,Np)
                    𝚽z12[p,1] += sources[p,6][ix,iy]
                    for is in range(1,Nm[3])
                        𝚽z12[p,is] += 𝚽z12⁻[p,is,ix,iy,2]
                    end
                end
            end

            # Sweep over z-axis
            for iz in z_sweep

                # Y-boundary initialization (sources + reflective/periodic)
                if (iy == 1 && sy > 0) || (iy == Ny && sy < 0)
                    if sy > 0
                        # Surface Y-
                        for p in range(1,Np)
                            𝚽y12[p,1,iz] += sources[p,3][ix,iz]
                            for is in range(1,Nm[2])
                                𝚽y12[p,is,iz] += 𝚽y12⁻[p,is,ix,iz,1]
                            end
                        end
                    else
                        # Surface Y+
                        for p in range(1,Np)
                            𝚽y12[p,1,iz] += sources[p,4][ix,iz]
                            for is in range(1,Nm[2])
                                𝚽y12[p,is,iz] += 𝚽y12⁻[p,is,ix,iz,2]
                            end
                        end
                    end
                end

                # X-boundary initialization (sources + reflective/periodic)
                if (ix == 1 && sx > 0) || (ix == Nx && sx < 0)
                    if sx > 0
                        # Surface X-
                        for p in range(1,Np)
                            𝚽x12[p,1,iy,iz] += sources[p,1][iy,iz]
                            for is in range(1,Nm[1])
                                𝚽x12[p,is,iy,iz] += 𝚽x12⁻[p,is,iy,iz,1]
                            end
                        end
                    else
                        # Surface X+
                        for p in range(1,Np)
                            𝚽x12[p,1,iy,iz] += sources[p,2][iy,iz]
                            for is in range(1,Nm[1])
                                𝚽x12[p,is,iy,iz] += 𝚽x12⁻[p,is,iy,iz,2]
                            end
                        end
                    end
                end

                # Flux calculation
                if ~is_CSD
                    𝚽l[:,:,ix,iy,iz],𝚽x12[:,:,iy,iz],𝚽y12[:,:,iz],𝚽z12 = gn_3D_BTE(sx,sy,sz,Σt[mat[ix,iy,iz]],Δx[ix],Δy[iy],Δz[iz],Ql[:,:,ix,iy,iz],𝚽x12[:,:,iy,iz],𝚽y12[:,:,iz],𝚽z12,𝒪x,𝒪y,𝒪z,Np,C,ω[1],ω[2],ω[3],𝒩x,𝒩y,𝒩z,isFC)
                else
                    𝚽l[:,:,ix,iy,iz],𝚽x12[:,:,iy,iz],𝚽y12[:,:,iz],𝚽z12,𝚽E12[:,:,ix,iy,iz] = gn_3D_BFP(sx,sy,sz,Σt[mat[ix,iy,iz]],S⁻[mat[ix,iy,iz]],S⁺[mat[ix,iy,iz]],S[mat[ix,iy,iz],:],Δx[ix],Δy[iy],Δz[iz],Ql[:,:,ix,iy,iz],𝚽x12[:,:,iy,iz],𝚽y12[:,:,iz],𝚽z12,𝚽E12[:,:,ix,iy,iz],Nm[1],Nm[2],Nm[3],Nm[4],Np,C,ω[1],ω[2],ω[3],ω[4],𝒩x,𝒩y,𝒩z,𝒲,isFC)
                end

                # Save boundary fluxes along x-axis (far boundary for this sweep)
                if (ix == Nx && sx > 0) || (ix == 1 && sx < 0)
                    for p in range(1,Np), is in range(1,Nm[1])
                        if sx > 0 # Surface X+
                            𝚽x12⁺[p,is,iy,iz,2] = 𝚽x12[p,is,iy,iz]
                        else # Surface X-
                            𝚽x12⁺[p,is,iy,iz,1] = 𝚽x12[p,is,iy,iz]
                        end
                    end
                end
            end

            # Save boundary fluxes along z-axis (far boundary for this sweep)
            for p in range(1,Np), is in range(1,Nm[3])
                if sz > 0 # Surface Z+
                    𝚽z12⁺[p,is,ix,iy,2] = 𝚽z12[p,is]
                else # Surface Z-
                    𝚽z12⁺[p,is,ix,iy,1] = 𝚽z12[p,is]
                end
            end
        end

        # Save boundary fluxes along y-axis (far boundary for this sweep)
        for p in range(1,Np), is in range(1,Nm[2]), iz in range(1,Nz)
            if sy > 0 # Surface Y+
                𝚽y12⁺[p,is,ix,iz,2] = 𝚽y12[p,is,iz]
            else # Surface Y-
                𝚽y12⁺[p,is,ix,iz,1] = 𝚽y12[p,is,iz]
            end
        end
    end

    return 𝚽l, 𝚽E12, 𝚽x12⁺, 𝚽y12⁺, 𝚽z12⁺
end
