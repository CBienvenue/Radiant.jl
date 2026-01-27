
function pn_sweep_3D(sx::Int64,sy::Int64,sz::Int64,𝚽l::Array{Float64,5},Ql::Array{Float64,5},Σt::Vector{Float64},mat::Array{Int64},Nx::Int64,Ny::Int64,Nz::Int64,Δx::Vector{Float64},Δy::Vector{Float64},Δz::Vector{Float64},Np::Int64,Np_source::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},C::Vector{Float64},ω::Vector{Array{Float64}},sources::Matrix{Union{Float64,Array{Float64}}},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},𝚽E12::Array{Float64},𝒲::Array{Float64},isFC::Bool,is_CSD::Bool,𝒩x::Matrix{Float64},𝒩y::Matrix{Float64},𝒩z::Matrix{Float64},𝚽x12⁺_in::Array{Float64,5},𝚽x12⁻_in::Array{Float64,5},𝚽y12⁺_in::Array{Float64,5},𝚽y12⁻_in::Array{Float64,5},𝚽z12⁺_in::Array{Float64,5},𝚽z12⁻_in::Array{Float64,5},boundary_conditions::Vector{Int64})

    # Initialization
    𝒪x = 𝒪[1]
    𝒪y = 𝒪[2]
    𝒪z = 𝒪[3]
    𝒪E = 𝒪[4]
    if (sx > 0) x_sweep = (1:Nx) else x_sweep = (Nx:-1:1) end
    if (sy > 0) y_sweep = (1:Ny) else y_sweep = (Ny:-1:1) end
    if (sz > 0) z_sweep = (1:Nz) else z_sweep = (Nz:-1:1) end
    𝚽x12_out = zeros(Np,Nm[1],2,Ny,Nz)
    𝚽y12_out = zeros(Np,Nm[2],2,Nx,Nz)
    𝚽z12_out = zeros(Np,Nm[3],2,Nx,Ny)

    # Sweep over x-axis
    𝚽x12 = zeros(Np,Nm[1],Ny,Nz)
    for ix in x_sweep
        𝚽y12 = zeros(Np,Nm[2],Nz)

        # Sweep over y-axis
        for iy in y_sweep
            𝚽z12 = zeros(Np,Nm[3])

            # Z-boundary initialization
            if sz > 0
                # Surface Z-
                for p in range(1,Np_source)
                    𝚽z12[p,1] += sources[p,5][ix,iy]
                end
                if boundary_conditions[5] != 0 # Not void
                    for p in range(1,Np), is in range(1,Nm[3])
                        if boundary_conditions[5] == 1 # Reflective at Z-
                            𝚽z12[p,is] += 𝚽z12⁻_in[p,is,1,ix,iy]
                        elseif boundary_conditions[5] == 2 # Periodic at Z-
                            𝚽z12[p,is] += 𝚽z12⁺_in[p,is,2,ix,iy]
                        end
                    end
                end
            else
                # Surface Z+
                for p in range(1,Np_source)
                    𝚽z12[p,1] += sources[p,6][ix,iy]
                end
                if boundary_conditions[6] != 0 # Not void
                    for p in range(1,Np), is in range(1,Nm[3])
                        if boundary_conditions[6] == 1 # Reflective at Z+
                            𝚽z12[p,is] += 𝚽z12⁻_in[p,is,2,ix,iy]
                        elseif boundary_conditions[6] == 2 # Periodic at Z+
                            𝚽z12[p,is] += 𝚽z12⁺_in[p,is,1,ix,iy]
                        end
                    end
                end
            end

            # Sweep over z-axis
            for iz in z_sweep

                # Y-boundary initialization
                if (iy == 1 && sy > 0) || (iy == Ny && sy < 0)
                    if sy > 0
                        # Surface Y-
                        for p in range(1,Np_source)
                            𝚽y12[p,1,iz] += sources[p,3][ix,iz]
                        end
                        if boundary_conditions[3] != 0 # Not void
                            for p in range(1,Np), is in range(1,Nm[2])
                                if boundary_conditions[3] == 1 # Reflective at Y-
                                    𝚽y12[p,is,iz] += 𝚽y12⁻_in[p,is,1,ix,iz]
                                elseif boundary_conditions[3] == 2 # Periodic at Y-
                                    𝚽y12[p,is,iz] += 𝚽y12⁺_in[p,is,2,ix,iz]
                                end
                            end
                        end
                    else
                        # Surface Y+
                        for p in range(1,Np_source)
                            𝚽y12[p,1,iz] += sources[p,4][ix,iz]
                        end
                        if boundary_conditions[4] != 0 # Not void
                            for p in range(1,Np), is in range(1,Nm[2])
                                if boundary_conditions[4] == 1 # Reflective at Y+
                                    𝚽y12[p,is,iz] += 𝚽y12⁻_in[p,is,2,ix,iz]
                                elseif boundary_conditions[4] == 2 # Periodic at Y+
                                    𝚽y12[p,is,iz] += 𝚽y12⁺_in[p,is,1,ix,iz]
                                end
                            end
                        end
                    end
                end

                # X-boundary initialization
                if (ix == 1 && sx > 0) || (ix == Nx && sx < 0)
                    if sx > 0
                        # Surface X-
                        for p in range(1,Np_source)
                            𝚽x12[p,1,iy,iz] += sources[p,1][iy,iz]
                        end
                        if boundary_conditions[1] != 0 # Not void
                            for p in range(1,Np), is in range(1,Nm[1])
                                if boundary_conditions[1] == 1 # Reflective at X-
                                    𝚽x12[p,is,iy,iz] += 𝚽x12⁻_in[p,is,1,iy,iz]
                                elseif boundary_conditions[1] == 2 # Periodic at X-
                                    𝚽x12[p,is,iy,iz] += 𝚽x12⁺_in[p,is,2,iy,iz]
                                end
                            end
                        end
                    else
                        # Surface X+
                        for p in range(1,Np_source)
                            𝚽x12[p,1,iy,iz] += sources[p,2][iy,iz]
                        end
                        if boundary_conditions[2] != 0 # Not void
                            for p in range(1,Np), is in range(1,Nm[1])
                                if boundary_conditions[2] == 1 # Reflective at X+
                                    𝚽x12[p,is,iy,iz] += 𝚽x12⁻_in[p,is,2,iy,iz]
                                elseif boundary_conditions[2] == 2 # Periodic at X+
                                    𝚽x12[p,is,iy,iz] += 𝚽x12⁺_in[p,is,1,iy,iz]
                                end
                            end
                        end
                    end
                end

                # Flux calculation
                if ~is_CSD
                    𝚽l[:,:,ix,iy,iz],𝚽x12[:,:,iy,iz],𝚽y12[:,:,iz],𝚽z12 = pn_3D_BTE(sx,sy,sz,Σt[mat[ix,iy,iz]],Δx[ix],Δy[iy],Δz[iz],Ql[:,:,ix,iy,iz],𝚽x12[:,:,iy,iz],𝚽y12[:,:,iz],𝚽z12,𝒪x,𝒪y,𝒪z,Np,C,copy(ω[1]),copy(ω[2]),copy(ω[3]),𝒩x,𝒩y,𝒩z,isFC)
                else
                    error("CSD method is not yet implemented for PN in 3D.")
                end

                # Save boundary fluxes along x-axis
                if (ix == Nx && sx > 0) || (ix == 1 && sx < 0)
                    for p in range(1,Np), is in range(1,Nm[1])
                        if sx > 0 # Surface X+
                            𝚽x12_out[p,is,2,iy,iz] = 𝚽x12[p,is,iy,iz]
                        else # Surface X-
                            𝚽x12_out[p,is,1,iy,iz] = 𝚽x12[p,is,iy,iz]
                        end
                    end
                end
            end

            # Save boundary fluxes along z-axis
            for p in range(1,Np), is in range(1,Nm[3])
                if sz > 0 # Surface Z+
                    𝚽z12_out[p,is,2,ix,iy] = 𝚽z12[p,is]
                else # Surface Z-
                    𝚽z12_out[p,is,1,ix,iy] = 𝚽z12[p,is]
                end
            end
        end

        # Save boundary fluxes along y-axis
        for iz in range(1,Nz)
            for p in range(1,Np), is in range(1,Nm[2])
                if sy > 0 # Surface Y+
                    𝚽y12_out[p,is,2,ix,iz] = 𝚽y12[p,is,iz]
                else # Surface Y-
                    𝚽y12_out[p,is,1,ix,iz] = 𝚽y12[p,is,iz]
                end
            end
        end
    end

    return 𝚽l, 𝚽E12, 𝚽x12_out, 𝚽y12_out, 𝚽z12_out
end