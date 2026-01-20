"""
    compute_sweep_2D(𝚽l::Array{Float64,4},Ql::Array{Float64,4},Σt::Vector{Float64},
    mat::Array{Int64,2},Ns::Vector{Int64},Δs::Vector{Vector{Float64}},Ω::Vector{Float64},
    Mn::Vector{Float64},Dn::Vector{Float64},P::Int64,Mnx⁻::Vector{Float64},
    Dnx⁻::Vector{Float64},Mny⁻::Vector{Float64},Dny⁻::Vector{Float64},Np_surf::Int64,
    𝒪::Vector{Int64},Nm::Vector{Int64},C::Vector{Float64},ω::Vector{Array{Float64}},
    sources::Matrix{Union{Float64,Array{Float64}}},isAdapt::Bool,isCSD::Bool,ΔE::Float64,
    𝚽E12::Array{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},
    𝒲::Array{Float64},isFC::Bool)

Compute the flux solution along one direction in 2D geometry.

# Input Argument(s)
- `𝚽l::Array{Float64,4}`: Legendre components of the in-cell flux.
- `Ql::Array{Float64,4}`: Legendre components of the in-cell source.
- `Σt::Vector{Float64}`: total cross-sections.
- `mat::Array{Int64,2}`: material identifier per voxel.
- `Ns::Vector{Int64}`: number of voxels along x- and y-axis.
- `Δs::Vector{Vector{Float64}}`: size of voxels along x- and y-axis.
- `Ω::Vector{Float64}`: direction cosines μ and η.
- `Mn::Vector{Float64}`: moment-to-discrete matrix.
- `Dn::Vector{Float64}`: discrete-to-moment matrix.
- `P::Int64`: number of angular interpolation basis.
- `Mnx⁻::Vector{Float64}`: moment-to-discrete matrix for surfaces along x-axis.
- `Dnx⁻::Vector{Float64}`: discrete-to-moment matrix for surfaces along x-axis.
- `Mny⁻::Vector{Float64}`: moment-to-discrete matrix for surfaces along y-axis.
- `Dny⁻::Vector{Float64}`: discrete-to-moment matrix for surfaces along y-axis.
- `Np_surf::Int64`: number of angular interpolation basis for surfaces.
- `𝒪::Vector{Int64}`: spatial and/or energy closure relation order.
- `Nm::Vector{Int64}`: number of spatial and/or energy moments.
- `C::Vector{Float64}`: constants related to the spatial and energy normalized
   Legendre expansion.
- `ω::Vector{Array{Float64}}`: weighting factors of the closure relations.
- `sources::Matrix{Union{Float64, Array{Float64}}}`: surface sources intensities.
- `isAdapt::Bool`: boolean for adaptive calculations.
- `isCSD::Bool`: boolean to indicate if continuous slowing-down term is treated in
   calculations.
- `ΔE::Float64`: energy group width.
- `𝚽E12::Array{Float64}`: incoming flux along the energy axis.
- `S⁻::Vector{Float64}`: stopping power at higher energy group boundary.
- `S⁺::Vector{Float64}`: stopping power at lower energy group boundary.
- `S::Array{Float64}`: stopping powers.
- `𝒲::Array{Float64}`: weighting constants.
- `isFC::Bool`: boolean indicating if the high-order incoming moments are fully coupled.

# Output Argument(s)
- `𝚽l::Array{Float64}`: Legendre components of the in-cell flux.
- `𝚽E12::Array{Float64}`: outgoing flux along the energy axis.

# Reference(s)
N/A

"""
function sn_sweep_2D(𝚽l::Array{Float64,4},Ql::Array{Float64,4},Σt::Vector{Float64},mat::Array{Int64,2},Ns::Vector{Int64},Δs::Vector{Vector{Float64}},Ω::Vector{Float64},Mn::Vector{Float64},Dn::Vector{Float64},P::Int64,Mnx⁻::Vector{Float64},Dnx⁻::Vector{Float64},Mny⁻::Vector{Float64},Dny⁻::Vector{Float64},Np_surf::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},C::Vector{Float64},ω::Vector{Array{Float64}},sources::Matrix{Union{Float64,Array{Float64}}},isAdapt::Bool,isCSD::Bool,ΔE::Float64,𝚽E12::Array{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},𝒲::Array{Float64},isFC::Bool,𝚽x12⁻,𝚽y12⁻,boundary_conditions,Np_source,pm_surf)

    # Initialization
    𝒪x = 𝒪[1]; 𝒪y = 𝒪[2]; 𝒪E = 𝒪[4]
    μ = Ω[1]; η = Ω[2]
    Δx = Δs[1]; Δy = Δs[2]
    Nx = Ns[1]; Ny = Ns[2]
    if (μ >= 0) x_sweep = (1:Nx) else x_sweep = (Nx:-1:1) end
    if (η >= 0) y_sweep = (1:Ny) else y_sweep = (Ny:-1:1) end
    𝚽x12⁺ = zeros(Np_surf,Nm[1],2,Ny)
    𝚽y12⁺ = zeros(Np_surf,Nm[2],2,Nx)

    # Sweep over x-axis
    𝚽x12 = zeros(Nm[1],Ny)
    for ix in x_sweep
        𝚽y12 = zeros(Nm[2])
        if η ≥ 0
            # Surface Y-
            for p in range(1,Np_source)
                𝚽y12[1] += Mny⁻[p] * sources[p,3][ix]
            end
            if boundary_conditions[3] != 0 # Not void
                for p in range(1,Np_surf), is in range(1,Nm[2])
                    if boundary_conditions[3] == 1 # Reflective
                        𝚽y12[is] += Mny⁻[p] * 𝚽y12⁻[p,is,1,ix] * (-1)^pm_surf[3][p]
                    elseif boundary_conditions[3] == 2 # Periodic
                        𝚽y12[is] += Mny⁻[p] * 𝚽y12⁻[p,is,2,ix] * (-1)^pm_surf[3][p]
                    end
                end
            end
        else
            # Surface Y+
            for p in range(1,Np_source)
                𝚽y12[1] += Mny⁻[p] * sources[p,4][ix]
            end
            if boundary_conditions[4] != 0 # Not void
                for p in range(1,Np_surf), is in range(1,Nm[2])
                    if boundary_conditions[4] == 1 # Reflective
                        𝚽y12[is] += Mny⁻[p] * 𝚽y12⁻[p,is,2,ix] * (-1)^pm_surf[4][p]
                    elseif boundary_conditions[4] == 2 # Periodic
                        𝚽y12[is] += Mny⁻[p] * 𝚽y12⁻[p,is,1,ix] * (-1)^pm_surf[4][p]
                    end
                end
            end
        end

        # Sweep over y-axis
        for iy in y_sweep
            if (ix == 1 && μ ≥ 0) || (ix == Nx && μ < 0 )
                if μ ≥ 0
                    # Surface X-
                    for p in range(1,Np_source)
                        𝚽x12[1,iy] += Mnx⁻[p] * sources[p,1][iy]
                    end
                    if boundary_conditions[1] != 0 # Not void
                        for p in range(1,Np_surf), is in range(1,Nm[1])
                            if boundary_conditions[1] == 1 # Reflective
                                𝚽x12[is,iy] += Mnx⁻[p] * 𝚽x12⁻[p,is,1,iy] * (-1)^pm_surf[1][p]
                            elseif boundary_conditions[1] == 2 # Periodic
                                𝚽x12[is,iy] += Mnx⁻[p] * 𝚽x12⁻[p,is,2,iy] * (-1)^pm_surf[1][p]
                            end
                        end
                    end
                else
                    # Surface X+
                    for p in range(1,Np_source)
                        𝚽x12[1,iy] += Mnx⁻[p] * sources[p,2][iy]
                    end
                    if boundary_conditions[2] != 0 # Not void
                        for p in range(1,Np_surf), is in range(1,Nm[1])
                            if boundary_conditions[2] == 1 # Reflective
                                𝚽x12[is,iy] += Mnx⁻[p] * 𝚽x12⁻[p,is,2,iy] * (-1)^pm_surf[2][p]
                            elseif boundary_conditions[2] == 2 # Periodic
                                𝚽x12[is,iy] += Mnx⁻[p] * 𝚽x12⁻[p,is,1,iy] * (-1)^pm_surf[2][p]
                            end
                        end
                    end
                end
            end

            # Source term
            Qn = zeros(Nm[5])
            for is in range(1,Nm[5]), p in range(1,P)
                Qn[is] += Mn[p] * Ql[p,is,ix,iy]
            end

            # Flux calculation
            if ~isCSD
                𝚽n,𝚽x12[:,iy],𝚽y12 = flux_2D_BTE(μ,η,Σt[mat[ix,iy]],Δx[ix],Δy[iy],Qn,𝚽x12[:,iy],𝚽y12,𝒪x,𝒪y,C,copy(ω[1]),copy(ω[2]),isAdapt,isFC)
            else
                𝚽n,𝚽x12[:,iy],𝚽y12,𝚽E12[:,ix,iy] = flux_2D_BFP(μ,η,Σt[mat[ix,iy]],S⁻[mat[ix,iy]],S⁺[mat[ix,iy]],S[mat[ix,iy],:],ΔE,Δx[ix],Δy[iy],Qn,𝚽x12[:,iy],𝚽y12,𝚽E12[:,ix,iy],𝒪E,𝒪x,𝒪y,C,copy(ω[1]),copy(ω[2]),copy(ω[3]),isAdapt,𝒲,isFC)
            end

            # Calculation of the Legendre components of the flux
            for is in range(1,Nm[5]), p in range(1,P)
                𝚽l[p,is,ix,iy] += Dn[p] * 𝚽n[is]
            end

            # Save boundary fluxes along x-axis
            if (ix == Nx && μ ≥ 0) || (ix == 1 && μ < 0 )
                for p in range(1,Np_surf)
                    for is in range(1,Nm[1])
                        # Surface X+
                        if μ ≥ 0
                            𝚽x12⁺[p,is,2,iy] += Dnx⁻[p] * 𝚽x12[is,iy]
                        # Surface X-
                        else
                            𝚽x12⁺[p,is,1,iy] += Dnx⁻[p] * 𝚽x12[is,iy]
                        end
                    end
                end
            end

            # Save boundary fluxes along y-axis
            if (iy == Ny && η ≥ 0) || (iy == 1 && η < 0 )
                for p in range(1,Np_surf)
                    for is in range(1,Nm[2])
                        # Surface Y+
                        if η ≥ 0
                            𝚽y12⁺[p,is,2,ix] += Dny⁻[p] * 𝚽y12[is]
                        # Surface Y-
                        else
                            𝚽y12⁺[p,is,1,ix] += Dny⁻[p] * 𝚽y12[is]
                        end
                    end
                end
            end
        end
    end
    return 𝚽l, 𝚽E12, 𝚽x12⁺, 𝚽y12⁺
end