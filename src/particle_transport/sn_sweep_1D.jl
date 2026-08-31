"""
    compute_sweep_1D(𝚽l::Array{Float64,3},Ql::Array{Float64,3},Σt::Vector{Float64},
    mat::Vector{Int64},Nx::Int64,Δx::Vector{Float64},μ::Float64,Mn::Vector{Float64},
    Dn::Vector{Float64},Np::Int64,Mnx⁻::Vector{Float64},Dnx⁻::Vector{Float64},
    Np_surf::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},C::Vector{Float64},
    ω::Vector{Array{Float64}},sources::Matrix{Union{Float64,Array{Float64}}},isAdapt::Bool,
    isCSD::Bool,ΔE::Float64,𝚽E12::Array{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64},
    S::Array{Float64},𝒲::Array{Float64},isFC::Bool)

Compute the flux solution along one direction in 1D geometry.

# Input Argument(s)
- `𝚽l::Array{Float64,3}`: Legendre components of the in-cell flux.
- `Ql::Array{Float64,3}`: Legendre components of the in-cell source.
- `Σt::Vector{Float64}`: total cross-sections.
- `mat::Vector{Int64}`: material identifier per voxel.
- `Nx::Int64`: number of voxels along x-axis.
- `Δx::Vector{Float64}`: size of voxels along x-axis.
- `μ::Float64`: direction cosines.
- `Mn::Vector{Float64}`: moment-to-discrete matrix.
- `Dn::Vector{Float64}`: discrete-to-moment matrix.
- `Np::Int64`: number of angular interpolation basis.
- `Mnx⁻::Vector{Float64}`: moment-to-discrete matrix for surfaces along x-axis.
- `Dnx⁻::Vector{Float64}`: discrete-to-moment matrix for surfaces along x-axis.
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
- `S⁻::Vector{Float64}`: stopping powers at higher energy group boundary.
- `S⁺::Vector{Float64}`: stopping powers at lower energy group boundary.
- `S::Array{Float64}`: stopping powers.
- `𝒲::Array{Float64}`: weighting constants.
- `isFC::Bool`: boolean indicating if the high-order incoming moments are fully coupled.

# Output Argument(s)
- `𝚽l::Array{Float64}`: Legendre components of the in-cell flux.
- `𝚽E12::Array{Float64}`: outgoing flux along the energy axis.

# Reference(s)
N/A

"""
function sn_sweep_1D(𝚽l::Array{Float64,3},Ql::Array{Float64,3},Σt::Vector{Float64},mat::Vector{Int64},Nx::Int64,Δx::Vector{Float64},μ::Float64,Mn::Vector{Float64},Dn::Vector{Float64},Np::Int64,Mnx⁻::Vector{Float64},Dnx⁻::Vector{Float64},Np_surf::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},C::Vector{Float64},ω::Vector{Array{Float64}},sources::Matrix{Union{Float64, Array{Float64}}},isAdapt::Bool,isCSD::Bool,ΔE::Float64,𝚽E12::Array{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},𝒲::Array{Float64},isFC::Bool,𝚽x12⁻::Array{Float64,3},boundary_conditions::Vector{Int64},Np_source)

    # Initialization
    𝒪x = 𝒪[1]
    𝒪E = 𝒪[4]
    𝚽x12 = zeros(Nm[1])
    𝚽x12⁺ = zeros(Np_surf,Nm[1],2)

    # Boundary conditions and sources
    if μ ≥ 0
        x_sweep = 1:Nx
        # Surface X-
        for p in range(1,Np_source)
            𝚽x12[1] += Mnx⁻[p] * sources[p,1]
        end
        if boundary_conditions[1] != 0 # Not void
            for p in range(1,Np_surf), is in range(1,Nm[1])
                if boundary_conditions[1] == 1 # Reflective
                    𝚽x12[is] += Mnx⁻[p] * 𝚽x12⁻[p,is,1]
                elseif boundary_conditions[1] == 2 # Periodic
                    𝚽x12[is] += Mnx⁻[p] * 𝚽x12⁻[p,is,2]
                end
            end
        end
    else
        x_sweep = Nx:-1:1
        # Surface X+
        for p in range(1,Np_source)
            𝚽x12[1] += Mnx⁻[p] * sources[p,2]
        end
        if boundary_conditions[2] != 0 # Not void
            for p in range(1,Np_surf), is in range(1,Nm[1])
                if boundary_conditions[2] == 1 # Reflective
                    𝚽x12[is] += Mnx⁻[p] * 𝚽x12⁻[p,is,2]
                elseif boundary_conditions[2] == 2 # Periodic
                    𝚽x12[is] += Mnx⁻[p] * 𝚽x12⁻[p,is,1]
                end
            end
        end
    end

    for ix in x_sweep

        # Source term
        Qn = zeros(Nm[5])
        for is in range(1,Nm[5]), p in range(1,Np)
            Qn[is] += Mn[p] * Ql[p,is,ix]
        end

        # Flux calculation
        if ~isCSD
            𝚽n,𝚽x12 = flux_1D_BTE(μ,Σt[mat[ix]],Δx[ix],Qn,𝚽x12[1],𝒪x,C,copy(ω[1]),isAdapt)
        else
            𝚽n,𝚽x12,𝚽E12[:,ix] = flux_1D_BFP(μ,Σt[mat[ix]],Δx[ix],Qn,𝚽x12,S⁻[mat[ix]],S⁺[mat[ix]],S[mat[ix],:],ΔE,𝚽E12[:,ix],𝒪E,𝒪x,C,copy(ω[1]),copy(ω[2]),isAdapt,𝒲,isFC)
        end

        # Calculation of the Legendre components of the flux
        for is in range(1,Nm[5]), p in range(1,Np)
            𝚽l[p,is,ix] += Dn[p] * 𝚽n[is]
        end
        
    end

    # Save boundary fluxes
    for p in range(1,Np_surf), is in range(1,Nm[1])
        if μ ≥ 0 # Surface X+
            𝚽x12⁺[p,is,2] += Dnx⁻[p] * 𝚽x12[is]
        else # Surface X-
            𝚽x12⁺[p,is,1] += Dnx⁻[p] * 𝚽x12[is]
        end
    end

    return 𝚽l, 𝚽E12, 𝚽x12⁺
end
"""
    sn_sweep_1D_fast!(𝚿,o𝚿,Ql,Mnn,Np,𝚽E12,𝚽E12o,oE,mat,Nx,srcx,outx,mom,cells,d,ws,Nmf,NmEf,
    isCSD,do_E,zero_E,save_out)

Optimized counterpart of `sn_sweep_1D`: sweep the spatial grid along one discrete ordinate.

Identical in structure to `sn_sweep_3D_fast!` with the y and z axes collapsed; see it for the
four changes to the reference.
"""
function sn_sweep_1D_fast!(𝚿::Vector{Float64},o𝚿::Int64,Ql::Vector{Float64},Mnn::Vector{Float64},Np::Int64,𝚽E12::Vector{Float64},𝚽E12o::Vector{Float64},oE::Int64,mat::Array{Int64,3},Nx::Int64,srcx::Vector{Float64},outx::Vector{Float64},mom::GNFastMoments{NMOM},cells::GNFastCells,d::SNFastDir{NMOM},ws::GNFastWorkspace,Nmf::Vector{Int64},NmEf::Int64,isCSD::Bool,do_E::Bool,zero_E::Bool,save_out::Bool) where {NMOM}

    Nm1 = Nmf[1]
    if d.sx > 0; x_sweep = 1:Nx else x_sweep = Nx:-1:1 end

    bufx = ws.bufx
    fill!(bufx, 0.0)
    @inbounds for i in 1:Nm1
        bufx[i] += srcx[i]
    end

    @inbounds for ix in x_sweep
        ikx = Int64(cells.kx[ix])
        icell = ix-1
        m = mat[ix,1,1]
        conf = gn_fast_conf(cells,m,ix,1,1)

        if ~isCSD
            sn_1D_BTE_fast!(𝚿,o𝚿+NMOM*icell,bufx,0,Ql,Np*NMOM*icell,Mnn,Np,mom,d,ws,conf,ikx)
        else
            sn_1D_BFP_fast!(𝚿,o𝚿+NMOM*icell,bufx,0,𝚽E12,𝚽E12o,oE+NmEf*icell,
                            Ql,Np*NMOM*icell,Mnn,Np,mom,d,ws,conf,ikx,m,do_E,zero_E)
        end
    end

    if save_out
        @inbounds for i in 1:Nm1
            outx[i] = bufx[i]
        end
    end

    return nothing
end
