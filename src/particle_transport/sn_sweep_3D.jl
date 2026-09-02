"""
    compute_sweep_3D(𝚽l::Array{Float64,5},Ql::Array{Float64,5},Σt::Vector{Float64},
    mat::Array{Int64,3},Ns::Vector{Int64},Δs::Vector{Vector{Float64}},Ω::Vector{Float64},
    Mn::Vector{Float64},Dn::Vector{Float64},P::Int64,Mnx⁻::Vector{Float64},
    Dnx⁻::Vector{Float64},Mny⁻::Vector{Float64},Dny⁻::Vector{Float64},Mnz⁻::Vector{Float64},
    Dnz⁻::Vector{Float64},Np_surf::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},
    C::Vector{Float64},ω::Vector{Array{Float64}},
    sources::Matrix{Union{Float64,Array{Float64}}},isAdapt::Bool,isCSD::Bool,ΔE::Float64,
    𝚽E12::Array{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},
    𝒲::Array{Float64},isFC::Bool)

Compute the flux solution along one direction in 3D geometry.

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
- `Mnz⁻::Vector{Float64}`: moment-to-discrete matrix for surfaces along z-axis.
- `Dnz⁻::Vector{Float64}`: discrete-to-moment matrix for surfaces along z-axis.
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
function sn_sweep_3D(𝚽l::Array{Float64,5},Ql::Array{Float64,5},Σt::Vector{Float64},mat::Array{Int64,3},Ns::Vector{Int64},Δs::Vector{Vector{Float64}},Ω::Vector{Float64},Mn::Vector{Float64},Dn::Vector{Float64},P::Int64,Mnx⁻::Vector{Float64},Dnx⁻::Vector{Float64},Mny⁻::Vector{Float64},Dny⁻::Vector{Float64},Mnz⁻::Vector{Float64},Dnz⁻::Vector{Float64},Np_surf::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},C::Vector{Float64},ω::Vector{Array{Float64}},sources::Matrix{Union{Float64,Array{Float64}}},isAdapt::Bool,isCSD::Bool,ΔE::Float64,𝚽E12::Array{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},𝒲::Array{Float64},isFC::Bool,𝚽x12⁻,𝚽y12⁻,𝚽z12⁻,boundary_conditions,Np_source)

    # Initialization
    𝒪x = 𝒪[1]; 𝒪y = 𝒪[2]; 𝒪z = 𝒪[3]; 𝒪E = 𝒪[4]
    μ = Ω[1]; η = Ω[2]; ξ = Ω[3]
    Δx = Δs[1]; Δy = Δs[2]; Δz = Δs[3]
    Nx = Ns[1]; Ny = Ns[2]; Nz = Ns[3]
    if (μ >= 0) x_sweep = (1:Nx) else x_sweep = (Nx:-1:1) end
    if (η >= 0) y_sweep = (1:Ny) else y_sweep = (Ny:-1:1) end
    if (ξ >= 0) z_sweep = (1:Nz) else z_sweep = (Nz:-1:1) end
    𝚽x12⁺ = zeros(Np_surf,Nm[1],2,Ny,Nz)
    𝚽y12⁺ = zeros(Np_surf,Nm[2],2,Nx,Nz)
    𝚽z12⁺ = zeros(Np_surf,Nm[3],2,Nx,Ny)

    # Sweeping over x-axis
    𝚽x12 = zeros(Nm[1],Ny,Nz)
    for ix in x_sweep

        # Sweeping over y-axis
        𝚽y12 = zeros(Nm[2],Nz)
        for iy in y_sweep
            𝚽z12 = zeros(Nm[3])
            if ξ ≥ 0
                # Surface Z-
                for p in range(1,Np_source)
                    𝚽z12[1] += Mnz⁻[p] * sources[p,5][ix,iy]  
                end
                if boundary_conditions[5] != 0 # Not void
                    for p in range(1,Np_surf), is in range(1,Nm[3])
                        if boundary_conditions[5] == 1 # Reflective
                            𝚽z12[is] += Mnz⁻[p] * 𝚽z12⁻[p,is,1,ix,iy]
                        elseif boundary_conditions[5] == 2 # Periodic
                            𝚽z12[is] += Mnz⁻[p] * 𝚽z12⁻[p,is,2,ix,iy]
                        end
                    end
                end
            else
                # Surface Z+
                for p in range(1,Np_source)
                    𝚽z12[1] += Mnz⁻[p] * sources[p,6][ix,iy]  
                end
                if boundary_conditions[6] != 0 # Not void
                    for p in range(1,Np_surf), is in range(1,Nm[3])
                        if boundary_conditions[6] == 1 # Reflective
                            𝚽z12[is] += Mnz⁻[p] * 𝚽z12⁻[p,is,2,ix,iy]
                        elseif boundary_conditions[6] == 2 # Periodic
                            𝚽z12[is] += Mnz⁻[p] * 𝚽z12⁻[p,is,1,ix,iy]
                        end
                    end
                end
            end

            # Sweeping over z-axis
            for iz in z_sweep
                if (iy == 1 && η ≥ 0) || (iy == Ny && η < 0 )
                    if η ≥ 0
                        # Surface Y-
                        for p in range(1,Np_source)
                            𝚽y12[1,iz] += Mny⁻[p] * sources[p,3][ix,iz]  
                        end
                        if boundary_conditions[3] != 0 # Not void
                            for p in range(1,Np_surf), is in range(1,Nm[2])
                                if boundary_conditions[3] == 1 # Reflective
                                    𝚽y12[is,iz] += Mny⁻[p] * 𝚽y12⁻[p,is,1,ix,iz]
                                elseif boundary_conditions[3] == 2 # Periodic
                                    𝚽y12[is,iz] += Mny⁻[p] * 𝚽y12⁻[p,is,2,ix,iz]
                                end
                            end
                        end
                    else
                        # Surface Y+
                        for p in range(1,Np_source)
                            𝚽y12[1,iz] += Mny⁻[p] * sources[p,4][ix,iz]  
                        end
                        if boundary_conditions[4] != 0 # Not void
                            for p in range(1,Np_surf), is in range(1,Nm[2])
                                if boundary_conditions[4] == 1 # Reflective
                                    𝚽y12[is,iz] += Mny⁻[p] * 𝚽y12⁻[p,is,2,ix,iz]
                                elseif boundary_conditions[4] == 2 # Periodic
                                    𝚽y12[is,iz] += Mny⁻[p] * 𝚽y12⁻[p,is,1,ix,iz]
                                end
                            end
                        end
                    end
                end
                if (ix == 1 && μ ≥ 0) || (ix == Nx && μ < 0 )
                    if μ ≥ 0
                        # Surface X-
                        for p in range(1,Np_source)
                            𝚽x12[1,iy,iz] += Mnx⁻[p] * sources[p,1][iy,iz]  
                        end
                        if boundary_conditions[1] != 0 # Not void
                            for p in range(1,Np_surf), is in range(1,Nm[1])
                                if boundary_conditions[1] == 1 # Reflective
                                    𝚽x12[is,iy,iz] += Mnx⁻[p] * 𝚽x12⁻[p,is,1,iy,iz]
                                elseif boundary_conditions[1] == 2 # Periodic
                                    𝚽x12[is,iy,iz] += Mnx⁻[p] * 𝚽x12⁻[p,is,2,iy,iz]
                                end
                            end
                        end
                    else
                        # Surface X+
                        for p in range(1,Np_source)
                            𝚽x12[1,iy,iz] += Mnx⁻[p] * sources[p,2][iy,iz]  
                        end
                        if boundary_conditions[2] != 0 # Not void
                            for p in range(1,Np_surf), is in range(1,Nm[1])
                                if boundary_conditions[2] == 1 # Reflective
                                    𝚽x12[is,iy,iz] += Mnx⁻[p] * 𝚽x12⁻[p,is,2,iy,iz]
                                elseif boundary_conditions[2] == 2 # Periodic
                                    𝚽x12[is,iy,iz] += Mnx⁻[p] * 𝚽x12⁻[p,is,1,iy,iz]
                                end
                            end
                        end
                    end
                end

                # Source term
                Qn = zeros(Nm[5])
                for is in range(1,Nm[5]), p in range(1,P)
                    Qn[is] += Mn[p] * Ql[p,is,ix,iy,iz]
                end

                # Flux calculation
                if ~isCSD
                    𝚽n,𝚽x12[:,iy,iz],𝚽y12[:,iz],𝚽z12 = flux_3D_BTE(μ,η,ξ,Σt[mat[ix,iy,iz]],Δx[ix],Δy[iy],Δz[iz],Qn,𝚽x12[:,iy,iz],𝚽y12[:,iz],𝚽z12,𝒪x,𝒪y,𝒪z,C,copy(ω[1]),copy(ω[2]),copy(ω[3]),isAdapt,isFC)
                else
                    𝚽n,𝚽x12[:,iy,iz],𝚽y12[:,iz],𝚽z12,𝚽E12[:,ix,iy,iz] = flux_3D_BFP(μ,η,ξ,Σt[mat[ix,iy,iz]],S⁻[mat[ix,iy,iz]],S⁺[mat[ix,iy,iz]],S[mat[ix,iy,iz],:],ΔE,Δx[ix],Δy[iy],Δz[iz],Qn,𝚽x12[:,iy,iz],𝚽y12[:,iz],𝚽z12,𝚽E12[:,ix,iy,iz],𝒪E,𝒪x,𝒪y,𝒪z,C,copy(ω[1]),copy(ω[2]),copy(ω[3]),copy(ω[4]),isAdapt,𝒲,isFC)
                end

                # Calculation of the Legendre components of the flux
                for is in range(1,Nm[5]), p in range(1,P)
                    𝚽l[p,is,ix,iy,iz] += Dn[p] * 𝚽n[is]
                end

                # Save boundary fluxes along x-axis
                if (ix == Nx && μ ≥ 0) || (ix == 1 && μ < 0 )
                    for p in range(1,Np_surf)
                        for is in range(1,Nm[1])
                            # Surface X+
                            if μ ≥ 0
                                𝚽x12⁺[p,is,2,iy,iz] += Dnx⁻[p] * 𝚽x12[is,iy,iz]
                            # Surface X-
                            else
                                𝚽x12⁺[p,is,1,iy,iz] += Dnx⁻[p] * 𝚽x12[is,iy,iz]
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
                                𝚽y12⁺[p,is,2,ix,iz] += Dny⁻[p] * 𝚽y12[is,iz]
                            # Surface Y-
                            else
                                𝚽y12⁺[p,is,1,ix,iz] += Dny⁻[p] * 𝚽y12[is,iz]
                            end
                        end
                    end
                end

                # Save boundary fluxes along z-axis
                if (iz == Nz && ξ ≥ 0) || (iz == 1 && ξ < 0 )
                    for p in range(1,Np_surf)
                        for is in range(1,Nm[3])
                            # Surface Z+
                            if ξ ≥ 0
                                𝚽z12⁺[p,is,2,ix,iy] += Dnz⁻[p] * 𝚽z12[is]
                            # Surface Z-
                            else
                                𝚽z12⁺[p,is,1,ix,iy] += Dnz⁻[p] * 𝚽z12[is]
                            end
                        end
                    end
                end
            end
        end
    end
    return 𝚽l, 𝚽E12, 𝚽x12⁺, 𝚽y12⁺, 𝚽z12⁺
end
"""
    sn_sweep_3D_fast!(𝚿,o𝚿,Ql,Mnn,Np,𝚽E12,𝚽E12o,oE,mat,Nx,Ny,Nz,srcx,srcy,srcz,outx,outy,outz,
    mom,cells,d,ws,Nmf,NmEf,isCSD,do_E,zero_E,save_out)

Optimized counterpart of `sn_sweep_3D`: sweep the spatial grid along one discrete ordinate.

Same sweep ordering and moving-boundary bookkeeping as the reference, with four changes.

The three `copy(ω[1])`, `copy(ω[2])`, `copy(ω[3])` the reference makes *in every cell* are
gone: they exist only because the adaptive scheme mutates the weights, and the optimized chain
does not serve that case (`sn_fast_applicable`).

The boundary slices `𝚽x12[:,iy,iz]` and `𝚽y12[:,iz]`, copied in and reassigned out at each
cell, are replaced by moving buffers addressed by offset — and the incoming half-range
transform is applied once per direction, ahead of the sweep, into `srcx`/`srcy`/`srcz` rather
than face by face inside it.

The outgoing face moments are written raw into `outx`/`outy`/`outz`; the reference's
`Dnx⁻[p]·𝚽x12[is]` accumulation, `Np_surf × Nm` scalar updates per boundary cell, becomes a
single rank-1 update per direction in the caller. Likewise the in-cell moments go to `𝚿`, one
value per moment, instead of `𝚽l[p,is,ix,iy,iz] += Dn[p]·𝚽n[is]` — `Np × Nm` scattered updates
per cell, in a 5-D array, for every direction of every pass.

The loop nesting is `iz` outer, `ix` inner, so the innermost index is the fastest one of the
arrays; the moving buffers follow suit — `bufz` spans the x-y plane, `bufy` the x line, `bufx`
is a single face. Either nesting satisfies the sweep's dependencies, since each cell is still
visited after its three upstream neighbours.
"""
function sn_sweep_3D_fast!(𝚿::Vector{Float64},o𝚿::Int64,Ql::Vector{Float64},Mnn::Vector{Float64},Np::Int64,𝚽E12::Vector{Float64},𝚽E12o::Vector{Float64},oE::Int64,mat::Array{Int64,3},Nx::Int64,Ny::Int64,Nz::Int64,srcx::Vector{Float64},srcy::Vector{Float64},srcz::Vector{Float64},outx::Vector{Float64},outy::Vector{Float64},outz::Vector{Float64},mom::GNFastMoments{NMOM},cells::GNFastCells,d::SNFastDir{NMOM},ws::GNFastWorkspace,Nmf::Vector{Int64},NmEf::Int64,isCSD::Bool,do_E::Bool,zero_E::Bool,save_out::Bool) where {NMOM}

    sx = d.sx; sy = d.sy; sz = d.sz
    Nm1 = Nmf[1]; Nm2 = Nmf[2]; Nm3 = Nmf[3]

    # Sweep ordering
    if sx > 0; x_sweep = 1:Nx else x_sweep = Nx:-1:1 end
    if sy > 0; y_sweep = 1:Ny else y_sweep = Ny:-1:1 end
    if sz > 0; z_sweep = 1:Nz else z_sweep = Nz:-1:1 end

    bufx = ws.bufx; bufy = ws.bufy; bufz = ws.bufz

    # Reset the moving z-boundary workspace (an x-y plane: z is the outer sweep axis)
    fill!(bufz, 0.0)

    @inbounds for iz in z_sweep
        fill!(bufy, 0.0)
        ikz = Int64(cells.kz[iz])
        is_z_entry = (iz == 1 && sz > 0) || (iz == Nz && sz < 0)
        is_z_exit  = save_out && ((iz == Nz && sz > 0) || (iz == 1 && sz < 0))

        # Z-boundary initialization (at the iz entrance plane only)
        if is_z_entry
            for i in 1:Nm3*Nx*Ny
                bufz[i] += srcz[i]
            end
        end

        # Y-boundary initialization (at the iy entrance line of this plane). Each ix slot is
        # read by the kernel of that same ix, so it can be filled ahead of the x sweep.
        for ix in 1:Nx
            b0 = Nm2*(ix-1); s0 = Nm2*((ix-1) + Nx*(iz-1))
            for i in 1:Nm2
                bufy[b0+i] += srcy[s0+i]
            end
        end

        for iy in y_sweep
            fill!(bufx, 0.0)
            iky = Int64(cells.ky[iy])
            is_y_exit = save_out && ((iy == Ny && sy > 0) || (iy == 1 && sy < 0))

            # X-boundary initialization (at the ix entrance face only)
            sxb = Nm1*((iy-1) + Ny*(iz-1))
            for i in 1:Nm1
                bufx[i] += srcx[sxb+i]
            end

            for ix in x_sweep
                ikx = Int64(cells.kx[ix])
                icell = (ix-1) + Nx*((iy-1) + Ny*(iz-1))
                oby = Nm2*(ix-1)
                obz = Nm3*((ix-1) + Nx*(iy-1))
                m = mat[ix,iy,iz]
                conf = gn_fast_conf(cells,m,ix,iy,iz)

                if ~isCSD
                    sn_3D_BTE_fast2!(𝚿,o𝚿+NMOM*icell,bufx,0,bufy,oby,bufz,obz,
                                     Ql,Np*NMOM*icell,Mnn,Np,mom,d,ws,conf,ikx,iky,ikz)
                else
                    sn_3D_BFP_fast!(𝚿,o𝚿+NMOM*icell,bufx,0,bufy,oby,bufz,obz,
                                    𝚽E12,𝚽E12o,oE+NmEf*icell,
                                    Ql,Np*NMOM*icell,Mnn,Np,mom,d,ws,conf,ikx,iky,ikz,m,do_E,zero_E)
                end

                # Save the z-outgoing boundary at the exit plane
                if is_z_exit
                    for i in 1:Nm3
                        outz[obz+i] = bufz[obz+i]
                    end
                end
            end

            # Save the y-outgoing boundary (end of the x sweep of this line)
            if is_y_exit
                for ix in 1:Nx
                    b0 = Nm2*(ix-1); s0 = Nm2*((ix-1) + Nx*(iz-1))
                    for i in 1:Nm2
                        outy[s0+i] = bufy[b0+i]
                    end
                end
            end

            # Save the x-outgoing boundary (end of the x sweep for this (iy,iz))
            if save_out
                for i in 1:Nm1
                    outx[sxb+i] = bufx[i]
                end
            end
        end
    end

    return nothing
end
