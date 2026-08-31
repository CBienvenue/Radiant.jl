"""
    sn_inner_pass!(...)

Apply one in-group source-iteration pass for the discrete-ordinates solver: build the source
from the current flux `𝚽l`, sweep all `Nd` ordinates, and update `𝚽l` and the incoming
boundary fluxes (`𝚽x12_in`, plus `𝚽y12_in`/`𝚽z12_in` in 2D/3D). This is the building block 
every in-group acceleration scheme is built on (see `sn_one_speed`).

The pass is the affine map `T(z) = A·z + c` over the state `z = (𝚽l, boundary fluxes)`:
- `homogeneous = false` gives `T(z)` (all fixed sources active);
- `homogeneous = true` gives the linear part `A·z`, by dropping the fixed sources 
   (`Qlout = 0`,`Np_source = 0`, zeroed incoming energy flux `𝚽E12`).

`𝚽E12_temp` receives the outgoing energy flux (when `isCSD`); `Ql` and `𝚽*12_temp` are
scratch. Remaining arguments mirror `sn_one_speed`.
"""
function sn_inner_pass!(𝚽l,𝚽x12_in,𝚽y12_in,𝚽z12_in,𝚽x12_temp,𝚽y12_temp,𝚽z12_temp,𝚽E12_temp,Ql,Qlout,𝚽E12_in,Σt,Σs,mat,ndims,Nd,Ns,Δs,Ω,Mn,Dn,Np,pl,Mn_surf,Dn_surf,Np_surf,n_to_n⁺,𝒪,Nm,isFC,C,ω,sources,isAdapt,isCSD,solver,ΔE,S⁻,S⁺,S,T,ℳ,is_EM,ℳ_EM,𝒲,boundary_conditions,Np_source;homogeneous::Bool)

    # Calculation of the Legendre components of the source (in-scattering)
    if homogeneous Ql .= 0.0 else Ql .= Qlout end
    if solver ∉ [4,5,6] Ql = scattering_source(Ql,𝚽l,Σs,mat,Np,pl,Nm[5],Ns) end

    # Finite element treatment of the angular Fokker-Planck term
    if solver ∈ [2,4] Ql = fokker_planck_source(Np,Nm[5],T,𝚽l,Ql,Ns,mat,ℳ) end

    # Electromagnetic source
    if is_EM
        for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3]), is in range(1,Nm[5])
            for p in range(1,Np), q in range(1,Np)
                Ql[p,is,ix,iy,iz] += ℳ_EM[p,q] * 𝚽l[q,is,ix,iy,iz]
            end
        end
    end

    # Fixed-source switches: drop the external surface source and incoming energy flux for the
    # homogeneous operator (out-of-group source already dropped above via Ql .= 0).
    Np_source_eff = homogeneous ? 0 : Np_source
    𝚽E12_eff = (isCSD && homogeneous) ? zero(𝚽E12_in) : 𝚽E12_in

    #----
    # Loop over all discrete ordinates
    #----
    𝚽l .= 0
    for n in range(1,Nd)
        if isCSD 𝚽E12_out = 𝚽E12_eff[n,:,:,:,:] else 𝚽E12_out = Array{Float64}(undef) end
        if ndims == 1
            nx⁻ = n_to_n⁺[1][n]
            nx⁺ = n_to_n⁺[2][n]
            if nx⁻ != 0
                Mnx⁻ = Mn_surf[1][nx⁻,:]
                Dnx⁻ = Dn_surf[1][:,nx⁻]
            elseif nx⁺ != 0
                Mnx⁻ = Mn_surf[2][nx⁺,:]
                Dnx⁻ = Dn_surf[2][:,nx⁺]
            else
                Mnx⁻ = zeros(Np)
                Dnx⁻ = zeros(Np)
            end
            𝚽l[:,:,:,1,1], 𝚽E12_out,𝚽x12_out = sn_sweep_1D(𝚽l[:,:,:,1,1],Ql[:,:,:,1,1],Σt,mat[:,1,1],Ns[1],Δs[1],Ω[1][n],Mn[n,:],Dn[:,n],Np,Mnx⁻,Dnx⁻,Np_surf,𝒪,Nm,C,ω,sources,isAdapt,isCSD,ΔE,𝚽E12_out,S⁻,S⁺,S,𝒲,isFC,𝚽x12_in,boundary_conditions,Np_source_eff)
        elseif ndims == 2
            nx⁻ = n_to_n⁺[1][n]
            nx⁺ = n_to_n⁺[2][n]
            ny⁻ = n_to_n⁺[3][n]
            ny⁺ = n_to_n⁺[4][n]
            if nx⁻ != 0
                Mnx⁻ = Mn_surf[1][nx⁻,:]
                Dnx⁻ = Dn_surf[1][:,nx⁻]
            elseif nx⁺ != 0
                Mnx⁻ = Mn_surf[2][nx⁺,:]
                Dnx⁻ = Dn_surf[2][:,nx⁺]
            else
                Mnx⁻ = zeros(Np)
                Dnx⁻ = zeros(Np)
            end
            if ny⁻ != 0
                Mny⁻ = Mn_surf[3][ny⁻,:]
                Dny⁻ = Dn_surf[3][:,ny⁻]
            elseif ny⁺ != 0
                Mny⁻ = Mn_surf[4][ny⁺,:]
                Dny⁻ = Dn_surf[4][:,ny⁺]
            else
                Mny⁻ = zeros(Np)
                Dny⁻ = zeros(Np)
            end
            𝚽l[:,:,:,:,1],𝚽E12_out,𝚽x12_out,𝚽y12_out = sn_sweep_2D(𝚽l[:,:,:,:,1],Ql[:,:,:,:,1],Σt,mat[:,:,1],Ns[1:2],Δs[1:2],[Ω[1][n],Ω[2][n]],Mn[n,:],Dn[:,n],Np,Mnx⁻,Dnx⁻,Mny⁻,Dny⁻,Np_surf,𝒪,Nm,C,ω,sources,isAdapt,isCSD,ΔE,𝚽E12_out,S⁻,S⁺,S,𝒲,isFC,𝚽x12_in,𝚽y12_in,boundary_conditions,Np_source_eff)
        elseif ndims == 3
            nx⁻ = n_to_n⁺[1][n]
            nx⁺ = n_to_n⁺[2][n]
            ny⁻ = n_to_n⁺[3][n]
            ny⁺ = n_to_n⁺[4][n]
            nz⁻ = n_to_n⁺[5][n]
            nz⁺ = n_to_n⁺[6][n]
            if nx⁻ != 0
                Mnx⁻ = Mn_surf[1][nx⁻,:]
                Dnx⁻ = Dn_surf[1][:,nx⁻]
            elseif nx⁺ != 0
                Mnx⁻ = Mn_surf[2][nx⁺,:]
                Dnx⁻ = Dn_surf[2][:,nx⁺]
            else
                Mnx⁻ = zeros(Np)
                Dnx⁻ = zeros(Np)
            end
            if ny⁻ != 0
                Mny⁻ = Mn_surf[3][ny⁻,:]
                Dny⁻ = Dn_surf[3][:,ny⁻]
            elseif ny⁺ != 0
                Mny⁻ = Mn_surf[4][ny⁺,:]
                Dny⁻ = Dn_surf[4][:,ny⁺]
            else
                Mny⁻ = zeros(Np)
                Dny⁻ = zeros(Np)
            end
            if nz⁻ != 0
                Mnz⁻ = Mn_surf[5][nz⁻,:]
                Dnz⁻ = Dn_surf[5][:,nz⁻]
            elseif nz⁺ != 0
                Mnz⁻ = Mn_surf[6][nz⁺,:]
                Dnz⁻ = Dn_surf[6][:,nz⁺]
            else
                Mnz⁻ = zeros(Np)
                Dnz⁻ = zeros(Np)
            end
            𝚽l,𝚽E12_out,𝚽x12_out,𝚽y12_out,𝚽z12_out = sn_sweep_3D(𝚽l,Ql,Σt,mat,Ns,Δs,[Ω[1][n],Ω[2][n],Ω[3][n]],Mn[n,:],Dn[:,n],Np,Mnx⁻,Dnx⁻,Mny⁻,Dny⁻,Mnz⁻,Dnz⁻,Np_surf,𝒪,Nm,C,ω,sources,isAdapt,isCSD,ΔE,𝚽E12_out,S⁻,S⁺,S,𝒲,isFC,𝚽x12_in,𝚽y12_in,𝚽z12_in,boundary_conditions,Np_source_eff)
        else
            error("Dimension is not 1, 2 or 3.")
        end
        𝚽x12_temp .+= 𝚽x12_out
        if ndims >= 2 𝚽y12_temp .+= 𝚽y12_out end
        if ndims >= 3 𝚽z12_temp .+= 𝚽z12_out end
        if isCSD 𝚽E12_temp[n,:,:,:,:] = 𝚽E12_out end
    end
    𝚽x12_in .= 𝚽x12_temp; 𝚽x12_temp .= 0.0
    if ndims >= 2 𝚽y12_in .= 𝚽y12_temp; 𝚽y12_temp .= 0.0 end
    if ndims >= 3 𝚽z12_in .= 𝚽z12_temp; 𝚽z12_temp .= 0.0 end

    return 𝚽l
end

"""
    _sn_fast_faces!(src,Mh,n,Np_surf,Nmfa,Nfca,𝚽12_in,slot,bc,sfc,Np_source,use_bc)

Build the incoming face moments of one axis for one direction.

The reference applies this half-range transform *inside* the sweep, face cell by face cell,
re-testing the boundary condition at every one of them. It is a per-direction quantity, so it
is done once here, ahead of the sweep. The arithmetic and the accumulation order are the
reference's.
"""
function _sn_fast_faces!(src::Vector{Float64},Mh::Matrix{Float64},n::Int64,Np_surf::Int64,Nmfa::Int64,Nfca::Int64,𝚽12_in::Array{Float64},slot::Int64,bc::Int64,sfc::Array{Float64,3},Np_source::Int64,use_bc::Bool)

    fill!(src,0.0)
    @inbounds begin
        if Np_source > 0
            S = reshape(sfc,size(sfc,1),Nfca)
            for c in 1:Nfca
                a = 0.0
                for p in 1:Np_source; a += Mh[p,n] * S[p,c] end
                src[1+Nmfa*(c-1)] += a
            end
        end
        if use_bc && bc != 0
            A = reshape(𝚽12_in,Np_surf,Nmfa,2,Nfca)
            for c in 1:Nfca, is in 1:Nmfa
                a = 0.0
                for p in 1:Np_surf; a += Mh[p,n] * A[p,is,slot,c] end
                src[is+Nmfa*(c-1)] += a
            end
        end
    end
    return nothing
end

"""
    _sn_fast_faces_out!(𝚽12_temp,Dh,n,Np_surf,Nmfa,Nfca,out,slot)

Accumulate one axis' outgoing face moments of one direction into the pass total.

The reference does this as `Np_surf × Nmfa` scalar updates inside the sweep, at every cell of
the exit face. Here the sweep leaves the raw face moments in `out` and the transform is a
rank-1 update — outside the threaded region, since the target is shared by every direction.
"""
function _sn_fast_faces_out!(𝚽12_temp::Array{Float64},Dh::Matrix{Float64},n::Int64,Np_surf::Int64,Nmfa::Int64,Nfca::Int64,out::Vector{Float64},slot::Int64)

    @inbounds begin
        A = reshape(𝚽12_temp,Np_surf,Nmfa,2,Nfca)
        for c in 1:Nfca, is in 1:Nmfa
            v = out[is+Nmfa*(c-1)]
            if v == 0.0 continue end
            for p in 1:Np_surf; A[p,is,slot,c] += Dh[p,n] * v end
        end
    end
    return nothing
end

"""
    sn_inner_pass_fast!(...)

Optimized counterpart of `sn_inner_pass!`: one in-group source-iteration pass of the
discrete-ordinates solver.

Same affine map `T(z) = A·z + c` over the same state, with the same `homogeneous` switch, so
it is interchangeable with the reference pass inside every acceleration scheme.

What changes is how the directions are traversed. The reference sweeps them one at a time,
slicing `Mn[n,:]`, `Dn[:,n]`, `Mn_surf[…]`, `Dn_surf[…]` and — when `isCSD` — copying the whole
incoming energy flux `𝚽E12_eff[n,:,:,:,:]` for each one, then accumulating the result into
`𝚽l` cell by cell. Here the directions are swept in blocks: the block's sweeps run
concurrently, each into its own stripe of `𝚿`, and the discrete-to-moment transform of the
whole block is a single `gemm`. The half-range transforms and the surface-source arrays come
precomputed from `sn_fast_scratch`.

`need_boundary_flux` is false when every boundary is void, in which case the boundary fluxes
are identically zero and neither the incoming transform nor the outgoing accumulation is
performed at all.

`reconstruct` marks the final pass on the converged state — the only one whose outgoing energy
flux is read — so the energy closure runs only there. The incoming energy flux is left
untouched throughout, which is what makes that possible; the reference has to work on a copy
because its kernel overwrites it.
"""
function sn_inner_pass_fast!(𝚽l::Array{Float64},𝚽x12_in::Array{Float64},𝚽y12_in::Array{Float64},𝚽z12_in::Array{Float64},𝚽x12_temp::Array{Float64},𝚽y12_temp::Array{Float64},𝚽z12_temp::Array{Float64},𝚽E12_temp::Array{Float64},Ql::Array{Float64},Qlout::Array{Float64},𝚽E12p::Vector{Float64},𝚽E12o::Vector{Float64},Σs::Array{Float64},mat::Array{Int64,3},ndims::Int64,Nd::Int64,Ns::Vector{Int64},Mn::Array{Float64,2},Dn::Array{Float64,2},Np::Int64,pl::Vector{Int64},Np_surf::Int64,Np_source::Int64,Nm::Vector{Int64},isCSD::Bool,solver::Int64,T::Vector{Float64},ℳ::Array{Float64},is_EM::Bool,ℳ_EM::Array{Float64},boundary_conditions::Vector{Int64},ctx::SNFastContext{NMOM},sc::SNFastScratch,need_boundary_flux::Bool;homogeneous::Bool,reconstruct::Bool) where {NMOM}

    Nx,Ny,Nz = Ns[1],Ns[2],Ns[3]
    Nvox = Nx*Ny*Nz
    NmEf = isCSD ? Nm[4] : 1

    # Calculation of the Legendre components of the source (in-scattering)
    if homogeneous Ql .= 0.0 else Ql .= Qlout end
    if solver ∉ [4,5,6] scattering_source_fast(Ql,𝚽l,Σs,mat,Np,pl,Nm[5],Ns) end

    # Finite element treatment of the angular Fokker-Planck term
    if solver ∈ [2,4] fokker_planck_source_fast(Np,Nm[5],T,𝚽l,Ql,Ns,mat,ℳ) end

    # Electromagnetic source
    if is_EM
        @inbounds for ix in range(1,Nx), iy in range(1,Ny), iz in range(1,Nz), is in range(1,Nm[5])
            for p in range(1,Np), q in range(1,Np)
                Ql[p,is,ix,iy,iz] += ℳ_EM[p,q] * 𝚽l[q,is,ix,iy,iz]
            end
        end
    end

    # Fixed-source switches, as in the reference pass.
    Np_source_eff = homogeneous ? 0 : Np_source
    zero_E = isCSD && homogeneous
    do_E   = isCSD && reconstruct
    if do_E copyto!(𝚽E12o,𝚽E12p) end

    Qlf = vec(Ql)
    𝚿m  = reshape(sc.𝚿, sc.M, sc.Nblk)
    𝚽lm = reshape(𝚽l, Np, NMOM*Nvox)
    𝚽l .= 0

    # Incoming/outgoing boundary slots per axis, from the direction's octant. The reference
    # picks the reflective slot on the axis' own side and the periodic one on the opposite.
    Nmf = sc.Nmf; Nfc = sc.Nfc

    n0 = 1
    while n0 ≤ Nd
        nb = min(sc.Nblk, Nd-n0+1)

        Threads.@threads for j in 1:nb
            n = n0+j-1
            d = ctx.dir[n]
            @inbounds for p in 1:Np; sc.Mnn[j][p] = Mn[n,p] end

            fin_x = d.sx > 0 ? 1 : 2
            bcx   = boundary_conditions[fin_x]
            slx   = bcx == 1 ? (d.sx > 0 ? 1 : 2) : (d.sx > 0 ? 2 : 1)
            _sn_fast_faces!(sc.srcx[j],sc.Mx,n,Np_surf,Nmf[1],Nfc[1],𝚽x12_in,slx,bcx,
                            sc.sfc[fin_x],Np_source_eff,need_boundary_flux)
            if ndims ≥ 2
                fin_y = d.sy > 0 ? 3 : 4
                bcy   = boundary_conditions[fin_y]
                sly   = bcy == 1 ? (d.sy > 0 ? 1 : 2) : (d.sy > 0 ? 2 : 1)
                _sn_fast_faces!(sc.srcy[j],sc.My,n,Np_surf,Nmf[2],Nfc[2],𝚽y12_in,sly,bcy,
                                sc.sfc[fin_y],Np_source_eff,need_boundary_flux)
            end
            if ndims ≥ 3
                fin_z = d.sz > 0 ? 5 : 6
                bcz   = boundary_conditions[fin_z]
                slz   = bcz == 1 ? (d.sz > 0 ? 1 : 2) : (d.sz > 0 ? 2 : 1)
                _sn_fast_faces!(sc.srcz[j],sc.Mz,n,Np_surf,Nmf[3],Nfc[3],𝚽z12_in,slz,bcz,
                                sc.sfc[fin_z],Np_source_eff,need_boundary_flux)
            end

            oE = isCSD ? NmEf*Nvox*(n-1) : 0
            if ndims == 3
                sn_sweep_3D_fast!(sc.𝚿,sc.M*(j-1),Qlf,sc.Mnn[j],Np,𝚽E12p,𝚽E12o,oE,mat,Nx,Ny,Nz,
                                  sc.srcx[j],sc.srcy[j],sc.srcz[j],sc.outx[j],sc.outy[j],sc.outz[j],
                                  ctx.mom,ctx.cells,ctx.dir[n],sc.ws[j],Nmf,NmEf,isCSD,do_E,zero_E,
                                  need_boundary_flux)
            elseif ndims == 2
                sn_sweep_2D_fast!(sc.𝚿,sc.M*(j-1),Qlf,sc.Mnn[j],Np,𝚽E12p,𝚽E12o,oE,mat,Nx,Ny,
                                  sc.srcx[j],sc.srcy[j],sc.outx[j],sc.outy[j],
                                  ctx.mom,ctx.cells,ctx.dir[n],sc.ws[j],Nmf,NmEf,isCSD,do_E,zero_E,
                                  need_boundary_flux)
            elseif ndims == 1
                sn_sweep_1D_fast!(sc.𝚿,sc.M*(j-1),Qlf,sc.Mnn[j],Np,𝚽E12p,𝚽E12o,oE,mat,Nx,
                                  sc.srcx[j],sc.outx[j],
                                  ctx.mom,ctx.cells,ctx.dir[n],sc.ws[j],Nmf,NmEf,isCSD,do_E,zero_E,
                                  need_boundary_flux)
            else
                error("Dimension is not 1, 2 or 3.")
            end
        end

        # Discrete-to-moment transform of the whole block, in one gemm
        @views mul!(𝚽lm, Dn[:,n0:n0+nb-1], transpose(𝚿m[:,1:nb]), 1.0, 1.0)

        # Outgoing boundary fluxes, accumulated outside the threaded region
        if need_boundary_flux
            for j in 1:nb
                n = n0+j-1; d = ctx.dir[n]
                _sn_fast_faces_out!(𝚽x12_temp,sc.Dx,n,Np_surf,Nmf[1],Nfc[1],sc.outx[j],d.sx > 0 ? 2 : 1)
                if ndims ≥ 2
                    _sn_fast_faces_out!(𝚽y12_temp,sc.Dy,n,Np_surf,Nmf[2],Nfc[2],sc.outy[j],d.sy > 0 ? 2 : 1)
                end
                if ndims ≥ 3
                    _sn_fast_faces_out!(𝚽z12_temp,sc.Dz,n,Np_surf,Nmf[3],Nfc[3],sc.outz[j],d.sz > 0 ? 2 : 1)
                end
            end
        end

        n0 += nb
    end

    𝚽x12_in .= 𝚽x12_temp; 𝚽x12_temp .= 0.0
    if ndims >= 2 𝚽y12_in .= 𝚽y12_temp; 𝚽y12_temp .= 0.0 end
    if ndims >= 3 𝚽z12_in .= 𝚽z12_temp; 𝚽z12_temp .= 0.0 end

    # The outgoing energy flux is stored per direction and per cell; the reference's
    # 𝚽E12_temp puts the direction first, so the final pass writes it back in that layout.
    if do_E
        @inbounds for n in 1:Nd, c in 1:Nvox, is in 1:NmEf
            𝚽E12_temp[n + Nd*((is-1) + NmEf*(c-1))] = 𝚽E12o[is + NmEf*((c-1) + Nvox*(n-1))]
        end
    end

    return 𝚽l
end
