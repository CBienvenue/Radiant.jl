"""
    flux_3D_BFP(μ::Float64,η::Float64,ξ::Float64,Σt::Float64,S⁻::Float64,S⁺::Float64,
    S::Vector{Float64},ΔE::Float64,Δx::Float64,Δy::Float64,Δz::Float64,Qn::Vector{Float64},
    𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},𝚽z12::Vector{Float64},
    𝚽E12::Vector{Float64},𝒪E::Int64,𝒪x::Int64,𝒪y::Int64,𝒪z::Int64,C::Vector{Float64},
    ωE::Array{Float64},ωx::Array{Float64},ωy::Array{Float64},ωz::Array{Float64},
    isAdapt::Bool,𝒲::Array{Float64},isFC::Bool)

Compute flux solution in a cell in 3D Cartesian geometry for the Boltzmann Fokker-Planck
equation.

# Input Argument(s)
- `μ::Float64`: direction cosine.
- `η::Float64`: direction cosine.
- `ξ::Float64`: direction cosine.
- `Σt::Float64`: total cross-sections.
- `Δx::Float64`: size of voxels along x-axis.
- `Δy::Float64`: size of voxels along y-axis.
- `Δz::Float64`: size of voxels along z-axis.
- `Qn::Vector{Float64}`: angular in-cell source.
- `𝚽x12::Vector{Float64}`: incoming angular flux along x-axis.
- `𝚽y12::Vector{Float64}`: incoming angular flux along y-axis.
- `𝚽z12::Vector{Float64}`: incoming angular flux along z-axis.
- `S⁻::Float64`: restricted stopping power at upper energy group boundary.
- `S⁺::Float64`: restricted stopping power at lower energy group boundary.
- `ΔE::Float64`: energy group width.
- `𝚽E12::Vector{Float64}`: incoming angular flux along E-axis.
- `𝒪E::Int64`: energy closure relation order.
- `𝒪x::Int64`: spatial closure relation order.
- `𝒪y::Int64`: spatial closure relation order.
- `𝒪z::Int64`: spatial closure relation order.
- `C::Vector{Float64}`: constants related to normalized Legendre.
- `ωE::Array{Float64}`: weighting factors of the E-axis scheme.
- `ωx::Array{Float64}`: weighting factors of the x-axis scheme.
- `ωy::Array{Float64}`: weighting factors of the y-axis scheme.
- `ωz::Array{Float64}`: weighting factors of the z-axis scheme.
- `isAdapt::Bool`: boolean for adaptive calculations.
- `𝒲::Array{Float64}` : weighting constants.
- `isFC::Bool`: boolean indicating if the high-order incoming moments are fully coupled.

# Output Argument(s)
- `𝚽n::Vector{Float64}`: angular in-cell flux.
- `𝚽x12::Vector{Float64}`: outgoing angular flux along x-axis.
- `𝚽y12::Vector{Float64}`: outgoing angular flux along y-axis.
- `𝚽z12::Vector{Float64}`: outgoing angular flux along z-axis.
- `𝚽E12::Vector{Float64}`: outgoing angular flux along E-axis.

# Reference(s)
N/A

"""
function flux_3D_BFP(μ::Float64,η::Float64,ξ::Float64,Σt::Float64,S⁻::Float64,S⁺::Float64,S::Vector{Float64},ΔE::Float64,Δx::Float64,Δy::Float64,Δz::Float64,Qn::Vector{Float64},𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},𝚽z12::Vector{Float64},𝚽E12::Vector{Float64},𝒪E::Int64,𝒪x::Int64,𝒪y::Int64,𝒪z::Int64,C::Vector{Float64},ωE::Array{Float64},ωx::Array{Float64},ωy::Array{Float64},ωz::Array{Float64},isAdapt::Bool,𝒲::Array{Float64},isFC::Bool)

# Initialization
sx = sign(μ)
sy = sign(η)
sz = sign(ξ)
hx = abs(μ)/Δx
hy = abs(η)/Δy
hz = abs(ξ)/Δz
if isFC Nm = 𝒪x*𝒪y*𝒪z*𝒪E else Nm = 𝒪E+𝒪x+𝒪y+𝒪z-3 end
𝒮 = zeros(Nm,Nm)
Q = zeros(Nm)
𝚽n = Q

# Adaptive weight calculations
if isAdapt ωx,ωy,ωz,ωE = adaptive(𝒪x,𝒪y,𝒪z,𝒪E,ωx,ωy,ωz,ωE,hx,hy,hz,1/ΔE,sx,sy,sz,-1,𝚽x12,𝚽y12,𝚽z12,𝚽E12,Qn,Σt,isFC) end

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,𝒪x), jx in range(1,𝒪x), iy in range(1,𝒪y), jy in range(1,𝒪y), iz in range(1,𝒪z), jz in range(1,𝒪z), iE in range(1,𝒪E), jE in range(1,𝒪E)
    if isFC
        i = 𝒪y*𝒪x*𝒪E*(iz-1) + 𝒪x*𝒪E * (iy-1) + 𝒪E * (ix-1) + iE
        j = 𝒪y*𝒪x*𝒪E*(jz-1) + 𝒪x*𝒪E * (jy-1) + 𝒪E * (jx-1) + jE
    else
        if count(>(1),(iE,ix,iy,iz)) ≥ 2 || count(>(1),(jE,jx,jy,jz)) ≥ 2 continue end
        i = 1 + (iE-1) + (ix-1) + (iy-1) + (iz-1)
        j = 1 + (jE-1) + (jx-1) + (jy-1) + (jz-1)
        if ix > 1 i += 𝒪E-1 end
        if iy > 1 i += 𝒪E-1 + 𝒪x-1 end
        if iz > 1 i += 𝒪E-1 + 𝒪x-1 + 𝒪y-1 end
        if jx > 1 j += 𝒪E-1 end
        if jy > 1 j += 𝒪E-1 + 𝒪x-1 end
        if jz > 1 j += 𝒪E-1 + 𝒪x-1 + 𝒪y-1 end
    end

    # Collision term
    if (i == j) 𝒮[i,j] += Σt end

    # Streaming term - x
    if iy == jy && iz == jz && iE == jE
        if (ix ≥ jx + 1) 𝒮[i,j] -= C[ix] * hx * sx * C[jx] * (1-(-1)^(ix-jx)) end
        𝒮[i,j] += C[ix] * hx * sx^(ix-1) * C[jx] * sx^(jx-1) * ωx[jx+1,jy,jz,jE]
    end

    # Streaming term - y
    if ix == jx && iz == jz && iE == jE
        if (iy ≥ jy + 1) 𝒮[i,j] -= C[iy] * hy * sy * C[jy] * (1-(-1)^(iy-jy)) end 
        𝒮[i,j] += C[iy] * hy * sy^(iy-1) * C[jy] * sy^(jy-1) * ωy[jy+1,jx,jz,jE]
    end

    # Streaming term - z
    if ix == jx && iy == jy && iE == jE
        if (iz ≥ jz + 1) 𝒮[i,j] -= C[iz] * hz * sz * C[jz] * (1-(-1)^(iz-jz)) end 
        𝒮[i,j] += C[iz] * hz * sz^(iz-1) * C[jz] * sz^(jz-1) * ωz[jz+1,jx,jy,jE]
    end

    # CSD term
    if ix == jx && iy == jy && iz == jz
        for kE in range(1,iE-1), wE in range(1,𝒪E)
            𝒮[i,j] += C[iE] * C[jE] * C[kE] * C[wE] * (1-(-1)^(iE-kE)) * S[wE] * 𝒲[jE,kE,wE]
        end
        𝒮[i,j] += C[iE] * S⁺ * (-1)^(iE-1) * C[jE] * (-1)^(jE-1) * ωE[jE+1,jx,jy,jz]
    end
end

# Source vector
for jx in range(1,𝒪x), jy in range(1,𝒪y), jz in range(1,𝒪z), jE in range(1,𝒪E)
    if isFC
        j = 𝒪y*𝒪x*𝒪E*(jz-1) + 𝒪x*𝒪E * (jy-1) + 𝒪E * (jx-1) + jE
        jEm = 𝒪y*𝒪x*(jz-1)+𝒪x*(jy-1)+jx
        jxm = 𝒪y*𝒪E*(jz-1)+𝒪E*(jy-1)+jE
        jym = 𝒪x*𝒪E*(jz-1)+𝒪E*(jx-1)+jE
        jzm = 𝒪x*𝒪E*(jy-1)+𝒪E*(jx-1)+jE
    else
        if count(>(1),(jE,jx,jy,jz)) ≥ 2 continue end
        j = 1 + (jE-1) + (jx-1) + (jy-1) + (jz-1)
        jEm = 1 + (jx-1) + (jy-1) + (jz-1)
        jxm = 1 + (jE-1) + (jy-1) + (jz-1)
        jym = 1 + (jE-1) + (jx-1) + (jz-1)
        jzm = 1 + (jE-1) + (jx-1) + (jy-1)
        if jx > 1 j += 𝒪E-1 end
        if jy > 1 j += 𝒪E-1 + 𝒪x-1 end
        if jz > 1 j += 𝒪E-1 + 𝒪x-1 + 𝒪y-1 end
        if jy > 1 jEm += 𝒪x-1 end
        if jz > 1 jEm += 𝒪x-1 + 𝒪y-1 end
        if jy > 1 jxm += 𝒪E-1 end
        if jz > 1 jxm += 𝒪E-1 + 𝒪y-1 end
        if jx > 1 jym += 𝒪E-1 end
        if jz > 1 jym += 𝒪E-1 + 𝒪x-1 end
        if jx > 1 jzm += 𝒪E-1 end
        if jy > 1 jzm += 𝒪E-1 + 𝒪x-1 end
    end
    Q[j] += Qn[j]
    Q[j] -= C[jx] * hx * (sx^(jx-1) * ωx[1,jy,jz,jE] - (-sx)^(jx-1)) * 𝚽x12[jxm] 
    Q[j] -= C[jy] * hy * (sy^(jy-1) * ωy[1,jx,jz,jE] - (-sy)^(jy-1)) * 𝚽y12[jym]
    Q[j] -= C[jz] * hz * (sz^(jz-1) * ωz[1,jx,jy,jE] - (-sz)^(jz-1)) * 𝚽z12[jzm]  
    Q[j] -= C[jE] * ((-1)^(jE-1)*S⁺*ωE[1,jx,jy,jz] - S⁻) * 𝚽E12[jEm]
end

# Solve the equation system
𝚽n = 𝒮\Q

# Closure relation
for jx in range(1,𝒪x), jy in range(1,𝒪y), jz in range(1,𝒪z), jE in range(1,𝒪E)
    if isFC
        j = 𝒪y*𝒪x*𝒪E*(jz-1) + 𝒪x*𝒪E * (jy-1) + 𝒪E * (jx-1) + jE
        jEm = 𝒪y*𝒪x*(jz-1)+𝒪x*(jy-1)+jx
        jxm = 𝒪y*𝒪E*(jz-1)+𝒪E*(jy-1)+jE
        jym = 𝒪x*𝒪E*(jz-1)+𝒪E*(jx-1)+jE
        jzm = 𝒪x*𝒪E*(jy-1)+𝒪E*(jx-1)+jE
    else
        if count(>(1),(jE,jx,jy,jz)) ≥ 2 continue end
        j = 1 + (jE-1) + (jx-1) + (jy-1) + (jz-1)
        jEm = 1 + (jx-1) + (jy-1) + (jz-1)
        jxm = 1 + (jE-1) + (jy-1) + (jz-1)
        jym = 1 + (jE-1) + (jx-1) + (jz-1)
        jzm = 1 + (jE-1) + (jx-1) + (jy-1)
        if jx > 1 j += 𝒪E-1 end
        if jy > 1 j += 𝒪E-1 + 𝒪x-1 end
        if jz > 1 j += 𝒪E-1 + 𝒪x-1 + 𝒪y-1 end
        if jy > 1 jEm += 𝒪x-1 end
        if jz > 1 jEm += 𝒪x-1 + 𝒪y-1 end
        if jy > 1 jxm += 𝒪E-1 end
        if jz > 1 jxm += 𝒪E-1 + 𝒪y-1 end
        if jx > 1 jym += 𝒪E-1 end
        if jz > 1 jym += 𝒪E-1 + 𝒪x-1 end
        if jx > 1 jzm += 𝒪E-1 end
        if jy > 1 jzm += 𝒪E-1 + 𝒪x-1 end
    end
    if (jE == 1) 𝚽E12[jEm] = ωE[1,jx,jy,jz] * 𝚽E12[jEm] end
    if (jx == 1) 𝚽x12[jxm] = ωx[1,jy,jz,jE] * 𝚽x12[jxm] end
    if (jy == 1) 𝚽y12[jym] = ωy[1,jx,jz,jE] * 𝚽y12[jym] end
    if (jz == 1) 𝚽z12[jzm] = ωz[1,jx,jy,jE] * 𝚽z12[jzm] end
    𝚽x12[jxm] += C[jx] * sx^(jx-1) * ωx[jx+1,jy,jz,jE] * 𝚽n[j]
    𝚽y12[jym] += C[jy] * sy^(jy-1) * ωy[jy+1,jx,jz,jE] * 𝚽n[j]
    𝚽z12[jzm] += C[jz] * sz^(jz-1) * ωz[jz+1,jx,jy,jE] * 𝚽n[j]
    𝚽E12[jEm] += C[jE] * (-1)^(jE-1) * ωE[jE+1,jx,jy,jz] * 𝚽n[j]
end

# Returning solutions
return 𝚽n, 𝚽x12, 𝚽y12, 𝚽z12, 𝚽E12
end
"""
    sn_3D_BFP_matrix!(𝒮,μ,η,ξ,Σt,S⁺,S,ΔE,Δx,Δy,Δz,𝒪E,𝒪x,𝒪y,𝒪z,C,ωE,ωx,ωy,ωz,𝒲,isFC)

Assemble the cell matrix `𝒮` of the 3D BFP discrete-ordinates kernel, for the optimized chain.

The body is the assembly loop of `flux_3D_BFP` transcribed unchanged, so the matrix is
bit-identical to the one the reference kernel builds. It is split out because the matrix does
not depend on the voxel beyond its widths and its material: `sn_fast_context_3D` calls this
once per (material, mesh-width triple, direction) and factorizes the result, where the
reference rebuilds it — and refactorizes it, through `𝒮\\Q` — in every voxel of every sweep.

Valid only when `~isAdapt`: the adaptive scheme recomputes `ω` per cell from the local flux,
which makes the matrix voxel-dependent. See `sn_fast_applicable`.
"""
function sn_3D_BFP_matrix!(𝒮::Matrix{Float64},μ::Float64,η::Float64,ξ::Float64,Σt::Float64,S⁺::Float64,S::Vector{Float64},ΔE::Float64,Δx::Float64,Δy::Float64,Δz::Float64,𝒪E::Int64,𝒪x::Int64,𝒪y::Int64,𝒪z::Int64,C::Vector{Float64},ωE::Array{Float64},ωx::Array{Float64},ωy::Array{Float64},ωz::Array{Float64},𝒲::Array{Float64},isFC::Bool)

fill!(𝒮,0.0)
sx = sign(μ)
sy = sign(η)
sz = sign(ξ)
hx = abs(μ)/Δx
hy = abs(η)/Δy
hz = abs(ξ)/Δz

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,𝒪x), jx in range(1,𝒪x), iy in range(1,𝒪y), jy in range(1,𝒪y), iz in range(1,𝒪z), jz in range(1,𝒪z), iE in range(1,𝒪E), jE in range(1,𝒪E)
    if isFC
        i = 𝒪y*𝒪x*𝒪E*(iz-1) + 𝒪x*𝒪E * (iy-1) + 𝒪E * (ix-1) + iE
        j = 𝒪y*𝒪x*𝒪E*(jz-1) + 𝒪x*𝒪E * (jy-1) + 𝒪E * (jx-1) + jE
    else
        if count(>(1),(iE,ix,iy,iz)) ≥ 2 || count(>(1),(jE,jx,jy,jz)) ≥ 2 continue end
        i = 1 + (iE-1) + (ix-1) + (iy-1) + (iz-1)
        j = 1 + (jE-1) + (jx-1) + (jy-1) + (jz-1)
        if ix > 1 i += 𝒪E-1 end
        if iy > 1 i += 𝒪E-1 + 𝒪x-1 end
        if iz > 1 i += 𝒪E-1 + 𝒪x-1 + 𝒪y-1 end
        if jx > 1 j += 𝒪E-1 end
        if jy > 1 j += 𝒪E-1 + 𝒪x-1 end
        if jz > 1 j += 𝒪E-1 + 𝒪x-1 + 𝒪y-1 end
    end

    # Collision term
    if (i == j) 𝒮[i,j] += Σt end

    # Streaming term - x
    if iy == jy && iz == jz && iE == jE
        if (ix ≥ jx + 1) 𝒮[i,j] -= C[ix] * hx * sx * C[jx] * (1-(-1)^(ix-jx)) end
        𝒮[i,j] += C[ix] * hx * sx^(ix-1) * C[jx] * sx^(jx-1) * ωx[jx+1,jy,jz,jE]
    end

    # Streaming term - y
    if ix == jx && iz == jz && iE == jE
        if (iy ≥ jy + 1) 𝒮[i,j] -= C[iy] * hy * sy * C[jy] * (1-(-1)^(iy-jy)) end 
        𝒮[i,j] += C[iy] * hy * sy^(iy-1) * C[jy] * sy^(jy-1) * ωy[jy+1,jx,jz,jE]
    end

    # Streaming term - z
    if ix == jx && iy == jy && iE == jE
        if (iz ≥ jz + 1) 𝒮[i,j] -= C[iz] * hz * sz * C[jz] * (1-(-1)^(iz-jz)) end 
        𝒮[i,j] += C[iz] * hz * sz^(iz-1) * C[jz] * sz^(jz-1) * ωz[jz+1,jx,jy,jE]
    end

    # CSD term
    if ix == jx && iy == jy && iz == jz
        for kE in range(1,iE-1), wE in range(1,𝒪E)
            𝒮[i,j] += C[iE] * C[jE] * C[kE] * C[wE] * (1-(-1)^(iE-kE)) * S[wE] * 𝒲[jE,kE,wE]
        end
        𝒮[i,j] += C[iE] * S⁺ * (-1)^(iE-1) * C[jE] * (-1)^(jE-1) * ωE[jE+1,jx,jy,jz]
    end
end

return nothing

end

"""
    sn_3D_BFP_fast!(𝚿,o𝚿,𝚽x12,ox,𝚽y12,oy,𝚽z12,oz,𝚽E12,𝚽E12o,oE,Ql,oQ,Mnn,Np,mom,d,ws,
    conf,ikx,iky,ikz,m,do_E,zero_E)

Optimized counterpart of `flux_3D_BFP`: solve one voxel of the 3D BFP equation.

Same arithmetic as the reference, with everything voxel-independent taken out. The reference
allocates `𝒮`, `Q` and the solution in every cell, re-derives the index maps and the assembly
coefficients there, and finishes with `𝚽n = 𝒮\\Q` — a full factorization per cell, per
direction, per pass. Here the arrays are flat with a per-voxel offset, `conf` selects a cached
factorization (`sn_fast_context`), and nothing is allocated at all.

The energy axis differs from the space axes in that its incoming flux is not a moving sweep
buffer but a per-cell array, and the reference *overwrites* it with the outgoing value — which
is why it works on a copy. `𝚽E12` (read) and `𝚽E12o` (written) are kept separate here so the
incoming flux survives the iteration and no copy is needed. `do_E` skips the outgoing energy
closure altogether on the ordinary passes, where nothing reads it; `zero_E` drops the incoming
energy flux, as the homogeneous operator requires.
"""
@inline function sn_3D_BFP_fast!(𝚿::Vector{Float64},o𝚿::Int64,𝚽x12::Vector{Float64},ox::Int64,𝚽y12::Vector{Float64},oy::Int64,𝚽z12::Vector{Float64},oz::Int64,𝚽E12::Vector{Float64},𝚽E12o::Vector{Float64},oE::Int64,Ql::Vector{Float64},oQ::Int64,Mnn::Vector{Float64},Np::Int64,mom::GNFastMoments{NMOM},d::SNFastDir{NMOM},ws::GNFastWorkspace,conf::Int64,ikx::Int64,iky::Int64,ikz::Int64,m::Int64,do_E::Bool,zero_E::Bool) where {NMOM}

    Q = ws.Q
    𝚽 = ws.𝚽

    # Source vector
    @inbounds for k in 1:NMOM
        col = mom.col[k]
        q = 0.0
        b = oQ + Np*(col-1)
        @simd for p in 1:Np
            q += Mnn[p] * Ql[b+p]
        end
        q += d.qcx[k,ikx] * 𝚽x12[ox+mom.cyz[k]] +
             d.qcy[k,iky] * 𝚽y12[oy+mom.cxz[k]] +
             d.qcz[k,ikz] * 𝚽z12[oz+mom.cxy[k]]
        if ~zero_E q += d.qcE[k,m] * 𝚽E12[oE+mom.cxyz[k]] end
        Q[col] = q
    end

    # Solve the equation system from the cached factorization
    gn_fast_solve!(𝚽,d.LU,d.ipiv,Q,NMOM,conf)

    # Closure relations
    if do_E gn_fast_closure!(𝚽E12o,oE,d.clE,𝚽,Val(NMOM),Val(1)) end
    gn_fast_closure!(𝚽x12,ox,d.clx,𝚽,Val(NMOM),Val(1))
    gn_fast_closure!(𝚽y12,oy,d.cly,𝚽,Val(NMOM),Val(1))
    gn_fast_closure!(𝚽z12,oz,d.clz,𝚽,Val(NMOM),Val(1))
    @inbounds for c in 1:NMOM
        𝚿[o𝚿+c] = 𝚽[c]
    end

    return nothing
end
