"""
    flux_3D_BTE(μ::Float64,η::Float64,ξ::Float64,Σt::Float64,Δx::Float64,Δy::Float64,
    Δz::Float64,Qn::Vector{Float64},𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},
    𝚽z12::Vector{Float64},𝒪x::Int64,𝒪y::Int64,𝒪z::Int64,C::Vector{Float64},
    ωx::Array{Float64},ωy::Array{Float64},ωz::Array{Float64},isAdapt::Bool,isFC::Bool)

Compute flux solution in a cell in 3D Cartesian geometry for the Boltzmann transport
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
- `𝒪x::Int64`: spatial closure relation order.
- `𝒪y::Int64`: spatial closure relation order.
- `𝒪z::Int64`: spatial closure relation order.
- `C::Vector{Float64}`: constants related to normalized Legendre.
- `ωx::Array{Float64}`: weighting factors of the x-axis scheme.
- `ωy::Array{Float64}`: weighting factors of the y-axis scheme.
- `ωz::Array{Float64}`: weighting factors of the z-axis scheme.
- `isAdapt::Bool`: boolean for adaptive calculations.
- `isFC::Bool`: boolean indicating if the high-order incoming moments are fully coupled.

# Output Argument(s)
- `𝚽n::Vector{Float64}`: angular in-cell flux.
- `𝚽x12::Vector{Float64}`: outgoing angular flux along x-axis.
- `𝚽y12::Vector{Float64}`: outgoing angular flux along y-axis.
- `𝚽z12::Vector{Float64}`: outgoing angular flux along z-axis.

# Reference(s)
N/A

"""
function flux_3D_BTE(μ::Float64,η::Float64,ξ::Float64,Σt::Float64,Δx::Float64,Δy::Float64,Δz::Float64,Qn::Vector{Float64},𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},𝚽z12::Vector{Float64},𝒪x::Int64,𝒪y::Int64,𝒪z::Int64,C::Vector{Float64},ωx::Array{Float64},ωy::Array{Float64},ωz::Array{Float64},isAdapt::Bool,isFC::Bool)

# Initialization
sx = sign(μ)
sy = sign(η)
sz = sign(ξ)
hx = abs(μ)/Δx
hy = abs(η)/Δy
hz = abs(ξ)/Δz
if isFC Nm = 𝒪x*𝒪y*𝒪z else Nm = 𝒪x+𝒪y+𝒪z-2 end
𝒮 = zeros(Nm,Nm)
Q = zeros(Nm)
𝚽n = Q

# Adaptive weight calculations
if isAdapt ωx,ωy,ωz = adaptive(𝒪x,𝒪y,𝒪z,ωx,ωy,ωz,hx,hy,hz,sx,sy,sz,𝚽x12,𝚽y12,𝚽z12,Qn,Σt,isFC) end

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,𝒪x), jx in range(1,𝒪x), iy in range(1,𝒪y), jy in range(1,𝒪y), iz in range(1,𝒪z), jz in range(1,𝒪z)
    if isFC
        i = 𝒪y*𝒪x*(iz-1) + 𝒪x * (iy-1) + ix
        j = 𝒪y*𝒪x*(jz-1) + 𝒪x * (jy-1) + jx
    else
        if count(>(1),(ix,iy,iz)) ≥ 2 || count(>(1),(jx,jy,jz)) ≥ 2 continue end
        i = 1 + (ix-1) + (iy-1) + (iz-1)
        j = 1 + (jx-1) + (jy-1) + (jz-1)
        if iy > 1 i += 𝒪x-1 end
        if iz > 1 i += 𝒪x-1 + 𝒪y-1 end
        if jy > 1 j += 𝒪x-1 end
        if jz > 1 j += 𝒪x-1 + 𝒪y-1 end
    end

    # Collision term
    if (i == j) 𝒮[i,j] += Σt end

    # Streaming term - x
    if iy == jy && iz == jz
        if (ix ≥ jx + 1) 𝒮[i,j] -= C[ix] * hx * sx * C[jx] * (1-(-1)^(ix-jx)) end
        𝒮[i,j] += C[ix] * hx * sx^(ix-1) * C[jx] * sx^(jx-1) * ωx[jx+1,jy,jz]
    end

    # Streaming term - y
    if ix == jx && iz == jz
        if (iy ≥ jy + 1) 𝒮[i,j] -= C[iy] * hy * sy * C[jy] * (1-(-1)^(iy-jy)) end 
        𝒮[i,j] += C[iy] * hy * sy^(iy-1) * C[jy] * sy^(jy-1) * ωy[jy+1,jx,jz]
    end

    # Streaming term - z
    if ix == jx && iy == jy
        if (iz ≥ jz + 1) 𝒮[i,j] -= C[iz] * hz * sz * C[jz] * (1-(-1)^(iz-jz)) end 
        𝒮[i,j] += C[iz] * hz * sz^(iz-1) * C[jz] * sz^(jz-1) * ωz[jz+1,jx,jy]
    end
end

# Source vector
for jx in range(1,𝒪x), jy in range(1,𝒪y), jz in range(1,𝒪z)
    if isFC
        j = 𝒪y*𝒪x*(jz-1) + 𝒪x * (jy-1) + jx
        jxm = 𝒪y*(jz-1) + jy
        jym = 𝒪x*(jz-1) + jx
        jzm = 𝒪x*(jy-1) + jx
    else
        if count(>(1),(jx,jy,jz)) ≥ 2 continue end
        j = 1 + (jx-1) + (jy-1) + (jz-1)
        jxm = 1 + (jy-1) + (jz-1)
        jym = 1 + (jx-1) + (jz-1)
        jzm = 1 + (jx-1) + (jy-1)
        if jy > 1 j += 𝒪x-1 end
        if jz > 1 j += 𝒪x-1 + 𝒪y-1 end
        if jz > 1 jxm += 𝒪y-1 end
        if jz > 1 jym += 𝒪x-1 end
        if jy > 1 jzm += 𝒪x-1 end
    end
    Q[j] = Qn[j]
    Q[j] -= C[jx] * hx * (sx^(jx-1) * ωx[1,jy,jz] - (-sx)^(jx-1)) * 𝚽x12[jxm] 
    Q[j] -= C[jy] * hy * (sy^(jy-1) * ωy[1,jx,jz] - (-sy)^(jy-1)) * 𝚽y12[jym] 
    Q[j] -= C[jz] * hz * (sz^(jz-1) * ωz[1,jx,jy] - (-sz)^(jz-1)) * 𝚽z12[jzm]
end

# Solve the equation system
𝚽n = 𝒮\Q

# Closure relation
for jx in range(1,𝒪x), jy in range(1,𝒪y), jz in range(1,𝒪z)
    if isFC
        j = 𝒪y*𝒪x*(jz-1) + 𝒪x * (jy-1) + jx
        jxm = 𝒪y*(jz-1) + jy
        jym = 𝒪x*(jz-1) + jx
        jzm = 𝒪x*(jy-1) + jx
    else
        if count(>(1),(jx,jy,jz)) ≥ 2 continue end
        j = 1 + (jx-1) + (jy-1) + (jz-1)
        jxm = 1 + (jy-1) + (jz-1)
        jym = 1 + (jx-1) + (jz-1)
        jzm = 1 + (jx-1) + (jy-1)
        if jy > 1 j += 𝒪x-1 end
        if jz > 1 j += 𝒪x-1 + 𝒪y-1 end
        if jz > 1 jxm += 𝒪y-1 end
        if jz > 1 jym += 𝒪x-1 end
        if jy > 1 jzm += 𝒪x-1 end
    end
    if (jx == 1) 𝚽x12[jxm] = ωx[1,jy,jz] * 𝚽x12[jxm] end
    if (jy == 1) 𝚽y12[jym] = ωy[1,jx,jz] * 𝚽y12[jym] end
    if (jz == 1) 𝚽z12[jzm] = ωz[1,jx,jy] * 𝚽z12[jzm] end
    𝚽x12[jxm] += C[jx] * sx^(jx-1) * ωx[jx+1,jy,jz] * 𝚽n[j]
    𝚽y12[jym] += C[jy] * sy^(jy-1) * ωy[jy+1,jx,jz] * 𝚽n[j]
    𝚽z12[jzm] += C[jz] * sz^(jz-1) * ωz[jz+1,jx,jy] * 𝚽n[j]
end

# Returning solutions
return 𝚽n, 𝚽x12, 𝚽y12, 𝚽z12
end
"""
    sn_3D_BTE_matrix!(𝒮,μ,η,ξ,Σt,Δx,Δy,Δz,𝒪x,𝒪y,𝒪z,C,ωx,ωy,ωz,isFC)

Assemble the cell matrix `𝒮` of the 3D SN BTE kernel, for the optimized solver chain.

The body is the assembly loop of `flux_3D_BTE` transcribed unchanged, so the matrix is
bit-identical to the one the reference kernel builds. It is split out because the matrix
does not depend on the voxel beyond its widths: it is a function of
`(Σt, |μ|/Δx, |η|/Δy, |ξ|/Δz, direction signs, C, ω)` only. The optimized chain therefore
assembles and factorizes it once per (material, mesh-width triple, direction) instead of
once per voxel per direction per pass — the reference builds *and* solves it in every cell,
allocating `𝒮`, `Q` and the solution each time.

**Only valid for the non-adaptive schemes.** With `isAdapt` the reference calls `adaptive()`
to recompute `ω` from the *local* flux, which makes the matrix voxel-dependent and the cache
wrong; the SN fast path falls back to the reference chain in that case.

See `set_fast_path`.
"""
function sn_3D_BTE_matrix!(𝒮::Matrix{Float64},μ::Float64,η::Float64,ξ::Float64,Σt::Float64,Δx::Float64,Δy::Float64,Δz::Float64,𝒪x::Int64,𝒪y::Int64,𝒪z::Int64,C::Vector{Float64},ωx::Array{Float64},ωy::Array{Float64},ωz::Array{Float64},isFC::Bool)

fill!(𝒮,0.0)
sx = sign(μ); sy = sign(η); sz = sign(ξ)
hx = abs(μ)/Δx; hy = abs(η)/Δy; hz = abs(ξ)/Δz

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,𝒪x), jx in range(1,𝒪x), iy in range(1,𝒪y), jy in range(1,𝒪y), iz in range(1,𝒪z), jz in range(1,𝒪z)
    if isFC
        i = 𝒪y*𝒪x*(iz-1) + 𝒪x * (iy-1) + ix
        j = 𝒪y*𝒪x*(jz-1) + 𝒪x * (jy-1) + jx
    else
        if count(>(1),(ix,iy,iz)) ≥ 2 || count(>(1),(jx,jy,jz)) ≥ 2 continue end
        i = 1 + (ix-1) + (iy-1) + (iz-1)
        j = 1 + (jx-1) + (jy-1) + (jz-1)
        if iy > 1 i += 𝒪x-1 end
        if iz > 1 i += 𝒪x-1 + 𝒪y-1 end
        if jy > 1 j += 𝒪x-1 end
        if jz > 1 j += 𝒪x-1 + 𝒪y-1 end
    end

    # Collision term
    if (i == j) 𝒮[i,j] += Σt end

    # Streaming term - x
    if iy == jy && iz == jz
        if (ix ≥ jx + 1) 𝒮[i,j] -= C[ix] * hx * sx * C[jx] * (1-(-1)^(ix-jx)) end
        𝒮[i,j] += C[ix] * hx * sx^(ix-1) * C[jx] * sx^(jx-1) * ωx[jx+1,jy,jz]
    end

    # Streaming term - y
    if ix == jx && iz == jz
        if (iy ≥ jy + 1) 𝒮[i,j] -= C[iy] * hy * sy * C[jy] * (1-(-1)^(iy-jy)) end
        𝒮[i,j] += C[iy] * hy * sy^(iy-1) * C[jy] * sy^(jy-1) * ωy[jy+1,jx,jz]
    end

    # Streaming term - z
    if ix == jx && iy == jy
        if (iz ≥ jz + 1) 𝒮[i,j] -= C[iz] * hz * sz * C[jz] * (1-(-1)^(iz-jz)) end
        𝒮[i,j] += C[iz] * hz * sz^(iz-1) * C[jz] * sz^(jz-1) * ωz[jz+1,jx,jy]
    end
end

return nothing

end

"""
    sn_3D_BTE_fast!(𝚽n,o𝚽n,𝚽x12,ox,𝚽y12,oy,𝚽z12,oz,Qn,oQ,mom,d,ws,conf,ikx,iky,ikz)

Optimized counterpart of `flux_3D_BTE`: solve one voxel of the 3D SN BTE.

Same arithmetic, with the per-cell work reduced to what actually depends on the voxel. The
reference allocates `𝒮`, `Q` and the solution in every cell and factorizes a fresh matrix
each time (`𝚽n = 𝒮\\Q`); here `conf` selects a cached factorization, the arrays are flat with
a per-voxel offset, and nothing is allocated at all.

SN solves one direction at a time, so the cell system is the moment system alone — there is
no angular block, and the face fluxes carry one value per face moment.

The right-hand side and the closures accumulate in the reference order, so the result is
bit-identical up to the triangular solves (`gn_fast_solve!`).
"""
@inline function sn_3D_BTE_fast!(𝚽n::Vector{Float64},o𝚽n::Int64,𝚽x12::Vector{Float64},ox::Int64,𝚽y12::Vector{Float64},oy::Int64,𝚽z12::Vector{Float64},oz::Int64,Qn::Vector{Float64},oQ::Int64,mom::GNFastMoments{NMOM},d::SNFastDir{NMOM},ws::GNFastWorkspace,conf::Int64,ikx::Int64,iky::Int64,ikz::Int64) where {NMOM}

    Q = ws.Q
    𝚽 = ws.𝚽

    # Source vector
    @inbounds for k in 1:NMOM
        col = mom.col[k]
        Q[col] = Qn[oQ+col] +
                 d.qcx[k,ikx] * 𝚽x12[ox+mom.cyz[k]] +
                 d.qcy[k,iky] * 𝚽y12[oy+mom.cxz[k]] +
                 d.qcz[k,ikz] * 𝚽z12[oz+mom.cxy[k]]
    end

    # Solve the equation system from the cached factorization
    gn_fast_solve!(𝚽,d.LU,d.ipiv,Q,NMOM,conf)

    # Closure relations
    gn_fast_closure!(𝚽x12,ox,d.clx,𝚽,Val(NMOM),Val(1))
    gn_fast_closure!(𝚽y12,oy,d.cly,𝚽,Val(NMOM),Val(1))
    gn_fast_closure!(𝚽z12,oz,d.clz,𝚽,Val(NMOM),Val(1))
    @inbounds for c in 1:NMOM
        𝚽n[o𝚽n+c] = 𝚽[c]
    end

    return nothing
end

"""
    sn_3D_BTE_fast2!(𝚿,o𝚿,𝚽x12,ox,𝚽y12,oy,𝚽z12,oz,Ql,oQ,Mnn,Np,mom,d,ws,conf,ikx,iky,ikz)

`sn_3D_BTE_fast!` with the moment-to-discrete transform of the source folded in.

The reference builds `Qn` in the sweep, one `zeros(Nm[5])` allocation per cell followed by
`Nm[5]×P` scalar updates through a 5-D index. Here `Ql` arrives flat, the transform is a small
contiguous gemv accumulated straight into the right-hand side, and nothing is allocated.
"""
@inline function sn_3D_BTE_fast2!(𝚿::Vector{Float64},o𝚿::Int64,𝚽x12::Vector{Float64},ox::Int64,𝚽y12::Vector{Float64},oy::Int64,𝚽z12::Vector{Float64},oz::Int64,Ql::Vector{Float64},oQ::Int64,Mnn::Vector{Float64},Np::Int64,mom::GNFastMoments{NMOM},d::SNFastDir{NMOM},ws::GNFastWorkspace,conf::Int64,ikx::Int64,iky::Int64,ikz::Int64) where {NMOM}

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
        Q[col] = q +
                 d.qcx[k,ikx] * 𝚽x12[ox+mom.cyz[k]] +
                 d.qcy[k,iky] * 𝚽y12[oy+mom.cxz[k]] +
                 d.qcz[k,ikz] * 𝚽z12[oz+mom.cxy[k]]
    end

    # Solve the equation system from the cached factorization
    gn_fast_solve!(𝚽,d.LU,d.ipiv,Q,NMOM,conf)

    # Closure relations
    gn_fast_closure!(𝚽x12,ox,d.clx,𝚽,Val(NMOM),Val(1))
    gn_fast_closure!(𝚽y12,oy,d.cly,𝚽,Val(NMOM),Val(1))
    gn_fast_closure!(𝚽z12,oz,d.clz,𝚽,Val(NMOM),Val(1))
    @inbounds for c in 1:NMOM
        𝚿[o𝚿+c] = 𝚽[c]
    end

    return nothing
end
