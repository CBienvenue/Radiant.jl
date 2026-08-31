"""
    flux_2D_BFP(μ::Float64,η::Float64,Σt::Float64,S⁻::Float64,S⁺::Float64,
    S::Vector{Float64},ΔE::Float64,Δx::Float64,Δy::Float64,Qn::Vector{Float64},
    𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},𝚽E12::Vector{Float64},𝒪E::Int64,𝒪x::Int64,
    𝒪y::Int64,C::Vector{Float64},ωE::Array{Float64},ωx::Array{Float64},ωy::Array{Float64},
    isAdapt::Bool,𝒲::Array{Float64},isFC::Bool)

Compute flux solution in a cell in 2D Cartesian geometry for the Boltzmann Fokker-Planck
equation.

# Input Argument(s)
- `μ::Float64`: direction cosine.
- `η::Float64`: direction cosine.
- `Σt::Float64`: total cross-sections.
- `Δx::Float64`: size of voxels along x-axis.
- `Δy::Float64`: size of voxels along y-axis.
- `Qn::Vector{Float64}`: angular in-cell source.
- `𝚽x12::Vector{Float64}`: incoming angular flux along x-axis.
- `𝚽y12::Vector{Float64}`: incoming angular flux along y-axis.
- `S⁻::Float64`: restricted stopping power at upper energy group boundary.
- `S⁺::Float64`: restricted stopping power at lower energy group boundary.
- `ΔE::Float64`: energy group width.
- `𝚽E12::Vector{Float64}`: incoming angular flux along E-axis.
- `𝒪E::Int64`: energy closure relation order.
- `𝒪x::Int64`: spatial closure relation order.
- `𝒪y::Int64`: spatial closure relation order.
- `C::Vector{Float64}`: constants related to normalized Legendre.
- `ωE::Array{Float64}`: weighting factors of the E-axis scheme.
- `ωx::Array{Float64}`: weighting factors of the x-axis scheme.
- `ωy::Array{Float64}`: weighting factors of the y-axis scheme.
- `isAdapt::Bool`: boolean for adaptive calculations.
- `𝒲::Array{Float64}` : weighting constants.
- `isFC::Bool`: boolean indicating if the high-order incoming moments are fully coupled.

# Output Argument(s)
- `𝚽n::Vector{Float64}`: angular in-cell flux.
- `𝚽x12::Vector{Float64}`: outgoing angular flux along x-axis.
- `𝚽y12::Vector{Float64}`: outgoing angular flux along y-axis.
- `𝚽E12::Vector{Float64}`: outgoing angular flux along E-axis.

# Reference(s)
N/A

"""
function flux_2D_BFP(μ::Float64,η::Float64,Σt::Float64,S⁻::Float64,S⁺::Float64,S::Vector{Float64},ΔE::Float64,Δx::Float64,Δy::Float64,Qn::Vector{Float64},𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},𝚽E12::Vector{Float64},𝒪E::Int64,𝒪x::Int64,𝒪y::Int64,C::Vector{Float64},ωE::Array{Float64},ωx::Array{Float64},ωy::Array{Float64},isAdapt::Bool,𝒲::Array{Float64},isFC::Bool)

# Initialization
sx = sign(μ)
sy = sign(η)
hx = abs(μ)/Δx
hy = abs(η)/Δy
if isFC Nm = 𝒪E*𝒪x*𝒪y else Nm = 𝒪E+𝒪x+𝒪y-2 end
𝒮 = zeros(Nm,Nm)
Q = zeros(Nm)
𝚽n = Q

# Adaptive weight calculations
if isAdapt ωx,ωy,ωE = adaptive(𝒪x,𝒪y,𝒪E,ωx,ωy,ωE,hx,hy,1/ΔE,sx,sy,-1,𝚽x12,𝚽y12,𝚽E12,Qn,Σt,isFC) end

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,𝒪x), jx in range(1,𝒪x), iy in range(1,𝒪y), jy in range(1,𝒪y), iE in range(1,𝒪E), jE in range(1,𝒪E)
    if isFC
        i = 𝒪x*𝒪E * (iy-1) + 𝒪E * (ix-1) + iE
        j = 𝒪x*𝒪E * (jy-1) + 𝒪E * (jx-1) + jE
    else
        if count(>(1),(ix,iy,iE)) ≥ 2 || count(>(1),(jx,jy,jE)) ≥ 2 continue end
        i = 1 + (iE-1) + (ix-1) + (iy-1)
        j = 1 + (jE-1) + (jx-1) + (jy-1)
        if ix > 1 i += 𝒪E-1 end
        if iy > 1 i += 𝒪E-1 + 𝒪x-1 end
        if jx > 1 j += 𝒪E-1 end
        if jy > 1 j += 𝒪E-1 + 𝒪x-1 end
    end

    # Collision term
    if (i == j) 𝒮[i,j] += Σt end

    # Streaming term - x
    if iy == jy && iE == jE
        if (ix ≥ jx + 1) 𝒮[i,j] -= C[ix] * hx * sx * C[jx] * (1-(-1)^(ix-jx)) end
        𝒮[i,j] += C[ix] * hx * sx^(ix-1) * C[jx] * sx^(jx-1) * ωx[jx+1,jy,jE]
    end

    # Streaming term - y
    if ix == jx && iE == jE
        if (iy ≥ jy + 1) 𝒮[i,j] -= C[iy] * hy * sy * C[jy] * (1-(-1)^(iy-jy)) end 
        𝒮[i,j] += C[iy] * hy * sy^(iy-1) * C[jy] * sy^(jy-1) * ωy[jy+1,jx,jE]
    end

    # CSD term
    if ix == jx && iy == jy
        for kE in range(1,iE-1), wE in range(1,𝒪E)
            𝒮[i,j] += C[iE] * C[jE] * C[kE] * C[wE] * (1-(-1)^(iE-kE)) * S[wE] * 𝒲[jE,kE,wE]
        end
        𝒮[i,j] += C[iE] * S⁺ * (-1)^(iE-1) * C[jE] * (-1)^(jE-1) * ωE[jE+1,jx,jy]
    end
end

# Source vector
for jx in range(1,𝒪x), jy in range(1,𝒪y), jE in range(1,𝒪E)
    if isFC
        j = 𝒪x*𝒪E * (jy-1) + 𝒪E * (jx-1) + jE
        jEm = 𝒪x*(jy-1) + jx
        jxm = 𝒪E*(jy-1) + jE
        jym = 𝒪E*(jx-1) + jE
    else
        if count(>(1),(jx,jy,jE)) ≥ 2 continue end
        j = 1 + (jE-1) + (jx-1) + (jy-1)
        jEm = 1 + (jx-1) + (jy-1)
        jxm = 1 + (jE-1) + (jy-1)
        jym = 1 + (jE-1) + (jx-1)
        if jx > 1 j += 𝒪E-1 end
        if jy > 1 j += 𝒪E-1 + 𝒪x-1 end
        if jy > 1 jEm += 𝒪x-1 end
        if jy > 1 jxm += 𝒪E-1 end
        if jx > 1 jym += 𝒪E-1 end
    end
    Q[j] += Qn[j]
    Q[j] -= C[jx] * hx * (sx^(jx-1) * ωx[1,jy,jE] - (-sx)^(jx-1)) * 𝚽x12[jxm] 
    Q[j] -= C[jy] * hy * (sy^(jy-1) * ωy[1,jx,jE] - (-sy)^(jy-1)) * 𝚽y12[jym] 
    Q[j] -= C[jE] * ((-1)^(jE-1)*S⁺*ωE[1,jx,jy] - S⁻) * 𝚽E12[jEm]
end

𝚽n = 𝒮\Q

# Closure relation
for jx in range(1,𝒪x), jy in range(1,𝒪y), jE in range(1,𝒪E)
    if isFC
        j = 𝒪x*𝒪E * (jy-1) + 𝒪E * (jx-1) + jE
        jEm = 𝒪x*(jy-1) + jx
        jxm = 𝒪E*(jy-1) + jE
        jym = 𝒪E*(jx-1) + jE
    else
        if count(>(1),(jx,jy,jE)) ≥ 2 continue end
        j = 1 + (jE-1) + (jx-1) + (jy-1)
        jEm = 1 + (jx-1) + (jy-1)
        jxm = 1 + (jE-1) + (jy-1)
        jym = 1 + (jE-1) + (jx-1)
        if jx > 1 j += 𝒪E-1 end
        if jy > 1 j += 𝒪E-1 + 𝒪x-1 end
        if jy > 1 jEm += 𝒪x-1 end
        if jy > 1 jxm += 𝒪E-1 end
        if jx > 1 jym += 𝒪E-1 end
    end
    if (jx == 1) 𝚽x12[jxm] = ωx[1,jy,jE] * 𝚽x12[jxm] end
    if (jy == 1) 𝚽y12[jym] = ωy[1,jx,jE] * 𝚽y12[jym] end
    if (jE == 1) 𝚽E12[jEm] = ωE[1,jx,jy] * 𝚽E12[jEm] end
    𝚽x12[jxm] += C[jx] * sx^(jx-1) * ωx[jx+1,jy,jE] * 𝚽n[j]
    𝚽y12[jym] += C[jy] * sy^(jy-1) * ωy[jy+1,jx,jE] * 𝚽n[j]
    𝚽E12[jEm] += C[jE] * (-1)^(jE-1) * ωE[jE+1,jx,jy] * 𝚽n[j]
end

# Returning solutions
return 𝚽n, 𝚽x12, 𝚽y12, 𝚽E12
end
"""
    sn_2D_BFP_matrix!(...)

Assemble the cell matrix `𝒮` of `flux_2D_BFP`, for the optimized solver chain.

The body is the reference kernel's assembly loop transcribed unchanged, so the matrix is
bit-identical to the one it builds. It is split out because the matrix does not depend on the
voxel beyond its widths and its material: `sn_fast_context` calls this once per (material,
mesh-width combination, direction) and factorizes the result, where the reference rebuilds and
refactorizes it in every cell of every sweep. Valid only when `~isAdapt`.
"""
function sn_2D_BFP_matrix!(𝒮::Matrix{Float64},μ::Float64,η::Float64,Σt::Float64,S⁺::Float64,S::Vector{Float64},ΔE::Float64,Δx::Float64,Δy::Float64,𝒪E::Int64,𝒪x::Int64,𝒪y::Int64,C::Vector{Float64},ωE::Array{Float64},ωx::Array{Float64},ωy::Array{Float64},𝒲::Array{Float64},isFC::Bool)

fill!(𝒮,0.0)
sx = sign(μ)
sy = sign(η)
hx = abs(μ)/Δx
hy = abs(η)/Δy

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,𝒪x), jx in range(1,𝒪x), iy in range(1,𝒪y), jy in range(1,𝒪y), iE in range(1,𝒪E), jE in range(1,𝒪E)
    if isFC
        i = 𝒪x*𝒪E * (iy-1) + 𝒪E * (ix-1) + iE
        j = 𝒪x*𝒪E * (jy-1) + 𝒪E * (jx-1) + jE
    else
        if count(>(1),(ix,iy,iE)) ≥ 2 || count(>(1),(jx,jy,jE)) ≥ 2 continue end
        i = 1 + (iE-1) + (ix-1) + (iy-1)
        j = 1 + (jE-1) + (jx-1) + (jy-1)
        if ix > 1 i += 𝒪E-1 end
        if iy > 1 i += 𝒪E-1 + 𝒪x-1 end
        if jx > 1 j += 𝒪E-1 end
        if jy > 1 j += 𝒪E-1 + 𝒪x-1 end
    end

    # Collision term
    if (i == j) 𝒮[i,j] += Σt end

    # Streaming term - x
    if iy == jy && iE == jE
        if (ix ≥ jx + 1) 𝒮[i,j] -= C[ix] * hx * sx * C[jx] * (1-(-1)^(ix-jx)) end
        𝒮[i,j] += C[ix] * hx * sx^(ix-1) * C[jx] * sx^(jx-1) * ωx[jx+1,jy,jE]
    end

    # Streaming term - y
    if ix == jx && iE == jE
        if (iy ≥ jy + 1) 𝒮[i,j] -= C[iy] * hy * sy * C[jy] * (1-(-1)^(iy-jy)) end 
        𝒮[i,j] += C[iy] * hy * sy^(iy-1) * C[jy] * sy^(jy-1) * ωy[jy+1,jx,jE]
    end

    # CSD term
    if ix == jx && iy == jy
        for kE in range(1,iE-1), wE in range(1,𝒪E)
            𝒮[i,j] += C[iE] * C[jE] * C[kE] * C[wE] * (1-(-1)^(iE-kE)) * S[wE] * 𝒲[jE,kE,wE]
        end
        𝒮[i,j] += C[iE] * S⁺ * (-1)^(iE-1) * C[jE] * (-1)^(jE-1) * ωE[jE+1,jx,jy]
    end
end

return nothing

end

"""
    sn_2D_BFP_fast!(𝚿,o𝚿,𝚽x12,ox,𝚽y12,oy,𝚽E12,𝚽E12o,oE,Ql,oQ,Mnn,Np,mom,d,ws,conf,ikx,iky,
    m,do_E,zero_E)

Optimized counterpart of `flux_2D_BFP`: solve one voxel of the 2D BFP equation.

Identical in structure to `sn_3D_BFP_fast!` minus the z axis; see it for the treatment of the
energy axis, whose incoming flux is left intact so no per-direction copy is needed.
"""
@inline function sn_2D_BFP_fast!(𝚿::Vector{Float64},o𝚿::Int64,𝚽x12::Vector{Float64},ox::Int64,𝚽y12::Vector{Float64},oy::Int64,𝚽E12::Vector{Float64},𝚽E12o::Vector{Float64},oE::Int64,Ql::Vector{Float64},oQ::Int64,Mnn::Vector{Float64},Np::Int64,mom::GNFastMoments{NMOM},d::SNFastDir{NMOM},ws::GNFastWorkspace,conf::Int64,ikx::Int64,iky::Int64,m::Int64,do_E::Bool,zero_E::Bool) where {NMOM}

    Q = ws.Q; 𝚽 = ws.𝚽
    @inbounds for k in 1:NMOM
        col = mom.col[k]
        q = 0.0; b = oQ + Np*(col-1)
        @simd for p in 1:Np
            q += Mnn[p] * Ql[b+p]
        end
        q += d.qcx[k,ikx] * 𝚽x12[ox+mom.cyz[k]] + d.qcy[k,iky] * 𝚽y12[oy+mom.cxz[k]]
        if ~zero_E q += d.qcE[k,m] * 𝚽E12[oE+mom.cxyz[k]] end
        Q[col] = q
    end
    gn_fast_solve!(𝚽,d.LU,d.ipiv,Q,NMOM,conf)
    if do_E gn_fast_closure!(𝚽E12o,oE,d.clE,𝚽,Val(NMOM),Val(1)) end
    gn_fast_closure!(𝚽x12,ox,d.clx,𝚽,Val(NMOM),Val(1))
    gn_fast_closure!(𝚽y12,oy,d.cly,𝚽,Val(NMOM),Val(1))
    @inbounds for c in 1:NMOM
        𝚿[o𝚿+c] = 𝚽[c]
    end
    return nothing
end
