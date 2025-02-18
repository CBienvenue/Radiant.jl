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
@inbounds for ix in range(1,𝒪x), jx in range(1,𝒪x), iy in range(1,𝒪y), jy in range(1,𝒪y), iz in range(1,𝒪z), jz in range(1,𝒪z)
    if isFC
        i = 𝒪y*𝒪x*(iz-1) + 𝒪x * (iy-1) + ix
        j = 𝒪y*𝒪x*(jz-1) + 𝒪x * (jy-1) + jx
    else
        if count(>(1),(ix,iy,iz)) ≥ 2 || count(>(1),(jx,jy,jz)) ≥ 2 continue end
        i = 1 + (ix-1) + (iy-1) + (iz-1)
        j = 1 + (jx-1) + (jy-1) + (iz-1)
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
@inbounds for jx in range(1,𝒪x), jy in range(1,𝒪y), jz in range(1,𝒪z)
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
@inbounds for jx in range(1,𝒪x), jy in range(1,𝒪y), jz in range(1,𝒪z)
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
    if (jz == 1) 𝚽z12[jzm] = ωz[1,jx,jy] * 𝚽y12[jzm] end
    𝚽x12[jxm] += C[jx] * sx^(jx-1) * ωx[jx+1,jy,jz] * 𝚽n[j]
    𝚽y12[jym] += C[jy] * sy^(jy-1) * ωy[jy+1,jx,jz] * 𝚽n[j]
    𝚽y12[jzm] += C[jz] * sz^(jz-1) * ωz[jz+1,jx,jy] * 𝚽n[j]
end

# Returning solutions
return 𝚽n, 𝚽x12, 𝚽y12, 𝚽z12
end