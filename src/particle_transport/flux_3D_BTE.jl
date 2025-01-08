"""
    flux_3D_BTE(μ::Float64,η::Float64,ξ::Float64,Σt::Float64,Δx::Float64,Δy::Float64,
    Δz::Float64,Qn::Vector{Float64},𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},
    𝚽z12::Vector{Float64},𝒪x::Int64,𝒪y::Int64,𝒪z::Int64,C::Vector{Float64},
    C::Vector{Float64},C::Vector{Float64},ωx::Array{Float64},ωy::Array{Float64},
    ωz::Array{Float64},isAdaptx::Bool,isAdapty::Bool,isAdaptz::Bool)

Compute flux solution in a cell in 3D Cartesian geometry for the Boltzmann transport
equation.

# Input Argument(s)
- 'μ::Float64': direction cosine.
- 'η::Float64': direction cosine.
- 'ξ::Float64': direction cosine.
- 'Σt::Float64': total cross-sections.
- 'Δx::Float64': size of voxels along x-axis.
- 'Δy::Float64': size of voxels along y-axis.
- 'Δz::Float64': size of voxels along z-axis.
- 'Qn::Vector{Float64}': angular in-cell source.
- '𝚽x12::Vector{Float64}': incoming angular flux along x-axis.
- '𝚽y12::Vector{Float64}': incoming angular flux along y-axis.
- '𝚽z12::Vector{Float64}': incoming angular flux along z-axis.
- '𝒪x::Int64': spatial closure relation order.
- '𝒪y::Int64': spatial closure relation order.
- '𝒪z::Int64': spatial closure relation order.
- 'C::Vector{Float64}': constants related to normalized Legendre.
- 'C::Vector{Float64}': constants related to normalized Legendre.
- 'C::Vector{Float64}': constants related to normalized Legendre.
- 'ωx::Array{Float64}': weighting factors of the x-axis scheme.
- 'ωy::Array{Float64}': weighting factors of the y-axis scheme.
- 'ωz::Array{Float64}': weighting factors of the z-axis scheme.
- 'isAdaptx::Bool': boolean for adaptive calculations.
- 'isAdapty::Bool': boolean for adaptive calculations.
- 'isAdaptz::Bool': boolean for adaptive calculations.

# Output Argument(s)
- '𝚽n::Vector{Float64}': angular in-cell flux.
- '𝚽x12::Vector{Float64}': outgoing angular flux along x-axis.
- '𝚽y12::Vector{Float64}': outgoing angular flux along y-axis.
- '𝚽z12::Vector{Float64}': outgoing angular flux along z-axis.

# Reference(s)
N/A

"""
function flux_3D_BTE(μ::Float64,η::Float64,ξ::Float64,Σt::Float64,Δx::Float64,Δy::Float64,Δz::Float64,Qn::Vector{Float64},𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},𝚽z12::Vector{Float64},𝒪x::Int64,𝒪y::Int64,𝒪z::Int64,C::Vector{Float64},ωx::Array{Float64},ωy::Array{Float64},ωz::Array{Float64},isAdapt::Bool)

# Initialization
sx = sign(μ)
sy = sign(η)
sz = sign(ξ)
hx = abs(μ)/Δx
hy = abs(η)/Δy
hz = abs(ξ)/Δz
Nm = 𝒪x*𝒪y*𝒪z
S = zeros(Nm,Nm)
Q = zeros(Nm)
𝚽n = Q

# Adaptive weight calculations
if isAdapt ωx,ωy,ωz = adaptive(𝒪x,𝒪y,𝒪z,ωx,ωy,ωz,hx,hy,hz,sx,sy,sz,𝚽x12,𝚽y12,𝚽z12,Qn,Σt) end

# Matrix of Legendre moment coefficients of the flux
@inbounds for ix in range(1,𝒪x), jx in range(1,𝒪x), iy in range(1,𝒪y), jy in range(1,𝒪y), iz in range(1,𝒪z), jz in range(1,𝒪z)
    i = 𝒪y*𝒪x*(iz-1) + 𝒪x * (iy-1) + ix
    j = 𝒪y*𝒪x*(jz-1) + 𝒪x * (jy-1) + jx
    if (i == j) S[i,j] += Σt end
    if iy == jy && iz == jz
        if (ix ≥ jx + 1) S[i,j] -= C[ix] * hx * sx * C[jx] * (1-(-1)^(ix-jx)) end
        S[i,j] += C[ix] * hx * sx^(ix-1) * C[jx] * sx^(jx-1) * ωx[jx+1,jy,jz]
    end
    if ix == jx && iz == jz
        if (iy ≥ jy + 1) S[i,j] -= C[iy] * hy * sy * C[jy] * (1-(-1)^(iy-jy)) end 
        S[i,j] += C[iy] * hy * sy^(iy-1) * C[jy] * sy^(jy-1) * ωy[jy+1,jx,jz]
    end
    if ix == jx && iy == jy
        if (iz ≥ jz + 1) S[i,j] -= C[iz] * hz * sz * C[jz] * (1-(-1)^(iz-jz)) end 
        S[i,j] += C[iz] * hz * sz^(iz-1) * C[jz] * sz^(jz-1) * ωz[jz+1,jx,jy]
    end
end

# Source vector
@inbounds for jx in range(1,𝒪x), jy in range(1,𝒪y), jz in range(1,𝒪z)
    j = 𝒪y*𝒪x*(jz-1) + 𝒪x * (jy-1) + jx
    jxm = 𝒪y*(jz-1) + jy
    jym = 𝒪x*(jz-1) + jx
    jzm = 𝒪x*(jy-1) + jx
    Q[j] = Qn[j]
    Q[j] -= C[jx] * hx * (sx^(jx-1) * ωx[1,jy,jz] - (-sx)^(jx-1)) * 𝚽x12[jxm] 
    Q[j] -= C[jy] * hy * (sy^(jy-1) * ωy[1,jx,jz] - (-sy)^(jy-1)) * 𝚽y12[jym] 
    Q[j] -= C[jz] * hz * (sz^(jz-1) * ωz[1,jx,jy] - (-sz)^(jz-1)) * 𝚽z12[jzm]
end

# Solve the equation system
𝚽n = S\Q

# Closure relation
@inbounds for jx in range(1,𝒪x), jy in range(1,𝒪y), jz in range(1,𝒪z)
    j = 𝒪x*(jy-1)+jx
    jxm = 𝒪y*(jz-1) + jy
    jym = 𝒪x*(jz-1) + jx
    jzm = 𝒪x*(jy-1) + jx
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