"""
    flux_2D_BTE(μ::Float64,η::Float64,Σt::Float64,Δx::Float64,Δy::Float64,
    Qn::Vector{Float64},𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},𝒪x::Int64,𝒪y::Int64,
    C::Vector{Float64},ωx::Array{Float64},ωy::Array{Float64},isAdapt::Bool)

Compute flux solution in a cell in 2D Cartesian geometry for the Boltzmann transport
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
- `𝒪x::Int64`: spatial closure relation order.
- `𝒪y::Int64`: spatial closure relation order.
- `C::Vector{Float64}`: constants related to normalized Legendre.
- `ωx::Array{Float64}`: weighting factors of the x-axis scheme.
- `ωy::Array{Float64}`: weighting factors of the y-axis scheme.
- `isAdapt::Bool`: boolean for adaptive calculations.
- `isFC::Bool`: boolean indicating if the high-order incoming moments are fully coupled.

# Output Argument(s)
- `𝚽n::Vector{Float64}`: angular in-cell flux.
- `𝚽x12::Vector{Float64}`: outgoing angular flux along x-axis.
- `𝚽y12::Vector{Float64}`: outgoing angular flux along y-axis.

# Reference(s)
N/A

"""
function flux_2D_BTE(μ::Float64,η::Float64,Σt::Float64,Δx::Float64,Δy::Float64,Qn::Vector{Float64},𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},𝒪x::Int64,𝒪y::Int64,C::Vector{Float64},ωx::Array{Float64},ωy::Array{Float64},isAdapt::Bool,isFC::Bool)

# Initialization
sx = sign(μ)
sy = sign(η)
hx = abs(μ)/Δx
hy = abs(η)/Δy
if isFC Nm = 𝒪x*𝒪y else Nm = 𝒪x+𝒪y-1 end
𝒮 = zeros(Nm,Nm)
Q = zeros(Nm)
𝚽n = Q

# Adaptive weight calculations
if isAdapt ωx,ωy = adaptive(𝒪x,𝒪y,ωx,ωy,hx,hy,sx,sy,𝚽x12,𝚽y12,Qn,Σt,isFC) end

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,𝒪x), jx in range(1,𝒪x), iy in range(1,𝒪y), jy in range(1,𝒪y)
    if isFC
        i = 𝒪x*(iy-1)+ix
        j = 𝒪x*(jy-1)+jx
    else
        if count(>(1),(ix,iy)) ≥ 2 || count(>(1),(jx,jy)) ≥ 2 continue end
        i = 1 + (iy-1) + (ix-1)
        j = 1 + (jy-1) + (jx-1)
        if iy > 1 i += 𝒪x-1 end
        if jy > 1 j += 𝒪x-1 end
    end

    # Collision term
    if (i == j) 𝒮[i,j] += Σt end

    # Streaming term - x
    if iy == jy
        if (ix ≥ jx + 1) 𝒮[i,j] -= C[ix] * hx * sx * C[jx] * (1-(-1)^(ix-jx)) end
    end
    𝒮[i,j] += C[ix] * hx * sx^(ix-1) * C[jx] * sx^(jx-1) * ωx[jx+1,jy,iy]

    # Streaming term - y
    if ix == jx
        if (iy ≥ jy + 1) 𝒮[i,j] -= C[iy] * hy * sy * C[jy] * (1-(-1)^(iy-jy)) end 
    end
    𝒮[i,j] += C[iy] * hy * sy^(iy-1) * C[jy] * sy^(jy-1) * ωy[jy+1,jx,ix]
end

# Source vector
for jx in range(1,𝒪x), jy in range(1,𝒪y)
    if isFC
        j = 𝒪x*(jy-1)+jx
    else
        if count(>(1),(jx,jy)) ≥ 2 continue end
        j = 1 + (jy-1) + (jx-1)
        if jy > 1 j += 𝒪x-1 end
    end
    Q[j] += Qn[j]
    Q[j] -= C[jx] * hx * (sx^(jx-1) * ωx[1,jy,jy] - (-sx)^(jx-1)) * 𝚽x12[jy] 
    Q[j] -= C[jy] * hy * (sy^(jy-1) * ωy[1,jx,jx] - (-sy)^(jy-1)) * 𝚽y12[jx] 
end

# Solve the equation system
𝚽n = 𝒮\Q

# Closure relations
for jx in range(1,𝒪x), jy in range(1,𝒪y)
    if isFC
        j = 𝒪x*(jy-1)+jx
    else
        if count(>(1),(jx,jy)) ≥ 2 continue end
        j = 1 + (jy-1) + (jx-1)
        if jy > 1 j += 𝒪x-1 end
    end
    if (jx == 1) 𝚽x12[jy] = ωx[1,jy,jy] * 𝚽x12[jy] end
    if (jy == 1) 𝚽y12[jx] = ωy[1,jx,jx] * 𝚽y12[jx] end
    for iy in range(1,𝒪y)
        𝚽x12[jy] += C[jx] * sx^(jx-1) * ωx[jx+1,jy,iy] * 𝚽n[j]
    end
    for ix in range(1,𝒪x)
        𝚽y12[jx] += C[jy] * sy^(jy-1) * ωy[jy+1,jx,ix] * 𝚽n[j]
    end
end

# Returning solutions
return 𝚽n, 𝚽x12, 𝚽y12
end