"""
    flux_2D_BFP(μ::Float64,η::Float64,Σt::Float64,S⁻::Float64,S⁺::Float64,
    S::Vector{Float64},ΔE::Float64,Δx::Float64,Δy::Float64,Qn::Vector{Float64},
    𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},𝚽E12::Vector{Float64},𝒪E::Int64,𝒪x::Int64,
    𝒪y::Int64,C::Vector{Float64},ωE::Array{Float64},ωx::Array{Float64},ωy::Array{Float64},
    isAdapt::Bool,𝒲::Array{Float64})

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

# Output Argument(s)
- `𝚽n::Vector{Float64}`: angular in-cell flux.
- `𝚽x12::Vector{Float64}`: outgoing angular flux along x-axis.
- `𝚽y12::Vector{Float64}`: outgoing angular flux along y-axis.
- `𝚽E12::Vector{Float64}`: outgoing angular flux along E-axis.

# Reference(s)
N/A

"""
function flux_2D_BFP(μ::Float64,η::Float64,Σt::Float64,S⁻::Float64,S⁺::Float64,S::Vector{Float64},ΔE::Float64,Δx::Float64,Δy::Float64,Qn::Vector{Float64},𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},𝚽E12::Vector{Float64},𝒪E::Int64,𝒪x::Int64,𝒪y::Int64,C::Vector{Float64},ωE::Array{Float64},ωx::Array{Float64},ωy::Array{Float64},isAdapt::Bool,𝒲::Array{Float64})

# Initialization
sx = sign(μ)
sy = sign(η)
hx = abs(μ)/Δx
hy = abs(η)/Δy
Nm = 𝒪E*𝒪x*𝒪y
𝒮 = zeros(Nm,Nm)
Q = zeros(Nm)
𝚽n = Q

# Adaptive weight calculations
if isAdapt ωx,ωy,ωE = adaptive(𝒪x,𝒪y,𝒪E,ωx,ωy,ωE,hx,hy,1/ΔE,sx,sy,-1,𝚽x12,𝚽y12,𝚽E12,Qn,Σt) end

# Matrix of Legendre moment coefficients of the flux
@inbounds for ix in range(1,𝒪x), jx in range(1,𝒪x), iy in range(1,𝒪y), jy in range(1,𝒪y), iE in range(1,𝒪E), jE in range(1,𝒪E)
    i = 𝒪x*𝒪E * (iy-1) + 𝒪E * (ix-1) + iE
    j = 𝒪x*𝒪E * (jy-1) + 𝒪E * (jx-1) + jE

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
@inbounds for jx in range(1,𝒪x), jy in range(1,𝒪y), jE in range(1,𝒪E)
    j = 𝒪x*𝒪E * (jy-1) + 𝒪E * (jx-1) + jE
    jEm = 𝒪x*(jy-1) + jx
    jxm = 𝒪E*(jy-1) + jE
    jym = 𝒪E*(jx-1) + jE
    Q[j] += Qn[j]
    Q[j] -= C[jx] * hx * (sx^(jx-1) * ωx[1,jy,jE] - (-sx)^(jx-1)) * 𝚽x12[jxm] 
    Q[j] -= C[jy] * hy * (sy^(jy-1) * ωy[1,jx,jE] - (-sy)^(jy-1)) * 𝚽y12[jym] 
    Q[j] -= C[jE] * ((-1)^(jE-1)*S⁺*ωE[1,jx,jy] - S⁻) * 𝚽E12[jEm]
end

𝚽n = 𝒮\Q

# Closure relation
@inbounds for jx in range(1,𝒪x), jy in range(1,𝒪y), jE in range(1,𝒪E)
    j = 𝒪x*𝒪E * (jy-1) + 𝒪E * (jx-1) + jE
    jEm = 𝒪x*(jy-1) + jx
    jxm = 𝒪E*(jy-1) + jE
    jym = 𝒪E*(jx-1) + jE
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