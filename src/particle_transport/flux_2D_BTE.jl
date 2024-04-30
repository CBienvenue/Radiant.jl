"""
    flux_2D_BTE(μ::Float64,η::Float64,Σt::Float64,Δx::Float64,Δy::Float64,
    Qn::Vector{Float64},𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},𝒪x::Int64,𝒪y::Int64,
    Cx::Vector{Float64},Cy::Vector{Float64},ωx::Array{Float64},ωy::Array{Float64},
    isAdaptx::Bool,isAdapty::Bool)

Compute flux solution in a cell in 2D Cartesian geometry for the Boltzmann transport
equation.

# Input Argument(s)
- 'μ::Float64': direction cosine.
- 'η::Float64': direction cosine.
- 'Σt::Float64': total cross-sections.
- 'Δx::Float64': size of voxels along x-axis.
- 'Δy::Float64': size of voxels along y-axis.
- 'Qn::Vector{Float64}': angular in-cell source.
- '𝚽x12::Vector{Float64}': incoming angular flux along x-axis.
- '𝚽y12::Vector{Float64}': incoming angular flux along y-axis.
- '𝒪x::Int64': spatial closure relation order.
- '𝒪y::Int64': spatial closure relation order.
- 'Cx::Vector{Float64}': constants related to normalized Legendre.
- 'Cy::Vector{Float64}': constants related to normalized Legendre.
- 'ωx::Array{Float64}': weighting factors of the x-axis scheme.
- 'ωy::Array{Float64}': weighting factors of the y-axis scheme.
- 'isAdaptx::Bool': boolean for adaptive calculations.
- 'isAdapty::Bool': boolean for adaptive calculations.

# Output Argument(s)
- '𝚽n::Vector{Float64}': angular in-cell flux.
- '𝚽x12::Vector{Float64}': outgoing angular flux along x-axis.
- '𝚽y12::Vector{Float64}': outgoing angular flux along y-axis.

# Reference(s)
N/A

"""
function flux_2D_BTE(μ::Float64,η::Float64,Σt::Float64,Δx::Float64,Δy::Float64,Qn::Vector{Float64},𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},𝒪x::Int64,𝒪y::Int64,Cx::Vector{Float64},Cy::Vector{Float64},ωx::Array{Float64},ωy::Array{Float64},isAdaptx::Bool,isAdapty::Bool)

# Initialization
μ = μ * Δy
η = η * Δx
Nm = 𝒪x*𝒪y
S = zeros(Nm,Nm)
Q = zeros(Nm)
𝚽n = Q
Tx = 0.0; Ty = 0.0

# Adaptive loop
isAdapt = isAdaptx && isAdapty
isFixed = false
while ~isFixed

# Matrix of Legendre moment coefficients of the flux
@inbounds for ix in range(1,𝒪x), jx in range(1,𝒪x), iy in range(1,𝒪y), jy in range(1,𝒪y)
    i = 𝒪x*(iy-1)+ix
    j = 𝒪y*(jy-1)+jx
    # Diagonal terms
    if i == j
        S[i,j] = Σt * Δx * Δy + Cx[ix]^2 * ωx[jx+1,jy] * abs(μ) + Cy[iy]^2 * ωx[jy+1,jx] * abs(η)
    # Upper diagonal terms
    elseif i < j
    # Space terms - x
    if iy == jy
    if mod(ix+jx,2) == 1
        S[i,j] = Cx[ix] * Cx[jx] * ωx[jx+1,jy] * μ
    else
        S[i,j] = Cx[ix] * Cx[jx] * ωx[jx+1,jy] * abs(μ)
    end
    # Space terms - y
    elseif ix == jx
    if mod(iy+jy,2) == 1
        S[i,j] = Cy[iy] * Cy[jy] * ωy[jy+1,jx] * η
    else
        S[i,j] = Cy[iy] * Cy[jy] * ωy[jy+1,jx] * abs(η)
    end
    end
    # Under diagonal terms
    else
    # Space terms - x
    if iy == jy
    if mod(ix+jx,2) == 1
        S[i,j] = Cx[ix] * Cx[jx] * (ωx[jx+1,jy]-2) * μ
    else
        S[i,j] = Cx[ix] * Cx[jx] * ωx[jx+1,jy] * abs(μ)
    end
    # Space terms - y
    elseif ix == jx
    if mod(iy+jy,2) == 1
        S[i,j] = Cy[iy] * Cy[jy] * (ωy[jy+1,jx]-2) * η
    else
        S[i,j] = Cy[iy] * Cy[jy] * ωy[jy+1,jx] * abs(η)
    end
    end
    end

    # T-factors
    if 𝒪x == 2 && 𝒪y == 2
        if i == 2 && j == 3
            S[i,j] += abs(η) * Ty
        elseif i == 3 && j == 2
            S[i,j] += abs(μ) * Tx
        elseif i == 4 && j == 2
            S[i,j] += sqrt(3) * sign(μ) * abs(μ) * Tx
        elseif i == 4 && j == 3
            S[i,j] += sqrt(3) * sign(η) * abs(η) * Ty
        end
    end

end

# Source vector
@inbounds for ix in range(1,𝒪x), iy in range(1,𝒪y)
    i = 𝒪x*(iy-1)+ix
    Q[i] = Qn[i] * Δx * Δy
    if mod(ix,2) == 1
        Q[i] += Cx[ix] * (1-ωx[1,iy]) * 𝚽x12[iy] * abs(μ)
    else
        Q[i] += -Cx[ix] * (1+ωx[1,iy]) * 𝚽x12[iy] * μ
    end
    if mod(iy,2) == 1
        Q[i] += Cy[iy] * (1-ωy[1,ix]) * 𝚽y12[ix] * abs(η)
    else
        Q[i] += -Cy[iy] * (1+ωy[1,ix]) * 𝚽y12[ix] * η
    end
end

𝚽n = S\Q

# Adaptive correction of weighting parameters
if isAdapt
    isFixed, ω, T = adaptive_2D([𝒪x,𝒪y],[ωx,ωy],𝚽n,[𝚽x12,𝚽y12],[sign(μ),sign(η)],[1.0,1.0],[Tx,Ty])
    ωx = ω[1]; ωy = ω[2]; Tx = T[1]; Ty = T[2];
else
    isFixed = true
end

end # End of adaptive loop

# Closure relations
@inbounds for ix in range(1,𝒪x), iy in range(1,𝒪y)
    if (ix == 1) 𝚽x12[iy] = ωx[1,iy] * 𝚽x12[iy] end
    if (iy == 1) 𝚽y12[ix] = ωy[1,ix] * 𝚽y12[ix] end
end
@inbounds for ix in range(1,𝒪x), iy in range(1,𝒪y)
    i = 𝒪x*(iy-1)+ix
    # Space terms - x
    if mod(ix,2) == 1
        𝚽x12[iy] += Cx[ix] * ωx[ix+1,iy] * 𝚽n[i]
    else
        𝚽x12[iy] += Cx[ix] * ωx[ix+1,iy] * 𝚽n[i] * sign(μ)
    end
    # Space terms - y
    if mod(iy,2) == 1
        𝚽y12[ix] += Cy[iy] * ωy[iy+1,ix] * 𝚽n[i]
    else
        𝚽y12[ix] += Cy[iy] * ωy[iy+1,ix] * 𝚽n[i] * sign(η)
    end
end
if 𝒪x == 2 && 𝒪y == 2
    𝚽x12[2] += Tx * 𝚽n[2]
    𝚽y12[2] += Ty * 𝚽n[3]
end

# Returning solutions
return 𝚽n, 𝚽x12, 𝚽y12

end