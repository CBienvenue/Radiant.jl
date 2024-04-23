"""
    flux_3D_BTE(μ::Float64,η::Float64,ξ::Float64,Σt::Float64,Δx::Float64,Δy::Float64,
    Δz::Float64,Qn::Vector{Float64},𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},
    𝚽z12::Vector{Float64},𝒪x::Int64,𝒪y::Int64,𝒪z::Int64,Cx::Vector{Float64},
    Cy::Vector{Float64},Cz::Vector{Float64},ωx::Array{Float64},ωy::Array{Float64},
    ωz::Array{Float64},isAdaptx::Bool,isAdapty::Bool,isAdaptz::Bool)

Compute flux solution in a cell in 3D Cartesian geometry for the Boltzmann transport
equation.

See also [`compute_sweep_3D`](@ref), [`flux_3D_BFP`](@ref), [`flux_1D_BTE`](@ref),
[`flux_2D_BTE`](@ref).

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
- 'Cx::Vector{Float64}': constants related to normalized Legendre.
- 'Cy::Vector{Float64}': constants related to normalized Legendre.
- 'Cz::Vector{Float64}': constants related to normalized Legendre.
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

# Author(s)
Charles Bienvenue

# Reference(s)
N/A

"""
function flux_3D_BTE(μ::Float64,η::Float64,ξ::Float64,Σt::Float64,Δx::Float64,Δy::Float64,Δz::Float64,Qn::Vector{Float64},𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},𝚽z12::Vector{Float64},𝒪x::Int64,𝒪y::Int64,𝒪z::Int64,Cx::Vector{Float64},Cy::Vector{Float64},Cz::Vector{Float64},ωx::Array{Float64},ωy::Array{Float64},ωz::Array{Float64},isAdaptx::Bool,isAdapty::Bool,isAdaptz::Bool)

# Initialization
μ = μ * Δy * Δz
η = η * Δx * Δz
ξ = ξ * Δx * Δy
Nm = 𝒪x*𝒪y*𝒪z
S = zeros(Nm,Nm)
Q = zeros(Nm)
𝚽n = Q

# Adaptive loop
isAdapt = isAdaptx && isAdapty && isAdaptz
isFixed = false
while ~isFixed

# Matrix of Legendre moment coefficients of the flux
@inbounds for ix in range(1,𝒪x), jx in range(1,𝒪x), iy in range(1,𝒪y), jy in range(1,𝒪y), iz in range(1,𝒪z), jz in range(1,𝒪z)
    i = 𝒪y*𝒪x*(iz-1) + 𝒪x * (iy-1) + ix
    j = 𝒪y*𝒪x*(jz-1) + 𝒪x * (jy-1) + jx
    # Diagonal terms
    if i == j
        S[i,j] = Σt * Δx * Δy * Δz + Cx[ix]^2 * ωx[jx+1,jy,jz] * abs(μ) + Cy[iy]^2 * ωy[jy+1,jx,jz] * abs(η) + Cz[iz]^2 * ωy[jz+1,jx,jy] * abs(ξ)
    # Upper diagonal terms
    elseif i < j
    if iz == jz
    if iy == jy
    # Space terms - x
    if mod(ix+jx,2) == 1
        S[i,j] = Cx[ix] * Cx[jx] * ωx[jx+1,jy,jz] * μ
    else
        S[i,j] = Cx[ix] * Cx[jx] * ωx[jx+1,jy,jz] * abs(μ)
    end
    elseif ix == jx
    # Space terms - y
    if mod(iy+jy,2) == 1
        S[i,j] = Cy[iy] * Cy[jy] * ωy[jy+1,jx,jz] * η
    else
        S[i,j] = Cy[iy] * Cy[jy] * ωy[jy+1,jx,jz] * abs(η)
    end
    end
    elseif iy == jy && ix == jx
    # Space terms - z
    if mod(iz+jz,2) == 1
        S[i,j] = Cz[iz] * Cz[jz] * ωz[jz+1,jx,jy] * ξ
    else
        S[i,j] = Cz[iz] * Cz[jz] * ωz[jz+1,jx,jy] * abs(ξ)
    end
    end
# Under diagonal terms
    else
    if iz == jz
    if iy == jy
    # Space terms - x
    if mod(ix+jx,2) == 1
        S[i,j] = Cx[ix] * Cx[jx] * (ωx[jx+1,jy,jz]-2) * μ
    else
        S[i,j] = Cx[ix] * Cx[jx] * ωx[jx+1,jy,jz] * abs(μ)
    end
    elseif ix == jx
    # Space terms - y
    if mod(iy+jy,2) == 1
        S[i,j] = Cy[iy] * Cy[jy] * (ωy[jy+1,jx,jz]-2) * η
    else
        S[i,j] = Cy[iy] * Cy[jy] * ωy[jy+1,jx,jz] * abs(η)
    end
    end
    elseif iy == jy && ix == jx
    # Space terms - z
    if mod(iz+jz,2) == 1
        S[i,j] = Cz[iz] * Cz[jz] * (ωz[jz+1,jx,jy]-2) * ξ
    else
        S[i,j] = Cz[iz] * Cz[jz] * ωz[jz+1,jx,jy] * abs(ξ)
    end
    end
    end
end

# Source vector
@inbounds for ix in range(1,𝒪x), iy in range(1,𝒪y), iz in range(1,𝒪z)
    i = 𝒪y*𝒪x*(iz-1) + 𝒪x * (iy-1) + ix
    ixm = 𝒪y*(iz-1) + iy
    iym = 𝒪x*(iz-1) + ix
    izm = 𝒪x*(iy-1) + ix
    Q[i] = Qn[i] * Δx * Δy * Δz
    # Space terms - x
    if mod(ix,2) == 1
        Q[i] += Cx[ix] * (1-ωx[1,iy,iz]) * 𝚽x12[ixm] * abs(μ)
    else
        Q[i] += -Cx[ix] * (1+ωx[1,iy,iz]) * 𝚽x12[ixm] * μ
    end
    # Space terms - y
    if mod(iy,2) == 1
        Q[i] += Cy[iy] * (1-ωy[1,ix,iz]) * 𝚽y12[iym] * abs(η)
    else
        Q[i] += -Cy[iy] * (1+ωy[1,ix,iz]) * 𝚽y12[iym] * η
    end
    # Space terms - z
    if mod(iz,2) == 1
        Q[i] += Cz[iz] * (1-ωz[1,ix,iy]) * 𝚽z12[izm] * abs(ξ)
    else
        Q[i] += -Cz[iz] * (1+ωz[1,ix,iy]) * 𝚽z12[izm] * ξ
    end
end

𝚽n = S\Q

# Adaptive correction of weighting parameters
if isAdapt
    error("Not implemented yet.")
    #isFixed, ω = adaptive(2,[𝒪x,𝒪y,𝒪z],[ωx,ωy,ωz],𝚽n,[𝚽x12,𝚽y12,𝚽z12],[sign(μ),sign(η),sign(ξ)],[1.0,1.0,1.0],[0.0,0.0,0.0])
    #ωx = ω[1]; ωy = ω[2]; ωz = ω[3];
else
    isFixed = true
end


end # End of adaptive loop

# Closure relation
@inbounds for ix in range(1,𝒪x), iy in range(1,𝒪y), iz in range(1,𝒪z)
    ixm = 𝒪y*(iz-1) + iy
    iym = 𝒪x*(iz-1) + ix
    izm = 𝒪x*(iy-1) + ix
    if (ix == 1) 𝚽x12[ixm] = ωx[1,iy,iz] * 𝚽x12[ixm] end
    if (iy == 1) 𝚽y12[iym] = ωy[1,ix,iz] * 𝚽y12[iym] end
    if (iz == 1) 𝚽z12[izm] = ωz[1,ix,iy] * 𝚽z12[izm] end
end
@inbounds for ix in range(1,𝒪x), iy in range(1,𝒪y), iz in range(1,𝒪z)
    i = 𝒪y*𝒪x*(iz-1) + 𝒪x * (iy-1) + ix
    ixm = 𝒪y*(iz-1) + iy
    iym = 𝒪x*(iz-1) + ix
    izm = 𝒪x*(iy-1) + ix
    # Space terms - x
    if mod(ix,2) == 1
        𝚽x12[ixm] += Cx[ix] * ωx[ix+1,iy,iz] * 𝚽n[i]
    else
        𝚽x12[ixm] += Cx[ix] * ωx[ix+1,iy,iz] * 𝚽n[i] * sign(μ)
    end
    # Space terms - y
    if mod(iy,2) == 1
        𝚽y12[iym] += Cy[iy] * ωy[iy+1,ix,iz] * 𝚽n[i]
    else
        𝚽y12[iym] += Cy[iy] * ωy[iy+1,ix,iz] * 𝚽n[i] * sign(η)
    end
    # Space terms - z
    if mod(iz,2) == 1
        𝚽z12[izm] += Cz[iz] * ωy[iz+1,ix,iy] * 𝚽n[i]
    else
        𝚽z12[izm] += Cz[iz] * ωy[iz+1,ix,iy] * 𝚽n[i] * sign(ξ)
    end
end

# Returning solutions
return 𝚽n, 𝚽x12, 𝚽y12, 𝚽z12

end