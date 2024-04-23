"""
    function flux_3D_BFP(μ::Float64,η::Float64,ξ::Float64,Σt::Float64,β⁻::Float64,
    β⁺::Float64,ΔE::Float64,Δx::Float64,Δy::Float64,Δz::Float64,Qn::Vector{Float64},
    𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},𝚽z12::Vector{Float64},
    𝚽E12::Vector{Float64},𝒪E::Int64,𝒪x::Int64,𝒪y::Int64,𝒪z::Int64,CE::Vector{Float64},
    Cx::Vector{Float64},Cy::Vector{Float64},Cz::Vector{Float64},ωE::Array{Float64},
    ωx::Array{Float64},ωy::Array{Float64},ωz::Array{Float64},isAdaptE::Bool,
    isAdaptx::Bool,isAdapty::Bool,isAdaptz::Bool)

Compute flux solution in a cell in 3D Cartesian geometry for the Boltzmann Fokker-Planck
equation.

See also [`compute_sweep_3D`](@ref), [`flux_3D_BTE`](@ref), [`flux_1D_BFP`](@ref),
[`flux_2D_BFP`](@ref).

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
- 'β⁻::Float64': restricted stopping power at upper energy group boundary.
- 'β⁺::Float64': restricted stopping power at lower energy group boundary.
- 'ΔE::Float64': energy group width.
- '𝚽E12::Vector{Float64}': incoming angular flux along E-axis.
- '𝒪E::Int64': energy closure relation order.
- '𝒪x::Int64': spatial closure relation order.
- '𝒪y::Int64': spatial closure relation order.
- '𝒪z::Int64': spatial closure relation order.
- 'CE::Vector{Float64}': constants related to normalized Legendre.
- 'Cx::Vector{Float64}': constants related to normalized Legendre.
- 'Cy::Vector{Float64}': constants related to normalized Legendre.
- 'Cz::Vector{Float64}': constants related to normalized Legendre.
- 'ωE::Array{Float64}': weighting factors of the E-axis scheme.
- 'ωx::Array{Float64}': weighting factors of the x-axis scheme.
- 'ωy::Array{Float64}': weighting factors of the y-axis scheme.
- 'ωz::Array{Float64}': weighting factors of the z-axis scheme.
- 'isAdaptE::Bool': boolean for adaptive calculations.
- 'isAdaptx::Bool': boolean for adaptive calculations.
- 'isAdapty::Bool': boolean for adaptive calculations.
- 'isAdaptz::Bool': boolean for adaptive calculations.

# Output Argument(s)
- '𝚽n::Vector{Float64}': angular in-cell flux.
- '𝚽x12::Vector{Float64}': outgoing angular flux along x-axis.
- '𝚽y12::Vector{Float64}': outgoing angular flux along y-axis.
- '𝚽z12::Vector{Float64}': outgoing angular flux along z-axis.
- '𝚽E12::Vector{Float64}': outgoing angular flux along E-axis.

# Author(s)
Charles Bienvenue

# Reference(s)
N/A

"""
function flux_3D_BFP(μ::Float64,η::Float64,ξ::Float64,Σt::Float64,β⁻::Float64,β⁺::Float64,ΔE::Float64,Δx::Float64,Δy::Float64,Δz::Float64,Qn::Vector{Float64},𝚽x12::Vector{Float64},𝚽y12::Vector{Float64},𝚽z12::Vector{Float64},𝚽E12::Vector{Float64},𝒪E::Int64,𝒪x::Int64,𝒪y::Int64,𝒪z::Int64,CE::Vector{Float64},Cx::Vector{Float64},Cy::Vector{Float64},Cz::Vector{Float64},ωE::Array{Float64},ωx::Array{Float64},ωy::Array{Float64},ωz::Array{Float64},isAdaptE::Bool,isAdaptx::Bool,isAdapty::Bool,isAdaptz::Bool)

# Initialization
μ = μ * Δy * Δz
η = η * Δx * Δz
ξ = ξ * Δx * Δy
Nm = 𝒪x*𝒪y*𝒪z*𝒪E
S = zeros(Nm,Nm)
Q = zeros(Nm)
𝚽n = Q

# Galerkin energy scheme weights
Λ = β⁻/β⁺
if abs(ωE[1,1,1,1]) > 0
    ωE[1,1,1,1] = ωE[1,1,1,1]*Λ
    ωE[2:𝒪E+1,1,1,1] = (ωE[2:𝒪E+1,1,1,1].-1).*Λ.+1
end

# Adaptive loop
isAdapt = isAdaptE && isAdaptx && isAdapty && isAdaptz
isFixed = false
while ~isFixed

# Matrix of Legendre moment coefficients of the flux
@inbounds for ix in range(1,𝒪x), jx in range(1,𝒪x), iy in range(1,𝒪y), jy in range(1,𝒪y), iz in range(1,𝒪z), jz in range(1,𝒪z), iE in range(1,𝒪E), jE in range(1,𝒪E)
    i = 𝒪y*𝒪x*𝒪E*(iz-1) + 𝒪x*𝒪E * (iy-1) + 𝒪E * (ix-1) + iE
    j = 𝒪y*𝒪x*𝒪E*(jz-1) + 𝒪x*𝒪E * (jy-1) + 𝒪E * (jx-1) + jE
    # Diagonal terms
    if i == j
        S[i,j] = (Σt + CE[iE]^2 * β⁺ * ωE[jE+1,jx,jy,jz] + (iE-1) * (β⁻-β⁺) ) * Δx * Δy * Δz + Cx[ix]^2 * ωx[jx+1,jy,jz,jE] * abs(μ) + Cy[iy]^2 * ωy[jy+1,jx,jz,jE] * abs(η) + Cz[iz]^2 * ωy[jz+1,jx,jy,jE] * abs(ξ)
    # Upper diagonal terms
    elseif i < j
    if iz == jz
    if iy == jy
    if ix == jx
    # Energy terms - E
    if mod(iE+jE,2) == 1
        S[i,j] = -CE[iE] * CE[jE] * β⁺ * ωE[jE+1,jx,jy,jz] * Δx * Δy * Δz
    else
        S[i,j] = CE[iE] * CE[jE] * β⁺ * ωE[jE+1,jx,jy,jz] * Δx * Δy * Δz
    end
    elseif iE == jE
    # Space terms - x
    if mod(ix+jx,2) == 1
        S[i,j] = Cx[ix] * Cx[jx] * ωx[jx+1,jy,jz,jE] * μ
    else
        S[i,j] = Cx[ix] * Cx[jx] * ωx[jx+1,jy,jz,jE] * abs(μ)
    end
    end
    elseif ix == jx && iE == jE
    # Space terms - y
    if mod(iy+jy,2) == 1
        S[i,j] = Cy[iy] * Cy[jy] * ωy[jy+1,jx,jz,jE] * η
    else
        S[i,j] = Cy[iy] * Cy[jy] * ωy[jy+1,jx,jz,jE] * abs(η)
    end
    end
    elseif ix == jx && iy == jy && iE == jE
    # Space terms - z
    if mod(iz+jz,2) == 1
        S[i,j] = Cz[iz] * Cz[jz] * ωz[jz+1,jx,jy,jE] * ξ
    else
        S[i,j] = Cz[iz] * Cz[jz] * ωz[jz+1,jx,jy,jE] * abs(ξ)
    end
    end
    # Under diagonal terms
    else
    if iz == jz
    if iy == jy
    if ix == jx
    # Energy terms - E
    if mod(iE+jE,2) == 1
        S[i,j] = -CE[iE] * CE[jE] * (β⁺*ωE[jE+1,jx,jy,jz]-β⁻-β⁺) * Δx * Δy * Δz
    else
        S[i,j] = CE[iE] * CE[jE] * (β⁺*ωE[jE+1,jx,jy,jz]+β⁻-β⁺) * Δx * Δy * Δz
    end
    elseif iE == jE
    # Space terms - x
    if mod(ix+jx,2) == 1
        S[i,j] = Cx[ix] * Cx[jx] * (ωx[jx+1,jy,jz,jE]-2) * μ
    else
        S[i,j] = Cx[ix] * Cx[jx] * ωx[jx+1,jy,jz,jE] * abs(μ)
    end
    end
    elseif ix == jx && iE == jE
    # Space terms - y
    if mod(iy+jy,2) == 1
        S[i,j] = Cy[iy] * Cy[jy] * (ωy[jy+1,jx,jz,jE]-2) * η
    else
        S[i,j] = Cy[iy] * Cy[jy] * ωy[jy+1,jx,jz,jE] * abs(η)
    end
    end
    elseif ix == jx && iy == jy && iE == jE
    # Space terms - z
    if mod(iz+jz,2) == 1
        S[i,j] = Cz[iz] * Cz[jz] * (ωz[jz+1,jx,jy,jE]-2) * ξ
    else
        S[i,j] = Cz[iz] * Cz[jz] * ωz[jz+1,jx,jy,jE] * abs(ξ)
    end
    end
    end
end

# Source vector
@inbounds for ix in range(1,𝒪x), iy in range(1,𝒪y), iz in range(1,𝒪z), iE in range(1,𝒪E)
    i = 𝒪y*𝒪x*𝒪E*(iz-1) + 𝒪x*𝒪E * (iy-1) + 𝒪E * (ix-1) + iE
    iEm = 𝒪y*𝒪x*(iz-1)+𝒪x*(iy-1)+ix
    ixm = 𝒪y*𝒪E*(iz-1)+𝒪E*(iy-1)+iE
    iym = 𝒪x*𝒪E*(iz-1)+𝒪E*(ix-1)+iE
    izm = 𝒪x*𝒪E*(iy-1)+𝒪E*(ix-1)+iE
    Q[i] = Qn[i] * Δx * Δy * Δz
    # Energy terms - E
    if mod(iE,2) == 1
        Q[i] += CE[iE] * (β⁻-β⁺*ωE[1,ix,iy,iz]) * 𝚽E12[iEm] * Δx * Δy * Δz
    else
        Q[i] += CE[iE] * (β⁻+β⁺*ωE[1,ix,iy,iz]) * 𝚽E12[iEm] * Δx * Δy * Δz
    end
    # Space terms - x
    if mod(ix,2) == 1
        Q[i] += Cx[ix] * (1-ωx[1,iy,iz,iE]) * 𝚽x12[ixm] * abs(μ)
    else
        Q[i] += -Cx[ix] * (1+ωx[1,iy,iz,iE]) * 𝚽x12[ixm] * μ
    end
    # Space terms - y
    if mod(iy,2) == 1
        Q[i] += Cy[iy] * (1-ωy[1,ix,iz,iE]) * 𝚽y12[iym] * abs(η)
    else
        Q[i] += -Cy[iy] * (1+ωy[1,ix,iz,iE]) * 𝚽y12[iym] * η
    end
    # Space terms - z
    if mod(iz,2) == 1
        Q[i] += Cz[iz] * (1-ωz[1,ix,iy,iE]) * 𝚽z12[izm] * abs(ξ)
    else
        Q[i] += -Cz[iz] * (1+ωz[1,ix,iy,iE]) * 𝚽z12[izm] * ξ
    end
end

𝚽n = S\Q

# Adaptive correction of weighting parameters
if isAdapt
    error("Not implemented yet.")
    #isFixed, ω = adaptive(3,[𝒪E,𝒪x,𝒪y,𝒪z],[ωE,ωx,ωy,ωz],𝚽n,[𝚽E12,𝚽x12,𝚽y12,𝚽z12],[-1.0,sign(μ),sign(η),sign(ξ)],[Λ,1.0,1.0,1.0],[0.0,0.0,0.0,0.0])
    #ωE = ω[1]; ωx = ω[2]; ωy = ω[3]; ωz = ω[4];
else
    isFixed = true
end

end # End of adaptive loop

# Closure relation
@inbounds for ix in range(1,𝒪x), iy in range(1,𝒪y), iz in range(1,𝒪z), iE in range(1,𝒪E)
    iEm = 𝒪y*𝒪x*(iz-1)+𝒪x*(iy-1)+ix
    ixm = 𝒪y*𝒪E*(iz-1)+𝒪E*(iy-1)+iE
    iym = 𝒪x*𝒪E*(iz-1)+𝒪E*(ix-1)+iE
    izm = 𝒪x*𝒪E*(iy-1)+𝒪E*(ix-1)+iE
    if (iE == 1) 𝚽E12[iEm] = ωE[1,ix,iy,iz] * 𝚽E12[iEm] end
    if (ix == 1) 𝚽x12[ixm] = ωx[1,iy,iz,iE] * 𝚽x12[ixm] end
    if (iy == 1) 𝚽y12[iym] = ωy[1,ix,iz,iE] * 𝚽y12[iym] end
    if (iz == 1) 𝚽z12[izm] = ωz[1,ix,iy,iE] * 𝚽z12[izm] end
end
@inbounds for ix in range(1,𝒪x), iy in range(1,𝒪y), iz in range(1,𝒪z), iE in range(1,𝒪E)
    i = 𝒪y*𝒪x*𝒪E*(iz-1) + 𝒪x*𝒪E * (iy-1) + 𝒪E * (ix-1) + iE
    iEm = 𝒪y*𝒪x*(iz-1)+𝒪x*(iy-1)+ix
    ixm = 𝒪y*𝒪E*(iz-1)+𝒪E*(iy-1)+iE
    iym = 𝒪x*𝒪E*(iz-1)+𝒪E*(ix-1)+iE
    izm = 𝒪x*𝒪E*(iy-1)+𝒪E*(ix-1)+iE
    # Energy terms - E
    if mod(iE,2) == 1
        𝚽E12[iEm] += CE[iE] * ωE[iE+1,ix,iy,iz] * 𝚽n[i]
    else
        𝚽E12[iEm] += -CE[iE] * ωE[iE+1,ix,iy,iz] * 𝚽n[i]
    end  
    # Space terms - x
    if mod(ix,2) == 1
        𝚽x12[ixm] += Cx[ix] * ωx[ix+1,iy,iz,iE] * 𝚽n[i]
    else
        𝚽x12[ixm] += Cx[ix] * ωx[ix+1,iy,iz,iE] * 𝚽n[i] * sign(μ)
    end
    # Space terms - y
    if mod(iy,2) == 1
        𝚽y12[iym] += Cy[iy] * ωy[iy+1,ix,iz,iE] * 𝚽n[i]
    else
        𝚽y12[iym] += Cy[iy] * ωy[iy+1,ix,iz,iE] * 𝚽n[i] * sign(η)
    end
    # Space terms - z
    if mod(iz,2) == 1
        𝚽z12[izm] += Cz[iz] * ωy[iz+1,ix,iy,iE] * 𝚽n[i]
    else
        𝚽z12[izm] += Cz[iz] * ωy[iz+1,ix,iy,iE] * 𝚽n[i] * sign(ξ)
    end
end
𝚽E12 .= 𝚽E12/ΔE

# Returning solutions
return 𝚽n, 𝚽x12, 𝚽y12, 𝚽z12, 𝚽E12

end