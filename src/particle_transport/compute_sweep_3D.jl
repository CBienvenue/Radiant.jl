"""
    compute_sweep_3D(𝚽ℓ::Array{Float64,5},Qℓ::Array{Float64,5},Σt::Vector{Float64},
    mat::Array{Int64,3},Ns::Vector{Int64},Δs::Vector{Vector{Float64}},Ω::Vector{Float64},
    Mn::Vector{Float64},Dn::Vector{Float64},P::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},
    C::Vector{Vector{Float64}},ω::Vector{Array{Float64}},
    S::Vector{Union{Float64,Array{Float64}}},isAdapt::Vector{Bool},isCSD::Bool,ΔE::Float64,
    𝚽E12::Array{Float64},β⁻::Vector{Float64},β⁺::Vector{Float64})

Compute the flux solution along one direction in 3D geometry.

See also [`compute_one_speed`](@ref), [`compute_sweep_1D`](@ref),
[`compute_sweep_2D`](@ref).

# Input Argument(s)
- '𝚽ℓ::Array{Float64,4}': Legendre components of the in-cell flux.
- 'Qℓ::Array{Float64,4}': Legendre components of the in-cell source.
- 'Σt::Vector{Float64}': total cross-sections.
- 'mat::Array{Int64,2}': material identifier per voxel.
- 'Ns::Vector{Int64}': number of voxels along x- and y-axis.
- 'Δs::Vector{Vector{Float64}}': size of voxels along x- and y-axis.
- 'Ω::Vector{Float64}': direction cosines μ and η.
- 'Mn::Vector{Float64}': moment-to-discrete matrix.
- 'Dn::Vector{Float64}': discrete-to-moment matrix.
- 'P::Int64': number of angular interpolation basis.
- '𝒪::Vector{Int64}': spatial and/or energy closure relation order.
- 'Nm::Vector{Int64}': number of spatial and/or energy moments.
- 'C::Vector{Vector{Float64}}': constants related to the spatial and energy normalized
   Legendre expansion.
- 'ω::Vector{Array{Float64}}': weighting factors of the closure relations.
- 'S::Vector{Union{Float64, Array{Float64}}}': surface sources intensities.
- 'isAdapt::Vector{Bool}': boolean for adaptive calculations.
- 'isCSD::Bool': boolean to indicate if continuous slowing-down term is treated in
   calculations.
- 'ΔE::Float64': energy group width.
- '𝚽E12::Array{Float64}': incoming flux along the energy axis.
- 'β⁻::Vector{Float64}': restricted stopping power at higher energy group boundary.
- 'β⁺::Vector{Float64}': restricted stopping power at lower energy group boundary.

# Output Argument(s)
- '𝚽ℓ::Array{Float64}': Legendre components of the in-cell flux.
- '𝚽E12::Array{Float64}': outgoing flux along the energy axis.

# Author(s)
Charles Bienvenue

# Reference(s)
N/A

"""
function compute_sweep_3D(𝚽ℓ::Array{Float64,5},Qℓ::Array{Float64,5},Σt::Vector{Float64},mat::Array{Int64,3},Ns::Vector{Int64},Δs::Vector{Vector{Float64}},Ω::Vector{Float64},Mn::Vector{Float64},Dn::Vector{Float64},P::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},C::Vector{Vector{Float64}},ω::Vector{Array{Float64}},S::Vector{Union{Float64,Array{Float64}}},isAdapt::Vector{Bool},isCSD::Bool,ΔE::Float64,𝚽E12::Array{Float64},β⁻::Vector{Float64},β⁺::Vector{Float64})

# Initialization
𝒪x = 𝒪[1]; 𝒪y = 𝒪[2]; 𝒪z = 𝒪[3]; 𝒪E = 𝒪[4]
μ = Ω[1]; η = Ω[2]; ξ = Ω[3]
Δx = Δs[1]; Δy = Δs[2]; Δz = Δs[3]
Nx = Ns[1]; Ny = Ns[2]; Nz = Ns[3]
if (μ >= 0) x_sweep = (1:Nx) else x_sweep = (Nx:-1:1) end
if (η >= 0) y_sweep = (1:Ny) else y_sweep = (Ny:-1:1) end
if (ξ >= 0) z_sweep = (1:Nz) else z_sweep = (Nz:-1:1) end

# Sweeping over x-axis
𝚽12x = zeros(Nm[1],Ny,Nz)
@inbounds for ix in x_sweep

# Sweeping over y-axis
𝚽12y = zeros(Nm[2],Nz)
@inbounds for iy in y_sweep
𝚽12z = zeros(Nm[3])
if ξ >= 0
    if S[5] != 0 # Surface Z-
        𝚽12z[1] += S[5][(iy-1)*Nx+ix]
    end
else
    if S[6] != 0  # Surface Z+
        𝚽12z[1] += S[6][(iy-1)*Nx+ix]
    end
end

# Sweeping over z-axis
@inbounds for iz in z_sweep
if (iy == 1 &&  η >= 0) || (iy == Ny && η < 0 )
    if η >= 0
        if S[3] != 0  # Surface Y-
            𝚽12y[1,iz] += S[3][(iz-1)*Nx+ix]
        end
    else
        if S[4] != 0  # Surface Y+
            𝚽12y[1,iz] += S[4][(iz-1)*Nx+ix]
        end
    end
end
if (ix == 1 && μ >= 0) || (ix == Nx && μ < 0 )
    if μ >= 0
        if S[1] != 0  # Surface X-
            𝚽12x[1,iy,iz] += S[1][(iz-1)*Ny+iy]
        end
    else
        if S[2] != 0  # Surface X+
            𝚽12x[1,iy,iz] += S[2][(iz-1)*Ny+iy]
        end
    end
end

# Source term
Qn = zeros(Nm[5])
@inbounds for is in range(1,Nm[5]), p in range(1,P)
    Qn[is] += Mn[p] * Qℓ[p,is,ix,iy,iz]
end

# Flux calculation
if ~isCSD
    𝚽n,𝚽12x[:,iy,iz],𝚽12y[:,iz],𝚽12z = flux_3D_BTE(μ,η,ξ,Σt[mat[ix,iy,iz]],Δx[ix],Δy[iy],Δz[iz],Qn,𝚽12x[:,iy,iz],𝚽12y[:,iz],𝚽12z,𝒪x,𝒪y,𝒪z,C[1],C[2],C[3],copy(ω[1][:,:,:,1]),copy(ω[2][:,:,:,1]),copy(ω[3][:,:,:,1]),isAdapt[1],isAdapt[2],isAdapt[3])
else
    𝚽n,𝚽12x[:,iy,iz],𝚽12y[:,iz],𝚽12z,𝚽E12[:,ix,iy,iz] = flux_3D_BFP(μ,η,ξ,Σt[mat[ix,iy,iz]],β⁻[mat[ix,iy,iz]],β⁺[mat[ix,iy,iz]],ΔE,Δx[ix],Δy[iy],Δz[iz],Qn,𝚽12x[:,iy,iz],𝚽12y[:,iz],𝚽12z,𝚽E12[:,ix,iy,iz],𝒪E,𝒪x,𝒪y,𝒪z,C[4],C[1],C[2],C[3],copy(ω[4]),copy(ω[1]),copy(ω[2]),copy(ω[3]),isAdapt[4],isAdapt[1],isAdapt[2],isAdapt[3])
end

# Calculation of the Legendre components of the flux
@inbounds for is in range(1,Nm[5]), p in range(1,P)
    𝚽ℓ[p,is,ix,iy,iz] += Dn[p] * 𝚽n[is]
end

end
end
end

# Return flux
return 𝚽ℓ, 𝚽E12

end