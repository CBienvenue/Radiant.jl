"""
    compute_sweep_2D(𝚽ℓ::Array{Float64,4},Qℓ::Array{Float64,4},Σt::Vector{Float64},
    mat::Array{Int64,2},Ns::Vector{Int64},Δs::Vector{Vector{Float64}},Ω::Vector{Float64},
    Mn::Vector{Float64},Dn::Vector{Float64},P::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},
    C::Vector{Vector{Float64}},ω::Vector{Array{Float64}},
    sources::Vector{Union{Float64,Array{Float64}}},isAdapt::Vector{Bool},isCSD::Bool,ΔE::Float64,
    𝚽E12::Array{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64})

Compute the flux solution along one direction in 2D geometry.

# Input Argument(s)
- `𝚽ℓ::Array{Float64,4}`: Legendre components of the in-cell flux.
- `Qℓ::Array{Float64,4}`: Legendre components of the in-cell source.
- `Σt::Vector{Float64}`: total cross-sections.
- `mat::Array{Int64,2}`: material identifier per voxel.
- `Ns::Vector{Int64}`: number of voxels along x- and y-axis.
- `Δs::Vector{Vector{Float64}}`: size of voxels along x- and y-axis.
- `Ω::Vector{Float64}`: direction cosines μ and η.
- `Mn::Vector{Float64}`: moment-to-discrete matrix.
- `Dn::Vector{Float64}`: discrete-to-moment matrix.
- `P::Int64`: number of angular interpolation basis.
- `𝒪::Vector{Int64}`: spatial and/or energy closure relation order.
- `Nm::Vector{Int64}`: number of spatial and/or energy moments.
- `C::Vector{Float64}`: constants related to the spatial and energy normalized
   Legendre expansion.
- `ω::Vector{Array{Float64}}`: weighting factors of the closure relations.
- `sources::Vector{Union{Float64, Array{Float64}}}`: surface sources intensities.
- `isAdapt::Bool`: boolean for adaptive calculations.
- `isCSD::Bool`: boolean to indicate if continuous slowing-down term is treated in
   calculations.
- `ΔE::Float64`: energy group width.
- `𝚽E12::Array{Float64}`: incoming flux along the energy axis.
- `S⁻::Vector{Float64}`: stopping power at higher energy group boundary.
- `S⁺::Vector{Float64}`: stopping power at lower energy group boundary.
- `S::Array{Float64}`: stopping powers.
- `𝒲::Array{Float64}`: weighting constants.

# Output Argument(s)
- `𝚽ℓ::Array{Float64}`: Legendre components of the in-cell flux.
- `𝚽E12::Array{Float64}`: outgoing flux along the energy axis.

# Reference(s)
N/A

"""
function compute_sweep_2D(𝚽ℓ::Array{Float64,4},Qℓ::Array{Float64,4},Σt::Vector{Float64},mat::Array{Int64,2},Ns::Vector{Int64},Δs::Vector{Vector{Float64}},Ω::Vector{Float64},Mn::Vector{Float64},Dn::Vector{Float64},P::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},C::Vector{Float64},ω::Vector{Array{Float64}},sources::Vector{Union{Float64,Array{Float64}}},isAdapt::Bool,isCSD::Bool,ΔE::Float64,𝚽E12::Array{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},𝒲::Array{Float64})

# Initialization
𝒪x = 𝒪[1]; 𝒪y = 𝒪[2]; 𝒪E = 𝒪[4]
μ = Ω[1]; η = Ω[2]
Δx = Δs[1]; Δy = Δs[2]
Nx = Ns[1]; Ny = Ns[2]
if (μ >= 0) x_sweep = (1:Nx) else x_sweep = (Nx:-1:1) end
if (η >= 0) y_sweep = (1:Ny) else y_sweep = (Ny:-1:1) end

# Sweep over x-axis
𝚽x12 = zeros(Nm[1],Ny)
@inbounds for ix in x_sweep
𝚽y12 = zeros(Nm[2])
if η >= 0
    if sources[3] != 0  # Surface Y-
        𝚽y12[1,1] += sources[3][ix]
    end
else
    if sources[4] != 0  # Surface Y+
        𝚽y12[1,1] += sources[4][ix]
    end
end

# Sweep over y-axis
for iy in y_sweep
if (ix == 1 &&  μ >= 0) || (ix == Nx && μ < 0 )
    if μ >= 0
        if sources[1] != 0  # Surface X-
            𝚽x12[1,iy] += sources[1][iy]
        end
    else
        if sources[2] != 0  # Surface X+
            𝚽x12[1,iy] += sources[2][iy]
        end
    end
end

# Source term
Qn = zeros(Nm[5])
for is in range(1,Nm[5]), p in range(1,P)
    Qn[is] += Mn[p] * Qℓ[p,is,ix,iy]
end

# Flux calculation
if ~isCSD
    𝚽n,𝚽x12[:,iy],𝚽y12 = flux_2D_BTE(μ,η,Σt[mat[ix,iy]],Δx[ix],Δy[iy],Qn,𝚽x12[:,iy],𝚽y12,𝒪x,𝒪y,C,copy(ω[1]),copy(ω[2]),isAdapt)
else
    𝚽n,𝚽x12[:,iy],𝚽y12,𝚽E12[:,ix,iy] = flux_2D_BFP(μ,η,Σt[mat[ix,iy]],S⁻[mat[ix,iy]],S⁺[mat[ix,iy]],S[mat[ix,iy],:],ΔE,Δx[ix],Δy[iy],Qn,𝚽x12[:,iy],𝚽y12,𝚽E12[:,ix,iy],𝒪E,𝒪x,𝒪y,C,copy(ω[1]),copy(ω[2]),copy(ω[3]),isAdapt,𝒲)
end

# Calculation of the Legendre components of the flux
@inbounds for is in range(1,Nm[5]), p in range(1,P)
    𝚽ℓ[p,is,ix,iy] += Dn[p] * 𝚽n[is]
end

end
end

# Return flux
return 𝚽ℓ, 𝚽E12

end