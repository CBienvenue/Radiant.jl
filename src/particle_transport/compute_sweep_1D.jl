"""
    compute_sweep_1D(𝚽ℓ::Array{Float64,3},Qℓ::Array{Float64,3},Σt::Vector{Float64},
    mat::Vector{Int64},Nx::Int64,Δx::Vector{Float64},μ::Float64,Mn::Vector{Float64},
    Dn::Vector{Float64},P::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},isFC::Bool
    C::Vector{Vector{Float64}},ω::Vector{Array{Float64}},
    sources::Vector{Union{Float64, Array{Float64}}},isAdapt::Vector{Bool},isCSD::Bool,
    ΔE::Float64,𝚽E12::Array{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64})

Compute the flux solution along one direction in 1D geometry.

# Input Argument(s)
- '𝚽ℓ::Array{Float64,3}': Legendre components of the in-cell flux.
- 'Qℓ::Array{Float64,3}': Legendre components of the in-cell source.
- 'Σt::Vector{Float64}': total cross-sections.
- 'mat::Vector{Int64}': material identifier per voxel.
- 'Nx::Int64': number of voxels along x-axis.
- 'Δx::Vector{Float64}': size of voxels along x-axis.
- 'μ::Float64': direction cosines.
- 'Mn::Vector{Float64}': moment-to-discrete matrix.
- 'Dn::Vector{Float64}': discrete-to-moment matrix.
- 'P::Int64': number of angular interpolation basis.
- '𝒪::Vector{Int64}': spatial and/or energy closure relation order.
- 'Nm::Vector{Int64}': number of spatial and/or energy moments.
- 'isFC::Bool': boolean to indicate if full coupling or not.
- 'C::Vector{Vector{Float64}}': constants related to the spatial and energy normalized
   Legendre expansion.
- 'ω::Vector{Array{Float64}}': weighting factors of the closure relations.
- 'sources::Vector{Union{Float64, Array{Float64}}}': surface sources intensities.
- 'isAdapt::Vector{Bool}': boolean for adaptive calculations.
- 'isCSD::Bool': boolean to indicate if continuous slowing-down term is treated in
   calculations.
- 'ΔE::Float64': energy group width.
- '𝚽E12::Array{Float64}': incoming flux along the energy axis.
- 'S⁻::Vector{Float64}': restricted stopping power at higher energy group boundary.
- 'S⁺::Vector{Float64}': restricted stopping power at lower energy group boundary.

# Output Argument(s)
- '𝚽ℓ::Array{Float64}': Legendre components of the in-cell flux.
- '𝚽E12::Array{Float64}': outgoing flux along the energy axis.

# Reference(s)
N/A

"""
function compute_sweep_1D(𝚽ℓ::Array{Float64,3},Qℓ::Array{Float64,3},Σt::Vector{Float64},mat::Vector{Int64},Nx::Int64,Δx::Vector{Float64},μ::Float64,Mn::Vector{Float64},Dn::Vector{Float64},P::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},isFC::Bool,C,ω,sources::Vector{Union{Float64, Array{Float64}}},isAdapt,isCSD::Bool,ΔE::Float64,𝚽E12::Array{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64},S)

# Initialization
𝒪x = 𝒪[1]
𝒪E = 𝒪[4]
𝚽x12 = zeros(Nm[1])

# Monodirectional boundary sources
if μ >= 0
    x_sweep = 1:Nx
    if sources[1] != 0   # Surface X-
        𝚽x12[1] += sources[1]
    end
else
    x_sweep = Nx:-1:1
    if sources[2] != 0   # Surface X+
        𝚽x12[1] += sources[2]
    end
end

@inbounds for ix in x_sweep

    # Source term
    Qn = zeros(Nm[5])
    for is in range(1,Nm[5]), p in range(1,P)
        Qn[is] += Mn[p] * Qℓ[p,is,ix]
    end

    # Flux calculation
    if ~isCSD
        𝚽n,𝚽x12 = flux_1D_BTE(μ,Σt[mat[ix]],Δx[ix],Qn,𝚽x12[1],𝒪x,C,copy(ω),isAdapt)
    else
        𝚽n,𝚽x12,𝚽E12[:,ix] = flux_1D_BFP(isFC,μ,Σt[mat[ix]],Δx[ix],Qn,𝚽x12,S⁻[mat[ix]],S⁺[mat[ix]],S[mat[ix],:],ΔE,𝚽E12[:,ix],𝒪E,𝒪x,C,copy(ω[1]),copy(ω[2]),isAdapt)
    end

    # Calculation of the Legendre components of the flux
    for is in range(1,Nm[5]), p in range(1,P)
        𝚽ℓ[p,is,ix] += Dn[p] * 𝚽n[is]
    end
    
end

return 𝚽ℓ, 𝚽E12

end