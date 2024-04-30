"""
    flux_1D_BTE(μ::Float64,Σt::Float64,Δx::Float64,Qn::Vector{Float64},𝚽x12::Float64,
    𝒪x::Int64,Cx::Vector{Float64},ωx::Vector{Float64},isAdaptx::Bool)

Compute flux solution in a cell in 1D Cartesian geometry for the Boltzmann transport
equation.

# Input Argument(s)
- 'μ::Float64': direction cosine.
- 'Σt::Float64': total cross-sections.
- 'Δx::Float64': size of voxels along x-axis.
- 'Qn::Vector{Float64}': angular in-cell source.
- '𝚽x12::Vector{Float64}': incoming angular flux along x-axis.
- '𝒪x::Int64': spatial closure relation order.
- 'Cx::Vector{Float64}': constants related to normalized Legendre.
- 'ωx::Vector{Float64}': weighting factors of the x-axis scheme.
- 'isAdaptx::Bool': boolean for adaptive calculations.

# Output Argument(s)
- '𝚽n::Vector{Float64}': angular in-cell flux.
- '𝚽x12::Float64': outgoing angular flux along x-axis.

# Reference(s)
N/A

"""
function flux_1D_BTE(μ::Float64,Σt::Float64,Δx::Float64,Qn::Vector{Float64},𝚽x12::Float64,𝒪x::Int64,Cx::Vector{Float64},ωx::Vector{Float64},isAdaptx::Bool)

# Initialization
S = zeros(𝒪x,𝒪x)
Q = zeros(𝒪x)
𝚽n = Q

# Adaptive loop
isFixed = false
while ~isFixed

# Matrix of Legendre moment coefficients of the flux
@inbounds for ix in range(1,𝒪x), jx in range(1,𝒪x)
    # Diagonal terms
    if ix == jx
        S[ix,jx] = Σt * Δx + Cx[ix]^2 * ωx[jx+1] * abs(μ)
    # Upper diagonal terms
    elseif ix < jx
    if mod(ix+jx,2) == 1
        S[ix,jx] = Cx[ix] * Cx[jx] * ωx[jx+1] * μ
    else
        S[ix,jx] = Cx[ix] * Cx[jx] * ωx[jx+1] * abs(μ)
    end
    # Under diagonal terms
    else
    if mod(ix+jx,2) == 1
        S[ix,jx] = Cx[ix] * Cx[jx] * (ωx[jx+1]-2) * μ
    else
        S[ix,jx] = Cx[ix] * Cx[jx] * ωx[jx+1] * abs(μ)
    end
    end
end

# Source vector
@inbounds for ix in range(1,𝒪x)
    Q[ix] = Qn[ix] * Δx
    if mod(ix,2) == 1
        Q[ix] += Cx[ix] * (1-ωx[1]) * 𝚽x12 * abs(μ)
    else
        Q[ix] += -Cx[ix] * (1+ωx[1]) * 𝚽x12 * μ
    end
end

𝚽n = S\Q

# Adaptive correction of weighting parameters
if isAdaptx
    isFixed, ωx = adaptive_1D(𝒪x,ωx,𝚽n,𝚽x12,sign(μ),1.0)
else
    isFixed = true
end

end # End of adaptive loop

# Closure relation
𝚽x12 = ωx[1] * 𝚽x12
@inbounds for ix in range(1,𝒪x)
    if mod(ix,2) == 1
        𝚽x12 += Cx[ix] * ωx[ix+1] * 𝚽n[ix]
    else
        𝚽x12 += Cx[ix] * ωx[ix+1] * 𝚽n[ix] * sign(μ)
    end
end

# Returning solutions
return 𝚽n, 𝚽x12

end