"""
    flux_1D_BTE(μ::Float64,Σt::Float64,Δx::Float64,Qn::Vector{Float64},𝚽x12::Float64,
    𝒪x::Int64,C::Vector{Float64},ωx::Vector{Float64},isAdapt::Bool)

Compute flux solution in a cell in 1D Cartesian geometry for the Boltzmann transport
equation.

# Input Argument(s)
- `μ::Float64`: direction cosine.
- `Σt::Float64`: total cross-sections.
- `Δx::Float64`: size of voxels along x-axis.
- `Qn::Vector{Float64}`: angular in-cell source.
- `𝚽x12::Vector{Float64}`: incoming angular flux along x-axis.
- `𝒪x::Int64`: spatial closure relation order.
- `C::Vector{Float64}`: constants related to normalized Legendre.
- `ωx::Vector{Float64}`: weighting factors of the x-axis scheme.
- `isAdapt::Bool`: boolean for adaptive calculations.

# Output Argument(s)
- `𝚽n::Vector{Float64}`: angular in-cell flux.
- `𝚽x12::Float64`: outgoing angular flux along x-axis.

# Reference(s)
N/A

"""
function flux_1D_BTE(μ::Float64,Σt::Float64,Δx::Float64,Qn::Vector{Float64},𝚽x12::Float64,𝒪x::Int64,C::Vector{Float64},ωx::Vector{Float64},isAdapt::Bool)

# Initialization
sx = sign(μ)
hx = abs(μ)/Δx
S = zeros(𝒪x,𝒪x)
Q = zeros(𝒪x)
𝚽n = Q

# Adaptive weight calculations
if isAdapt ωx = adaptive(𝒪x,ωx,hx,sx,𝚽x12,Qn,Σt) end

# Matrix of Legendre moment coefficients of the flux
@inbounds for ix in range(1,𝒪x), jx in range(1,𝒪x)
    if (ix == jx) S[ix,jx] += Σt end
    if (ix ≥ jx + 1) S[ix,jx] -= C[ix] * hx * sx * C[jx] * (1-(-1)^(ix-jx)) end
    S[ix,jx] += C[ix] * hx * sx^(ix-1) * C[jx] * sx^(jx-1) * ωx[jx+1]
end

# Source vector
@inbounds for jx in range(1,𝒪x)
    Q[jx] += Qn[jx]
    Q[jx] -= C[jx] * hx * (sx^(jx-1) * ωx[1] - (-sx)^(jx-1)) * 𝚽x12
end

# Solve the equation system
𝚽n = S\Q

# Closure relation
𝚽x12 = ωx[1] * 𝚽x12
@inbounds for jx in range(1,𝒪x)
    𝚽x12 += C[jx] * sx^(jx-1) * ωx[jx+1] * 𝚽n[jx]
end

# Returning solutions
return 𝚽n, 𝚽x12

end