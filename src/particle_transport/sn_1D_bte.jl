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
isTangential = (abs(μ) ≤ 1e-10)
if isTangential sx = 1.0; μ = 0.0 end
hx = abs(μ)/Δx
𝒮 = zeros(𝒪x,𝒪x)
Q = zeros(𝒪x)
𝚽n = Q

# Adaptive weight calculations
if isAdapt && !isTangential ωx = adaptive(𝒪x,ωx,hx,sx,𝚽x12,Qn,Σt) end

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,𝒪x), jx in range(1,𝒪x)
    if (ix == jx) 𝒮[ix,jx] += Σt end
    if (ix ≥ jx + 1) 𝒮[ix,jx] -= C[ix] * hx * sx * C[jx] * (1-(-1)^(ix-jx)) end
    𝒮[ix,jx] += C[ix] * hx * sx^(ix-1) * C[jx] * sx^(jx-1) * ωx[jx+1]
end

# Source vector
for jx in range(1,𝒪x)
    Q[jx] += Qn[jx]
    Q[jx] -= C[jx] * hx * (sx^(jx-1) * ωx[1] - (-sx)^(jx-1)) * 𝚽x12
end

# Solve the equation system
𝚽n = 𝒮\Q
if !isTangential
    𝚽x12 = ωx[1] * 𝚽x12
    for jx in range(1,𝒪x)
        𝚽x12 += C[jx] * sx^(jx-1) * ωx[jx+1] * 𝚽n[jx]
    end
end

# Returning solutions
return 𝚽n, 𝚽x12

end
"""
    sn_1D_BTE_matrix!(...)

Assemble the cell matrix `𝒮` of `flux_1D_BTE`, for the optimized solver chain.

The body is the reference kernel's assembly loop transcribed unchanged, so the matrix is
bit-identical to the one it builds. It is split out because the matrix does not depend on the
voxel beyond its widths and its material: `sn_fast_context` calls this once per (material,
mesh-width combination, direction) and factorizes the result, where the reference rebuilds and
refactorizes it in every cell of every sweep. Valid only when `~isAdapt`.
"""
function sn_1D_BTE_matrix!(𝒮::Matrix{Float64},μ::Float64,Σt::Float64,Δx::Float64,𝒪x::Int64,C::Vector{Float64},ωx::Array{Float64},isFC::Bool)

fill!(𝒮,0.0)
sx = sign(μ)
isTangential = (abs(μ) ≤ 1e-10)
if isTangential sx = 1.0; μ = 0.0 end
hx = abs(μ)/Δx

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,𝒪x), jx in range(1,𝒪x)
    if (ix == jx) 𝒮[ix,jx] += Σt end
    if (ix ≥ jx + 1) 𝒮[ix,jx] -= C[ix] * hx * sx * C[jx] * (1-(-1)^(ix-jx)) end
    𝒮[ix,jx] += C[ix] * hx * sx^(ix-1) * C[jx] * sx^(jx-1) * ωx[jx+1]
end

return nothing

end

"""
    sn_1D_BTE_fast!(𝚿,o𝚿,𝚽x12,ox,Ql,oQ,Mnn,Np,mom,d,ws,conf,ikx)

Optimized counterpart of `flux_1D_BTE`: solve one voxel of the 1D BTE.

Identical in structure to `sn_3D_BTE_fast2!` with the y and z axes collapsed. The reference's
tangential case (`|μ| ≤ 1e-10`, which skips the closure) is resolved once per direction when the
context is built, and appears here as an empty closure relation.
"""
@inline function sn_1D_BTE_fast!(𝚿::Vector{Float64},o𝚿::Int64,𝚽x12::Vector{Float64},ox::Int64,Ql::Vector{Float64},oQ::Int64,Mnn::Vector{Float64},Np::Int64,mom::GNFastMoments{NMOM},d::SNFastDir{NMOM},ws::GNFastWorkspace,conf::Int64,ikx::Int64) where {NMOM}

    Q = ws.Q; 𝚽 = ws.𝚽
    @inbounds for k in 1:NMOM
        col = mom.col[k]
        q = 0.0; b = oQ + Np*(col-1)
        @simd for p in 1:Np
            q += Mnn[p] * Ql[b+p]
        end
        Q[col] = q + d.qcx[k,ikx] * 𝚽x12[ox+mom.cyz[k]]
    end
    gn_fast_solve!(𝚽,d.LU,d.ipiv,Q,NMOM,conf)
    gn_fast_closure!(𝚽x12,ox,d.clx,𝚽,Val(NMOM),Val(1))
    @inbounds for c in 1:NMOM
        𝚿[o𝚿+c] = 𝚽[c]
    end
    return nothing
end
