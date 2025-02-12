"""
    flux_1D_BFP(isFC::Bool,μ::Float64,Σt::Float64,Δx::Float64,Qn::Vector{Float64},
    𝚽x12::Vector{Float64},S⁻::Float64,S⁺::Float64,S::Vector{Float64},ΔE::Float64,
    𝚽E12::Vector{Float64},𝒪E::Int64,𝒪x::Int64,C::Vector{Float64},ωE::Array{Float64},
    ωx::Array{Float64},isAdapt::Bool,𝒲::Array{Float64})

Compute flux solution in a cell in 1D Cartesian geometry for the Boltzmann Fokker-Planck
equation.

# Input Argument(s)
- 'isFC::Bool': boolean to indicate if full coupling or not.
- 'μ::Float64': direction cosine.
- 'Σt::Float64': total cross-sections.
- 'Δx::Float64': size of voxels along x-axis.
- 'Qn::Vector{Float64}': angular in-cell source.
- '𝚽x12::Vector{Float64}': incoming angular flux along x-axis.
- 'S⁻::Float64': stopping powers at upper energy group boundary.
- 'S⁺::Float64': stopping powers at lower energy group boundary.
- 'S::Vector{Float64}': stopping powers.
- 'ΔE::Float64': energy group width.
- '𝚽E12::Vector{Float64}': incoming angular flux along E-axis.
- '𝒪E::Int64': energy closure relation order.
- '𝒪x::Int64': spatial closure relation order.
- 'C::Vector{Float64}': constants related to normalized Legendre.
- 'ωE::Array{Float64}': weighting factors of the E-axis scheme.
- 'ωx::Array{Float64}': weighting factors of the x-axis scheme.
- 'isAdapt::Bool': boolean for adaptive calculations.
- '𝒲::Array{Float64}' : weighting constants.

# Output Argument(s)
- '𝚽n::Vector{Float64}': angular in-cell flux.
- '𝚽x12::Vector{Float64}': outgoing angular flux along x-axis.
- '𝚽E12::Vector{Float64}': outgoing angular flux along E-axis.

# Reference(s)
N/A

"""
function flux_1D_BFP(isFC::Bool,μ::Float64,Σt::Float64,Δx::Float64,Qn::Vector{Float64},𝚽x12::Vector{Float64},S⁻::Float64,S⁺::Float64,S::Vector{Float64},ΔE::Float64,𝚽E12::Vector{Float64},𝒪E::Int64,𝒪x::Int64,C::Vector{Float64},ωE::Array{Float64},ωx::Array{Float64},isAdapt::Bool,𝒲::Array{Float64})

# Initialization
sx = sign(μ)
hx = abs(μ)/Δx
Nm = 𝒪x*𝒪E
𝒮 = zeros(Nm,Nm)
Q = zeros(Nm)
𝚽n = Q

# Adaptive weight calculations
if isAdapt ωx,ωE = adaptive(𝒪x,𝒪E,ωx,ωE,hx,1/ΔE,sx,-1,𝚽x12,𝚽E12,Qn,Σt) end

# Matrix of Legendre moment coefficients of the flux
@inbounds for ix in range(1,𝒪x), jx in range(1,𝒪x), iE in range(1,𝒪E), jE in range(1,𝒪E)
    i = 𝒪E*(ix-1)+iE
    j = 𝒪E*(jx-1)+jE

    # Collision term
    if (i == j) 𝒮[i,j] += Σt end

    # Streaming term - x
    if iE == jE
        if (ix ≥ jx + 1) 𝒮[i,j] -= C[ix] * hx * sx * C[jx] * (1-(-1)^(ix-jx)) end
    end
    𝒮[i,j] += C[ix] * hx * sx^(ix-1) * C[jx] * sx^(jx-1) * ωx[jx+1,jE,iE]

    # CSD term
    if ix == jx
        for kE in range(1,iE-1), wE in range(1,𝒪E)
            𝒮[i,j] += C[iE] * C[jE] * C[kE] * C[wE] * (1-(-1)^(iE-kE)) * S[wE] * 𝒲[jE,kE,wE]
        end
    end
    𝒮[i,j] += C[iE] * S⁺ * (-1)^(iE-1) * C[jE] * (-1)^(jE-1) * ωE[jE+1,jx,ix]

end

# Source vector
@inbounds for jx in range(1,𝒪x), jE in range(1,𝒪E)
    j = 𝒪E*(jx-1)+jE
    Q[j] += Qn[j]
    Q[j] -= C[jx] * hx * (sx^(jx-1) * ωx[1,jE,jE] - (-sx)^(jx-1)) * 𝚽x12[jE] 
    Q[j] -= C[jE] * ((-1)^(jE-1)*S⁺*ωE[1,jx,jx] - S⁻) * 𝚽E12[jx]
end

# Solve the equation system
𝚽n = 𝒮\Q

# Closure relations
@inbounds for jx in range(1,𝒪x), jE in range(1,𝒪E)
    j = 𝒪E*(jx-1)+jE
    if (jx == 1) 𝚽x12[jE] = ωx[1,jE,jE] * 𝚽x12[jE] end
    if (jE == 1) 𝚽E12[jx] = ωE[1,jx,jx] * 𝚽E12[jx] end
    for iE in range(1,𝒪E)
        𝚽x12[jE] += C[jx] * sx^(jx-1) * ωx[jx+1,jE,iE] * 𝚽n[j]
    end
    for ix in range(1,𝒪x)
        𝚽E12[jx] += C[jE] * (-1)^(jE-1) * ωE[jE+1,jx,ix] * 𝚽n[j]
    end
end

# Returning solutions
return 𝚽n, 𝚽x12, 𝚽E12
end