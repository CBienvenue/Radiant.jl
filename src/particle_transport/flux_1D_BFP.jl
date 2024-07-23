"""
    flux_1D_BFP(isFC::Bool,μ::Float64,Σt::Float64,Δx::Float64,Qn::Vector{Float64},
    𝚽x12::Vector{Float64},β⁻::Float64,β⁺::Float64,ΔE::Float64,𝚽E12::Vector{Float64},
    𝒪E::Int64,𝒪x::Int64,CE::Vector{Float64},Cx::Vector{Float64},ωE::Array{Float64},
    ωx::Array{Float64},isAdaptE::Bool,isAdaptx::Bool)

Compute flux solution in a cell in 1D Cartesian geometry for the Boltzmann Fokker-Planck
equation.

# Input Argument(s)
- 'isFC::Bool': boolean to indicate if full coupling or not.
- 'μ::Float64': direction cosine.
- 'Σt::Float64': total cross-sections.
- 'Δx::Float64': size of voxels along x-axis.
- 'Qn::Vector{Float64}': angular in-cell source.
- '𝚽x12::Vector{Float64}': incoming angular flux along x-axis.
- 'β⁻::Float64': restricted stopping power at upper energy group boundary.
- 'β⁺::Float64': restricted stopping power at lower energy group boundary.
- 'ΔE::Float64': energy group width.
- '𝚽E12::Vector{Float64}': incoming angular flux along E-axis.
- '𝒪E::Int64': energy closure relation order.
- '𝒪x::Int64': spatial closure relation order.
- 'CE::Vector{Float64}': constants related to normalized Legendre.
- 'Cx::Vector{Float64}': constants related to normalized Legendre.
- 'ωE::Array{Float64}': weighting factors of the E-axis scheme.
- 'ωx::Array{Float64}': weighting factors of the x-axis scheme.
- 'isAdaptE::Bool': boolean for adaptive calculations.
- 'isAdaptx::Bool': boolean for adaptive calculations.

# Output Argument(s)
- '𝚽n::Vector{Float64}': angular in-cell flux.
- '𝚽x12::Vector{Float64}': outgoing angular flux along x-axis.
- '𝚽E12::Vector{Float64}': outgoing angular flux along E-axis.

# Reference(s)
N/A

"""
function flux_1D_BFP(isFC::Bool,μ::Float64,Σt::Float64,Δx::Float64,Qn::Vector{Float64},𝚽x12::Vector{Float64},β⁻::Float64,β⁺::Float64,ΔE::Float64,𝚽E12::Vector{Float64},𝒪E::Int64,𝒪x::Int64,CE::Vector{Float64},Cx::Vector{Float64},ωE::Array{Float64},ωx::Array{Float64},isAdaptE::Bool,isAdaptx::Bool)

#----
# Full coupling
#----
if isFC

# Initialization
S = zeros(𝒪x*𝒪E,𝒪x*𝒪E)
Q = zeros(𝒪x*𝒪E)
𝚽n = Q
TE = 0.0; Tx = 0.0

# Galerkin energy scheme weights
Λ = β⁻/β⁺
if abs(ωE[1,1]) > 0
    ωE[1,1] = ωE[1,1]*Λ
    ωE[2:𝒪E+1,1] = (ωE[2:𝒪E+1,1].-1).*Λ.+1
end

# Adaptive loop
isAdapt = isAdaptE && isAdaptx
isFixed = false
while ~isFixed

# Matrix of Legendre moment coefficients of the flux
@inbounds for ix in range(1,𝒪x), jx in range(1,𝒪x), iE in range(1,𝒪E), jE in range(1,𝒪E)
    i = 𝒪E*(ix-1)+iE
    j = 𝒪E*(jx-1)+jE
    # Diagonal terms
    if i == j
        S[i,j] = ( Σt + CE[iE]^2 * β⁺ * ωE[jE+1,jx] + (iE-1) * (β⁻-β⁺) ) * Δx + Cx[ix]^2 * ωx[jx+1,jE] * abs(μ)
    # Upper diagonal terms
    elseif i < j
    # Energy terms - E
    if ix == jx
    if mod(iE+jE,2) == 1
        S[i,j] = -CE[iE] * CE[jE] * β⁺ * ωE[jE+1,jx] * Δx
    else
        S[i,j] = CE[iE] * CE[jE] * β⁺ * ωE[jE+1,jx] * Δx
    end
    # Space terms - x
    elseif  iE == jE 
    if mod(ix+jx,2) == 1
        S[i,j] = Cx[ix] * Cx[jx] * ωx[jx+1,jE] * μ
    else
        S[i,j] = Cx[ix] * Cx[jx] * ωx[jx+1,jE] * abs(μ)
    end
    end
    # Under diagonal terms
    else
    # Energy terms - E
    if ix == jx
    if mod(iE+jE,2) == 1
        S[i,j] = -CE[iE] * CE[jE] * (β⁺*ωE[jE+1,jx]-β⁻-β⁺) * Δx
    else
        S[i,j] = CE[iE] * CE[jE] * (β⁺*ωE[jE+1,jx]+β⁻-β⁺) * Δx
    end
    # Space terms - x
    elseif  iE == jE 
    if mod(ix+jx,2) == 1
        S[i,j] = Cx[ix] * Cx[jx] * (ωx[jx+1,jE]-2) * μ
    else
        S[i,j] = Cx[ix] * Cx[jx] * ωx[jx+1,jE] * abs(μ)
    end
    end
    end

    # T-factors
    if 𝒪x == 2 && 𝒪E == 2
        if i == 2 && j == 3
            S[i,j] += abs(μ) * Tx
        elseif i == 3 && j == 2
            S[i,j] += β⁺ * Δx * TE
        elseif i == 4 && j == 2
            S[i,j] += sqrt(3) * sign(-1.0) * β⁺ * Δx * TE
        elseif i == 4 && j == 3
            S[i,j] += sqrt(3) * sign(μ) * abs(μ) * Tx
        end
    end

end

# Source vector
@inbounds for ix in range(1,𝒪x), iE in range(1,𝒪E)
    i = 𝒪E*(ix-1)+iE
    Q[i] = Qn[i] * Δx
    # Energy terms - E
    if mod(iE,2) == 1
        Q[i] += CE[iE] * (β⁻-β⁺*ωE[1,ix]) * 𝚽E12[ix] * Δx
    else
        Q[i] += CE[iE] * (β⁻+β⁺*ωE[1,ix]) * 𝚽E12[ix] * Δx
    end
    # Space terms - x
    if mod(ix,2) == 1
        Q[i] += Cx[ix] * (1-ωx[1,iE]) * 𝚽x12[iE] * abs(μ)
    else
        Q[i] += -Cx[ix] * (1+ωx[1,iE]) * 𝚽x12[iE] * μ
    end
end

𝚽n = S\Q

# Adaptive correction of weighting parameters
if isAdapt
    isFixed, ω, T = adaptive_2D([𝒪E,𝒪x],[ωE,ωx],𝚽n,[𝚽E12,𝚽x12],[-1.0,sign(μ)],[Λ,1.0],[TE,Tx],Qn,[1/(ΔE*Σt),abs(μ)/(Δx*Σt)],Q)
    ωE = ω[1]; ωx = ω[2]; TE = T[1]; Tx = T[2];
else
    isFixed = true
end

end # End of adaptive loop

# Closure relation
@inbounds for ix in range(1,𝒪x), iE in range(1,𝒪E)
    if (iE == 1) 𝚽E12[ix] = ωE[1,ix] * 𝚽E12[ix] end
    if (ix == 1) 𝚽x12[iE] = ωx[1,iE] * 𝚽x12[iE] end
end
@inbounds for ix in range(1,𝒪x), iE in range(1,𝒪E)
    i = 𝒪E*(ix-1)+iE
    # Energy terms - E
    if mod(iE,2) == 1
        𝚽E12[ix] += CE[iE] * ωE[iE+1,ix] * 𝚽n[i]
    else
        𝚽E12[ix] += -CE[iE] * ωE[iE+1,ix] * 𝚽n[i]
    end  
    # Space terms - x
    if mod(ix,2) == 1
        𝚽x12[iE] += Cx[ix] * ωx[ix+1,iE] * 𝚽n[i]
    else
        𝚽x12[iE] += Cx[ix] * ωx[ix+1,iE] * 𝚽n[i] * sign(μ)
    end
end
if 𝒪x == 2 && 𝒪E == 2
    𝚽E12[2] += TE * 𝚽n[2]
    𝚽x12[2] += Tx * 𝚽n[3]
end
𝚽E12 .= 𝚽E12/ΔE

# Returning solutions
return 𝚽n, 𝚽x12, 𝚽E12

#----
# Partial coupling
#----
else

# Initialization
S = zeros(𝒪x+𝒪E-1,𝒪x+𝒪E-1)
Q = zeros(𝒪x+𝒪E-1)
𝚽n = Q
TE = 0.0; Tx = 0.0

# Galerkin energy scheme weights
Λ = β⁻/β⁺
if abs(ωE[1,1]) > 0
    ωE[1,1] = ωE[1,1]*Λ
    ωE[2:𝒪E+1,1] = (ωE[2:𝒪E+1,1].-1).*Λ.+1
end

# Adaptive loop
isAdapt = isAdaptE && isAdaptx
isFixed = false
while ~isFixed

# Matrix of Legendre moment coefficients of the flux
S[1,1] = ( Σt + CE[1]^2 * β⁺ * ωE[2,1] + (β⁻-β⁺) ) * Δx + Cx[1]^2 * ωx[2,1] * abs(μ)
for iE in range(1,𝒪E), jE in range(1,𝒪E)
    i = iE
    j = jE
    # Diagonal terms
    if i == j && i != 1
        S[i,j] = ( Σt + CE[iE]^2 * β⁺ * ωE[jE+1,1] + (iE-1) * (β⁻-β⁺) ) * Δx
        S[i,j] += abs(μ)
    # Upper diagonal terms
    elseif i < j
    if mod(iE+jE,2) == 1
        S[i,j] = -CE[iE] * CE[jE] * β⁺ * ωE[jE+1,1] * Δx
    else
        S[i,j] = CE[iE] * CE[jE] * β⁺ * ωE[jE+1,1] * Δx
    end
    # Under diagonal terms
    elseif i > j 
    if mod(iE+jE,2) == 1
        S[i,j] = -CE[iE] * CE[jE] * (β⁺*ωE[jE+1,1]-β⁻-β⁺) * Δx
    else
        S[i,j] = CE[iE] * CE[jE] * (β⁺*ωE[jE+1,1]+β⁻-β⁺) * Δx
    end
    end
end
for ix in range(1,𝒪x), jx in range(1,𝒪x)
    i = ix
    j = jx
    if ix > 1 i += 𝒪E-1 end
    if jx > 1 j += 𝒪E-1 end
    # Diagonal terms
    if i == j && i != 1
        S[i,j] = Σt * Δx + Cx[ix]^2 * ωx[jx+1,1] * abs(μ)
        S[i,j] += β⁻ * Δx
    # Upper diagonal terms
    elseif i < j
    if mod(ix+jx,2) == 1
        S[i,j] = Cx[ix] * Cx[jx] * ωx[jx+1,1] * μ
    else
        S[i,j] = Cx[ix] * Cx[jx] * ωx[jx+1,1] * abs(μ)
    end
    # Under diagonal terms
    elseif i > j
    if mod(ix+jx,2) == 1
        S[i,j] = Cx[ix] * Cx[jx] * (ωx[jx+1,1]-2) * μ
    else
        S[i,j] = Cx[ix] * Cx[jx] * ωx[jx+1,1] * abs(μ)
    end
    end
end
# T-factors
if 𝒪x == 2 && 𝒪E == 2
    S[2,3] += abs(μ) * Tx
    S[3,2] += β⁺ * Δx * TE
end

# Source vector
for i in range(1,𝒪x+𝒪E-1)
    if i == 1
        iE = 1; ix = 1;
    elseif i > 1 && i <= 1 + (𝒪E-1)
        iE = i; ix = 1;
    else
        iE = 1; ix = i - (𝒪E-1);
    end
    Q[i] = Qn[i] * Δx
    # Energy terms - E
    if ix == 1
        if mod(iE,2) == 1
            Q[i] += CE[iE] * (β⁻-β⁺*ωE[1,1]) * 𝚽E12[ix] * Δx
        else
            Q[i] += CE[iE] * (β⁻+β⁺*ωE[1,1]) * 𝚽E12[ix] * Δx
        end
    else
       Q[i] += β⁻ * 𝚽E12[ix] * Δx
    end
    # Space terms - x
    if iE == 1
        if mod(ix,2) == 1
            Q[i] += Cx[ix] * (1-ωx[1,1]) * 𝚽x12[iE] * abs(μ)
        else
            Q[i] += -Cx[ix] * (1+ωx[1,1]) * 𝚽x12[iE] * μ
        end
    else
        Q[i] += 𝚽x12[iE] * abs(μ)
    end
end

𝚽n = S\Q

# Adaptive correction of weighting parameters
if isAdapt
    error("Not implemented yet.")
else
    isFixed = true
end

end # End of adaptive loop

# Closure relation
𝚽E12[1] = ωE[1,1] * 𝚽E12[1]
𝚽x12[1] = ωx[1,1] * 𝚽x12[1]
𝚽E12[2:end] .= 0.0
𝚽x12[2:end] .= 0.0
for i in range(1,𝒪x+𝒪E-1)
    if i == 1
        iE = 1; ix = 1;
    elseif i > 1 && i <= 1 + (𝒪E-1)
        iE = i; ix = 1;
    else
        iE = 1; ix = i - (𝒪E-1);
    end
    # Energy terms - E
    if ix == 1
        if mod(iE,2) == 1
            𝚽E12[ix] += CE[iE] * ωE[iE+1,1] * 𝚽n[i]
        else
            𝚽E12[ix] += -CE[iE] * ωE[iE+1,1] * 𝚽n[i]
        end  
    else
        𝚽E12[ix] += 𝚽n[i]
    end
    # Space terms - x
    if iE == 1
        if mod(ix,2) == 1
            𝚽x12[iE] += Cx[ix] * ωx[ix+1,1] * 𝚽n[i]
        else
            𝚽x12[iE] += Cx[ix] * ωx[ix+1,1] * 𝚽n[i] * sign(μ)
        end
    else
        𝚽x12[iE] += 𝚽n[i]
    end
end
if 𝒪x == 2 && 𝒪E == 2
    𝚽E12[2] += TE * 𝚽n[2]
    𝚽x12[2] += Tx * 𝚽n[3]
end
𝚽E12 .= 𝚽E12/ΔE

# Returning solutions
return 𝚽n, 𝚽x12, 𝚽E12

end 
end

#=
#----
# Full coupling
#----
if isFC

# Initialization
S = zeros(𝒪x*𝒪E,𝒪x*𝒪E)
Q = zeros(𝒪x*𝒪E)
𝚽n = Q
TE = 0.0; Tx = 0.0

# Galerkin energy scheme weights
Λ = 1 #β⁻/β⁺
if abs(ωE[1,1]) > 0
    ωE[1,1] = ωE[1,1]*Λ
    ωE[2:𝒪E+1,1] = (ωE[2:𝒪E+1,1].-1).*Λ.+1
end

hx = abs(μ/(Σt*Δx))
hE = abs(-(β⁺+β⁻)/(2*Σt))
sx = sign(μ)
sE = -1

# Adaptive loop
isAdapt = isAdaptE && isAdaptx
isFixed = false
while ~isFixed

# Matrix of Legendre moment coefficients of the flux
@inbounds for ix in range(1,𝒪x), jx in range(1,𝒪x), iE in range(1,𝒪E), jE in range(1,𝒪E)
    i = 𝒪E*(ix-1)+iE
    j = 𝒪E*(jx-1)+jE
    # Diagonal terms
    if i == j
        S[i,j] = 1 + CE[iE]^2 * hE * ωE[jE+1,jx] + Cx[ix]^2 * hx * ωx[jx+1,jE]
    # Upper diagonal terms
    elseif i < j
    # Energy terms - E
    if ix == jx
    if mod(iE+jE,2) == 1
        S[i,j] = CE[iE] * CE[jE] * sE * hE * ωE[jE+1,jx]
    else
        S[i,j] = CE[iE] * CE[jE] * hE * ωE[jE+1,jx]
    end
    # Space terms - x
    elseif  iE == jE 
    if mod(ix+jx,2) == 1
        S[i,j] = Cx[ix] * Cx[jx] * ωx[jx+1,jE] * sx * hx
    else
        S[i,j] = Cx[ix] * Cx[jx] * ωx[jx+1,jE] * hx
    end
    end
    # Under diagonal terms
    else
    # Energy terms - E
    if ix == jx
    if mod(iE+jE,2) == 1
        S[i,j] = CE[iE] * CE[jE] * (ωE[jE+1,jx]-2) * sE * hE
    else
        S[i,j] = CE[iE] * CE[jE] * ωE[jE+1,jx] * hE
    end
    # Space terms - x
    elseif  iE == jE 
    if mod(ix+jx,2) == 1
        S[i,j] = Cx[ix] * Cx[jx] * (ωx[jx+1,jE]-2) * sx * hx
    else
        S[i,j] = Cx[ix] * Cx[jx] * ωx[jx+1,jE] * hx
    end
    end
    end

    # T-factors
    if 𝒪x == 2 && 𝒪E == 2
        if i == 2 && j == 3
            S[i,j] += hx * sx * sE * Tx
        elseif i == 3 && j == 2
            S[i,j] += hE * sx * sE * TE
        elseif i == 4 && j == 2
            S[i,j] += sqrt(3) * sx * hE * TE
        elseif i == 4 && j == 3
            S[i,j] += sqrt(3) * sE * hx * Tx
        end
    end

end

# Source vector
@inbounds for ix in range(1,𝒪x), iE in range(1,𝒪E)
    i = 𝒪E*(ix-1)+iE
    Q[i] = Qn[i]/Σt
    # Energy terms - E
    if mod(iE,2) == 1
        Q[i] += CE[iE] * (1-ωE[1,ix]) * 𝚽E12[ix] * hE
    else
        Q[i] += -CE[iE] * (1+ωE[1,ix]) * 𝚽E12[ix] * sE * hE
    end
    # Space terms - x
    if mod(ix,2) == 1
        Q[i] += Cx[ix] * (1-ωx[1,iE]) * 𝚽x12[iE] * hx
    else
        Q[i] += -Cx[ix] * (1+ωx[1,iE]) * 𝚽x12[iE] * sx * hx
    end
end

𝚽n = S\Q

# Adaptive correction of weighting parameters
if isAdapt
    isFixed, ω, T = adaptive_2D([𝒪E,𝒪x],[ωE,ωx],𝚽n,[𝚽E12,𝚽x12],[-1.0,sign(μ)],[Λ,1.0],[TE,Tx],Q)
    ωE = ω[1]; ωx = ω[2]; TE = T[1]; Tx = T[2];
else
    isFixed = true
end

end # End of adaptive loop

# Closure relation
@inbounds for ix in range(1,𝒪x), iE in range(1,𝒪E)
    if (iE == 1) 𝚽E12[ix] = ωE[1,ix] * 𝚽E12[ix] end
    if (ix == 1) 𝚽x12[iE] = ωx[1,iE] * 𝚽x12[iE] end
end
@inbounds for ix in range(1,𝒪x), iE in range(1,𝒪E)
    i = 𝒪E*(ix-1)+iE
    # Energy terms - E
    if mod(iE,2) == 1
        𝚽E12[ix] += CE[iE] * ωE[iE+1,ix] * 𝚽n[i]
    else
        𝚽E12[ix] += CE[iE] * ωE[iE+1,ix] * 𝚽n[i] * sE
    end  
    # Space terms - x
    if mod(ix,2) == 1
        𝚽x12[iE] += Cx[ix] * ωx[ix+1,iE] * 𝚽n[i]
    else
        𝚽x12[iE] += Cx[ix] * ωx[ix+1,iE] * 𝚽n[i] * sx
    end
end
if 𝒪x == 2 && 𝒪E == 2
    𝚽E12[2] += sx * sE * TE * 𝚽n[2]
    𝚽x12[2] += sx * sE * Tx * 𝚽n[3]
end
𝚽E12 .= 𝚽E12/ΔE

# Returning solutions
return 𝚽n, 𝚽x12, 𝚽E12
=#