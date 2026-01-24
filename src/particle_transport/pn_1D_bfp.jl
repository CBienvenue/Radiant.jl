function pn_1D_BFP(sx::Int64,Σt::Float64,Δx::Float64,Qn::Array{Float64},𝚽x12::Array{Float64},S⁻::Float64,S⁺::Float64,S::Vector{Float64},𝚽E12::Array{Float64},𝒪E::Int64,𝒪x::Int64,Np::Int64,C::Vector{Float64},ωE::Array{Float64},ωx::Array{Float64},𝒲::Array{Float64},isFC::Bool,𝒩::Matrix{Float64})

# Initialization
𝒮 = zeros(𝒪x*𝒪E*Np,𝒪x*𝒪E*Np)
Q = zeros(𝒪x*𝒪E*Np)
𝚽 = Q
𝚽n = copy(Qn)
g(n,sx) = (1+sx)/2 - (1-sx)/2 * (-1)^n
jpm(jx,jE) = 𝒪E*(jx-1)+jE

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,𝒪x), jx in range(1,𝒪x), iE in range(1,𝒪E), jE in range(1,𝒪E)
    factor = C[ix]*sx/Δx * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1,jE,iE] - (jx ≤ ix-1)*(1-(-1)^(ix-jx)))
    for ip in range(1,Np), jp in range(1,Np)
        i = 𝒪x*𝒪E*(ip-1)+𝒪E*(ix-1)+iE
        j = 𝒪x*𝒪E*(jp-1)+𝒪E*(jx-1)+jE
        
        # Collision term
        if (i == j) 𝒮[i,j] += Σt end

        # Streaming term - x
        if iE == jE
            𝒮[i,j] += factor * 𝒩[ip,jp]
        end

        # CSD term
        if ip == jp
            if ix == jx
                for kE in range(1,iE-1), wE in range(1,𝒪E)
                    𝒮[i,j] += C[iE] * C[jE] * C[kE] * C[wE] * (1-(-1)^(iE-kE)) * S[wE] * 𝒲[jE,kE,wE]
                end
            end
            𝒮[i,j] += C[iE] * S⁺ * (-1)^(iE-1) * C[jE] * (-1)^(jE-1) * ωE[jE+1,jx,ix]
        end
    end
end

# Source vector
for jx in range(1,𝒪x), jE in range(1,𝒪E)
    factor = - C[jx]*sx/Δx * (g(jx-1,sx)*ωx[1,jE,jE]+g(jx-1,-sx))
    for jp in range(1,Np)
        j = 𝒪x*𝒪E*(jp-1)+𝒪E*(jx-1)+jE

        # Volume sources
        Q[j] += Qn[jp,jpm(jx,jE)]

        # Incoming boundary sources
        for ip in range(1,Np)
            Q[j] += factor * 𝒩[ip,jp] * 𝚽x12[ip,jE]
        end

        # CSD incoming sources
        Q[j] -= C[jE] * ((-1)^(jE-1)*S⁺*ωE[1,jx,jx] - S⁻) * 𝚽E12[jp,jx]
    end
end

# Solve the equation system
𝚽 = 𝒮\Q

# Closure relations
for jp in range(1,Np), jE in range(1,𝒪E), jx in range(1,𝒪x)
    j = 𝒪x*𝒪E*(jp-1)+𝒪E*(jx-1)+jE
    if (jx == 1) 𝚽x12[jp,jE] = ωx[1,jE,jE] * 𝚽x12[jp,jE] end
    if (jE == 1) 𝚽E12[jp,jx] = ωE[1,jx,jx] * 𝚽E12[jp,jx] end
    for iE in range(1,𝒪E)
        𝚽x12[jp,jE] += C[jx] * sx^(jx-1) * ωx[jx+1,jE,iE] * 𝚽[j]
    end
    for ix in range(1,𝒪x)
        𝚽E12[jp,jx] += C[jE] * (-1)^(jE-1) * ωE[jE+1,jx,ix] * 𝚽[j]
    end
    𝚽n[jp,jpm(jx,jE)] = 𝚽[j]
end

# Returning solutions
return 𝚽n, 𝚽x12, 𝚽E12

end
