function pn_1D_BTE(sx::Int64,Σt::Float64,Δx::Float64,Qn::Array{Float64},𝚽x12::Vector{Float64},𝒪x::Int64,Np::Int64,C::Vector{Float64},ωx::Vector{Float64},𝒩::Matrix{Float64})

# Initialization
𝒮 = zeros(𝒪x*Np,𝒪x*Np)
Q = zeros(𝒪x*Np)
𝚽 = Q
𝚽n = copy(Qn)
g(n,sx) = (1+sx)/2 - (1-sx)/2 * (-1)^n
jxp(jx,jy) = 𝒪x*(jy-1)+jx

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,𝒪x), jx in range(1,𝒪x)
    factor = C[ix]*sx/Δx * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx)))
    for ip in range(1,Np), jp in range(1,Np)
        i = jxp(ix,ip)
        j = jxp(jx,jp)

        # Collision term
        if (ip == jp) && (ix == jx) 𝒮[i,j] += Σt end

        # Streaming term
        𝒮[i,j] += factor * 𝒩[ip,jp]
    end
end

# Source vector
for jx in range(1,𝒪x)
    factor = -C[jx]*sx/Δx * (g(jx-1,sx)*ωx[1]+g(jx-1,-sx))
    for jp in range(1,Np)
        j = jxp(jx,jp)

        # Volume sources
        Q[j] += Qn[jp,jx]

        # Incoming boundary sources
        for ip in range(1,Np)
            Q[j] += factor * 𝒩[ip,jp] * 𝚽x12[ip]
        end
    end
end

# Solve the equation system
𝚽 = 𝒮\Q

# Closure relations
for jp in range(1,Np)
    𝚽x12[jp] = ωx[1] * 𝚽x12[jp]
    for jx in range(1,𝒪x)
        j = jxp(jx,jp)
        𝚽x12[jp] += C[jx] * sx^(jx-1) * ωx[jx+1] * 𝚽[j]
        𝚽n[jp,jx] = 𝚽[j]
    end
end

# Returning solutions
return 𝚽n, 𝚽x12

end