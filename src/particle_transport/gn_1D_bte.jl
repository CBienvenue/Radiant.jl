function gn_1D_BTE(sx::Int64,Σt::Float64,Δx::Float64,Qn::Array{Float64},𝚽x12::Vector{Float64},Nmx::Int64,Np::Int64,C::Vector{Float64},ωx::Vector{Float64},𝒩x::Matrix{Float64})

# Initialization
Nm = Nmx*Np
𝒮 = zeros(Nm,Nm)
Q = zeros(Nm)
𝚽 = Q
𝚽n = copy(Qn)
g(n,sx) = if sx > 0 return 1 else return -(-1)^(n-1) end
index_xp(ix,ip) = Nmx*(ip-1)+ix

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,Nmx), jx in range(1,Nmx)
    factor = C[ix]/Δx * C[jx] * (g(ix,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix)*(1-(-1)^(ix-jx)))
    for ip in range(1,Np), jp in range(1,Np)
        i = index_xp(ix,ip)
        j = index_xp(jx,jp)

        # Collision term
        if (i == j) 𝒮[i,j] += Σt end

        # Streaming term
        𝒮[i,j] += factor * 𝒩x[ip,jp]
    end
end

# Source vector
for ix in range(1,Nmx)
    factor = -C[ix]/Δx * (g(ix,sx)*ωx[1]+g(ix,-sx))
    for ip in range(1,Np)
        i = index_xp(ix,ip)

        # Volume sources
        Q[i] += Qn[ip,ix]

        # Incoming boundary sources
        for jp in range(1,Np)
            Q[i] += factor * 𝒩x[ip,jp] * 𝚽x12[jp]
        end
    end
end

# Solve the equation system
𝚽 = 𝒮\Q

# Closure relations
for ip in range(1,Np)
    𝚽x12[ip] = ωx[1] * 𝚽x12[ip]
    for ix in range(1,Nmx)
        i = index_xp(ix,ip)
        𝚽x12[ip] += C[ix] * sx^(ix-1) * ωx[ix+1] * 𝚽[i]
        𝚽n[ip,ix] = 𝚽[i]
    end
end

# Returning solutions
return 𝚽n, 𝚽x12

end