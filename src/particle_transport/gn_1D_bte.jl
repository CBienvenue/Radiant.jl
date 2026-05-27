function gn_1D_BTE!(𝚽n::AbstractArray{Float64,2},𝚽x12::AbstractVector{Float64},sx::Int64,Σt::Float64,Δx::Float64,Qn::AbstractArray{Float64,2},𝒮::Matrix{Float64},Q::Vector{Float64},𝚽::Vector{Float64},Nmx::Int64,Np::Int64,C::Vector{Float64},ωx::Vector{Float64},𝒩x::AbstractMatrix{Float64})

# Initialization
Nm = Nmx*Np
@inbounds for j in 1:Nm
    Q[j] = 0.0
    for i in 1:Nm
        𝒮[i,j] = 0.0
    end
end
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

# Solve the equation system (in place)
F = lu!(𝒮)
ldiv!(𝚽, F, Q)

# Closure relations
for ip in range(1,Np)
    𝚽x12[ip] = ωx[1] * 𝚽x12[ip]
    for ix in range(1,Nmx)
        i = index_xp(ix,ip)
        𝚽x12[ip] += C[ix] * sx^(ix-1) * ωx[ix+1] * 𝚽[i]
        𝚽n[ip,ix] = 𝚽[i]
    end
end

return nothing

end
