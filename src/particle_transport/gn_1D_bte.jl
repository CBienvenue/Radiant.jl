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

"""
    gn_1D_BTE_matrix!(𝒮,sx,Σt,Δx,Nmx,Np,C,ωx,𝒩x)

Assemble the cell matrix `𝒮` of the 1D BTE kernel, for the optimized solver chain.

The body is the assembly loop of `gn_1D_BTE!` transcribed unchanged — including its
`(jx ≤ ix)` guard, which differs in form from the `(jx ≤ ix-1)` of 2D/3D but not in value
(the `(1-(-1)^(ix-jx))` factor vanishes at `jx == ix`). It is split out because the matrix
does not depend on the voxel beyond its width `Δx`: `gn_fast_context` calls this once per
(material, mesh width, angular patch) and factorizes the result.

See `set_fast_path`.
"""
function gn_1D_BTE_matrix!(𝒮::Matrix{Float64},sx::Int64,Σt::Float64,Δx::Float64,Nmx::Int64,Np::Int64,C::Vector{Float64},ωx::Vector{Float64},𝒩x::Matrix{Float64})

fill!(𝒮,0.0)
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

return nothing

end

"""
    gn_1D_BTE_fast!(𝚽n,o𝚽n,𝚽x12,ox,Qn,oQ,mom,pat,ws,conf,ikx)

Optimized counterpart of `gn_1D_BTE!`: solve one voxel of the 1D BTE.

Identical in structure to `gn_2D_BTE_fast!` minus the y axis. The moment tables are the 3D
ones evaluated at `Nmy = Nmz = NmE = 1`, which reproduce the reference 1D index maps
exactly; the x face then carries a single moment.
"""
@inline function gn_1D_BTE_fast!(𝚽n::Vector{Float64},o𝚽n::Int64,𝚽x12::Vector{Float64},ox::Int64,Qn::Vector{Float64},oQ::Int64,mom::GNFastMoments{NMOM},pat::GNFastPatch{NQ},ws::GNFastWorkspace,conf::Int64,ikx::Int64) where {NMOM,NQ}

    Nq   = NQ
    Nmom = NMOM
    Q    = ws.Q
    𝚽    = ws.𝚽
    𝒩x   = pat.𝒩x

    # Source vector
    @inbounds for k in 1:NMOM
        col = mom.col[k]
        bx = ox + Nq*(mom.cyz[k]-1); bQ = oQ + Nq*(col-1)
        fx = pat.qcx[k,ikx]
        for ip in 1:Nq
            q = Qn[bQ+ip]
            for jp in 1:Nq; q += fx * 𝒩x[ip,jp] * 𝚽x12[bx+jp] end
            Q[col + (ip-1)*Nmom] = q
        end
    end

    # Solve the equation system from the cached factorization
    gn_fast_solve!(𝚽,pat.LU,pat.ipiv,Q,pat.Nm,conf)

    # Closure relations
    gn_fast_closure!(𝚽x12,ox,pat.clx,𝚽,Val(NMOM),Val(NQ))
    @inbounds for c in 1:NMOM
        b = o𝚽n + Nq*(c-1); r = c - Nmom
        for ip in 1:Nq
            𝚽n[b+ip] = 𝚽[r + ip*Nmom]
        end
    end

    return nothing
end
