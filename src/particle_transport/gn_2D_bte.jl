function gn_2D_BTE!(𝚽n::AbstractArray{Float64,2},𝚽x12::AbstractArray{Float64,2},𝚽y12::AbstractArray{Float64,2},sx::Int64,sy::Int64,Σt::Float64,Δx::Float64,Δy::Float64,Qn::AbstractArray{Float64,2},𝒮::Matrix{Float64},Q::Vector{Float64},𝚽::Vector{Float64},Nmx::Int64,Nmy::Int64,Np::Int64,C::Vector{Float64},ωx::Array{Float64},ωy::Array{Float64},𝒩x::AbstractMatrix{Float64},𝒩y::AbstractMatrix{Float64},isFC::Bool)

# Initialization
Nm = isFC ? Nmx*Nmy*Np : (Nmx+Nmy-1)*Np
@inbounds for j in 1:Nm
    Q[j] = 0.0
    for i in 1:Nm
        𝒮[i,j] = 0.0
    end
end
g(n,sx) = if sx > 0 return 1 else return -(-1)^(n-1) end
function index_xy(ix,iy)
    if isFC
        return Nmx*(iy-1)+ix
    else
        i = 1 + (iy-1) + (ix-1)
        if iy > 1 i += Nmx-1 end
        return i
    end
end
function index_xyp(ix,iy,ip)
    if isFC
        i = Nmx*Nmy*(ip-1) + Nmx*(iy-1) + ix
    else
        i = 1 + (iy-1) + (ix-1)
        if iy > 1 i += Nmx-1 end
        i += (ip-1)*(Nmx+Nmy-1)
    end
end

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,Nmx), jx in range(1,Nmx), iy in range(1,Nmy), jy in range(1,Nmy)
    if (~isFC) && (count(>(1),(ix,iy)) ≥ 2 || count(>(1),(jx,jy)) ≥ 2) continue end
    fx = C[ix]/Δx * C[jx] * (g(ix,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx)))
    fy = C[iy]/Δy * C[jy] * (g(iy,sy)*sy^(jy-1)*ωy[jy+1] - (jy ≤ iy-1)*(1-(-1)^(iy-jy)))
    for ip in range(1,Np), jp in range(1,Np)
        i = index_xyp(ix,iy,ip)
        j = index_xyp(jx,jy,jp)

        # Collision term
        if (i == j) 𝒮[i,j] += Σt end

        # Streaming term - x
        if (iy == jy)
            𝒮[i,j] += fx * 𝒩x[ip,jp]
        end

        # Streaming term - y
        if (ix == jx)
            𝒮[i,j] += fy * 𝒩y[ip,jp]
        end
    end
end

# Source vector
for ix in range(1,Nmx), iy in range(1,Nmy)
    if (~isFC) && (count(>(1),(ix,iy)) ≥ 2) continue end
    fx = -C[ix]/Δx * (g(ix,sx)*ωx[1] + g(ix,-sx))
    fy = -C[iy]/Δy * (g(iy,sy)*ωy[1] + g(iy,-sy))
    for ip in range(1,Np)
        i = index_xyp(ix,iy,ip)

        # Volume sources
        Q[i] += Qn[ip,index_xy(ix,iy)]

        # Incoming boundary sources - x
        for jp in range(1,Np)
            Q[i] += fx * 𝒩x[ip,jp] * 𝚽x12[jp,iy]
        end

        # Incoming boundary sources - y
        for jp in range(1,Np)
            Q[i] += fy * 𝒩y[ip,jp] * 𝚽y12[jp,ix]
        end
    end
end

# Solve the equation system (in place)
F = lu!(𝒮)
ldiv!(𝚽, F, Q)

# Closure relations
for ip in 1:Np
    for iy in 1:Nmy
        𝚽x12[ip,iy] = ωx[1] * 𝚽x12[ip,iy]
        for ix in 1:Nmx
            if (~isFC) && (count(>(1),(ix,iy)) ≥ 2) continue end
            i = index_xyp(ix,iy,ip)
            𝚽x12[ip,iy] += C[ix] * sx^(ix-1) * ωx[ix+1] * 𝚽[i]
        end
    end
    for ix in 1:Nmx
        𝚽y12[ip,ix] = ωy[1] * 𝚽y12[ip,ix]
        for iy in 1:Nmy
            if (~isFC) && (count(>(1),(ix,iy)) ≥ 2) continue end
            i = index_xyp(ix,iy,ip)
            𝚽y12[ip,ix] += C[iy] * sy^(iy-1) * ωy[iy+1] * 𝚽[i]
        end
    end
    for ix in 1:Nmx, iy in 1:Nmy
        if (~isFC) && (count(>(1),(ix,iy)) ≥ 2) continue end
        i = index_xyp(ix,iy,ip)
        𝚽n[ip,index_xy(ix,iy)] = 𝚽[i]
    end
end

return nothing

end

"""
    gn_2D_BTE_matrix!(𝒮,sx,sy,Σt,Δx,Δy,Nmx,Nmy,Np,C,ωx,ωy,𝒩x,𝒩y,isFC)

Assemble the cell matrix `𝒮` of the 2D BTE kernel, for the optimized solver chain.

The body is the assembly loop of `gn_2D_BTE!` transcribed unchanged, so the matrix is
bit-identical to the one the reference kernel builds. It is split out because the matrix
does not depend on the voxel beyond its widths `(Δx,Δy)`: `gn_fast_context` calls this once
per (material, mesh-width pair, angular patch) and factorizes the result, instead of
rebuilding and refactorizing it in every voxel of every sweep.

See `set_fast_path`.
"""
function gn_2D_BTE_matrix!(𝒮::Matrix{Float64},sx::Int64,sy::Int64,Σt::Float64,Δx::Float64,Δy::Float64,Nmx::Int64,Nmy::Int64,Np::Int64,C::Vector{Float64},ωx::Vector{Float64},ωy::Vector{Float64},𝒩x::Matrix{Float64},𝒩y::Matrix{Float64},isFC::Bool)

fill!(𝒮,0.0)
g(n,sx) = if sx > 0 return 1 else return -(-1)^(n-1) end
function index_xyp(ix,iy,ip)
    if isFC
        i = Nmx*Nmy*(ip-1) + Nmx*(iy-1) + ix
    else
        i = 1 + (iy-1) + (ix-1)
        if iy > 1 i += Nmx-1 end
        i += (ip-1)*(Nmx+Nmy-1)
    end
end

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,Nmx), jx in range(1,Nmx), iy in range(1,Nmy), jy in range(1,Nmy)
    if (~isFC) && (count(>(1),(ix,iy)) ≥ 2 || count(>(1),(jx,jy)) ≥ 2) continue end
    fx = C[ix]/Δx * C[jx] * (g(ix,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx)))
    fy = C[iy]/Δy * C[jy] * (g(iy,sy)*sy^(jy-1)*ωy[jy+1] - (jy ≤ iy-1)*(1-(-1)^(iy-jy)))
    for ip in range(1,Np), jp in range(1,Np)
        i = index_xyp(ix,iy,ip)
        j = index_xyp(jx,jy,jp)

        # Collision term
        if (i == j) 𝒮[i,j] += Σt end

        # Streaming term - x
        if (iy == jy)
            𝒮[i,j] += fx * 𝒩x[ip,jp]
        end

        # Streaming term - y
        if (ix == jx)
            𝒮[i,j] += fy * 𝒩y[ip,jp]
        end
    end
end

return nothing

end

"""
    gn_2D_BTE_fast!(𝚽n,o𝚽n,𝚽x12,ox,𝚽y12,oy,Qn,oQ,mom,pat,ws,conf,ikx,iky)

Optimized counterpart of `gn_2D_BTE!`: solve one voxel of the 2D BTE.

Identical in structure to `gn_3D_BTE_fast!` minus the z axis; see it for why the arrays come
as flat vectors plus a per-voxel offset, and why the cell matrix is not rebuilt. The moment
tables are the 3D ones evaluated at `Nmz = NmE = 1`, which reproduce the reference 2D index
maps exactly.
"""
@inline function gn_2D_BTE_fast!(𝚽n::Vector{Float64},o𝚽n::Int64,𝚽x12::Vector{Float64},ox::Int64,𝚽y12::Vector{Float64},oy::Int64,Qn::Vector{Float64},oQ::Int64,mom::GNFastMoments{NMOM},pat::GNFastPatch{NQ},ws::GNFastWorkspace,conf::Int64,ikx::Int64,iky::Int64) where {NMOM,NQ}

    Nq   = NQ
    Nmom = NMOM
    Q    = ws.Q
    𝚽    = ws.𝚽
    𝒩x   = pat.𝒩x; 𝒩y = pat.𝒩y

    # Source vector
    @inbounds for k in 1:NMOM
        col = mom.col[k]
        bx = ox + Nq*(mom.cyz[k]-1); by = oy + Nq*(mom.cxz[k]-1)
        bQ = oQ + Nq*(col-1)
        fx = pat.qcx[k,ikx]; fy = pat.qcy[k,iky]
        for ip in 1:Nq
            q = Qn[bQ+ip]
            for jp in 1:Nq; q += fx * 𝒩x[ip,jp] * 𝚽x12[bx+jp] end
            for jp in 1:Nq; q += fy * 𝒩y[ip,jp] * 𝚽y12[by+jp] end
            Q[col + (ip-1)*Nmom] = q
        end
    end

    # Solve the equation system from the cached factorization
    gn_fast_solve!(𝚽,pat.LU,pat.ipiv,Q,pat.Nm,conf)

    # Closure relations
    gn_fast_closure!(𝚽x12,ox,pat.clx,𝚽,Val(NMOM),Val(NQ))
    gn_fast_closure!(𝚽y12,oy,pat.cly,𝚽,Val(NMOM),Val(NQ))
    @inbounds for c in 1:NMOM
        b = o𝚽n + Nq*(c-1); r = c - Nmom
        for ip in 1:Nq
            𝚽n[b+ip] = 𝚽[r + ip*Nmom]
        end
    end

    return nothing
end
