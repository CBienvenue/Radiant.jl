function gn_3D_BTE!(𝚽n::AbstractArray{Float64,2},𝚽x12::AbstractArray{Float64,2},𝚽y12::AbstractArray{Float64,2},𝚽z12::AbstractArray{Float64,2},sx::Int64,sy::Int64,sz::Int64,Σt::Float64,Δx::Float64,Δy::Float64,Δz::Float64,Qn::AbstractArray{Float64,2},𝒮::Matrix{Float64},Q::Vector{Float64},𝚽::Vector{Float64},Nmx::Int64,Nmy::Int64,Nmz::Int64,Np::Int64,C::Vector{Float64},ωx::Array{Float64},ωy::Array{Float64},ωz::Array{Float64},𝒩x::AbstractMatrix{Float64},𝒩y::AbstractMatrix{Float64},𝒩z::AbstractMatrix{Float64},isFC::Bool)

# Initialization
Nm = isFC ? Nmx*Nmy*Nmz*Np : (Nmx+Nmy+Nmz-2)*Np
@inbounds for j in 1:Nm
    Q[j] = 0.0
    for i in 1:Nm
        𝒮[i,j] = 0.0
    end
end
g(n,sx) = if sx > 0 return 1 else return -(-1)^(n-1) end
function index_xy(ix,iy)
    if isFC
        return Nmx*(iy-1) + ix
    else
        i = 1 + (ix-1) + (iy-1)
        if iy > 1 i += Nmx-1 end
        return i
    end
end
function index_yz(iy,iz)
    if isFC
        return Nmy*(iz-1) + iy
    else
        i = 1 + (iy-1) + (iz-1)
        if iz > 1 i += Nmy-1 end
        return i
    end
end
function index_xz(ix,iz)
    if isFC
        return Nmx*(iz-1) + ix
    else
        i = 1 + (ix-1) + (iz-1)
        if iz > 1 i += Nmx-1 end
        return i
    end
end
function index_xyz(ix,iy,iz)
    if isFC
        return Nmy*Nmx*(iz-1) + Nmx*(iy-1) + ix
    else
        i = 1 + (ix-1) + (iy-1) + (iz-1)
        if iy > 1 i += Nmx-1 end
        if iz > 1 i += Nmx-1 + Nmy-1 end
        return i
    end
end
function index_xyzp(ix,iy,iz,jp)
    if isFC
        i = Nmx*Nmy*Nmz*(jp-1) + Nmx*Nmy*(iz-1) + Nmx*(iy-1) + ix
    else
        i = 1 + (ix-1) + (iy-1) + (iz-1)
        if iy > 1 i += Nmx-1 end
        if iz > 1 i += Nmx-1 + Nmy-1 end
        i += (jp-1)*(Nmx+Nmy+Nmz-2)
    end
end

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,Nmx), jx in range(1,Nmx), iy in range(1,Nmy), jy in range(1,Nmy), iz in range(1,Nmz), jz in range(1,Nmz)
    fx = C[ix]/Δx * C[jx] * (g(ix,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx)))
    fy = C[iy]/Δy * C[jy] * (g(iy,sy)*sy^(jy-1)*ωy[jy+1] - (jy ≤ iy-1)*(1-(-1)^(iy-jy)))
    fz = C[iz]/Δz * C[jz] * (g(iz,sz)*sz^(jz-1)*ωz[jz+1] - (jz ≤ iz-1)*(1-(-1)^(iz-jz)))
    for ip in range(1,Np), jp in range(1,Np)
        if (~isFC) && (count(>(1),(ix,iy,iz)) ≥ 2 || count(>(1),(jx,jy,jz)) ≥ 2) continue end
        i = index_xyzp(ix,iy,iz,ip)
        j = index_xyzp(jx,jy,jz,jp)

        # Collision term
        if (i == j) 𝒮[i,j] += Σt end

        # Streaming term - x
        if (iy == jy) && (iz == jz)
            𝒮[i,j] += fx * 𝒩x[ip,jp]
        end

        # Streaming term - y
        if (ix == jx) && (iz == jz)
            𝒮[i,j] += fy * 𝒩y[ip,jp]
        end

        # Streaming term - z
        if (ix == jx) && (iy == jy)
            𝒮[i,j] += fz * 𝒩z[ip,jp]
        end
    end
end

# Source vector
for ix in range(1,Nmx), iy in range(1,Nmy), iz in range(1,Nmz)
    fx = -C[ix]/Δx * (g(ix,sx)*ωx[1]+g(ix,-sx))
    fy = -C[iy]/Δy * (g(iy,sy)*ωy[1]+g(iy,-sy))
    fz = -C[iz]/Δz * (g(iz,sz)*ωz[1]+g(iz,-sz))
    for ip in range(1,Np)
        if (~isFC) && (count(>(1),(ix,iy,iz)) ≥ 2) continue end
        j = index_xyzp(ix,iy,iz,ip)

        # Volume sources
        Q[j] += Qn[ip,index_xyz(ix,iy,iz)]

        # Incoming boundary sources - x
        for jp in range(1,Np)
            Q[j] += fx * 𝒩x[ip,jp] * 𝚽x12[jp,index_yz(iy,iz)]
        end

        # Incoming boundary sources - y
        for jp in range(1,Np)
            Q[j] += fy * 𝒩y[ip,jp] * 𝚽y12[jp,index_xz(ix,iz)]
        end

        # Incoming boundary sources - z
        for jp in range(1,Np)
            Q[j] += fz * 𝒩z[ip,jp] * 𝚽z12[jp,index_xy(ix,iy)]
        end
    end
end

# Solve the equation system (in place: lu! mutates 𝒮; ldiv! writes solution into 𝚽)
F = lu!(𝒮)
ldiv!(𝚽, F, Q)

# Closure relations
for ip in 1:Np
    for iy in 1:Nmy, iz in 1:Nmz
        if (~isFC) && (count(>(1),(iy,iz)) ≥ 2) continue end
        𝚽x12[ip,index_yz(iy,iz)] = ωx[1] * 𝚽x12[ip,index_yz(iy,iz)]
        for ix in 1:Nmx
            if (~isFC) && (count(>(1),(ix,iy,iz)) ≥ 2) continue end
            j = index_xyzp(ix,iy,iz,ip)
            𝚽x12[ip,index_yz(iy,iz)] += C[ix] * sx^(ix-1) * ωx[ix+1] * 𝚽[j]
        end
    end
    for ix in 1:Nmx, iz in 1:Nmz
        if (~isFC) && (count(>(1),(ix,iz)) ≥ 2) continue end
        𝚽y12[ip,index_xz(ix,iz)] = ωy[1] * 𝚽y12[ip,index_xz(ix,iz)]
        for iy in 1:Nmy
            if (~isFC) && (count(>(1),(ix,iy,iz)) ≥ 2) continue end
            j = index_xyzp(ix,iy,iz,ip)
            𝚽y12[ip,index_xz(ix,iz)] += C[iy] * sy^(iy-1) * ωy[iy+1] * 𝚽[j]
        end
    end
    for ix in 1:Nmx, iy in 1:Nmy
        if (~isFC) && (count(>(1),(ix,iy)) ≥ 2) continue end
        𝚽z12[ip,index_xy(ix,iy)] = ωz[1] * 𝚽z12[ip,index_xy(ix,iy)]
        for iz in 1:Nmz
            if (~isFC) && (count(>(1),(ix,iy,iz)) ≥ 2) continue end
            j = index_xyzp(ix,iy,iz,ip)
            𝚽z12[ip,index_xy(ix,iy)] += C[iz] * sz^(iz-1) * ωz[iz+1] * 𝚽[j]
        end
    end
    for ix in 1:Nmx, iy in 1:Nmy, iz in 1:Nmz
        if (~isFC) && (count(>(1),(ix,iy,iz)) ≥ 2) continue end
        j = index_xyzp(ix,iy,iz,ip)
        𝚽n[ip,index_xyz(ix,iy,iz)] = 𝚽[j]
    end
end

return nothing

end

"""
    gn_3D_BTE_matrix!(𝒮,sx,sy,sz,Σt,Δx,Δy,Δz,Nmx,Nmy,Nmz,Np,C,ωx,ωy,ωz,𝒩x,𝒩y,𝒩z,isFC)

Assemble the cell matrix `𝒮` of the 3D BTE kernel, for the optimized solver chain.

The body is the assembly loop of `gn_3D_BTE!` transcribed unchanged, so the matrix is
bit-identical to the one the reference kernel builds. It is split out because the matrix
does not depend on the voxel beyond its widths `(Δx,Δy,Δz)`: `gn_fast_context` calls this
once per (material, mesh-width triple, angular patch) and factorizes the result, instead of
rebuilding and refactorizing it in every voxel of every sweep.

See `set_fast_path`.
"""
function gn_3D_BTE_matrix!(𝒮::Matrix{Float64},sx::Int64,sy::Int64,sz::Int64,Σt::Float64,Δx::Float64,Δy::Float64,Δz::Float64,Nmx::Int64,Nmy::Int64,Nmz::Int64,Np::Int64,C::Vector{Float64},ωx::Vector{Float64},ωy::Vector{Float64},ωz::Vector{Float64},𝒩x::Matrix{Float64},𝒩y::Matrix{Float64},𝒩z::Matrix{Float64},isFC::Bool)

fill!(𝒮,0.0)
g(n,sx) = if sx > 0 return 1 else return -(-1)^(n-1) end
function index_xyzp(ix,iy,iz,jp)
    if isFC
        i = Nmx*Nmy*Nmz*(jp-1) + Nmx*Nmy*(iz-1) + Nmx*(iy-1) + ix
    else
        i = 1 + (ix-1) + (iy-1) + (iz-1)
        if iy > 1 i += Nmx-1 end
        if iz > 1 i += Nmx-1 + Nmy-1 end
        i += (jp-1)*(Nmx+Nmy+Nmz-2)
    end
end

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,Nmx), jx in range(1,Nmx), iy in range(1,Nmy), jy in range(1,Nmy), iz in range(1,Nmz), jz in range(1,Nmz)
    fx = C[ix]/Δx * C[jx] * (g(ix,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx)))
    fy = C[iy]/Δy * C[jy] * (g(iy,sy)*sy^(jy-1)*ωy[jy+1] - (jy ≤ iy-1)*(1-(-1)^(iy-jy)))
    fz = C[iz]/Δz * C[jz] * (g(iz,sz)*sz^(jz-1)*ωz[jz+1] - (jz ≤ iz-1)*(1-(-1)^(iz-jz)))
    for ip in range(1,Np), jp in range(1,Np)
        if (~isFC) && (count(>(1),(ix,iy,iz)) ≥ 2 || count(>(1),(jx,jy,jz)) ≥ 2) continue end
        i = index_xyzp(ix,iy,iz,ip)
        j = index_xyzp(jx,jy,jz,jp)

        # Collision term
        if (i == j) 𝒮[i,j] += Σt end

        # Streaming term - x
        if (iy == jy) && (iz == jz)
            𝒮[i,j] += fx * 𝒩x[ip,jp]
        end

        # Streaming term - y
        if (ix == jx) && (iz == jz)
            𝒮[i,j] += fy * 𝒩y[ip,jp]
        end

        # Streaming term - z
        if (ix == jx) && (iy == jy)
            𝒮[i,j] += fz * 𝒩z[ip,jp]
        end
    end
end

return nothing

end

"""
    gn_3D_BTE_fast!(𝚽n,o𝚽n,𝚽x12,ox,𝚽y12,oy,𝚽z12,oz,Qn,oQ,mom,pat,ws,conf,ikx,iky,ikz)

Optimized counterpart of `gn_3D_BTE!`: solve one voxel of the 3D BTE.

Same arithmetic, reorganized around the precomputed `GNFastMoments`/`GNFastPatch` tables
(see `gn_fast_context.jl`). The index closures, the activity tests and the assembly
coefficients are gone — they are voxel-independent — and the cell matrix is not rebuilt:
`conf` selects its cached factorization, leaving only the right-hand side, two triangular
solves and the closure relations.

Every array arrives as a **flat vector plus the voxel's base offset**, not as a 2-D slice;
see `gn_3D_BFP_fast!` for why.

The right-hand side and the closures accumulate in the reference kernel's order, so the
result is bit-identical to `gn_3D_BTE!` up to the triangular solves (`gn_fast_solve!`).

`ikx`, `iky`, `ikz` are the distinct-mesh-width indices of the voxel, `conf` its cell-system
configuration; both come from the `GNFastCells` map.
"""
@inline function gn_3D_BTE_fast!(𝚽n::Vector{Float64},o𝚽n::Int64,𝚽x12::Vector{Float64},ox::Int64,𝚽y12::Vector{Float64},oy::Int64,𝚽z12::Vector{Float64},oz::Int64,Qn::Vector{Float64},oQ::Int64,mom::GNFastMoments{NMOM},pat::GNFastPatch{NQ},ws::GNFastWorkspace,conf::Int64,ikx::Int64,iky::Int64,ikz::Int64) where {NMOM,NQ}

    Nq   = NQ
    Nmom = NMOM
    Q    = ws.Q
    𝚽    = ws.𝚽
    𝒩x   = pat.𝒩x; 𝒩y = pat.𝒩y; 𝒩z = pat.𝒩z

    # Source vector
    @inbounds for k in 1:NMOM
        col = mom.col[k]
        bx = ox + Nq*(mom.cyz[k]-1); by = oy + Nq*(mom.cxz[k]-1)
        bz = oz + Nq*(mom.cxy[k]-1); bQ = oQ + Nq*(col-1)
        fx = pat.qcx[k,ikx]; fy = pat.qcy[k,iky]; fz = pat.qcz[k,ikz]
        for ip in 1:Nq
            q = Qn[bQ+ip]
            for jp in 1:Nq; q += fx * 𝒩x[ip,jp] * 𝚽x12[bx+jp] end
            for jp in 1:Nq; q += fy * 𝒩y[ip,jp] * 𝚽y12[by+jp] end
            for jp in 1:Nq; q += fz * 𝒩z[ip,jp] * 𝚽z12[bz+jp] end
            Q[col + (ip-1)*Nmom] = q
        end
    end

    # Solve the equation system from the cached factorization
    gn_fast_solve!(𝚽,pat.LU,pat.ipiv,Q,pat.Nm,conf)

    # Closure relations
    gn_fast_closure!(𝚽x12,ox,pat.clx,𝚽,Val(NMOM),Val(NQ))
    gn_fast_closure!(𝚽y12,oy,pat.cly,𝚽,Val(NMOM),Val(NQ))
    gn_fast_closure!(𝚽z12,oz,pat.clz,𝚽,Val(NMOM),Val(NQ))
    @inbounds for c in 1:NMOM
        b = o𝚽n + Nq*(c-1); r = c - Nmom
        for ip in 1:Nq
            𝚽n[b+ip] = 𝚽[r + ip*Nmom]
        end
    end

    return nothing
end
