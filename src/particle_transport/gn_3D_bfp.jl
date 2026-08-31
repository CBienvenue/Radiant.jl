function gn_3D_BFP!(𝚽n::AbstractArray{Float64,2},𝚽x12::AbstractArray{Float64,2},𝚽y12::AbstractArray{Float64,2},𝚽z12::AbstractArray{Float64,2},𝚽E12::AbstractArray{Float64,2},sx::Int64,sy::Int64,sz::Int64,Σt::Float64,S⁻::Float64,S⁺::Float64,S::AbstractVector{Float64},Δx::Float64,Δy::Float64,Δz::Float64,Qn::AbstractArray{Float64,2},𝒮::Matrix{Float64},Q::Vector{Float64},𝚽::Vector{Float64},Nmx::Int64,Nmy::Int64,Nmz::Int64,NmE::Int64,Np::Int64,C::Vector{Float64},ωx::Array{Float64},ωy::Array{Float64},ωz::Array{Float64},ωE::Array{Float64},𝒩x::AbstractMatrix{Float64},𝒩y::AbstractMatrix{Float64},𝒩z::AbstractMatrix{Float64},𝒲::Array{Float64},isFC::Bool)

# Initialization
Nm = isFC ? Nmx*Nmy*Nmz*NmE*Np : (Nmx+Nmy+Nmz+NmE-3)*Np
@inbounds for j in 1:Nm
    Q[j] = 0.0
    for i in 1:Nm
        𝒮[i,j] = 0.0
    end
end
g(n,sx) = if sx > 0 return 1 else return -(-1)^(n-1) end
function index_Exy(iE,ix,iy)
    if isFC
        i = Nmx*NmE*(iy-1) + NmE*(ix-1) + iE
    else
        i = 1 + (iy-1) + (ix-1) + (iE-1)
        if ix > 1 i += NmE-1 end
        if iy > 1 i += NmE-1 + Nmx-1 end
        return i
    end
end
function index_Exz(iE,ix,iz)
    if isFC
        i = Nmx*NmE*(iz-1) + NmE*(ix-1) + iE
    else
        i = 1 + (iz-1) + (ix-1) + (iE-1)
        if ix > 1 i += NmE-1 end
        if iz > 1 i += NmE-1 + Nmx-1 end
        return i
    end
end
function index_Eyz(iE,iy,iz)
    if isFC
        i = Nmy*NmE*(iz-1) + NmE*(iy-1) + iE
    else
        i = 1 + (iz-1) + (iy-1) + (iE-1)
        if iy > 1 i += NmE-1 end
        if iz > 1 i += NmE-1 + Nmy-1 end
        return i
    end
end
function index_xyz(ix,iy,iz)
    if isFC
        i = Nmy*Nmx*(iz-1) + Nmx*(iy-1) + ix
    else
        i = 1 + (iz-1) + (iy-1) + (ix-1)
        if iy > 1 i += Nmx-1 end
        if iz > 1 i += Nmx-1 + Nmy-1 end
        return i
    end
end
function index_Exyz(iE,ix,iy,iz)
    if isFC
        i = Nmx*Nmy*NmE*(iz-1) + Nmx*NmE*(iy-1) + NmE*(ix-1) + iE
    else
        i = 1 + (iz-1) + (iy-1) + (ix-1) + (iE-1)
        if ix > 1 i += NmE-1 end
        if iy > 1 i += NmE-1 + Nmx-1 end
        if iz > 1 i += NmE-1 + Nmx-1 + Nmy-1 end
        return i
    end
end
function index_Exyzp(iE,ix,iy,iz,jp)
    if isFC
        i = Nmx*Nmy*Nmz*NmE*(jp-1) + Nmx*Nmy*NmE*(iz-1) + Nmx*NmE*(iy-1) + NmE*(ix-1) + iE
    else
        i = 1 + (iz-1) + (iy-1) + (ix-1) + (iE-1)
        if ix > 1 i += NmE-1 end
        if iy > 1 i += NmE-1 + Nmx-1 end
        if iz > 1 i += NmE-1 + Nmx-1 + Nmy-1 end
        i += (jp-1)*(NmE+Nmx+Nmy+Nmz-3)
        return i
    end
end

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,Nmx), jx in range(1,Nmx), iy in range(1,Nmy), jy in range(1,Nmy), iz in range(1,Nmz), jz in range(1,Nmz), iE in range(1,NmE), jE in range(1,NmE)
    if (~isFC) && (count(>(1),(iE,ix,iy,iz)) ≥ 2 || count(>(1),(jE,jx,jy,jz)) ≥ 2) continue end
    fx = C[ix]/Δx * C[jx] * (g(ix,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx)))
    fy = C[iy]/Δy * C[jy] * (g(iy,sy)*sy^(jy-1)*ωy[jy+1] - (jy ≤ iy-1)*(1-(-1)^(iy-jy)))
    fz = C[iz]/Δz * C[jz] * (g(iz,sz)*sz^(jz-1)*ωz[jz+1] - (jz ≤ iz-1)*(1-(-1)^(iz-jz)))
    for ip in range(1,Np), jp in range(1,Np)
        i = index_Exyzp(iE,ix,iy,iz,ip)
        j = index_Exyzp(jE,jx,jy,jz,jp)

        # Collision term
        if (i == j) 𝒮[i,j] += Σt end

        # Streaming term - x
        if (iE == jE) && (iy == jy) && (iz == jz)
            𝒮[i,j] += fx * 𝒩x[ip,jp]
        end

        # Streaming term - y
        if (iE == jE) && (ix == jx) && (iz == jz)
            𝒮[i,j] += fy * 𝒩y[ip,jp]
        end

        # Streaming term - z
        if (iE == jE) && (ix == jx) && (iy == jy)
            𝒮[i,j] += fz * 𝒩z[ip,jp]
        end

        # CSD term
        if ip == jp && ix == jx && iy == jy && iz == jz
            for kE in range(1,iE-1), wE in range(1,NmE)
                𝒮[i,j] += C[iE] * C[jE] * C[kE] * C[wE] * (1-(-1)^(iE-kE)) * S[wE] * 𝒲[jE,kE,wE]
            end
            𝒮[i,j] += C[iE] * S⁺ * (-1)^(iE-1) * C[jE] * (-1)^(jE-1) * ωE[jE+1]
        end
    end
end

# Source vector
for ix in range(1,Nmx), iy in range(1,Nmy), iz in range(1,Nmz), iE in range(1,NmE)
    if (~isFC) && (count(>(1),(iE,ix,iy,iz)) ≥ 2) continue end
    fx = -C[ix]/Δx * (g(ix,sx)*ωx[1] + g(ix,-sx))
    fy = -C[iy]/Δy * (g(iy,sy)*ωy[1] + g(iy,-sy))
    fz = -C[iz]/Δz * (g(iz,sz)*ωz[1] + g(iz,-sz))
    for ip in range(1,Np)
        i = index_Exyzp(iE,ix,iy,iz,ip)

        # Volume sources
        Q[i] += Qn[ip,index_Exyz(iE,ix,iy,iz)]

        # Incoming boundary sources - x
        for jp in range(1,Np)
            Q[i] += fx * 𝒩x[ip,jp] * 𝚽x12[jp,index_Eyz(iE,iy,iz)]
        end

        # Incoming boundary sources - y
        for jp in range(1,Np)
            Q[i] += fy * 𝒩y[ip,jp] * 𝚽y12[jp,index_Exz(iE,ix,iz)]
        end

        # Incoming boundary sources - z
        for jp in range(1,Np)
            Q[i] += fz * 𝒩z[ip,jp] * 𝚽z12[jp,index_Exy(iE,ix,iy)]
        end

        # CSD incoming sources
        Q[i] += C[iE] * ((-1)^iE*S⁺*ωE[1] + S⁻) * 𝚽E12[ip,index_xyz(ix,iy,iz)]
    end
end

# Solve the equation system (in place)
F = lu!(𝒮)
ldiv!(𝚽, F, Q)

# Closure relations
for ip in 1:Np
    for ix in 1:Nmx, iy in 1:Nmy, iz in 1:Nmz
        if (~isFC) && (count(>(1),(ix,iy,iz)) ≥ 2) continue end
        𝚽E12[ip,index_xyz(ix,iy,iz)] = ωE[1] * 𝚽E12[ip,index_xyz(ix,iy,iz)]
        for iE in 1:NmE
            if (~isFC) && (count(>(1),(iE,ix,iy,iz)) ≥ 2) continue end
            i = index_Exyzp(iE,ix,iy,iz,ip)
            𝚽E12[ip,index_xyz(ix,iy,iz)] += C[iE] * (-1)^(iE-1) * ωE[iE+1] * 𝚽[i]
        end
    end
    for iE in 1:NmE, iy in 1:Nmy, iz in 1:Nmz
        if (~isFC) && (count(>(1),(iE,iy,iz)) ≥ 2) continue end
        𝚽x12[ip,index_Eyz(iE,iy,iz)] = ωx[1] * 𝚽x12[ip,index_Eyz(iE,iy,iz)]
        for ix in 1:Nmx
            if (~isFC) && (count(>(1),(iE,ix,iy,iz)) ≥ 2) continue end
            i = index_Exyzp(iE,ix,iy,iz,ip)
            𝚽x12[ip,index_Eyz(iE,iy,iz)] += C[ix] * sx^(ix-1) * ωx[ix+1] * 𝚽[i]
        end
    end
    for iE in 1:NmE, ix in 1:Nmx, iz in 1:Nmz
        if (~isFC) && (count(>(1),(iE,ix,iz)) ≥ 2) continue end
        𝚽y12[ip,index_Exz(iE,ix,iz)] = ωy[1] * 𝚽y12[ip,index_Exz(iE,ix,iz)]
        for iy in 1:Nmy
            if (~isFC) && (count(>(1),(iE,ix,iy,iz)) ≥ 2) continue end
            i = index_Exyzp(iE,ix,iy,iz,ip)
            𝚽y12[ip,index_Exz(iE,ix,iz)] += C[iy] * sy^(iy-1) * ωy[iy+1] * 𝚽[i]
        end
    end
    for iE in 1:NmE, ix in 1:Nmx, iy in 1:Nmy
        if (~isFC) && (count(>(1),(iE,ix,iy)) ≥ 2) continue end
        𝚽z12[ip,index_Exy(iE,ix,iy)] = ωz[1] * 𝚽z12[ip,index_Exy(iE,ix,iy)]
        for iz in 1:Nmz
            if (~isFC) && (count(>(1),(iE,ix,iy,iz)) ≥ 2) continue end
            i = index_Exyzp(iE,ix,iy,iz,ip)
            𝚽z12[ip,index_Exy(iE,ix,iy)] += C[iz] * sz^(iz-1) * ωz[iz+1] * 𝚽[i]
        end
    end
    for iE in 1:NmE, ix in 1:Nmx, iy in 1:Nmy, iz in 1:Nmz
        if (~isFC) && (count(>(1),(iE,ix,iy,iz)) ≥ 2) continue end
        i = index_Exyzp(iE,ix,iy,iz,ip)
        𝚽n[ip,index_Exyz(iE,ix,iy,iz)] = 𝚽[i]
    end
end

return nothing

end

"""
    gn_3D_BFP_matrix!(𝒮,sx,sy,sz,Σt,S⁺,S,Δx,Δy,Δz,Nmx,Nmy,Nmz,NmE,Np,C,ωx,ωy,ωz,ωE,
    𝒩x,𝒩y,𝒩z,𝒲,isFC)

Assemble the cell matrix `𝒮` of the 3D BFP kernel, for the optimized solver chain.

The body is the assembly loop of `gn_3D_BFP!` transcribed unchanged, so the matrix is
bit-identical to the one the reference kernel builds. It is split out because the matrix
does not depend on the voxel beyond its widths `(Δx,Δy,Δz)`: `gn_fast_context` calls this
once per (material, mesh-width triple, angular patch, energy group) and factorizes the
result, instead of rebuilding and refactorizing it in every voxel of every sweep.

See `set_fast_path`.
"""
function gn_3D_BFP_matrix!(𝒮::Matrix{Float64},sx::Int64,sy::Int64,sz::Int64,Σt::Float64,S⁺::Float64,S::AbstractVector{Float64},Δx::Float64,Δy::Float64,Δz::Float64,Nmx::Int64,Nmy::Int64,Nmz::Int64,NmE::Int64,Np::Int64,C::Vector{Float64},ωx::Vector{Float64},ωy::Vector{Float64},ωz::Vector{Float64},ωE::Vector{Float64},𝒩x::Matrix{Float64},𝒩y::Matrix{Float64},𝒩z::Matrix{Float64},𝒲::Array{Float64},isFC::Bool)

fill!(𝒮,0.0)
g(n,sx) = if sx > 0 return 1 else return -(-1)^(n-1) end
function index_Exyzp(iE,ix,iy,iz,jp)
    if isFC
        i = Nmx*Nmy*Nmz*NmE*(jp-1) + Nmx*Nmy*NmE*(iz-1) + Nmx*NmE*(iy-1) + NmE*(ix-1) + iE
    else
        i = 1 + (iz-1) + (iy-1) + (ix-1) + (iE-1)
        if ix > 1 i += NmE-1 end
        if iy > 1 i += NmE-1 + Nmx-1 end
        if iz > 1 i += NmE-1 + Nmx-1 + Nmy-1 end
        i += (jp-1)*(NmE+Nmx+Nmy+Nmz-3)
        return i
    end
end

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,Nmx), jx in range(1,Nmx), iy in range(1,Nmy), jy in range(1,Nmy), iz in range(1,Nmz), jz in range(1,Nmz), iE in range(1,NmE), jE in range(1,NmE)
    if (~isFC) && (count(>(1),(iE,ix,iy,iz)) ≥ 2 || count(>(1),(jE,jx,jy,jz)) ≥ 2) continue end
    fx = C[ix]/Δx * C[jx] * (g(ix,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx)))
    fy = C[iy]/Δy * C[jy] * (g(iy,sy)*sy^(jy-1)*ωy[jy+1] - (jy ≤ iy-1)*(1-(-1)^(iy-jy)))
    fz = C[iz]/Δz * C[jz] * (g(iz,sz)*sz^(jz-1)*ωz[jz+1] - (jz ≤ iz-1)*(1-(-1)^(iz-jz)))
    for ip in range(1,Np), jp in range(1,Np)
        i = index_Exyzp(iE,ix,iy,iz,ip)
        j = index_Exyzp(jE,jx,jy,jz,jp)

        # Collision term
        if (i == j) 𝒮[i,j] += Σt end

        # Streaming term - x
        if (iE == jE) && (iy == jy) && (iz == jz)
            𝒮[i,j] += fx * 𝒩x[ip,jp]
        end

        # Streaming term - y
        if (iE == jE) && (ix == jx) && (iz == jz)
            𝒮[i,j] += fy * 𝒩y[ip,jp]
        end

        # Streaming term - z
        if (iE == jE) && (ix == jx) && (iy == jy)
            𝒮[i,j] += fz * 𝒩z[ip,jp]
        end

        # CSD term
        if ip == jp && ix == jx && iy == jy && iz == jz
            for kE in range(1,iE-1), wE in range(1,NmE)
                𝒮[i,j] += C[iE] * C[jE] * C[kE] * C[wE] * (1-(-1)^(iE-kE)) * S[wE] * 𝒲[jE,kE,wE]
            end
            𝒮[i,j] += C[iE] * S⁺ * (-1)^(iE-1) * C[jE] * (-1)^(jE-1) * ωE[jE+1]
        end
    end
end

return nothing

end

"""
    gn_3D_BFP_fast!(𝚽n,o𝚽n,𝚽x12,ox,𝚽y12,oy,𝚽z12,oz,𝚽E12,oE,Qn,oQ,mom,pat,ws,conf,
    ikx,iky,ikz,m)

Optimized counterpart of `gn_3D_BFP!`: solve one voxel of the 3D BFP equation.

Same arithmetic, reorganized around the precomputed `GNFastMoments`/`GNFastPatch` tables
(see `gn_fast_context.jl`). The index closures, the activity tests and the assembly
coefficients are gone — they are voxel-independent — and the cell matrix is not rebuilt:
`conf` selects its cached factorization, leaving only the right-hand side, two triangular
solves and the closure relations.

Every array arrives as a **flat vector plus the voxel's base offset**, not as a 2-D slice.
All these arrays carry the angular index first, so element `(ip, c)` of one of them is at
`off + ip + Nq*(c-1)`. The reference form took `AbstractArray{Float64,2}` slices, which made
the sweep build five `SubArray`s per voxel — 6,7 % of the whole run at the call site alone,
and ~16 % once the generic `getindex`/`size`/bounds machinery they drag in is counted.

The right-hand side and the closures accumulate in the reference kernel's order, so the
result is bit-identical to `gn_3D_BFP!` up to the triangular solves (`gn_fast_solve!`).

`ikx`, `iky`, `ikz` are the distinct-mesh-width indices of the voxel, `m` its material and
`conf` its cell-system configuration; all come from the `GNFastCells` map.
"""
@inline function gn_3D_BFP_fast!(𝚽n::Vector{Float64},o𝚽n::Int64,𝚽x12::Vector{Float64},ox::Int64,𝚽y12::Vector{Float64},oy::Int64,𝚽z12::Vector{Float64},oz::Int64,𝚽E12::Vector{Float64},oE::Int64,Qn::Vector{Float64},oQ::Int64,mom::GNFastMoments{NMOM},pat::GNFastPatch{NQ},ws::GNFastWorkspace,conf::Int64,ikx::Int64,iky::Int64,ikz::Int64,m::Int64,do_E::Bool) where {NMOM,NQ}

    Nq   = NQ
    Nmom = NMOM
    Q    = ws.Q
    𝚽    = ws.𝚽
    𝒩x   = pat.𝒩x; 𝒩y = pat.𝒩y; 𝒩z = pat.𝒩z

    # Source vector
    @inbounds for k in 1:NMOM
        col = mom.col[k]
        bx = ox + Nq*(mom.cyz[k]-1); by = oy + Nq*(mom.cxz[k]-1)
        bz = oz + Nq*(mom.cxy[k]-1); bE = oE + Nq*(mom.cxyz[k]-1)
        bQ = oQ + Nq*(col-1)
        fx = pat.qcx[k,ikx]; fy = pat.qcy[k,iky]; fz = pat.qcz[k,ikz]; fE = pat.qcE[k,m]
        for ip in 1:Nq
            q = Qn[bQ+ip]
            for jp in 1:Nq; q += fx * 𝒩x[ip,jp] * 𝚽x12[bx+jp] end
            for jp in 1:Nq; q += fy * 𝒩y[ip,jp] * 𝚽y12[by+jp] end
            for jp in 1:Nq; q += fz * 𝒩z[ip,jp] * 𝚽z12[bz+jp] end
            q += fE * 𝚽E12[bE+ip]
            Q[col + (ip-1)*Nmom] = q
        end
    end

    # Solve the equation system from the cached factorization
    gn_fast_solve!(𝚽,pat.LU,pat.ipiv,Q,pat.Nm,conf)

    # Closure relations. The energy closure overwrites the voxel's incoming energy flux with
    # the outgoing one, and that outgoing value only ever feeds 𝚽E12_temp — which is read
    # after convergence, not during the iteration. Skipping it on the ordinary passes
    # therefore changes nothing, and it leaves 𝚽E12 pristine so the caller need not
    # re-transform it at every pass (see gn_inner_pass_fast!).
    if do_E gn_fast_closure!(𝚽E12,oE,pat.clE,𝚽,Val(NMOM),Val(NQ)) end
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
