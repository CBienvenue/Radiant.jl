function gn_1D_BFP!(𝚽n::AbstractArray{Float64,2},𝚽x12::AbstractArray{Float64,2},𝚽E12::AbstractArray{Float64,2},sx::Int64,Σt::Float64,S⁻::Float64,S⁺::Float64,S::AbstractVector{Float64},Δx::Float64,Qn::AbstractArray{Float64,2},𝒮::Matrix{Float64},Q::Vector{Float64},𝚽::Vector{Float64},Nmx::Int64,NmE::Int64,Np::Int64,C::Vector{Float64},ωx::Vector{Float64},ωE::Vector{Float64},𝒩x::AbstractMatrix{Float64},𝒲::Array{Float64},isFC::Bool)

# Initialization
Nm = isFC ? Nmx*NmE*Np : (Nmx+NmE-1)*Np
@inbounds for j in 1:Nm
    Q[j] = 0.0
    for i in 1:Nm
        𝒮[i,j] = 0.0
    end
end
g(n,sx) = if sx > 0 return 1 else return -(-1)^(n-1) end
function index_Ex(iE,ix)
    if isFC
        return NmE*(ix-1)+iE
    else
        i = 1 + (ix-1) + (iE-1)
        if ix > 1 i += NmE-1 end
        return i
    end
end
function index_Exp(iE,ix,ip)
    if isFC
        i = Nmx*NmE*(ip-1) + NmE*(ix-1) + iE
    else
        i = 1 + (ix-1) + (iE-1)
        if ix > 1 i += NmE-1 end
        i += (ip-1)*(NmE+Nmx-1)
    end
end

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,Nmx), jx in range(1,Nmx), iE in range(1,NmE), jE in range(1,NmE)
    factor = C[ix]/Δx * C[jx] * (g(ix,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix)*(1-(-1)^(ix-jx)))
    if (~isFC) && (count(>(1),(ix,iE)) ≥ 2 || count(>(1),(jx,jE)) ≥ 2) continue end
    for ip in range(1,Np), jp in range(1,Np)
        i = index_Exp(iE,ix,ip)
        j = index_Exp(jE,jx,jp)

        # Collision term
        if (i == j) 𝒮[i,j] += Σt end

        # Streaming term - x
        if iE == jE
            𝒮[i,j] += factor * 𝒩x[ip,jp]
        end

        # CSD term
        if ip == jp && ix == jx
            for kE in range(1,iE-1), wE in range(1,NmE)
                𝒮[i,j] += C[iE] * C[jE] * C[kE] * C[wE] * (1-(-1)^(iE-kE)) * S[wE] * 𝒲[jE,kE,wE]
            end
            𝒮[i,j] += C[iE] * S⁺ * (-1)^(iE-1) * C[jE] * (-1)^(jE-1) * ωE[jE+1]
        end
    end
end

# Source vector
for ix in range(1,Nmx), iE in range(1,NmE)
    if (~isFC) && (count(>(1),(ix,iE)) ≥ 2) continue end
    factor = -C[ix]/Δx * (g(ix,sx)*ωx[1]+g(ix,-sx))
    for ip in range(1,Np)
        i = index_Exp(iE,ix,ip)

        # Volume sources
        Q[i] += Qn[ip,index_Ex(iE,ix)]

        # Incoming boundary sources
        for jp in range(1,Np)
            Q[i] += factor * 𝒩x[ip,jp] * 𝚽x12[jp,iE]
        end

        # CSD incoming sources
        Q[i] += C[iE] * ((-1)^iE*S⁺*ωE[1] + S⁻) * 𝚽E12[ip,ix]
    end
end

# Solve the equation system (in place)
F = lu!(𝒮)
ldiv!(𝚽, F, Q)

# Closure relations
for ip in 1:Np
    for iE in 1:NmE
        𝚽x12[ip,iE] = ωx[1] * 𝚽x12[ip,iE]
        for ix in 1:Nmx
            if (~isFC) && (count(>(1),(ix,iE)) ≥ 2) continue end
            i = index_Exp(iE,ix,ip)
            𝚽x12[ip,iE] += C[ix] * sx^(ix-1) * ωx[ix+1] * 𝚽[i]
        end
    end
    for ix in 1:Nmx
        𝚽E12[ip,ix] = ωE[1] * 𝚽E12[ip,ix]
        for iE in 1:NmE
            if (~isFC) && (count(>(1),(ix,iE)) ≥ 2) continue end
            i = index_Exp(iE,ix,ip)
            𝚽E12[ip,ix] += C[iE] * (-1)^(iE-1) * ωE[iE+1] * 𝚽[i]
        end
    end
    for ix in 1:Nmx, iE in 1:NmE
        if (~isFC) && (count(>(1),(ix,iE)) ≥ 2) continue end
        i = index_Exp(iE,ix,ip)
        𝚽n[ip,index_Ex(iE,ix)] = 𝚽[i]
    end
end

return nothing

end

"""
    gn_1D_BFP_matrix!(𝒮,sx,Σt,S⁺,S,Δx,Nmx,NmE,Np,C,ωx,ωE,𝒩x,𝒲,isFC)

Assemble the cell matrix `𝒮` of the 1D BFP kernel, for the optimized solver chain.

The body is the assembly loop of `gn_1D_BFP!` transcribed unchanged, so the matrix is
bit-identical to the one the reference kernel builds. See `gn_1D_BTE_matrix!`.
"""
function gn_1D_BFP_matrix!(𝒮::Matrix{Float64},sx::Int64,Σt::Float64,S⁺::Float64,S::AbstractVector{Float64},Δx::Float64,Nmx::Int64,NmE::Int64,Np::Int64,C::Vector{Float64},ωx::Vector{Float64},ωE::Vector{Float64},𝒩x::Matrix{Float64},𝒲::Array{Float64},isFC::Bool)

fill!(𝒮,0.0)
g(n,sx) = if sx > 0 return 1 else return -(-1)^(n-1) end
function index_Exp(iE,ix,ip)
    if isFC
        i = Nmx*NmE*(ip-1) + NmE*(ix-1) + iE
    else
        i = 1 + (ix-1) + (iE-1)
        if ix > 1 i += NmE-1 end
        i += (ip-1)*(NmE+Nmx-1)
    end
end

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,Nmx), jx in range(1,Nmx), iE in range(1,NmE), jE in range(1,NmE)
    factor = C[ix]/Δx * C[jx] * (g(ix,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix)*(1-(-1)^(ix-jx)))
    if (~isFC) && (count(>(1),(ix,iE)) ≥ 2 || count(>(1),(jx,jE)) ≥ 2) continue end
    for ip in range(1,Np), jp in range(1,Np)
        i = index_Exp(iE,ix,ip)
        j = index_Exp(jE,jx,jp)

        # Collision term
        if (i == j) 𝒮[i,j] += Σt end

        # Streaming term - x
        if iE == jE
            𝒮[i,j] += factor * 𝒩x[ip,jp]
        end

        # CSD term
        if ip == jp && ix == jx
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
    gn_1D_BFP_fast!(𝚽n,o𝚽n,𝚽x12,ox,𝚽E12,oE,Qn,oQ,mom,pat,ws,conf,ikx,m,do_E)

Optimized counterpart of `gn_1D_BFP!`: solve one voxel of the 1D BFP equation.

Identical in structure to `gn_2D_BFP_fast!` minus the y axis; `do_E` skips the energy
closure on the ordinary passes for the same reason (see `gn_3D_BFP_fast!`).
"""
@inline function gn_1D_BFP_fast!(𝚽n::Vector{Float64},o𝚽n::Int64,𝚽x12::Vector{Float64},ox::Int64,𝚽E12::Vector{Float64},oE::Int64,Qn::Vector{Float64},oQ::Int64,mom::GNFastMoments{NMOM},pat::GNFastPatch{NQ},ws::GNFastWorkspace,conf::Int64,ikx::Int64,m::Int64,do_E::Bool) where {NMOM,NQ}

    Nq   = NQ
    Nmom = NMOM
    Q    = ws.Q
    𝚽    = ws.𝚽
    𝒩x   = pat.𝒩x

    # Source vector
    @inbounds for k in 1:NMOM
        col = mom.col[k]
        bx = ox + Nq*(mom.cyz[k]-1)
        bE = oE + Nq*(mom.cxyz[k]-1); bQ = oQ + Nq*(col-1)
        fx = pat.qcx[k,ikx]; fE = pat.qcE[k,m]
        for ip in 1:Nq
            q = Qn[bQ+ip]
            for jp in 1:Nq; q += fx * 𝒩x[ip,jp] * 𝚽x12[bx+jp] end
            q += fE * 𝚽E12[bE+ip]
            Q[col + (ip-1)*Nmom] = q
        end
    end

    # Solve the equation system from the cached factorization
    gn_fast_solve!(𝚽,pat.LU,pat.ipiv,Q,pat.Nm,conf)

    # Closure relations
    if do_E gn_fast_closure!(𝚽E12,oE,pat.clE,𝚽,Val(NMOM),Val(NQ)) end
    gn_fast_closure!(𝚽x12,ox,pat.clx,𝚽,Val(NMOM),Val(NQ))
    @inbounds for c in 1:NMOM
        b = o𝚽n + Nq*(c-1); r = c - Nmom
        for ip in 1:Nq
            𝚽n[b+ip] = 𝚽[r + ip*Nmom]
        end
    end

    return nothing
end
