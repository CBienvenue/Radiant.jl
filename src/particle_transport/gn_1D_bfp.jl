function gn_1D_BFP(sx::Int64,Σt::Float64,S⁻::Float64,S⁺::Float64,S::Vector{Float64},Δx::Float64,Qn::Array{Float64},𝚽x12::Array{Float64},𝚽E12::Array{Float64},Nmx::Int64,NmE::Int64,Np::Int64,C::Vector{Float64},ωx::Vector{Float64},ωE::Vector{Float64},𝒩x::Matrix{Float64},𝒲::Array{Float64},isFC::Bool)

# Initialization
Nm = Nmx*NmE*Np
𝒮 = zeros(Nm,Nm)
Q = zeros(Nm)
𝚽 = Q
𝚽n = copy(Qn)
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

# Solve the equation system
𝚽 = 𝒮\Q

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

# Returning solutions
return 𝚽n, 𝚽x12, 𝚽E12

end