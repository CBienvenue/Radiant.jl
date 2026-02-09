function gn_2D_BFP(sx::Int64,sy::Int64,Σt::Float64,S⁻::Float64,S⁺::Float64,S::Vector{Float64},Δx::Float64,Δy::Float64,Qn::Array{Float64},𝚽x12::Array{Float64},𝚽y12::Array{Float64},𝚽E12::Array{Float64},Nmx::Int64,Nmy::Int64,NmE::Int64,Np::Int64,C::Vector{Float64},ωx::Array{Float64},ωy::Array{Float64},ωE::Array{Float64},𝒩x::Matrix{Float64},𝒩y::Matrix{Float64},𝒲::Array{Float64},isFC::Bool)

# Initialization
if isFC Nm = Nmx*Nmy*NmE*Np else Nm = (Nmx+Nmy+NmE-2)*Np end
𝒮 = zeros(Nm,Nm)
Q = zeros(Nm)
𝚽 = Q
𝚽n = copy(Qn)
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
function index_Ex(iE,ix)
    if isFC 
        return NmE*(ix-1)+iE 
    else  
        i = 1 + (ix-1) + (iE-1)
        if ix > 1 i += NmE-1 end
        return i
    end
end
function index_Ey(iE,iy)
    if isFC 
        return NmE*(iy-1)+iE 
    else  
        i = 1 + (iy-1) + (iE-1)
        if iy > 1 i += NmE-1 end
        return i
    end
end
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
function index_Exyp(iE,ix,iy,jp)
    if isFC
        i = NmE*Nmx*Nmy*(jp-1) + NmE*Nmx*(iy-1) + NmE*(ix-1) + iE
    else
        i = 1 + (iE-1) + (ix-1) + (iy-1)
        if ix > 1 i += NmE-1 end
        if iy > 1 i += NmE-1 + Nmx-1 end
        i += (jp-1)*(NmE+Nmx+Nmy-2)
    end
end

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,Nmx), jx in range(1,Nmx), iy in range(1,Nmy), jy in range(1,Nmy), iE in range(1,NmE), jE in range(1,NmE)
    if (~isFC) && (count(>(1),(iE,ix,iy)) ≥ 2 || count(>(1),(jE,jx,jy)) ≥ 2) continue end
    fx = C[ix]/Δx * C[jx] * (g(ix,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx)))
    fy = C[iy]/Δy * C[jy] * (g(iy,sy)*sy^(jy-1)*ωy[jy+1] - (jy ≤ iy-1)*(1-(-1)^(iy-jy)))
    for ip in range(1,Np), jp in range(1,Np)
        i = index_Exyp(iE,ix,iy,ip)
        j = index_Exyp(jE,jx,jy,jp)

        # Collision term
        if (i == j) 𝒮[i,j] += Σt end

        # Streaming term - x
        if (iE == jE) && (iy == jy)
            𝒮[i,j] += fx * 𝒩x[ip,jp]
        end

        # Streaming term - y
        if (iE == jE) && (ix == jx)
            𝒮[i,j] += fy * 𝒩y[ip,jp]
        end

        # CSD term
        if ip == jp && ix == jx && iy == jy
            for kE in range(1,iE-1), wE in range(1,NmE)
                𝒮[i,j] += C[iE] * C[jE] * C[kE] * C[wE] * (1-(-1)^(iE-kE)) * S[wE] * 𝒲[jE,kE,wE]
            end
            𝒮[i,j] += C[iE] * S⁺ * (-1)^(iE-1) * C[jE] * (-1)^(jE-1) * ωE[jE+1]
        end 
    end
end

# Source vector
for ix in range(1,Nmx), iy in range(1,Nmy), iE in range(1,NmE)
    if (~isFC) && (count(>(1),(iE,ix,iy)) ≥ 2) continue end
    fx = -C[ix]/Δx * (g(ix,sx)*ωx[1] + g(ix,-sx))
    fy = -C[iy]/Δy * (g(iy,sy)*ωy[1] + g(iy,-sy))
    for ip in range(1,Np)
        i = index_Exyp(iE,ix,iy,ip)

        # Volume sources
        Q[i] += Qn[ip,index_Exy(iE,ix,iy)]

        # Incoming boundary sources - x
        for jp in range(1,Np)
            Q[i] += fx * 𝒩x[ip,jp] * 𝚽x12[jp,index_Ey(iE,iy)]
        end

        # Incoming boundary sources - y
        for jp in range(1,Np)
            Q[i] += fy * 𝒩y[ip,jp] * 𝚽y12[jp,index_Ex(iE,ix)]
        end

        # CSD incoming sources
        Q[i] += C[iE] * ((-1)^iE*S⁺*ωE[1] + S⁻) * 𝚽E12[ip,index_xy(ix,iy)]
    end
end

# Solve the equation system
𝚽 = 𝒮\Q

# Closure relations
for ip in 1:Np
    for ix in 1:Nmx, iy in 1:Nmy
        if (~isFC) && (count(>(1),(ix,iy)) ≥ 2) continue end
        𝚽E12[ip,index_xy(ix,iy)] = ωE[1] * 𝚽E12[ip,index_xy(ix,iy)]
        for iE in 1:NmE
            if (~isFC) && (count(>(1),(iE,ix,iy)) ≥ 2) continue end
            i = index_Exyp(iE,ix,iy,ip)
            𝚽E12[ip,index_xy(ix,iy)] += C[iE] * (-1)^(iE-1) * ωE[iE+1] * 𝚽[i]
        end
    end
    for iE in 1:NmE, iy in 1:Nmy
        if (~isFC) && (count(>(1),(iE,iy)) ≥ 2) continue end
        𝚽x12[ip,index_Ey(iE,iy)] = ωx[1] * 𝚽x12[ip,index_Ey(iE,iy)]
        for ix in 1:Nmx
            if (~isFC) && (count(>(1),(iE,ix,iy)) ≥ 2) continue end
            i = index_Exyp(iE,ix,iy,ip)
            𝚽x12[ip,index_Ey(iE,iy)] += C[ix] * sx^(ix-1) * ωx[ix+1] * 𝚽[i]
        end
    end
    for iE in 1:NmE, ix in 1:Nmx
        if (~isFC) && (count(>(1),(iE,ix)) ≥ 2) continue end
        𝚽y12[ip,index_Ex(iE,ix)] = ωy[1] * 𝚽y12[ip,index_Ex(iE,ix)]
        for iy in 1:Nmy
            if (~isFC) && (count(>(1),(iE,ix,iy)) ≥ 2) continue end
            i = index_Exyp(iE,ix,iy,ip)
            𝚽y12[ip,index_Ex(iE,ix)] += C[iy] * sy^(iy-1) * ωy[iy+1] * 𝚽[i]
        end
    end
    for iE in 1:NmE, ix in 1:Nmx, iy in 1:Nmy
        if (~isFC) && (count(>(1),(iE,ix,iy)) ≥ 2) continue end
        i = index_Exyp(iE,ix,iy,ip)
        𝚽n[ip,index_Exy(iE,ix,iy)] = 𝚽[i]
    end
end

# Returning solutions
return 𝚽n, 𝚽x12, 𝚽y12, 𝚽E12

end