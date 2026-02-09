function gn_3D_BTE(sx::Int64,sy::Int64,sz::Int64,Σt::Float64,Δx::Float64,Δy::Float64,Δz::Float64,Qn::Array{Float64},𝚽x12::Array{Float64},𝚽y12::Array{Float64},𝚽z12::Array{Float64},Nmx::Int64,Nmy::Int64,Nmz::Int64,Np::Int64,C::Vector{Float64},ωx::Array{Float64},ωy::Array{Float64},ωz::Array{Float64},𝒩x::Matrix{Float64},𝒩y::Matrix{Float64},𝒩z::Matrix{Float64},isFC::Bool)

# Initialization
if isFC Nm = Nmx*Nmy*Nmz*Np else Nm = (Nmx+Nmy+Nmz-2)*Np end
𝒮 = zeros(Nm,Nm)
Q = zeros(Nm)
𝚽 = Q
𝚽n = copy(Qn)
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
    fx = -C[ix]/Δx * (g(ix,sx)*ωx[1,iy,iz]+g(ix,-sx))
    fy = -C[iy]/Δy * (g(iy,sy)*ωy[1,ix,iz]+g(iy,-sy))
    fz = -C[iz]/Δz * (g(iz,sz)*ωz[1,ix,iy]+g(iz,-sz))
    for ip in range(1,Np)
        if (~isFC) && (count(>(1),(ix,iy,iz)) ≥ 2) continue end
        j = index_xyzp(ix,iy,iz,ip)

        # Volume sources
        Q[j] += Qn[ip,index_xyz(ix,iy,iz)]

        # Incoming boundary sources - x
        for jp in range(1,Np)
            Q[j] += fx * 𝒩x[ip,jp] * 𝚽x12[jp,index_xy(iy,iz)]
        end

        # Incoming boundary sources - y
        for jp in range(1,Np)
            Q[j] += fy * 𝒩y[ip,jp] * 𝚽y12[jp,index_xy(ix,iz)]
        end

        # Incoming boundary sources - z
        for jp in range(1,Np)
            Q[j] += fz * 𝒩z[ip,jp] * 𝚽z12[jp,index_xy(ix,iy)]
        end
    end
end

# Solve the equation system
𝚽 = 𝒮\Q

# Closure relations
for ip in 1:Np
    for iy in 1:Nmy, iz in 1:Nmz
        if (~isFC) && (count(>(1),(iy,iz)) ≥ 2) continue end
        𝚽x12[ip,index_xy(iy,iz)] = ωx[1,iy,iz] * 𝚽x12[ip,index_xy(iy,iz)]
        for ix in 1:Nmx
            if (~isFC) && (count(>(1),(ix,iy,iz)) ≥ 2) continue end
            j = index_xyzp(ix,iy,iz,ip)
            𝚽x12[ip,index_xy(iy,iz)] += C[ix] * sx^(ix-1) * ωx[ix+1,iy,iz] * 𝚽[j]
        end
    end
    for ix in 1:Nmx, iz in 1:Nmz
        if (~isFC) && (count(>(1),(ix,iz)) ≥ 2) continue end
        𝚽y12[ip,index_xy(ix,iz)] = ωy[1,ix,iz] * 𝚽y12[ip,index_xy(ix,iz)]
        for iy in 1:Nmy
            if (~isFC) && (count(>(1),(ix,iy,iz)) ≥ 2) continue end
            j = index_xyzp(ix,iy,iz,ip)
            𝚽y12[ip,index_xy(ix,iz)] += C[iy] * sy^(iy-1) * ωy[iy+1,ix,iz] * 𝚽[j]
        end
    end
    for ix in 1:Nmx, iy in 1:Nmy
        if (~isFC) && (count(>(1),(ix,iy)) ≥ 2) continue end
        𝚽z12[ip,index_xy(ix,iy)] = ωz[1,ix,iy] * 𝚽z12[ip,index_xy(ix,iy)]
        for iz in 1:Nmz
            if (~isFC) && (count(>(1),(ix,iy,iz)) ≥ 2) continue end
            j = index_xyzp(ix,iy,iz,ip)
            𝚽z12[ip,index_xy(ix,iy)] += C[iz] * sz^(iz-1) * ωz[iz+1,ix,iy] * 𝚽[j]
        end
    end
    for ix in 1:Nmx, iy in 1:Nmy, iz in 1:Nmz
        if (~isFC) && (count(>(1),(ix,iy,iz)) ≥ 2) continue end
        j = index_xyzp(ix,iy,iz,ip)
        𝚽n[ip,index_xyz(ix,iy,iz)] = 𝚽[j]
    end
end

# Returning solutions
return 𝚽n, 𝚽x12, 𝚽y12, 𝚽z12

end