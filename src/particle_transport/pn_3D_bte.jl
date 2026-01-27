function pn_3D_BTE(sx::Int64,sy::Int64,sz::Int64,Σt::Float64,Δx::Float64,Δy::Float64,Δz::Float64,Qn::Array{Float64},𝚽x12::Array{Float64},𝚽y12::Array{Float64},𝚽z12::Array{Float64},𝒪x::Int64,𝒪y::Int64,𝒪z::Int64,Np::Int64,C::Vector{Float64},ωx::Array{Float64},ωy::Array{Float64},ωz::Array{Float64},𝒩x::Matrix{Float64},𝒩y::Matrix{Float64},𝒩z::Matrix{Float64},isFC::Bool)

# Initialization
if isFC Nm = 𝒪x*𝒪y*𝒪z*Np else Nm = (𝒪x+𝒪y+𝒪z-2)*Np end
𝒮 = zeros(Nm,Nm)
Q = zeros(Nm)
𝚽 = Q
𝚽n = copy(Qn)
g(n,sx) = (1+sx)/2 - (1-sx)/2 * (-1)^n
function jyz(jy,jz)
    if isFC 
        return 𝒪y*(jz-1) + jy
    else  
        j = 1 + (jy-1) + (jz-1)
        if jz > 1 j += 𝒪y-1 end
        return j
    end
end
function jxyz(jx,jy,jz)
    if isFC 
        return 𝒪y*𝒪x*(jz-1) + 𝒪x * (jy-1) + jx
    else  
        j = 1 + (jx-1) + (jy-1) + (jz-1)
        if jy > 1 j += 𝒪x-1 end
        if jz > 1 j += 𝒪x-1 + 𝒪y-1 end
        return j
    end
end
function jxyzp(jx,jy,jz,jp)
    if isFC
        j = 𝒪x*𝒪y*𝒪z*(jp-1) + 𝒪x*𝒪y*(jz-1) + 𝒪x*(jy-1) + jx
    else
        j = 1 + (jx-1) + (jy-1) + (jz-1)
        if jy > 1 j += 𝒪x-1 end
        if jz > 1 j += 𝒪x-1 + 𝒪y-1 end
        j += (jp-1)*(𝒪x+𝒪y+𝒪z-2)
    end
end

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,𝒪x), jx in range(1,𝒪x), iy in range(1,𝒪y), jy in range(1,𝒪y), iz in range(1,𝒪z), jz in range(1,𝒪z)
    fx = C[ix]*sx/Δx * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1,jy,jz] - (jx ≤ ix-1)*(1-(-1)^(ix-jx)))
    fy = C[iy]*sy/Δy * C[jy] * (g(iy-1,sy)*sy^(jy-1)*ωy[jy+1,jx,jz] - (jy ≤ iy-1)*(1-(-1)^(iy-jy)))
    fz = C[iz]*sz/Δz * C[jz] * (g(iz-1,sz)*sz^(jz-1)*ωz[jz+1,jx,jy] - (jz ≤ iz-1)*(1-(-1)^(iz-jz)))
    for ip in range(1,Np), jp in range(1,Np)
        if (~isFC) && (count(>(1),(ix,iy,iz)) ≥ 2 || count(>(1),(jx,jy,jz)) ≥ 2) continue end
        i = jxyzp(ix,iy,iz,ip)
        j = jxyzp(jx,jy,jz,jp)

        # Collision term
        if (i == j) 𝒮[i,j] += Σt end

        # Streaming term - x
        if (iy == jy) && (iz == jz)
            𝒮[i,j] += fx * 𝒩x[ip,jp]
        end

        # Streaming term - y
        if (ix == jx) && (iz == jz)
            if sy == sz
                𝒮[i,j] += fy * 𝒩y[ip,jp]
            else
                𝒮[i,j] += fy * 𝒩z[ip,jp]
            end
        end

        # Streaming term - z
        if (ix == jx) && (iy == jy)
            if sy == sz
                𝒮[i,j] += fz * 𝒩z[ip,jp]
            else
                𝒮[i,j] += fz * 𝒩y[ip,jp]
            end
        end
    end
end

# Source vector
for jx in range(1,𝒪x), jy in range(1,𝒪y), jz in range(1,𝒪z)
    fx = -C[jx]*sx/Δx * (g(jx-1,sx)*ωx[1,jy,jz]+g(jx-1,-sx))
    fy = -C[jy]*sy/Δy * (g(jy-1,sy)*ωy[1,jx,jz]+g(jy-1,-sy))
    fz = -C[jz]*sz/Δz * (g(jz-1,sz)*ωz[1,jx,jy]+g(jz-1,-sz))
    for jp in range(1,Np)
        if (~isFC) && (count(>(1),(jx,jy,jz)) ≥ 2) continue end
        j = jxyzp(jx,jy,jz,jp)

        # Volume sources
        Q[j] += Qn[jp,jxyz(jx,jy,jz)]

        # Incoming boundary sources - x
        for ip in range(1,Np)
            Q[j] += fx * 𝒩x[ip,jp] * 𝚽x12[ip,jyz(jy,jz)]
        end

        # Incoming boundary sources - y
        for ip in range(1,Np)
            if sy == sz
                Q[j] += fy * 𝒩y[ip,jp] * 𝚽y12[ip,jyz(jx,jz)]
            else
                Q[j] += fy * 𝒩z[ip,jp] * 𝚽y12[ip,jyz(jx,jz)]
            end
        end

        # Incoming boundary sources - z
        for ip in range(1,Np)
            if sy == sz
                Q[j] += fz * 𝒩z[ip,jp] * 𝚽z12[ip,jyz(jx,jy)]
            else
                Q[j] += fz * 𝒩y[ip,jp] * 𝚽z12[ip,jyz(jx,jy)]
            end
        end
    end
end

# Solve the equation system
𝚽 = 𝒮\Q

# Closure relations
for jp in 1:Np
    for jy in 1:𝒪y, jz in 1:𝒪z
        if (~isFC) && (count(>(1),(jy,jz)) ≥ 2) continue end
        𝚽x12[jp,jyz(jy,jz)] = ωx[1,jy,jz] * 𝚽x12[jp,jyz(jy,jz)]
        for jx in 1:𝒪x
            if (~isFC) && (count(>(1),(jx,jy,jz)) ≥ 2) continue end
            j = jxyzp(jx,jy,jz,jp)
            𝚽x12[jp,jyz(jy,jz)] += C[jx] * sx^(jx-1) * ωx[jx+1,jy,jz] * 𝚽[j]
        end
    end
    for jx in 1:𝒪x, jz in 1:𝒪z
        if (~isFC) && (count(>(1),(jx,jz)) ≥ 2) continue end
        𝚽y12[jp,jyz(jx,jz)] = ωy[1,jx,jz] * 𝚽y12[jp,jyz(jx,jz)]
        for jy in 1:𝒪y
            if (~isFC) && (count(>(1),(jx,jy,jz)) ≥ 2) continue end
            j = jxyzp(jx,jy,jz,jp)
            𝚽y12[jp,jyz(jx,jz)] += C[jy] * sy^(jy-1) * ωy[jy+1,jx,jz] * 𝚽[j]
        end
    end
    for jx in 1:𝒪x, jy in 1:𝒪y
        if (~isFC) && (count(>(1),(jx,jy)) ≥ 2) continue end
        𝚽z12[jp,jyz(jx,jy)] = ωz[1,jx,jy] * 𝚽z12[jp,jyz(jx,jy)]
        for jz in 1:𝒪z
            if (~isFC) && (count(>(1),(jx,jy,jz)) ≥ 2) continue end
            j = jxyzp(jx,jy,jz,jp)
            𝚽z12[jp,jyz(jx,jy)] += C[jz] * sz^(jz-1) * ωz[jz+1,jx,jy] * 𝚽[j]
        end
    end
    for jx in 1:𝒪x, jy in 1:𝒪y, jz in 1:𝒪z
        if (~isFC) && (count(>(1),(jx,jy,jz)) ≥ 2) continue end
        j = jxyzp(jx,jy,jz,jp)
        𝚽n[jp,jxyz(jx,jy,jz)] = 𝚽[j]
    end
end

# Returning solutions
return 𝚽n, 𝚽x12, 𝚽y12, 𝚽z12

end