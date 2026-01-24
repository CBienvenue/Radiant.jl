function pn_2D_BTE(sx::Int64,sy::Int64,Σt::Float64,Δx::Float64,Δy::Float64,Qn::Array{Float64},𝚽x12::Array{Float64},𝚽y12::Array{Float64},𝒪x::Int64,𝒪y::Int64,Np::Int64,C::Vector{Float64},ωx::Array{Float64},ωy::Array{Float64},𝒩x::Matrix{Float64},𝒩y::Matrix{Float64},isFC::Bool)

# Initialization
if isFC Nm = 𝒪x*𝒪y*Np else Nm = (𝒪x+𝒪y-1)*Np end
𝒮 = zeros(Nm,Nm)
Q = zeros(Nm)
𝚽 = Q
𝚽n = copy(Qn)
g(n,sx) = (1+sx)/2 - (1-sx)/2 * (-1)^n
function jxy(jx,jy)
    if isFC 
        return 𝒪x*(jy-1)+jx 
    else  
        idx = 1 + (jy-1) + (jx-1)
        if jy > 1 idx += 𝒪x-1 end
        return idx
    end
end
function jxyp(jx,jy,jp)
    if isFC
        j = 𝒪x*𝒪y*(jp-1) + 𝒪x*(jy-1) + jx
    else
        j = 1 + (jy-1) + (jx-1)
        if jy > 1 j += 𝒪x-1 end
        j += (jp-1)*(𝒪x+𝒪y-1)
    end
end

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,𝒪x), jx in range(1,𝒪x), iy in range(1,𝒪y), jy in range(1,𝒪y)
    fx = C[ix]*sx/Δx * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1,jy,iy] - (jx ≤ ix-1)*(1-(-1)^(ix-jx)))
    fy = C[iy]*sy/Δy * C[jy] * (g(iy-1,sy)*sy^(jy-1)*ωy[jy+1,jx,ix] - (jy ≤ iy-1)*(1-(-1)^(iy-jy)))
    for ip in range(1,Np), jp in range(1,Np)
        if (~isFC) && (count(>(1),(ix,iy)) ≥ 2 || count(>(1),(jx,jy)) ≥ 2) continue end
        i = jxyp(ix,iy,ip)
        j = jxyp(jx,jy,jp)

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
for jx in range(1,𝒪x), jy in range(1,𝒪y)
    fx = -C[jx]*sx/Δx * (g(jx-1,sx)*ωx[1,jy,jy]+g(jx-1,-sx))
    fy = -C[jy]*sy/Δy * (g(jy-1,sy)*ωy[1,jx,jx]+g(jy-1,-sy))
    for jp in range(1,Np)
        if (~isFC) && (count(>(1),(jx,jy)) ≥ 2) continue end
        j = jxyp(jx,jy,jp)

        # Volume sources
        Q[j] += Qn[jp,jxy(jx,jy)]

        # Incoming boundary sources - x
        for ip in range(1,Np)
            Q[j] += fx * 𝒩x[ip,jp] * 𝚽x12[ip,jy]
        end

        # Incoming boundary sources - y
        for ip in range(1,Np)
            Q[j] += fy * 𝒩y[ip,jp] * 𝚽y12[ip,jx]
        end
    end
end

# Solve the equation system
𝚽 = 𝒮\Q

# Closure relations
for jp in 1:Np
    for jy in 1:𝒪y
        𝚽x12[jp,jy] = ωx[1,jy,jy] * 𝚽x12[jp,jy]
        for jx in 1:𝒪x
            if (~isFC) && (count(>(1),(jx,jy)) ≥ 2) continue end
            j = jxyp(jx,jy,jp)
            𝚽x12[jp,jy] += C[jx] * sx^(jx-1) * ωx[jx+1,jy,1] * 𝚽[j]
        end
    end
    for jx in 1:𝒪x
        𝚽y12[jp,jx] = ωy[1,jx,jx] * 𝚽y12[jp,jx]
        for jy in 1:𝒪y
            if (~isFC) && (count(>(1),(jx,jy)) ≥ 2) continue end
            j = jxyp(jx,jy,jp)
            𝚽y12[jp,jx] += C[jy] * sy^(jy-1) * ωy[jy+1,jx,1] * 𝚽[j]
        end
    end
    for jx in 1:𝒪x, jy in 1:𝒪y
            if (~isFC) && (count(>(1),(jx,jy)) ≥ 2) continue end
            j = jxyp(jx,jy,jp)
        𝚽n[jp,jxy(jx,jy)] = 𝚽[j]
    end
end

# Returning solutions
return 𝚽n, 𝚽x12, 𝚽y12

end