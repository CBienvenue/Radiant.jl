function pn_1D_BTE(sx::Int64,Σt::Float64,Δx::Float64,Qn::Array{Float64},𝚽x12::Vector{Float64},𝒪x::Int64,Np::Int64,C::Vector{Float64},ωx::Vector{Float64},is_SPH::Bool,pl::Vector{Int64},pm::Vector{Int64})

# Initialization
𝒮 = zeros(𝒪x*Np,𝒪x*Np)
Q = zeros(𝒪x*Np)
𝚽 = Q
𝚽n = copy(Qn)
g(n,sx) = (1+sx)/2 - (1-sx)/2 * (-1)^n
p(l,m) = l^2 + l + m + 1

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,𝒪x), jx in range(1,𝒪x), ip in range(1,Np), jp in range(1,Np)
    i = 𝒪x*(ip-1)+ix
    j = 𝒪x*(jp-1)+jx
    il = pl[ip]
    im = pm[ip]
    jl = pl[jp]
    jm = pm[jp]

    # Collision term
    if (ip == jp) && (ix == jx) 𝒮[i,j] += Σt end

    # Streaming term - x
    if is_SPH
        if il == jl && im == jm
            𝒮[i,j] += C[ix]*sx/Δx * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx))) * (1/2)
        elseif il == jl - 1 && im == jm
            𝒮[i,j] += C[ix]*sx/Δx * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx))) * sqrt((il-abs(im)+1)*(il+abs(im)+1))/(2*sqrt(2*il+1)*sqrt(2*il+3))
        elseif il == jl + 1 && im == jm
            𝒮[i,j] += C[ix]*sx/Δx * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx))) * sqrt((il-abs(im))*(il+abs(im)))/(2*sqrt(2*il-1)*sqrt(2*il+1))
        end
    else
        if ip == jp
            𝒮[i,j] += C[ix]*sx/Δx * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx))) * (1/2)
        elseif ip == jp - 1
            𝒮[i,j] += C[ix]*sx/Δx * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx))) * (il+1)/(2*sqrt(2*il+1)*sqrt(2*il+3))
        elseif ip == jp + 1
            𝒮[i,j] += C[ix]*sx/Δx * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx))) * il/(2*sqrt(2*il-1)*sqrt(2*il+1))
        end
    end
end

# Source vector
for jx in range(1,𝒪x), jp in range(1,Np)
    j = 𝒪x*(jp-1)+jx
    jl = pl[jp]
    jm = pm[jp]
    Q[j] += Qn[jp,jx]
    if is_SPH
        if (jl > 0 && jl-1 ≥ abs(jm)) Q[j] -= C[jx]*sx/Δx * (g(jx-1,sx)*ωx[1]+g(jx-1,-sx)) * sqrt((jl-abs(jm))*(jl+abs(jm)))/(2*sqrt(2*jl-1)*sqrt(2*jl+1)) * 𝚽x12[p(jl-1,jm)] end
        if (jl < sqrt(Np)-1) Q[j] -= C[jx]*sx/Δx * (g(jx-1,sx)*ωx[1]+g(jx-1,-sx)) * sqrt((jl-abs(jm)+1)*(jl+abs(jm)+1))/(2*sqrt(2*jl+1)*sqrt(2*jl+3)) * 𝚽x12[p(jl+1,jm)]  end
        Q[j] -= C[jx]*sx/Δx * (g(jx-1,sx)*ωx[1]+g(jx-1,-sx)) * (1/2) * 𝚽x12[p(jl,jm)]
    else
        if (jp != 1) Q[j] -= C[jx]*sx/Δx * (g(jx-1,sx)*ωx[1]+g(jx-1,-sx)) * jl/(2*sqrt(2*jl-1)*sqrt(2*jl+1)) * 𝚽x12[jp-1] end
        if (jp != Np) Q[j] -= C[jx]*sx/Δx * (g(jx-1,sx)*ωx[1]+g(jx-1,-sx)) * (jl+1)/(2*sqrt(2*jl+1)*sqrt(2*jl+3)) * 𝚽x12[jp+1]  end
        Q[j] -= C[jx]*sx/Δx * (g(jx-1,sx)*ωx[1]+g(jx-1,-sx)) * (1/2) * 𝚽x12[jp]
    end
end

# Solve the equation system
𝚽 = 𝒮\Q

# Closure relations
for jp in range(1,Np)
    𝚽x12[jp] = ωx[1] * 𝚽x12[jp]
    for jx in range(1,𝒪x)
        j = 𝒪x*(jp-1)+jx
        𝚽x12[jp] += C[jx] * sx^(jx-1) * ωx[jx+1] * 𝚽[j]
        𝚽n[jp,jx] = 𝚽[j]
    end
end

# Returning solutions
return 𝚽n, 𝚽x12

end
