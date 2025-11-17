function pn_1D_BTE(sx::Int64,Σt::Float64,Δx::Float64,Qn::Array{Float64},𝚽x12::Vector{Float64},𝒪x::Int64,Np::Int64,C::Vector{Float64},ωx::Vector{Float64})

# Initialization
𝒮 = zeros(𝒪x*Np,𝒪x*Np)
Q = zeros(𝒪x*Np)
𝚽 = Q
𝚽n = copy(Qn)
g(n,sx) = (1+sx)/2 - (1-sx)/2 * (-1)^n

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,𝒪x), jx in range(1,𝒪x), ip in range(1,Np), jp in range(1,Np)
    i = 𝒪x*(ip-1)+ix
    j = 𝒪x*(jp-1)+jx
    il = ip - 1
    if ip == jp
        if (ix == jx) 𝒮[i,j] += Σt end
        𝒮[i,j] += C[ix]*sx/(2*Δx) * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx)))
    elseif ip == jp - 1
        𝒮[i,j] += C[ix]*sx/(2*Δx) * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx))) * (il+1)/(sqrt(2*il+1)*sqrt(2*il+3))
    elseif ip == jp + 1
        𝒮[i,j] += C[ix]*sx/(2*Δx) * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx))) * il/(sqrt(2*il-1)*sqrt(2*il+1))
    end
end

# Source vector
for jx in range(1,𝒪x), jp in range(1,Np)
    j = 𝒪x*(jp-1)+jx
    jl = jp - 1
    Q[j] += Qn[jp,jx]
    if (jp != 1) Q[j] -= C[jx]*sx/(2*Δx) * (g(jx-1,sx)*ωx[1]+g(jx-1,-sx)) * jl/(sqrt(2*jl-1)*sqrt(2*jl+1)) * 𝚽x12[jp-1] end
    if (jp != Np) Q[j] -= C[jx]*sx/(2*Δx) * (g(jx-1,sx)*ωx[1]+g(jx-1,-sx)) * (jl+1)/(sqrt(2*jl+1)*sqrt(2*jl+3)) * 𝚽x12[jp+1]  end
    Q[j] -= C[jx]*sx/(2*Δx) * (g(jx-1,sx)*ωx[1]+g(jx-1,-sx)) * 𝚽x12[jp] 
end

# Solve the equation system
𝚽 = 𝒮\Q

# Closure relation
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
