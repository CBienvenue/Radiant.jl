function pn_1D_BTE(sx::Int64,Σt::Float64,Δx::Float64,Qn::Array{Float64},𝚽x12::Vector{Float64},𝒪x::Int64,Np::Int64,C::Vector{Float64},ωx::Vector{Float64},is_SPH::Bool,pl::Vector{Int64},pm::Vector{Int64},PN_model::Int64,pa,pb,pc,𝒩⁻,𝒩,𝒩⁺)

# Initialization
𝒮 = zeros(𝒪x*Np,𝒪x*Np)
Q = zeros(𝒪x*Np)
𝚽 = Q
𝚽n = copy(Qn)
g(n,sx) = (1+sx)/2 - (1-sx)/2 * (-1)^n
p(l,m) = l^2 + l + m + 1
p(l,a,b) = div(l*(l+1)*(l+2),6) + (l+1)*a - div(a*(a-1),2) + b + 1

# Matrix of Legendre moment coefficients of the flux
for ix in range(1,𝒪x), jx in range(1,𝒪x), ip in range(1,Np), jp in range(1,Np)
    i = 𝒪x*(ip-1)+ix
    j = 𝒪x*(jp-1)+jx
    il = pl[ip]
    jl = pl[jp]
    if PN_model == 2
        im = pm[ip]
        jm = pm[jp]
    elseif PN_model == 3
        ia = pa[ip]
        ib = pb[ip]
        ic = pc[ip]
        ja = pa[jp]
        jb = pb[jp]
        jc = pc[jp]
    end

    # Collision term
    if (ip == jp) && (ix == jx) 𝒮[i,j] += Σt end

    # Streaming term - x
    if PN_model == 1 # Legendre polynomials
        if ip == jp
            𝒮[i,j] += C[ix]*sx/Δx * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx))) * 𝒩[ip]
        elseif ip == jp - 1
            𝒮[i,j] += C[ix]*sx/Δx * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx))) * 𝒩⁺[ip]
        elseif ip == jp + 1
            𝒮[i,j] += C[ix]*sx/Δx * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx))) * 𝒩⁻[ip]
        end
    elseif PN_model == 2 # Spherical harmonics
        if il == jl && im == jm
            𝒮[i,j] += C[ix]*sx/Δx * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx))) * 𝒩[ip]
        elseif il == jl - 1 && im == jm
            𝒮[i,j] += C[ix]*sx/Δx * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx))) * 𝒩⁺[ip]
        elseif il == jl + 1 && im == jm
            𝒮[i,j] += C[ix]*sx/Δx * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx))) * 𝒩⁻[ip]
        end
    elseif PN_model == 3 # Cartesian harmonics
        if ib == jb && ic == jc
            if il == jl && ia == ja
                𝒮[i,j] += C[ix]*sx/Δx * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx))) * 𝒩[ip]
            elseif il == jl - 1 && ia == ja - 1
                𝒮[i,j] += C[ix]*sx/Δx * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx))) * 𝒩⁺[ip]
            elseif il == jl + 1 && ia == ja + 1
                𝒮[i,j] += C[ix]*sx/Δx * C[jx] * (g(ix-1,sx)*sx^(jx-1)*ωx[jx+1] - (jx ≤ ix-1)*(1-(-1)^(ix-jx))) * 𝒩⁻[ip]
            end
        end
    end
end

# Source vector
for jx in range(1,𝒪x), jp in range(1,Np)
    j = 𝒪x*(jp-1)+jx
    jl = pl[jp]
    if PN_model == 2
        jm = pm[jp]
    else
        ja = pa[jp]
        jb = pb[jp]
        jc = pc[jp]
    end
    Q[j] += Qn[jp,jx]
    if PN_model == 1 # Legendre polynomials
        if (jp != 1) Q[j] -= C[jx]*sx/Δx * (g(jx-1,sx)*ωx[1]+g(jx-1,-sx)) * 𝒩⁻[jp] * 𝚽x12[jp-1] end
        if (jp != Np) Q[j] -= C[jx]*sx/Δx * (g(jx-1,sx)*ωx[1]+g(jx-1,-sx)) * 𝒩⁺[jp] * 𝚽x12[jp+1]  end
        Q[j] -= C[jx]*sx/Δx * (g(jx-1,sx)*ωx[1]+g(jx-1,-sx)) * 𝒩[jp] * 𝚽x12[jp]
    elseif PN_model == 2 # Spherical harmonics
        if (jl > 0 && jl-1 ≥ abs(jm)) Q[j] -= C[jx]*sx/Δx * (g(jx-1,sx)*ωx[1]+g(jx-1,-sx)) * 𝒩⁻[jp] * 𝚽x12[p(jl-1,jm)] end
        if (jl < sqrt(Np)-1) Q[j] -= C[jx]*sx/Δx * (g(jx-1,sx)*ωx[1]+g(jx-1,-sx)) * 𝒩⁺[jp] * 𝚽x12[p(jl+1,jm)]  end
        Q[j] -= C[jx]*sx/Δx * (g(jx-1,sx)*ωx[1]+g(jx-1,-sx)) * 𝒩[jp] * 𝚽x12[p(jl,jm)]
    elseif PN_model == 3 # Cartesian harmonics
        if (jl > 0 && ja > 0) Q[j] -= C[jx]*sx/Δx * (g(jx-1,sx)*ωx[1]+g(jx-1,-sx)) * 𝒩⁻[jp] * 𝚽x12[p(jl-1,ja-1,jb)] end
        if p(jl+1,ja+1,jb) ≤ length(𝚽x12) Q[j] -= C[jx]*sx/Δx * (g(jx-1,sx)*ωx[1]+g(jx-1,-sx)) * 𝒩⁺[jp] * 𝚽x12[p(jl+1,ja+1,jb)] end
        Q[j] -= C[jx]*sx/Δx * (g(jx-1,sx)*ωx[1]+g(jx-1,-sx)) * 𝒩[jp] * 𝚽x12[p(jl,ja,jb)]
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