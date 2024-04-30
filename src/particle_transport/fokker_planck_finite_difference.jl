"""
    fokker_planck_finite_difference(N::Int64,quadrature_type::String,Ndims::Int64)

Calculate the Fokker-Planck scattering matrix using finite-difference scheme.

# Input Argument(s)
- 'N::Int64': number of directions.
- 'quadrature_type::String': type of quadrature.
- 'Ndims::Int64': dimension of the geometry.

# Output Argument(s)
- 'ℳ::Array{Float64}': Fokker-Planck scattering matrix.
- 'λ₀::Float64': correction factor for total cross-section.

# Reference(s)
- Morel (1985) : An Improved Fokker-Planck Angular Differencing Scheme.
- Landesman (1988) : Angular Fokker-Pianck Decomposition and Representation Techniques.
- Larsen (2010) : Advances in Discrete-Ordinates Methodology.
- Morel (2007) : A Discretization Scheme for the Three-Dimensional Angular Fokker-Planck
  Operator.

"""
function fokker_planck_finite_difference(N::Int64,quadrature_type::String,Ndims::Int64,pℓ::Vector{Int64},pm::Vector{Int64},P::Int64,Nd::Int64,Mn::Array{Float64},Dn::Array{Float64})

Mn_FP = Mn
Dn_FP = Dn
N_FP = Nd

if Ndims == 1 

    # Extract quadrature informations
    μ,w = quadrature(N,quadrature_type,Ndims)

    # Reorder evaluation points and weights
    index = sortperm(μ)
    μn = μ[index]
    wn = w[index]

    # Compute the scattering matrix
    β12 = zeros(N+1)
    @inbounds for n in range(2,N)
        β12[n] = β12[n-1] - 2*wn[n-1]*μn[n-1]
    end

    Cn⁻ = zeros(N)
    @inbounds for n in range(2,N)
        Cn⁻[n] = β12[n]/(wn[n]*(μn[n]-μn[n-1]))
    end

    Cn⁺ = zeros(N)
    @inbounds for n in range(1,N-1)
        Cn⁺[n] = β12[n+1]/(wn[n]*(μn[n+1]-μn[n]))
    end
    
    Cn = Cn⁻ + Cn⁺
    λ₀ = 0
    ℳ_temp = zeros(N,N)
    @inbounds for n in range(1,N), m in range(1,N)
        if n == m
            ℳ_temp[n,m] = λ₀ - Cn[n]
        elseif m == n - 1
            ℳ_temp[n,m] = Cn⁻[n]
        elseif m == n + 1
            ℳ_temp[n,m] = Cn⁺[n]
        end
    end

    # Reorder the scattering matrix
    ℳ = zeros(N,N)
    @inbounds for n in range(1,N), m in range(1,N)
        ℳ[index[n],index[m]] = ℳ_temp[n,m]
    end

    #=
    Pℓ = zeros(N,P)
    @inbounds for n in range(1,N)
        Pℓ[n,:] = legendre_polynomials(P-1,μ[n])
    end

    Mn_FP = zeros(N,P)
    Dn_FP = zeros(P,N)
    for p in range(1,P),n in range(1,N)
        Mn_FP[n,p] = (2*pℓ[p]+1)/2 * Pℓ[n,p]
        Dn_FP[p,n] = w[n] * Pℓ[n,p]
    end
    Mn_FP = pinv(Dn_FP)
    N_FP = N
    =#

elseif Ndims == 2

    Ω_GC,w_GC = quadrature(N,"gauss-legendre-chebychev",Ndims)
    μ_GC = Ω_GC[1]; η_GC = Ω_GC[2]; ξ_GC = Ω_GC[3];
    N_GC = length(w_GC)
    Rℓm = zeros(N_GC,P)
    @inbounds for n in range(1,N_GC)
        if ξ_GC[n] != 0
            if η_GC[n] > 0
                ϕ = sign(ξ_GC[n]) * atan(abs(ξ_GC[n]/η_GC[n]))
            elseif η_GC[n] == 0
                ϕ = sign(ξ_GC[n]) * π/2
            elseif η_GC[n] < 0
                ϕ = sign(ξ_GC[n]) * (π - atan(abs(ξ_GC[n]/η_GC[n])))
            end
        else
            if η_GC[n] >= 0
                ϕ = 0.0
            elseif η_GC[n] < 0
                ϕ = 1.0 * π
            end
        end
        for p in range(1,P)
            Rℓm[n,p] = real_spherical_harmonics(pℓ[p],pm[p],μ_GC[n],ϕ)
        end
    end

    Mn_FP = zeros(N_GC,P)
    Dn_FP = zeros(P,N_GC)
    for p in range(1,P),n in range(1,N_GC)
        Mn_FP[n,p] = (2*pℓ[p]+1)/(4*π) * Rℓm[n,p]
        Dn_FP[p,n] = w_GC[n] * Rℓm[n,p]
    end
    N_FP = N_GC

    # Initialization and quadrature parameters
    λ₀ = 0.0
    μᵢ,wᵢ = gauss_legendre(N)
    index = sortperm(μᵢ); μᵢ = μᵢ[index]; wᵢ = wᵢ[index]
    Nᵢ = length(μᵢ)
    wⱼ = π/Nᵢ
    ωⱼ = zeros(1,Nᵢ)
    for n in range(1,Nᵢ) ωⱼ[n] = π*(n-0.5)/Nᵢ end

    # Compute the parameters
    β12 = zeros(Nᵢ+1); d12 = zeros(Nᵢ+1); γn = zeros(1,Nᵢ)
    for n in range(2,Nᵢ)
        β12[n] = β12[n-1] - 2*wᵢ[n-1]*μᵢ[n-1]
        d12[n] = (sqrt(1-μᵢ[n]^2)-sqrt(1-μᵢ[n-1]^2))/(μᵢ[n]-μᵢ[n-1])
    end
    for n in range(1,Nᵢ)
        cn = (β12[n+1]*d12[n+1] - β12[n]*d12[n])/wᵢ[n]
        Kn = 2*(1-μᵢ[n]^2)+cn*sqrt(1-μᵢ[n]^2)
        γn[n] = (π^2*Kn)/(2*Nᵢ^2*(1-cos(π/Nᵢ)))
    end

    # Compute the scattering matrix
    ℳ = zeros(Nᵢ^2,Nᵢ^2)
    for n in range(1,Nᵢ), i in range(1,Nᵢ), m in range(1,Nᵢ), j in range(1,Nᵢ)

        ii = i+Nᵢ*(n-1)
        jj = j+Nᵢ*(m-1)

        if m == n && j == i
            if n == 1
                ℳ[ii,jj] -= 1/wᵢ[n] * β12[n+1]/(μᵢ[n+1]-μᵢ[n])
            elseif n == Nᵢ
                ℳ[ii,jj] -= 1/wᵢ[n] * β12[n]/(μᵢ[n]-μᵢ[n-1])
            else
                ℳ[ii,jj] -= 1/wᵢ[n] * ( β12[n+1]/(μᵢ[n+1]-μᵢ[n]) + β12[n]/(μᵢ[n]-μᵢ[n-1]) )
            end
            if i == 1
                ℳ[ii,jj] -= 1/(1-μᵢ[n]^2)*γn[n]/wⱼ/(ωⱼ[i+1]-ωⱼ[i])
            elseif i == Nᵢ
                ℳ[ii,jj] -= 1/(1-μᵢ[n]^2)*γn[n]/wⱼ/(ωⱼ[i]-ωⱼ[i-1])
            else
                ℳ[ii,jj] -= 1/(1-μᵢ[n]^2)*γn[n]/wⱼ/(ωⱼ[i]-ωⱼ[i-1]) + 1/(1-μᵢ[n]^2)*γn[n]/wⱼ/(ωⱼ[i+1]-ωⱼ[i])
            end

        elseif m == n - 1 && j == i
            if n != 1
                ℳ[ii,jj] += 1/wᵢ[n] * β12[n]/(μᵢ[n]-μᵢ[n-1])
            end
            
        elseif m == n + 1 && j == i
            if n != Nᵢ
                ℳ[ii,jj] += 1/wᵢ[n] * β12[n+1]/(μᵢ[n+1]-μᵢ[n])
            end

        elseif m == n && (j == i - 1 || i == 1 && j == Nᵢ)
            if i != 1
                ℳ[ii,jj] += 1/(1-μᵢ[n]^2)*γn[n]/wⱼ/(ωⱼ[i]-ωⱼ[i-1])
            end

        elseif m == n && (j == i + 1 || i == Nᵢ && j == 1)
            if i != Nᵢ
                ℳ[ii,jj] += 1/(1-μᵢ[n]^2)*γn[n]/wⱼ/(ωⱼ[i+1]-ωⱼ[i])
            end
        end
    end

elseif Ndims == 3

    Ω_GC,w_GC = quadrature(N,"gauss-legendre-chebychev",Ndims)
    μ_GC = Ω_GC[1]; η_GC = Ω_GC[2]; ξ_GC = Ω_GC[3];
    N_GC = length(w_GC)
    Rℓm = zeros(N_GC,P)
    @inbounds for n in range(1,N_GC)
        if ξ_GC[n] != 0
            if η_GC[n] > 0
                ϕ = sign(ξ_GC[n]) * atan(abs(ξ_GC[n]/η_GC[n]))
            elseif η_GC[n] == 0
                ϕ = sign(ξ_GC[n]) * π/2
            elseif η_GC[n] < 0
                ϕ = sign(ξ_GC[n]) * (π - atan(abs(ξ_GC[n]/η_GC[n])))
            end
        else
            if η_GC[n] >= 0
                ϕ = 0.0
            elseif η_GC[n] < 0
                ϕ = 1.0 * π
            end
        end
        for p in range(1,P)
            Rℓm[n,p] = real_spherical_harmonics(pℓ[p],pm[p],μ_GC[n],ϕ)
        end
    end

    Mn_FP = zeros(N_GC,P)
    Dn_FP = zeros(P,N_GC)
    for p in range(1,P),n in range(1,N_GC)
        Mn_FP[n,p] = (2*pℓ[p]+1)/(4*π) * Rℓm[n,p]
        Dn_FP[p,n] = w_GC[n] * Rℓm[n,p]
    end
    N_FP = N_GC

    # Initialization and quadrature parameters
    λ₀ = 0.0
    μᵢ,wᵢ = gauss_legendre(N)
    index = sortperm(μᵢ); μᵢ = μᵢ[index]; wᵢ = wᵢ[index]
    Nᵢ = length(μᵢ)
    wⱼ = π/Nᵢ
    ωⱼ = zeros(1,2*Nᵢ)
    for n in range(1,2*Nᵢ) ωⱼ[n] = π*(n-0.5)/Nᵢ end

    # Compute the parameters
    β12 = zeros(Nᵢ+1); d12 = zeros(Nᵢ+1); γn = zeros(1,Nᵢ)
    for n in range(2,Nᵢ)
        β12[n] = β12[n-1] - 2*wᵢ[n-1]*μᵢ[n-1]
        d12[n] = (sqrt(1-μᵢ[n]^2)-sqrt(1-μᵢ[n-1]^2))/(μᵢ[n]-μᵢ[n-1])
    end
    for n in range(1,Nᵢ)
        cn = (β12[n+1]*d12[n+1] - β12[n]*d12[n])/wᵢ[n]
        Kn = 2*(1-μᵢ[n]^2)+cn*sqrt(1-μᵢ[n]^2)
        γn[n] = (π^2*Kn)/(2*Nᵢ^2*(1-cos(π/Nᵢ)))
    end

    # Compute the scattering matrix
    ℳ = zeros(2*Nᵢ^2,2*Nᵢ^2)
    for n in range(1,Nᵢ), i in range(1,2*Nᵢ), m in range(1,Nᵢ), j in range(1,2*Nᵢ)

        ii = i+2*Nᵢ*(n-1)
        jj = j+2*Nᵢ*(m-1)

        if m == n && j == i
            if n == 1
                ℳ[ii,jj] -= 1/wᵢ[n] * β12[n+1]/(μᵢ[n+1]-μᵢ[n])
            elseif n == Nᵢ
                ℳ[ii,jj] -= 1/wᵢ[n] * β12[n]/(μᵢ[n]-μᵢ[n-1])
            else
                ℳ[ii,jj] -= 1/wᵢ[n] * ( β12[n+1]/(μᵢ[n+1]-μᵢ[n]) + β12[n]/(μᵢ[n]-μᵢ[n-1]) )
            end
            if i == 1
                ℳ[ii,jj] -= 1/(1-μᵢ[n]^2)*γn[n]/wⱼ/(ωⱼ[i]+2*π-ωⱼ[end]) + 1/(1-μᵢ[n]^2)*γn[n]/wⱼ/(ωⱼ[i+1]-ωⱼ[i])
            elseif i == 2*Nᵢ
                ℳ[ii,jj] -= 1/(1-μᵢ[n]^2)*γn[n]/wⱼ/(ωⱼ[i]-ωⱼ[i-1]) + 1/(1-μᵢ[n]^2)*γn[n]/wⱼ/(ωⱼ[1]+2*π-ωⱼ[i])
            else
                ℳ[ii,jj] -= 1/(1-μᵢ[n]^2)*γn[n]/wⱼ/(ωⱼ[i]-ωⱼ[i-1]) + 1/(1-μᵢ[n]^2)*γn[n]/wⱼ/(ωⱼ[i+1]-ωⱼ[i])
            end

        elseif m == n - 1 && j == i
            if n != 1
                ℳ[ii,jj] += 1/wᵢ[n] * β12[n]/(μᵢ[n]-μᵢ[n-1])
            end
            
        elseif m == n + 1 && j == i
            if n != Nᵢ
                ℳ[ii,jj] += 1/wᵢ[n] * β12[n+1]/(μᵢ[n+1]-μᵢ[n])
            end

        elseif m == n && (j == i - 1 || i == 1 && j == 2*Nᵢ)
            if i != 1
                ℳ[ii,jj] += 1/(1-μᵢ[n]^2)*γn[n]/wⱼ/(ωⱼ[i]-ωⱼ[i-1])
            else
                ℳ[ii,jj] += 1/(1-μᵢ[n]^2)*γn[n]/wⱼ/(ωⱼ[i]+2*π-ωⱼ[end])
            end

        elseif m == n && (j == i + 1 || i == 2*Nᵢ && j == 1)
            if i != 2*Nᵢ
                ℳ[ii,jj] += 1/(1-μᵢ[n]^2)*γn[n]/wⱼ/(ωⱼ[i+1]-ωⱼ[i])
            else
                ℳ[ii,jj] += 1/(1-μᵢ[n]^2)*γn[n]/wⱼ/(ωⱼ[1]+2*π-ωⱼ[i])
            end
        end
    end

end

return ℳ, λ₀, Mn_FP, Dn_FP, N_FP

end