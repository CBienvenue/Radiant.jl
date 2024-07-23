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
λ₀ = 0.0

if Ndims == 1 

    #N = 6

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
    
    
    Σ = zeros(P,P)
    for p in range(1,P)
        ℓ = pℓ[p]
        Σ[p,p] = -ℓ*(ℓ+1)
    end
    ℳ2 = zeros(N,N)
    ℳ2 = Mn*Σ*Dn

    display(Mn*Dn_FP*ℳ*Mn_FP*Dn)
    display(ℳ2)

    error(maximum(abs.((ℳ2.-Mn*Dn_FP*ℳ*Mn_FP*Dn)./ℳ2)))
    =#

elseif Ndims == 2

    # 3D voronoi data
    Ω_3D,w_3D = quadrature(N,quadrature_type,3)
    voronoi_data = voronoi_sphere(Ω_3D)
    Nd_3D = length(w_3D)

    Ω,w = quadrature(N,quadrature_type,2)
    Nd = length(w)

    γ, edge_list = fokker_planck_weights_2D(Ω,w,Ω_3D,w_3D,voronoi_data)

    # 2D mapping
    map_3D_to_2D = Vector{Int64}(undef,Nd_3D)
    map_2D_to_3D = Vector{Int64}(undef,Nd)
    for n_3D in range(1,Nd_3D)
        ϵ_ξ = ones(Nd)*Inf
        for j in range(1,Nd)
            ϵ_ξ[j] = abs(Ω_3D[1][n_3D] - Ω[1][j])^2 + abs(Ω_3D[2][n_3D] - Ω[2][j])^2 + abs(Ω_3D[3][n_3D] - sign(Ω_3D[3][n_3D]) * Ω[3][j])^2
        end
        map_3D_to_2D[n_3D] = argmin(ϵ_ξ)
        if Ω_3D[3][n_3D] ≥ 0
            map_2D_to_3D[argmin(ϵ_ξ)] = n_3D
        end
    end

    # Build matrix
    ℳ = zeros(Nd,Nd)
    for n in range(1,Nd)
        voronoi_data_n = voronoi_data[map_2D_to_3D[n]]
        Nn_triangle = voronoi_data_n["Ni_triangle"]
        for j in range(1,Nn_triangle)
            voronoi_data_ij = voronoi_data_n["data_triangle_ij"][j]
            m_3D = voronoi_data_ij["index"]
            if Ω_3D[3][m_3D] < 0 continue end
            m = map_3D_to_2D[m_3D]
            if n == m continue end
            index_γ = findfirst(x -> x == (min(n,m),max(n,m)), edge_list)
            ℳ[n,m] += γ[index_γ]/w[n]
            ℳ[n,n] -= γ[index_γ]/w[n]
        end
    end
    
    #display(ℳ)

    ϕ = 10*rand(Nd)
    m0 = zeros(Nd)
    m1 = zeros(Nd)
    m2 = zeros(Nd)
    m3 = zeros(Nd)
    ℳn = zeros(Nd)
    for m in range(1,Nd)
        ℳn[m] = sum(ℳ[m,:].*ϕ)
    end
    m0 = sum(ℳn.*w)
    m1 = sum(ℳn.*Ω[1].*w) .+ 2 .* sum(Ω[1] .* w .*ϕ)
    m2 = sum(ℳn.*Ω[2].*w) .+ 2 .* sum(Ω[2] .* w .*ϕ)
    m3 = sum(ℳn.*Ω[3].*w) .+ 2 .* sum(Ω[3] .* w .*ϕ)

    #println([m0,m1,m2,m3])
    #error()
    #error([m0,m1,m2,m3])
    #=
    Rℓm = zeros(Nd,P)
    μ = Ω[1]; η = Ω[2]; ξ = Ω[3];
    @inbounds for n in range(1,Nd)
        if ξ[n] != 0
            if η[n] > 0
                ϕ = sign(ξ[n]) * atan(abs(ξ[n]/η[n]))
            elseif η[n] == 0
                ϕ = sign(ξ[n]) * π/2
            elseif η[n] < 0
                ϕ = sign(ξ[n]) * (π - atan(abs(ξ[n]/η[n])))
            end
        else
            if η[n] >= 0
                ϕ = 0.0
            elseif η[n] < 0
                ϕ = 1.0 * π
            end
        end
        for p in range(1,P)
            Rℓm[n,p] = real_spherical_harmonics(pℓ[p],pm[p],μ[n],ϕ)
        end
    end
    Mn_FP = zeros(Nd,P)
    Dn_FP = zeros(P,Nd)
    for p in range(1,P),n in range(1,Nd)
        Mn_FP[n,p] = (2*pℓ[p]+1)/(4*π) * Rℓm[n,p]
        Dn_FP[p,n] = w[n] * Rℓm[n,p]
    end
    Mn_FP = pinv(Dn_FP)
    =#
    #=
    ℳ = zeros(Nd,Nd)
    for n_3D in range(1,Nd_3D)
        if Ω_3D[3][n_3D] < 0 continue end
        voronoi_data_n = voronoi_data[n_3D]
        Nn_triangle = voronoi_data_n["Ni_triangle"]
        xn = voronoi_data_n["xi"]
        Sn = voronoi_data_n["Si"]
        n = map2D[n_3D]
        for j in range(1,Nn_triangle)
            voronoi_data_ij = voronoi_data_n["data_triangle_ij"][j]
            m_3D = voronoi_data_ij["index"]
            if Ω_3D[3][m_3D] < 0 continue end
            if abs(Ω_3D[3][n_3D]) < 1e-7 r = 0.5 else r = 1 end
            m = map2D[m_3D]
            Δx = voronoi_data_ij["Δxij"]
            ℓ = voronoi_data_ij["ℓij"]
            index_γ = findfirst(x -> x == (min(n_3D,m_3D),max(n_3D,m_3D)), edge_list)
            ℳ[n,m] += γ[index_γ]/(w[n]/2)
            ℳ[n,n] -= γ[index_γ]/(w[n]/2)
        end
    end

    ℳ_3D = zeros(Nd_3D,Nd_3D)
    for n_3D in range(1,Nd_3D)
        voronoi_data_n = voronoi_data[n_3D]
        Nn_triangle = voronoi_data_n["Ni_triangle"]
        for j in range(1,Nn_triangle)
            voronoi_data_ij = voronoi_data_n["data_triangle_ij"][j]
            m_3D = voronoi_data_ij["index"]
            index_γ = findfirst(x -> x == (min(n_3D,m_3D),max(n_3D,m_3D)), edge_list)
            ℳ_3D[n_3D,m_3D] += γ[index_γ]/w_3D[n_3D]
            ℳ_3D[n_3D,n_3D] -= γ[index_γ]/w_3D[n_3D]
        end
    end

    display(ℳ)
    display(ℳ_3D)

    ϕ = 10*rand(Nd_3D) #zeros(Nd_3D)
    #ϕ[1] = 1
    m0 = zeros(Nd_3D)
    m1 = zeros(Nd_3D)
    m2 = zeros(Nd_3D)
    m3 = zeros(Nd_3D)
    ℳn = zeros(Nd_3D)
    for m in range(1,Nd_3D)
        ℳn[m] = sum(ℳ_3D[m,:].*ϕ)
    end
    m0 = sum(ℳn.*w_3D)
    m1 = sum(ℳn.*Ω_3D[1].*w_3D) .+ 2 .* sum(Ω_3D[1] .* w_3D .*ϕ)
    m2 = sum(ℳn.*Ω_3D[2].*w_3D) .+ 2 .* sum(Ω_3D[2] .* w_3D .*ϕ)
    m3 = sum(ℳn.*Ω_3D[3].*w_3D) .+ 2 .* sum(Ω_3D[3] .* w_3D .*ϕ)
    println(m0)
    println(m1)
    println(m2)
    println(m3)
    =#
    #error()
    #=
    ϕ = 1:Nd
    m0 = zeros(Nd)
    m1 = zeros(Nd)
    m2 = zeros(Nd)
    m3 = zeros(Nd)
    for i in range(1,Nd)
        m0[i] = sum(ℳ[i,:].*w[i])
        m1[i] = sum(ℳ[i,:].*Ω[1][i].*w[i]/2 .*ϕ .+ 2 .* Ω[1][i] .* w[i]/2 .*ϕ[i])
    end

    println(m1)
    display(ℳ)
    error()
    =#
    #error()
    #for n in range(1,Nd)
    #    println(w[n],[Ω[1][n],Ω[2][n],Ω[3][n]])
    #end
    #println(Nd_3D)
    #println((count(!iszero, ℳ_3D)-Nd_3D)/2)
    #error()
    
    #=
    # Inversion of M-matrix
    Σ = zeros(P,P)
    for p in range(1,P)
        ℓ = pℓ[p]
        Σ[p,p] = -ℓ*(ℓ+1)
    end
    ℳ = zeros(N,N)
    ℳ = Mn*Σ*Dn

    display(ℳ)
    =#
    #=
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
    
    Mn_FP = pinv(Dn_FP)
    #display(Mn*Dn)
    #display(Dn_FP*Mn_FP)
    #display(Mn*Dn_FP*Mn_FP*Dn)
    #println(cond(Mn*Dn)," ",cond(Dn_FP*Mn_FP)," ",eigvals(Dn_FP*Mn_FP)," ",rank(Dn_FP*Mn_FP))
    #error()

    #Mn_FP = pinv(Dn_FP)
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

        elseif m == n && j == i - 1
            ℳ[ii,jj] += 1/(1-μᵢ[n]^2)*γn[n]/wⱼ/(ωⱼ[i]-ωⱼ[i-1])

        elseif m == n && j == i + 1
            ℳ[ii,jj] += 1/(1-μᵢ[n]^2)*γn[n]/wⱼ/(ωⱼ[i+1]-ωⱼ[i])
        end

    end
    
    display(Mn*Dn_FP*ℳ*Mn_FP*Dn)
    
    error()
    =#

elseif Ndims == 3

    Ω,w = quadrature(N,quadrature_type,3)
    voronoi_data = voronoi_sphere(Ω)
    Nd = length(w)
    γ, edge_list = fokker_planck_weights_3D(Ω,w,voronoi_data)

    # Build matrix
    ℳ = zeros(Nd,Nd)
    for n in range(1,Nd)
        voronoi_data_n = voronoi_data[n]
        Nn_triangle = voronoi_data_n["Ni_triangle"]
        for j in range(1,Nn_triangle)
            voronoi_data_ij = voronoi_data_n["data_triangle_ij"][j]
            m = voronoi_data_ij["index"]
            index_γ = findfirst(x -> x == (min(n,m),max(n,m)), edge_list)
            ℳ[n,m] += γ[index_γ]/w[n]
            ℳ[n,n] -= γ[index_γ]/w[n]
        end
    end

    #display(ℳ)

    ϕ = 10*rand(Nd)
    m0 = zeros(Nd)
    m1 = zeros(Nd)
    m2 = zeros(Nd)
    m3 = zeros(Nd)
    ℳn = zeros(Nd)
    for m in range(1,Nd)
        ℳn[m] = sum(ℳ[m,:].*ϕ)
    end
    m0 = sum(ℳn.*w)
    m1 = sum(ℳn.*Ω[1].*w) .+ 2 .* sum(Ω[1] .* w .*ϕ)
    m2 = sum(ℳn.*Ω[2].*w) .+ 2 .* sum(Ω[2] .* w .*ϕ)
    m3 = sum(ℳn.*Ω[3].*w) .+ 2 .* sum(Ω[3] .* w .*ϕ)

    #println([m0,m1,m2,m3])
    #error()

    #=
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
    Mn_FP = pinv(Dn_FP) # Moore-Penrose inverse
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
    =#

end

#println(maximum(eigvals(ℳ))," ",eigvals(ℳ))
#println(eigvals(Mn*Dn_FP*ℳ*Mn_FP*Dn))
#display(ℳ)
#display(Dn_FP*ℳ*Mn_FP)
#display(Mn*Dn_FP*ℳ*Mn_FP*Dn)
#println(eigvals(-ℳ))

return ℳ, λ₀, Mn_FP, Dn_FP, N_FP

end