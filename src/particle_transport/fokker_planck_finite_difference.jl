"""
    fokker_planck_finite_difference(N::Int64,quadrature_type::String,Ndims::Int64)

Calculate the Fokker-Planck scattering matrix using finite-difference scheme.

# Input Argument(s)
- `N::Int64`: quadrature order.
- `quadrature_type::String`: type of quadrature.
- `Ndims::Int64`: dimension of the geometry.
- `Nd::Int64` : number of directions.
- `Mn::Array{Float64}` : moment-to-discrete matrix.
- `Dn::Array{Float64}` : discrete-to-moment matrix.
- `Qdims::Int64` : quadrature dimension.

# Output Argument(s)
- `ℳ::Array{Float64}`: Fokker-Planck scattering matrix.
- `λ₀::Float64`: correction factor for total cross-section.

# Reference(s)
- Morel (1985) : An Improved Fokker-Planck Angular Differencing Scheme.
- Landesman and Morel (1988) : Angular Fokker-Pianck Decomposition and Representation
  Techniques.
- Larsen and Morel (2010) : Advances in Discrete-Ordinates Methodology.
- Morel et al. (2007) : A Discretization Scheme for the Three-Dimensional Angular
  Fokker-Planck Operator.

"""
function fokker_planck_finite_difference(N::Int64,quadrature_type::String,Ndims::Int64,Nd::Int64,Mn::Array{Float64},Dn::Array{Float64},Qdims::Int64)

λ₀ = 0.0

if Qdims == 1 

    # Extract quadrature informations
    μ,w = quadrature(N,quadrature_type,Ndims)

    # Reorder evaluation points and weights
    index = sortperm(μ)
    μn = μ[index]
    wn = w[index]

    # Compute the scattering matrix
    β12 = zeros(N+1)
    for n in range(2,N)
        β12[n] = β12[n-1] - 2*wn[n-1]*μn[n-1]
    end

    Cn⁻ = zeros(N)
    for n in range(2,N)
        Cn⁻[n] = β12[n]/(wn[n]*(μn[n]-μn[n-1]))
    end

    Cn⁺ = zeros(N)
    for n in range(1,N-1)
        Cn⁺[n] = β12[n+1]/(wn[n]*(μn[n+1]-μn[n]))
    end
    
    Cn = Cn⁻ + Cn⁺
    λ₀ = maximum(Cn)
    ℳ_temp = zeros(N,N)
    for n in range(1,N), m in range(1,N)
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
    for n in range(1,N), m in range(1,N)
        ℳ[index[n],index[m]] = ℳ_temp[n,m]
    end

elseif Qdims == 2

    if quadrature_type == "gauss-legendre-chebychev"

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
    else

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
                ϵ_ξ[j] = abs(Ω_3D[1][n_3D] - Ω[1][j])^2 + abs(Ω_3D[2][n_3D] - Ω[2][j])^2 + abs(sign(Ω_3D[3][n_3D]) * Ω_3D[3][n_3D] - Ω[3][j])^2
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
    end

    diag_ℳ = [abs(ℳ[i,i]) for i in range(1,Nd)]
    λ₀ = maximum(diag_ℳ)
    for i in range(1,Nd)
        ℳ[i,i] += λ₀
    end

elseif Qdims == 3

    if quadrature_type == "gauss-legendre-chebychev"

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
    else

        # Voronoi tesselation of the unit sphere and compute weighting parameters
        Ω,w = quadrature(N,quadrature_type,3)
        voronoi_data = voronoi_sphere(Ω)
        Nd = length(w)
        γ, edge_list = fokker_planck_weights_3D(Ω,w,voronoi_data)

        # Compute the scattering matrix
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
    end

    diag_ℳ = [abs(ℳ[i,i]) for i in range(1,Nd)]
    λ₀ = maximum(diag_ℳ)
    for i in range(1,Nd)
        ℳ[i,i] += λ₀
    end

end

return ℳ, λ₀
end