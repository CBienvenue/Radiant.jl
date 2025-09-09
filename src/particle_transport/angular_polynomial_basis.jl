"""
    angular_polynomial_basis(Ndims::Int64,Ω::Union{Vector{Vector{Float64}},
    Vector{Float64}},w::Vector{Float64},L::Int64,N::Int64,type::String)

Compute the polynomial interpolation basis, based on a choice of quadrature and
discretization, and produce both discrete-to-moments and moments-to-discrete matrices.  

# Input Argument(s)
- `Ω::Union{Vector{Vector{Float64}}` : director cosines.
- `w::Vector{Float64}` : quadrature weights.
- `L::Int64` : legendre order.
- `type::String` : type of scattering source treatment.
- `Qdims::Int64` : quadrature dimension.

# Output Argument(s)
- `Np::Int64` : number of interpolation basis.
- `Mn::Array{Float64}` : discrete-to-moment matrix.
- `Dn::Array{Float64}` : moment-to-discrete matrix.
- `pl::Vector{Int64}` : legendre order associated with each interpolation basis.
- `pm::Vector{Int64}` : spherical harmonics order associated with each interpolation basis.

# Reference(s)
- Lewis (1984), Computational Methods of Neutron Transport.
- Morel (1988), A Hybrid Collocation-Galerkin-Sn Method for Solving the Boltzmann
  Transport Equation.
- Sanchez and Ragusa (2011), On the Construction of Galerkin Angular Quadratures.
- Drumm et al. (2011), Least squares finite elements algorithms in the SCEPTRE radiation
  transport code.
- Morel et al. (2017), Comparison of Two Galerkin Quadrature Methods.
- Shands et al. (2024), A New Galerkin Quadrature Method Not Requiring a Matrix Inverse.

"""
function angular_polynomial_basis(Ω::Vector{Vector{Float64}},w::Vector{Float64},L::Int64,type::String,Qdims::Int64)

    #----
    # Compute Legendre or real spherical harmonics
    #----

    # Legendre polynomials
    if Qdims == 1
        μ = Ω[1]
        Nd = length(μ)
        Pl = [zeros(L+1) for n in range(1,Nd)]
        for n in range(1,Nd)
            Pl[n] = legendre_polynomials_up_to_L(L,μ[n])
        end
    # Real spherical harmonics
    else
        μ = Ω[1]; η = Ω[2]; ξ = Ω[3]; ϕ = atan.(ξ,η)
        Nd = length(w)
        Rlm = [[zeros(2*l+1) for l in 0:L] for n in range(1,Nd)]
        for n in range(1,Nd)
            Rlm[n] = real_spherical_harmonics_up_to_L(L,μ[n],ϕ[n])
        end
    end

    #----
    # Produce discrete-to-moments (D) and moments-to-discrete matrix (M)
    #----

    # 1D case
    if Qdims == 1

        # Standard SN
        if type == "standard"
            Np = L+1
            pl = zeros(Int64,Np)
            pm = zeros(Int64,Np)
            Mn = zeros(Nd,Np)
            Dn = zeros(Np,Nd)
            for l in range(0,L)
                pl[l+1] = l
                for n in range(1,Nd)
                    Mn[n,l+1] = (2*l+1)/2 * Pl[n][l+1]
                    Dn[l+1,n] = w[n] * Pl[n][l+1]
                end
            end

        # Galerkin method (inversion of M-matrix)
        elseif type == "galerkin-m"
            Np = Nd
            pl = zeros(Int64,Np)
            pm = zeros(Int64,Np)
            Mn = zeros(Nd,Nd)
            for l in range(0,Nd-1)
                pl[l+1] = l
                for n in range(1,Nd)
                    Mn[n,l+1] = (2*l+1)/2 * Pl[n][l+1]
                end
            end
            Dn = inv(Mn)

        # Galerkin method (inversion of D-matrix)
        elseif type == "galerkin-d"
            Np = Nd
            pl = zeros(Int64,Np)
            pm = zeros(Int64,Np)
            Dn = zeros(Nd,Nd)
            for l in range(0,Nd-1)
                pl[l+1] = l
                for n in range(1,Nd)
                    Dn[l+1,n] = w[n] * Pl[n][l+1]
                end
            end
            Mn = inv(Dn)
        else
            error("Unknown method.")
        end

    # 2D and 3D case
    elseif Qdims <= 3

        # Standard SN
        if type == "standard"
            if Qdims == 2
                Np = div((L+1)*(L+2),2)
            else
                Np = (L+1)^2
            end
            pl = zeros(Int64,Np)
            pm = zeros(Int64,Np)
            Mn = zeros(Nd,Np)
            Dn = zeros(Np,Nd)
            li = 1
            for l in range(0,L), m in range(-l,l)
                if Qdims == 2 && mod(l+m,2) == 1 continue end
                pl[li] = l
                pm[li] = m
                for n in range(1,Nd)
                    Mn[n,li] = (2*l+1)/(4*π) * Rlm[n][l+1][m+l+1]
                    Dn[li,n] = w[n] * Rlm[n][l+1][m+l+1]
                end
                li += 1
            end
        elseif type == "galerkin-d"
            Np,Mn,Dn,pl,pm = angular_matrix_gram_schmidt(Nd,L,Rlm,w,Qdims,1)
        elseif type == "galerkin-m"
            Np,Mn,Dn,pl,pm = angular_matrix_gram_schmidt(Nd,L,Rlm,w,Qdims,2)
        else
            error("Unknown method.")
        end
    else
        error("Domain dimension should be either 1, 2 or 3.")
    end
    return Np,Mn,Dn,pl,pm
end

"""
    surface_angular_polynomial_basis(Ω::Vector{Vector{Float64}},w::Vector{Float64},
    L::Int64,type::String,Qdims::Int64,Ndims::Int64,geo_type::String)

Compute the polynomial interpolation basis for each surface of the geometry, based on a
choice of quadrature and discretization, and produce both discrete-to-moments and 
moments-to-discrete matrices.

# Input Argument(s)
- `Ω::Union{Vector{Vector{Float64}}` : director cosines.
- `w::Vector{Float64}` : quadrature weights.
- `L::Int64` : legendre order.
- `type::String` : type of scattering source treatment.
- `Qdims::Int64` : quadrature dimension.
- `Ndims::Int64` : geometry dimension.
- `geo_type::String` : geometry type.

# Output Argument(s)
- `Np::Int64` : number of interpolation basis.
- `Mn::Array{Float64}` : discrete-to-moment matrix.
- `Dn::Array{Float64}` : moment-to-discrete matrix.
- `n⁺_to_n::Vector{Vector{Int64}}` : mapping from half-range to full-range indices.
- `n_to_n⁺::Vector{Vector{Int64}}` : mapping from full-range to half-range indices.
- `pl::Vector{Int64}` : legendre order associated with each interpolation basis.
- `pm::Vector{Int64}` : spherical harmonics order associated with each interpolation basis.

"""
function surface_angular_polynomial_basis(Ω::Vector{Vector{Float64}},w::Vector{Float64},L::Int64,type::String,Qdims::Int64,Ndims::Int64,geo_type::String)
    if geo_type == "cartesian"
        if Ndims == 1
            surfaces = ["X-","X+"]
        elseif Ndims == 2
            surfaces = ["X-","X+","Y-","Y+"]
        elseif Ndims == 3
            surfaces = ["X-","X+","Y-","Y+","Z-","Z+"]
        else
            error("Geometry dimension should be either 1, 2 or 3.")
        end
        Np = 0
        Mn = Vector{Array{Float64}}()
        Dn = Vector{Array{Float64}}()
        n⁺_to_n = Vector{Vector{Int64}}()
        n_to_n⁺ = Vector{Vector{Int64}}()
        pl = Vector{Vector{Int64}}()
        pm = Vector{Vector{Int64}}()
        for surf_i in surfaces
            Np,Mn_i,Dn_i,n⁺_to_n_i,n_to_n⁺_i,pl_i,pm_i = surface_angular_polynomial_basis(Ω,w,L,type,Qdims,surf_i)
            push!(Mn,Mn_i); push!(Dn,Dn_i); push!(n⁺_to_n,n⁺_to_n_i); push!(n_to_n⁺,n_to_n⁺_i); push!(pl,pl_i); push!(pm,pm_i)
        end
        return Np,Mn,Dn,n⁺_to_n,n_to_n⁺,pl,pm
    else
        error("Not implemented yet.")
    end
end

"""
    surface_angular_polynomial_basis(Ω::Vector{Vector{Float64}},w::Vector{Float64},L::Int64,
    type::String,Qdims::Int64,surface::String)

Compute the polynomial interpolation basis for a given surface, based on a
choice of quadrature and discretization, and produce both discrete-to-moments and 
moments-to-discrete matrices.

# Input Argument(s)
- `Ω::Union{Vector{Vector{Float64}}` : director cosines.
- `w::Vector{Float64}` : quadrature weights.
- `L::Int64` : legendre order.
- `type::String` : type of scattering source treatment.
- `Qdims::Int64` : quadrature dimension.
- `surface::String` : surface identifier.

# Output Argument(s)
- `Np::Int64` : number of interpolation basis.
- `Mn::Array{Float64}` : discrete-to-moment matrix.
- `Dn::Array{Float64}` : moment-to-discrete matrix.
- `n⁺_to_n::Vector{Vector{Int64}}` : mapping from half-range to full-range indices.
- `n_to_n⁺::Vector{Vector{Int64}}` : mapping from full-range to half-range indices.
- `pl::Vector{Int64}` : legendre order associated with each interpolation basis.
- `pm::Vector{Int64}` : spherical harmonics order associated with each interpolation basis.

"""
function surface_angular_polynomial_basis(Ω::Vector{Vector{Float64}},w::Vector{Float64},L::Int64,type::String,Qdims::Int64,surface::String)

    #----
    # Compute half-range Legendre polynomials or half-range real spherical harmonics
    #----

    # Half-range Legendre polynomials
    if Qdims == 1
        μ = Ω[1]
        if surface == "X-"
            n⁺_to_n = findall(x -> x > 0, μ)
            μ⁺ = μ[n⁺_to_n]
        elseif surface == "X+"
            n⁺_to_n = findall(x -> x < 0, μ)
            μ⁺ = -μ[n⁺_to_n]
        else
            error("Unknown surface for 1D angular domain.")
        end
        Nd⁺ = length(n⁺_to_n)
        w⁺ = w[n⁺_to_n]
        Pl = [zeros(L+1) for n⁺ in range(1,Nd⁺)]
        for n⁺ in range(1,Nd⁺)
            Pl[n⁺] = half_range_legendre_polynomials_up_to_L(L,μ⁺[n⁺])
        end

    # Half-range real spherical harmonics
    else
        μ = Ω[1]
        η = Ω[2]
        ξ = Ω[3]
        if surface == "X-"
            n⁺_to_n = findall(x -> x > 0, μ)
            μ⁺ = μ[n⁺_to_n]
            ϕ = atan.(ξ,η)
        elseif surface == "X+"
            n⁺_to_n = findall(x -> x < 0, μ)
            μ⁺ = -μ[n⁺_to_n]
            ϕ = atan.(ξ,η)
        elseif surface == "Y-"
            n⁺_to_n = findall(x -> x > 0, η)
            μ⁺ = η[n⁺_to_n]
            ϕ = atan.(μ,ξ)
        elseif surface == "Y+"
            n⁺_to_n = findall(x -> x < 0, η)
            μ⁺ = -η[n⁺_to_n]
            ϕ = atan.(μ,ξ)
        elseif surface == "Z-"
            n⁺_to_n = findall(x -> x > 0, ξ)
            μ⁺ = ξ[n⁺_to_n]
            ϕ = atan.(η,μ)
        elseif surface == "Z+"
            n⁺_to_n = findall(x -> x < 0, ξ)
            μ⁺ = -ξ[n⁺_to_n]
            ϕ = atan.(η,μ)
        else
            error("Unknown surface for 3D angular domain.")
        end
        Nd⁺ = length(n⁺_to_n)
        ϕ⁺ = ϕ[n⁺_to_n]
        w⁺ = w[n⁺_to_n]
        ψlm = [[zeros(2*l+1) for l in 0:L] for n⁺ in range(1,Nd⁺)]
        for n⁺ in range(1,Nd⁺)
            ψlm[n⁺] = real_half_range_spherical_harmonics_up_to_L(L,abs(μ⁺[n⁺]),ϕ⁺[n⁺])
        end
    end

    # Index mapping from n⁺ to n
    Nd = length(μ)
    n_to_n⁺ = zeros(Int, Nd)
    for (n⁺, n) in enumerate(n⁺_to_n)
        n_to_n⁺[n] = n⁺
    end

    #----
    # Produce discrete-to-moments (D) and moments-to-discrete matrix (M)
    #----

    # 1D case
    if Qdims == 1

        # Standard SN
        if type == "standard"
            Np = L+1
            pl = zeros(Int64,Np)
            pm = zeros(Int64,Np)
            Mn = zeros(Nd⁺,Np)
            Dn = zeros(Np,Nd⁺)
            for l in range(0,L)
                pl[l+1] = l
                for n⁺ in range(1,Nd⁺)
                    Mn[n⁺,l+1] = Pl[n⁺][l+1]
                    Dn[l+1,n⁺] = w⁺[n⁺] * Pl[n⁺][l+1]
                end
            end

        # Galerkin method (inversion of M-matrix)
        elseif type == "galerkin-m"
            Np = Nd⁺
            pl = zeros(Int64,Np)
            pm = zeros(Int64,Np)
            Mn = zeros(Nd⁺,Nd⁺)
            for l in range(0,Nd⁺-1)
                pl[l+1] = l
                for n⁺ in range(1,Nd⁺)
                    Mn[n⁺,l+1] = Pl[n⁺][l+1]
                end
            end
            Dn = inv(Mn)

        # Galerkin method (inversion of D-matrix)
        elseif type == "galerkin-d"
            Np = Nd⁺
            pl = zeros(Int64,Np)
            pm = zeros(Int64,Np)
            Dn = zeros(Nd⁺,Nd⁺)
            for l in range(0,Nd⁺-1)
                pl[l+1] = l
                for n⁺ in range(1,Nd⁺)
                    Dn[l+1,n⁺] = w⁺[n⁺] * Pl[n⁺][l+1]
                end
            end
            Mn = inv(Dn)
        else
            error("Unknown method.")
        end

    # 3D case
    elseif Qdims <= 3

        # Standard SN
        if type == "standard"
            if Qdims == 2
                Np = div((L+1)*(L+2),2)
            else
                Np = (L+1)^2
            end
            pl = zeros(Int64,Np)
            pm = zeros(Int64,Np)
            Mn = zeros(Nd⁺,Np)
            Dn = zeros(Np,Nd⁺)
            li = 1
            for l in range(0,L), m in range(-l,l)
                if Qdims == 2 && mod(l+m,2) == 1 continue end
                pl[li] = l
                pm[li] = m
                for n⁺ in range(1,Nd⁺)
                    Mn[n⁺,li] = ψlm[n⁺][l+1][m+l+1]
                    Dn[li,n⁺] = w[n⁺] * ψlm[n⁺][l+1][m+l+1]
                end
                li += 1
            end
        elseif type == "galerkin-d"
            Np,Mn,Dn,pl,pm = angular_matrix_gram_schmidt(Nd⁺,L,ψlm,w,Qdims,1)
        elseif type == "galerkin-m"
            Np,Mn,Dn,pl,pm = angular_matrix_gram_schmidt(Nd⁺,L,ψlm,w,Qdims,3)
        else
            error("Unknown method.")
        end
    else
        error("Domain dimension should be either 1 or 3.")
    end

    return Np,Mn,Dn,n⁺_to_n,n_to_n⁺,pl,pm
end

"""
    angular_matrix_gram_schmidt(Nd::Int64,L::Int64,Rlm::Array{Float64},w::Vector{Float64},
    Qdims::Int64,g_type::Int64=1)

Choose a suitable interpolation basis based on a Gram-Schmidt process and produce both
discrete-to-moments and moments-to-discrete matrices.

# Input Argument(s)
- `Nd::Int64` : number of directions in the quadrature set.
- `L::Int64` : legendre truncation order.
- `Rlm::Array{Float64}` : spherical harmonics.
- `w::Vector{Float64}` : quadrature weights.
- `Qdims::Int64` : quadrature dimension.
- `g_type::Int64` : type of Galerkin method.

# Output Argument(s)
- `Np::Int64` : number of interpolation basis.
- `Mn::Array{Float64}` : discrete-to-moment matrix.
- `Dn::Array{Float64}` : moment-to-discrete matrix.
- `pl::Vector{Int64}` : legendre order associated with each interpolation basis.
- `pm::Vector{INt64}` : spherical harmonics order associated with each interpolation basis.

# Reference(s)
- Drumm et al. (2011) : Least squares finite elements algorithms in the SCEPTRE radiation
  transport code.

"""
function angular_matrix_gram_schmidt(Nd::Int64,L::Int64,Rlm::Vector{Vector{Vector{Float64}}},w::Vector{Float64},Qdims::Int64,g_type::Int64=1)

    norm(a) = sqrt(sum(a.*a))
    inner_product(a,b) = sum(a.*b)

    Np = Nd
    pl = zeros(Int64,Nd)
    pm = zeros(Int64,Nd)
    Mn = zeros(Nd,Nd)
    Dn = zeros(Nd,Nd)

    # Gram-Schmidt procedure to find strongly independant set of spherical harmonics

    # Initialisation
    pl[1] = 0
    pm[1] = 0
    u = zeros(Nd,Nd)
    for n in range(1,Nd)
        if g_type == 1
            Dn[1,n] = w[n] * Rlm[n][1][1]
        elseif g_type == 2
            Mn[n,1] = 1/(4*π) * Rlm[n][1][1]
        elseif g_type == 3
            Mn[n,1] = Rlm[n][1][1]
        end
        vec_Rlm = zeros(Nd)
        for n2 in range(1,Nd)
            vec_Rlm[n2] = Rlm[n2][1][1]
        end
        u[n,1] = Rlm[n][1][1]/norm(vec_Rlm)
    end
    i = 1
    # Iterate over the spherical harmonics
    for l in range(1,L)
        if Qdims == 2 m_min = 0 else m_min = -l end
        for m in range(m_min,l)
            S = zeros(Nd)
            vec_Rlm = zeros(Nd)
            for n2 in range(1,Nd)
                vec_Rlm[n2] = Rlm[n2][l+1][m+l+1]
            end
            for n in range(1,i), n2 in range(1,Nd)
                S[n2] += inner_product(u[:,n],vec_Rlm) * u[n2,n]
            end
            vk = zeros(Nd)
            for n in range(1,Nd)
                vk[n] = Rlm[n][l+1][m+l+1] - S[n]
            end

            # If the set is strongly independant, keep it, otherwise go to next
            if norm(vk) > 1e-5
                i += 1 
                for n in range(1,Nd)
                    u[n,i] = vk[n]/norm(vk)
                    if g_type == 1
                        Dn[i,n] = w[n] * Rlm[n][l+1][m+l+1]
                    elseif g_type == 2
                        Mn[n,i] = (2*l+1)/(4*π) * Rlm[n][l+1][m+l+1]
                    elseif g_type == 3
                        Mn[n,i] = Rlm[n][l+1][m+l+1]
                    end
                end
                pl[i] = l
                pm[i] = m
            end
            if (Nd == i) break end
            if (l == L && m == L) error(string("The Gram-Schmidt procedure to find a suitable interpolation basis of spherical harmonics requires more of them (L should be > ",L,").")) end
        end
    end

    if g_type == 1
        Mn = inv(Dn)
    else
        Dn = inv(Mn)
    end

    return Np,Mn,Dn,pl,pm
end