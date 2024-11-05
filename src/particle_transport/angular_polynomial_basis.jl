"""
    angular_polynomial_basis(Ndims::Int64,Ω::Union{Vector{Vector{Float64}},
    Vector{Float64}},w::Vector{Float64},L::Int64,N::Int64,type::String)

Compute the polynomial interpolation basis, based on a choice of quadrature and
discretization, and produce both discrete-to-moments and moments-to-discrete matrices.  

# Input Argument(s)
- 'Ndims::Int64': geometry dimension.
- 'Ω::Union{Vector{Vector{Float64}}': director cosines.
- 'w::Vector{Float64}': quadrature weights.
- 'L::Int64': legendre order.
- 'N::Int64': quadrature order.
- 'type::String': type of scattering source treatment.

# Output Argument(s)
- 'P::Int64': number of interpolation basis.
- 'Mn::Array{Float64}': discrete-to-moment matrix.
- 'Dn::Array{Float64}': moment-to-discrete matrix.
- 'pℓ::Vector{Int64}': legendre order associated with each interpolation basis.

# Reference(s)
- Morel (1988) : A Hybrid Collocation-Galerkin-Sn Method for Solving the Boltzmann
  Transport Equation.
- Morel (2017) : Comparison of Two Galerkin Quadrature Methods.
- Sanchez (2011) : On the Construction of Galerkin Angular Quadratures.
- Lewis (1984) : Computational Methods of Neutron Transport.
- Drumm (2011) : Least squares finite elements algorithms in the SCEPTRE radiation
  transport code.

"""
function angular_polynomial_basis(Ndims::Int64,Ω::Vector{Vector{Float64}},w::Vector{Float64},L::Int64,N::Int64,type::String,Qdims::Int64)

#----
# Compute Legendre or real spherical harmonics
#----

if Qdims == 1
    μ = Ω[1]
    Pℓ = zeros(N,L+1,1)
    @inbounds for n in range(1,N)
        Pℓ[n,:] = legendre_polynomials(L,μ[n])
    end
else
    μ = Ω[1]; η = Ω[2]; ξ = Ω[3];
    N = length(w)
    Rℓm = zeros(N,L+1,2*L+1)
    @inbounds for n in range(1,N)
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
        for ℓ in range(0,L), m in range(-ℓ,ℓ)
            Rℓm[n,ℓ+1,m+ℓ+1] = real_spherical_harmonics(ℓ,m,μ[n],ϕ)
        end
    end
end

#----
# Produce discrete-to-moments (D) and moments-to-discrete matrix (M)
#----

# 1D case
if Qdims == 1

    # Standard SN
    if type == "standard"
        P = L+1
        pℓ = zeros(Int64,P)
        pm = zeros(Int64,P)
        Mn = zeros(N,P)
        Dn = zeros(P,N)
        for ℓ in range(0,L)
            pℓ[ℓ+1] = ℓ
            for n in range(1,N)
                Mn[n,ℓ+1] = (2*ℓ+1)/2 * Pℓ[n,ℓ+1]
                Dn[ℓ+1,n] = w[n] * Pℓ[n,ℓ+1]
            end
        end

    # Galerkin method (inversion of M-matrix)
    elseif type == "galerkin-m"
        P = N
        pℓ = zeros(Int64,P)
        pm = zeros(Int64,P)
        Mn = zeros(N,N)
        for ℓ in range(0,N-1)
            pℓ[ℓ+1] = ℓ
            for n in range(1,N)
                Mn[n,ℓ+1] = (2*ℓ+1)/2 * Pℓ[n,ℓ+1]
            end
        end
        Dn = inv(Mn)

    # Galerkin method (inversion of D-matrix)
    elseif type == "galerkin-d"
        P = N
        pℓ = zeros(Int64,P)
        pm = zeros(Int64,P)
        Dn = zeros(N,N)
        for ℓ in range(0,N-1)
            pℓ[ℓ+1] = ℓ
            for n in range(1,N)
                Dn[ℓ+1,n] = w[n] * Pℓ[n,ℓ+1]
            end
        end
        Mn = inv(Dn)
    else
        error("Unknown method.")
    end

# 2D case
elseif Qdims == 2

    # Standard SN
    if type == "standard"
        P = div((L+1)*(L+2),2)
        pℓ = zeros(Int64,P)
        pm = zeros(Int64,P)
        Mn = zeros(N,P)
        Dn = zeros(P,N)
        ℓi = 1
        for ℓ in range(0,L), m in range(-ℓ,ℓ)
            if mod(ℓ+m,2) == 1 continue end
            pℓ[ℓi] = ℓ
            pm[ℓi] = m
            for n in range(1,N)
                Mn[n,ℓi] = (2*ℓ+1)/(4*π) * Rℓm[n,ℓ+1,m+ℓ+1]
                Dn[ℓi,n] = w[n] * Rℓm[n,ℓ+1,m+ℓ+1]
            end
            ℓi += 1
        end
    elseif type == "galerkin-d"
        P,Mn,Dn,pℓ,pm = angular_matrix_gram_schmidt(N,L,Rℓm,w,Qdims)
    elseif type == "galerkin-m"
        P,Mn,Dn,pℓ,pm = angular_matrix_gram_schmidt(N,L,Rℓm,w,Qdims,2)
    else
        error("Unknown method.")
    end

# 3D case
elseif Qdims == 3

    # Standard SN
    if type == "standard"
        P = (L+1)^2
        pℓ = zeros(Int64,P)
        pm = zeros(Int64,P)
        Mn = zeros(N,P)
        Dn = zeros(P,N)
        ℓi = 1
        for ℓ in range(0,L), m in range(-ℓ,ℓ)
            pℓ[ℓi] = ℓ
            pm[ℓi] = m
            for n in range(1,N)
                Mn[n,ℓi] = (2*ℓ+1)/(4*π) * Rℓm[n,ℓ+1,m+ℓ+1]
                Dn[ℓi,n] = w[n] * Rℓm[n,ℓ+1,m+ℓ+1]
            end
            ℓi += 1
        end
    elseif type == "galerkin-d"
        P,Mn,Dn,pℓ,pm = angular_matrix_gram_schmidt(N,L,Rℓm,w,Qdims)
    elseif type == "galerkin-m"
        P,Mn,Dn,pℓ,pm = angular_matrix_gram_schmidt(N,L,Rℓm,w,Qdims,2)
    else
        error("Unknown method.")
    end
else
    error("Dimension should be either 1, 2 or 3.")
end

return P,Mn,Dn,pℓ,pm
end

function angular_matrix_gram_schmidt(N,L,Rℓm,w,Qdims,g_type=1)

    norm(a) = sqrt(sum(a.*a))
    inner_product(a,b) = sum(a.*b)

    P = N
    pℓ = zeros(Int64,N)
    pm = zeros(Int64,N)
    Mn = zeros(N,N)
    Dn = zeros(N,N)

    # Gram-Schmidt procedure to find strongly independant
    # set of spherical harmonics

    # Initialisation
    pℓ[1] = 0
    pm[1] = 0
    u = zeros(N,N)
    for n in range(1,N)
        if g_type == 1
            Dn[1,n] = w[n] * Rℓm[n,1,1]
        else
            Mn[n,1] = 1/(4*π) * Rℓm[n,1,1]
        end
        u[n,1] = Rℓm[n,1,1]/norm(Rℓm[:,1,1])
    end
    i = 1
    # Iterate over the spherical harmonics
    for ℓ in range(1,L)
        m_min = -ℓ
        for m in range(m_min,ℓ)
            S = zeros(N)
            for n in range(1,i), n2 in range(1,N)
                S[n2] += inner_product(u[:,n],Rℓm[:,ℓ+1,m+ℓ+1]) * u[n2,n]
            end
            vk = zeros(N)
            for n in range(1,N)
                vk[n] = Rℓm[n,ℓ+1,m+ℓ+1] - S[n]
            end

            # If the set is strongly independant, keep it, otherwise go to next
            if norm(vk) > 1e-5
                i += 1 
                for n in range(1,N)
                    u[n,i] = vk[n]/norm(vk)
                    if g_type == 1
                        Dn[i,n] = w[n] * Rℓm[n,ℓ+1,m+ℓ+1]
                    else
                        Mn[n,i] = (2*ℓ+1)/(4*π) * Rℓm[n,ℓ+1,m+ℓ+1]
                    end
                end
                pℓ[i] = ℓ
                pm[i] = m
            end
            if (N == i) break end
            if (ℓ == L && m == L) error(string("The Gram-Schmidt procedure to find a suitable interpolation basis of spherical harmonics requires more of them (L should be > ",L,").")) end
        end
    end

    if g_type == 1
        Mn = inv(Dn)
    else
        Dn = inv(Mn)
    end

    return P,Mn,Dn,pℓ,pm
end