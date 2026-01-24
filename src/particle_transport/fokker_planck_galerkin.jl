"""
    fokker_planck_galerkin(N::Int64,Mn::Array{Float64},Dn::Array{Float64},
    pl::Vector{Int64},P::Int64)

Compute the angular Fokker-Planck scattering matrix using Galerkin quadrature.

# Input Argument(s)
- `N::Int64`: number of directions.
- `Mn::Array{Float64}`: moment-to-discrete matrix.
- `Dn::Array{Float64}`: discrete-to-moment matrix.
- `pl::Vector{Int64}`: legendre order associated with each interpolation basis. 
- `P::Int64`: number of angular interpolation basis.

# Output Argument(s)
- `ℳ_p::Array{Float64}`: Fokker-Planck scattering matrix.
- `λ₀::Float64`: correction factor for total cross-section.

# Reference(s)
- Morel et al. (1988) : A Hybrid Collocation-Galerkin-Sn Method for Solving the Boltzmann
  Transport Equation.
- Landesman and Morel (1988) : Angular Fokker-Pianck Decomposition and Representation
  Techniques.

"""
function fokker_planck_galerkin(N::Int64,Mn::Array{Float64},Dn::Array{Float64},pl::Vector{Int64},P::Int64)

    # Initialisation
    λ₀ = 0
    ℳ_p = fokker_planck_galerkin(P,pl)

    # Making the diagonal positive
    ℳ_n = zeros(N,N)
    ℳ_n = Mn*ℳ_p*Dn
    for n in range(1,N)
        if (ℳ_n[n,n] < -λ₀) λ₀ = -ℳ_n[n,n] end
    end
    for n in range(1,N)
        ℳ_n[n,n] += λ₀
    end
    ℳ_p = Dn*ℳ_n*Mn
    return ℳ_p, λ₀

end

function fokker_planck_galerkin(P::Int64,pl::Vector{Int64})
    # Angular Fokker-Planck moments
    ℳ_p = zeros(P,P)
    for p in range(1,P)
        l = pl[p]
        ℳ_p[p,p] = -l*(l+1)
    end
    return ℳ_p
end