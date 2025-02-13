"""
    fokker_planck_galerkin(N::Int64,Mn::Array{Float64},Dn::Array{Float64},
    pℓ::Vector{Int64},P::Int64)

Compute the angular Fokker-Planck scattering matrix using Galerkin quadrature.

# Input Argument(s)
- `N::Int64`: number of directions.
- `Mn::Array{Float64}`: moment-to-discrete matrix.
- `Dn::Array{Float64}`: discrete-to-moment matrix.
- `pℓ::Vector{Int64}`: legendre order associated with each interpolation basis. 
- `P::Int64`: number of angular interpolation basis.

# Output Argument(s)
- `ℳ::Array{Float64}`: Fokker-Planck scattering matrix.
- `λ₀::Float64`: correction factor for total cross-section.

# Reference(s)
- Morel et al. (1988) : A Hybrid Collocation-Galerkin-Sn Method for Solving the Boltzmann
  Transport Equation.
- Landesman and Morel (1988) : Angular Fokker-Pianck Decomposition and Representation
  Techniques.

"""
function fokker_planck_galerkin(N::Int64,Mn::Array{Float64},Dn::Array{Float64},pℓ::Vector{Int64},P::Int64)

# Inversion of M-matrix
Σ = zeros(P,P)
for p in range(1,P)
    ℓ = pℓ[p]
    Σ[p,p] = -ℓ*(ℓ+1)
end
ℳ = zeros(N,N)
ℳ = Mn*Σ*Dn

# Making the diagonal positive
λ₀ = 0
for n in range(1,N)
    if (ℳ[n,n] < -λ₀) λ₀ = -ℳ[n,n] end
end
for n in range(1,N)
    ℳ[n,n] += λ₀
end

return ℳ, λ₀

end