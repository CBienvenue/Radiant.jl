"""
    gauss_legendre_chebychev(N::Int64)

Computing the product of Gauss-Legendre and Gauss-Chebyshev quadratures on the unit sphere.

See also [`quadrature`](@ref).

# Input Argument(s)
- 'N::Int64': quadrature order.

# Output Argument(s)
- 'μ::Vector{Float64}': vector with μᵢ points.
- 'η::Vector{Float64}': vector with ηᵢ points.
- 'ξ::Vector{Float64}': vector with ξᵢ points.
- 'w::Vector{Float64}': vector with wᵢ weights.

# Author(s)
Charles Bienvenue

# Reference(s)
N/A

"""
function gauss_legendre_chebychev(N::Int64)

if (mod(N,2)!=0) error("Gauss-Chebychev product-type quadrature order in not even.") end

μ = zeros(2*N^2); η = zeros(2*N^2); ξ = zeros(2*N^2); w = zeros(2*N^2)

μn,wn = gauss_legendre(N)
index = sortperm(μn)
μn = μn[index]; wn = wn[index]

i = 1
for n in range(1,N)
    for m in range(1,2*N) 

        ϕ = π/N * (m-0.5)
        μ[i] = μn[n]
        η[i] = sqrt(1-μn[n]^2) * cos(ϕ)
        ξ[i] = sqrt(1-μn[n]^2) * sin(ϕ)
        w[i] = wn[n] * π/N

        i += 1

    end
end

return [μ,η,ξ],w

end