"""
    legendre_polynomials(L::Int64,x::Float64)

Calculate the Legendre polynomials Pℓ(x) values for ℓ=0,L.

# Input Argument(s)
- 'L::Int64': truncation order.
- 'x::Float64': evaluation points.

# Output Argument(s)
- 'Pℓ::Vector{Float64}': Legendre polynomials.

# Reference(s)
- Weisstein (2023) : Legendre Polynomial [MathWorld - A Wolfram Web Resource].

"""
function legendre_polynomials(L::Int64,x::Float64)

# Verification of input paramters
if L < 0 error("Legendre order is greater or equal to zero.") end
if abs(x) > 1 error("Invalid evaluation point.") end

# Calculations of Legendre polynomials using Bonnet's recursion formula
Pℓ = ones(L+1)
if L ≥ 1 Pℓ[2] = x end
for ℓ in range(2,L)
    Pℓ[ℓ+1] = ((2*ℓ-1)*x*Pℓ[ℓ] - (ℓ-1)*Pℓ[ℓ-1])/ℓ
end
return Pℓ
end