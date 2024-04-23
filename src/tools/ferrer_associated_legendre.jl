"""
    ferrer_associated_legendre(ℓ::Int64,m::Int64,x::Float64)

Calculate the associated Legendre polynomial Pℓm(x), using Ferrer definition
(without (-1)ᵐ factor). 

See also [`real_spherical_harmonics`](@ref).

# Input Argument(s)
- 'ℓ::Int64': order of the associated Legendre polynomial.
- 'm::Int64': degree of the associated Legendre polynomial.
- 'x::Real': evaluation point.

# Output Argument(s)
- 'Pℓm::Float64': associated Legendre polynomials of order ℓ and degree m evaluated at
   point x.

# Author(s)
Charles Bienvenue

# Reference(s)
- Hébert (2016): Applied Reactor Physics

"""
function ferrer_associated_legendre(ℓ::Int64,m::Int64,x::Real)

# Verification of input paramters
if ℓ < 0 error("Legendre order is greater or equal to zero.") end
if m < 0 error("Negative degree m.") end
if m > ℓ error("Legendre order greater than m-order (m > ℓ).") end
if abs(x) > 1 error("Invalid evaluation point.") end

# Calculation of the associated Legendre polynomial (Ferrer definition)
Pℓm = 1
Pmm = 1; Pmm⁺ = 1; Pmm⁺⁺ = 1
if m > 0 Pmm = prod(1:2:2*m-1)*sqrt(1-x^2)^m end
if ℓ == m
    Pℓm = Pmm
else
    Pmm⁺ = (2*m+1) * x * Pmm
    if ℓ == m + 1
        Pℓm = Pmm⁺
    else
        for ℓi in range(m+2,ℓ)
            Pmm⁺⁺ = ((2*ℓi-1)*x*Pmm⁺ - (ℓi+m-1)*Pmm)/(ℓi-m)
            Pmm = Pmm⁺
            Pmm⁺ = Pmm⁺⁺
        end
        Pℓm = Pmm⁺⁺
    end
end

return Pℓm
end