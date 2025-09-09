"""
    ferrer_associated_legendre(l::Int64,m::Int64,x::Float64)

Calculate the associated Legendre polynomial Plm(x), using Ferrer definition
(without (-1)ᵐ factor). 

# Input Argument(s)
- `l::Int64`: order of the associated Legendre polynomial.
- `m::Int64`: degree of the associated Legendre polynomial.
- `x::Real`: evaluation point.

# Output Argument(s)
- `Plm::Float64`: associated Legendre polynomials of order l and degree m evaluated at
   point x.

# Reference(s)
- Hébert (2016), Applied Reactor Physics.

"""
function ferrer_associated_legendre(l::Int64,m::Int64,x::Real)

# Verification of input paramters
if l < 0 error("Legendre order is greater or equal to zero.") end
if m < 0 error("Negative degree m.") end
if m > l error("Legendre order greater than m-order (m > l).") end
if abs(x) > 1 error("Invalid evaluation point.") end

# Calculation of the associated Legendre polynomial (Ferrer definition)
Plm = 1
Pmm = 1; Pmm⁺ = 1; Pmm⁺⁺ = 1
if m > 0 Pmm = prod(1:2:2*m-1)*sqrt(1-x^2)^m end
if l == m
    Plm = Pmm
else
    Pmm⁺ = (2*m+1) * x * Pmm
    if l == m + 1
        Plm = Pmm⁺
    else
        for li in range(m+2,l)
            Pmm⁺⁺ = ((2*li-1)*x*Pmm⁺ - (li+m-1)*Pmm)/(li-m)
            Pmm = Pmm⁺
            Pmm⁺ = Pmm⁺⁺
        end
        Plm = Pmm⁺⁺
    end
end

return Plm
end