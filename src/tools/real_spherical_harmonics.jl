"""
    real_spherical_harmonics(ℓ::Int64,m::Int64,μ::Float64,ϕ::Float64)

Calculate the real spherical harmonics components Rℓm(μ,ϕ). 

See also [`ferrer_associated_legendre`](@ref).

# Input Argument(s)
- 'ℓ::Int64': order of the associated Legendre polynomial.
- 'm::Int64': degree of the associated Legendre polynomial.
- 'μ::Float64': direction cosine.
- 'ϕ::Float64': azimuthal angle.

# Output Argument(s)
- 'Rℓm::Float64': real spherical harmonics of order ℓ and degree m evaluated at μ and ϕ.

# Author(s)
Charles Bienvenue

# Reference(s)
- Hébert (2016): Applied Reactor Physics

"""
function real_spherical_harmonics(ℓ::Int64,m::Int64,μ::Float64,ϕ::Float64)

# Verification of input paramters
if ℓ < 0 error("Legendre order is greater or equal to zero.") end
if abs(m) > ℓ error("Legendre order smaller than m-order (|m| > ℓ).") end
if ~(-1 ≤ μ ≤ 1) error("Invalid direction cosine") end

# Compute the real spherical harmonics
if (m ≥ 0) 𝓣m = cos(m*ϕ) else 𝓣m = sin(abs(m)*ϕ) end
Pℓm = ferrer_associated_legendre(ℓ,abs(m),μ)
Cℓm = sqrt((2-(m == 0)) * exp( sum(log.(1:ℓ-abs(m))) - sum(log.(1:ℓ+abs(m))) ))
Rℓm = Cℓm * Pℓm * 𝓣m

return Rℓm
end