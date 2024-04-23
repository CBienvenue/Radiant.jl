"""
    newton_bissection(f::Function,dfdx::Function,x⁻::Real,x⁺::Real)
⁻
Mixed Newton-bissection iterative root-finding method.

# Input Argument(s)
- 'f::Function': function f(x).
- 'dfdx::Function': derivative of f(x).
- 'x⁻::Real': x-value such as f(x⁻) < 0.
- 'x⁺::Real': x-value such as f(x⁺) > 0.

# Output Argument(s)
- 'x::Real': root of f(x).

# Author(s)
Charles Bienvenue

# Reference(s)

"""
function newton_bissection(f::Function,dfdx::Function,x⁻::Real,x⁺::Real,ϵ₀::Float64=eps(),imax::Int64=500)

# Verification of input paramters
if (f(x⁻) > 0 || f(x⁺) < 0) error("Badly choosen interval for bissection method.") end

# Choose the point closest to x = 0
if abs(f(x⁻)) < abs(f(x⁺)) x₀ = x⁻ else x₀ = x⁺ end

# Root-finding iteration process
x = 0.0
ϵ = Inf
i = 1
while ϵ > ϵ₀

    # Newton method
    x = x₀ - f(x₀)/dfdx(x₀)

    # If out of interval, then apply bissection method
    if ~(x⁻ ≤ x ≤ x⁺) x = (x⁻+x⁺)/2 end

    # Squeeze the boundaries around f(x) = 0
    if f(x) < 0 x⁻ = x else x⁺ = x end

    # Convergence criterion
    ϵ = abs(f(x))
    i += 1
    if i > imax || isnan(x) || isinf(abs(x)) error("Unable to converge Newton-bissection method.") end
    x₀ = x

end

return x
end