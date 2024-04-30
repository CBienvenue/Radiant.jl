"""
    gamma(z::Number)

Compute the gamma function Γ(x).

# Input Argument(s)
- 'z::Number': input value of the function.

# Output Argument(s)
- 'Γ::Number': gamma function.

# Reference(s)
- Wikipedia contributors (2023) : Gamma function.

"""
function gamma(z::Number)

# Weierstrass definition of the Gamma function (except for non-positive integer)
if ~(imag(z) == 0 && round(-real(z)) == abs(real(z)))
    
    γ = 0.5772156649 # Euler–Mascheroni constant
    P = 1/(1+z) * exp(z)
    ϵ = Inf
    n = 2
    while ϵ > 1e-3
        P⁻ = copy(P)
        P *= 1/(1+z/n) * exp(z/n)
        if imag(P) != 0
            ϵ = max(abs((real(P)-real(P⁻))/real(P)),abs((imag(P)-imag(P⁻))/imag(P)))
        else
            ϵ = abs((real(P)-real(P⁻))/real(P))
        end
        n += 1
        if n >= 10000 || isnan(P) || isinf(P) error("Unable to compute the Gamma function.") end
    end

    Γ = P  * exp(-γ*z)/z
else
    error("Gamma function for non-positive integer not implemented yet.")
end

return Γ
end