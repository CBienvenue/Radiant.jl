"""
    jacobi_polynomials_up_to_L(L::Int64,α::Int64,β::Int64,x::Real)

Calculate the Jacobi polynomials up to order L.

# Input Argument(s)
- `L::Int64`: maximum Legendre order.
- `α::Int64`: alpha parameter.
- `β::Int64`: beta parameter.
- `x::Float64`: direction cosine.

# Output Argument(s)
- `Plαβ::Vector{Float64}`: Jacobi polynomials up to L, evaluated at μ and ϕ.

"""
function jacobi_polynomials_up_to_L(L::Int64,α::Int64,β::Int64,x::Real)
    Plαβ = zeros(L+1)

    # Endpoint cases
    if x == 1
        for l in 0:L
            Plαβ[l+1] = factorial_factor([l+α],[l,α])
        end
    elseif x == -1
        for l in 0:L
            Plαβ[l+1] = (-1)^l * factorial_factor([l+β],[l,β])
        end

    # General formula
    else
        Plαβ[1] = 1
        if (L ≥ 1) Plαβ[2] = (α+1) + 0.5*(α+β+2)*(x-1) end
        for l in 2:L

            # Recurrence formula
            if l < 125
                Plαβ[l+1] = ((2*l+α+β-1)*((2*l+α+β)*(2*l+α+β-2)*x+α^2-β^2)*Plαβ[l] - (2*(l+α-1)*(l+β-1)*(2*l+α+β))*Plαβ[l-1])/(2*l*(l+α+β)*(2*l+α+β-2))

            # Darboux asymptotic formula
            else
                θ = acos(x)
                kθ = 1/(sqrt(π)*sin(0.5*θ)^(α+0.5)*cos(0.5*θ)^(β+0.5)) 
                N = l + 0.5 * (α + β + 1)
                γ = -0.5*π * (α+0.5)
                Plαβ[l+1] = 1/sqrt(l) * kθ * cos(N*θ+γ)
            end
        end
    end
    return Plαβ
end