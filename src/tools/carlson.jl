"""
    carlson(N::Int64)

Computing the Carlson quadrature on the unit sphere.

# Input Argument(s)
- `N::Int64`: quadrature order.

# Output Argument(s)
- `μ::Vector{Float64}`: vector with μᵢ points.
- `η::Vector{Float64}`: vector with ηᵢ points.
- `ξ::Vector{Float64}`: vector with ξᵢ points.
- `w::Vector{Float64}`: vector with wᵢ weights.

# Reference(s)
- Carlson (1976), A Method of Characteristics and Other Improvements in Solution Methods
  for the Transport Equation.

"""
function carlson(N::Int64)

# Initialization
if N < 2 || (mod(N,2) != 0) error("Carlson quadrature available for even N ≥ 2.") end
μ = zeros(N*(N+2)); η = zeros(N*(N+2)); ξ = zeros(N*(N+2)); w = zeros(N*(N+2))

# Prepare weights and points
wn = zeros(div(N,2)); μn⁻ = zeros(div(N,2)); μn0 = zeros(div(N,2))
for n in range(1,div(N,2))
    wn[n] = 4*(N-2*n+2)/(N*(N+2))
    μn⁻[n] = 1 - (N-2*n+2)*(N-2*n+4)/(N*(N+2))
    μn0[n] = 1 - (N-2*n+2)^2/(N*(N+2))
end

if N == 2

    μn = zeros(1)
    μn[1] = 1/sqrt(3)

else
    # Newton/bissection method to compute x-values
    f(x) = 1/3 - sum(wn.*(μn0.+x.*μn⁻).^2)
    dfdx(x) = - sum(2 .*wn.*μn⁻.*(μn0.+x.*μn⁻))
    x = newton_bisection(f,dfdx,1,-1)

    # Compute μ-values
    μn = zeros(div(N,2))
    for n in range(1,div(N,2))
        μn[n] =  μn0[n] + x * μn⁻[n]
    end

end

# Compute A-values
w0 = 1/(N*(N+2))
if N == 2
    A = 1.0
else
    # Newton/bissection method to compute A-values
    function g(A)
        f = 0.0
        for n in range(1,div(N,2))
            f += π/2 * wn[n] * μn[n] - 4*π * w0/sqrt(2) * sin(π*A/4) * sqrt(1-μn[n]^2) / sin(π*A/(2*N-4*n+4))
        end
        return f
    end
    function dgdA(A)
        f = 0.0
        for n in range(1,div(N,2))
            f += -4*π*w0*sqrt(2)*π/8 * ( cos(π*A/4) * sqrt(1-μn[n]^2) / sin(π*A/(2*N-4*n+4)) - 2 * sin(π*A/4) * sqrt(1-μn[n]^2) * cos(π*A/(2*N-4*n+4)) / sin(π*A/(2*N-4*n+4))^2 / (N-2*n+2) )
        end
        return f
    end
    A = newton_bisection(g,dgdA,1/2,2)
end

# Compute the values for all octants
i = 1
s = [1 1 1; 1 1 -1; 1 -1 1; 1 -1 -1; -1 1 1; -1 1 -1; -1 -1 1; -1 -1 -1]

for t in range(1,8)
    for n in range(1,div(N,2))
        for m in range(1,div(N,2)-n+1)
            ϕ = π/2 * ((2*m-1)/(N-2*n+2)*A + (1-A)/2)
            μ[i] = s[t,1] * μn[n]
            η[i] = s[t,2] * sqrt(1-μ[i]^2) * cos(ϕ)
            ξ[i] = s[t,3] * sqrt(1-μ[i]^2) * sin(ϕ)
            w[i] = w0
            i += 1
        end
    end
end

return [μ,η,ξ] , 4*π*w

end