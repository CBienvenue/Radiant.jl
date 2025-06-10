"""
    gauss_legendre(N::Int64)

Compute the Gauss-Legendre quadrature nodes and weight.

# Input Argument(s)
- `N::Int64`: quadrature order.

# Output Argument(s)
- `μ::Vector{Float64}`: vector with μᵢ points.
- `w::Vector{Float64}`: vector with wᵢ weights.

# Reference(s)
- Hale and Townsend (2013), Fast and accurate computation of Gauss-Legendre and
  Gauss-Jacobi quadrature nodes and weights.

"""
function gauss_legendre(N::Int64)

# Initial guess (Tricomi)
x = zeros(N)
for i in range(1,div(N,2))
    ϕ = π/2 * (4*i-1)/(2*N+1)
    x[i] = (1 - (N-1)/(8*N^3) - 1/(384*N^4)*(39-28/(sin(ϕ)^2))) * cos(ϕ)
    x[N-i+1] = -x[i]
end

# Compute quadrature parameters
w = zeros(N)
dPdx = 0
for i in range(1,N)

    # Initial comparaison point outside the domain for x
    xold = 2

    # Compute the zeroth of the Legendre polynomial of order N at point xᵢ
    while abs(x[i]-xold) > eps()

        xold = x[i]

        # Legendre polynomial (Bonnet's recursion formula)
        P = [1,x[i]]
        for n in range(2,N)
           Pn = ((2*n-1)*x[i]*P[2]-(n-1)*P[1])/n
           P[1] = P[2]
           P[2] = Pn
        end

        # Legendre polynomial derivative of order N
        dPdx = N*(P[1]-x[i]*P[2])/(1-x[i]^2)

        # Newton-Raphson method
        x[i] -= P[2]/dPdx

    end

    # Compute the weight associated with point xᵢ
    w[i] = 2/((1-x[i]^2)*dPdx^2)

end

# Verify if all roots are unique
seen = Vector{Float64}()
for i in x
    if i in seen
        error("Unable to build Lobatto quadrature at order N = ",N,".")
    else
        push!(seen,i)
    end
end

# Output
return [x,w]
end
