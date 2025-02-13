"""
    gauss_lobatto(N::Int64)

Computing the Gauss-Lobatto quadrature.

# Input Argument(s)
- `N::Int64`: quadrature order.

# Output Argument(s)
- `μ::Vector{Float64}`: vector with μᵢ points.
- `w::Vector{Float64}`: vector with wᵢ weights.

# Reference(s)
- Hale and Townsend (2013), Fast and accurate computation of Gauss-Legendre and
  Gauss-Jacobi quadrature nodes and weights.

"""
function gauss_lobatto(N::Int64)

# Initial guess
x = zeros(N)
x = cos.(π*(0:(N-1))/(N-1))

# Compute quadrature parameters
w = zeros(N)
dPdx = 0
d²Pdx² = 0
P = [0.0,0.0]
w[1] = 2/(N*(N-1))
w[end] = w[1]
for i in range(2,N-1)

    # Initial comparaison point outside the domain for x
    xold = 2

    # Compute the zeroth of the Legendre polynomial derivative of order N-1 at point xᵢ
    while abs(x[i]-xold) > eps()

        xold = x[i]

        # Legendre polynomial (Bonnet's recursion formula)
        P = [1,x[i]]
        for n in range(2,N-1)
           Pn = ((2*n-1)*x[i]*P[2]-(n-1)*P[1])/n
           P[1] = P[2]
           P[2] = Pn
        end

        # Legendre polynomial first and second derivative of order N
        dPdx = N*(P[1]-x[i]*P[2])/(1-x[i]^2)
        d²Pdx² = (2*x[i]*dPdx-N*(N+1)*P[2])/(1-x[i]^2)

        # Newton-Raphson method
        x[i] -= dPdx/d²Pdx²

    end

    # Compute the weight associated with point xᵢ
    w[i] = 2/(N*(N-1)*P[2]^2)
    
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

return [x,w]
end