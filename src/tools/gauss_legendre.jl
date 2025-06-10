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

"""
    gauss_legendre(N::Int64,f::Function,a::Real,b::Real)

Compute the Gauss-Legendre quadrature of f(x) on x ∈ [a,b].

# Input Argument(s)
- `N::Int64`: quadrature order.
- `f::Function` : function to integrate.
- `a::Real` : lower bound of integration.
- `b:Real` : upper bound of integration.

# Output Argument(s)
- `F::Float64`: quadrature of f(x) over x ∈ [a,b].

"""
function gauss_legendre(N::Int64,f::Function,a::Real,b::Real)
    u,w = gauss_legendre(N)
    Δx = b-a
    xm = (a+b)/2
    F = 0
    for n in range(1,N)
        x = u[n]*Δx/2 + xm
        F += w[n] * Δx/2 * f(x)
    end
    return F
end

"""
    gauss_lobatto(N::Int64,f::Function,a::Real,b::Real)

Compute the Gauss-Lobatto quadrature of f(x) on x ∈ [a,b].

# Input Argument(s)
- `N::Int64`: quadrature order.
- `f::Function` : function to integrate.
- `a::Real` : lower bound of integration.
- `b:Real` : upper bound of integration.

# Output Argument(s)
- `F::Float64`: quadrature of f(x) over x ∈ [a,b].

"""
function gauss_lobatto(N::Int64,f::Function,a::Real,b::Real)
    u,w = gauss_lobatto(N)
    Δx = b-a
    xm = (a+b)/2
    F = 0
    for n in range(1,N)
        x = u[n]*Δx/2 + xm
        F += w[n] * Δx/2 * f(x)
    end
    return F
end

"""
    adaptive_gauss_legendre(N::Int64,f::Function,a::Real,b::Real,tol::Real=1e-8,depth::Int64=0,max_depth::Int64=10)

Compute the Gauss-Legendre quadrature of f(x) adaptively on x ∈ [a,b].

# Input Argument(s)
- `N::Int64`: quadrature order.
- `f::Function` : function to integrate.
- `a::Real` : lower bound of integration.
- `b:Real` : upper bound of integration.
- `tol::Real` : tolerance for the quadrature integration.
- `depth::Int64` : actual depth of the adaptive calculations.
- `max_depth` : maximum depth of the adaptive calculations.

# Output Argument(s)
- `F::Float64`: quadrature of f(x) over x ∈ [a,b].

"""
function adaptive_quadrature(N::Int64,f::Function,a::Real,b::Real;tol::Real=1e-8,depth::Int64=0,max_depth::Int64=50)
    if depth ≥ max_depth error("Unable to converge the quadrature at the given tolerance.") end
    xm = (a+b)/2
    I₁ = gauss_legendre(N,f,a,b)
    I₂ = gauss_legendre(N,f,a,xm) + gauss_legendre(N,f,xm,b)
    ϵ = abs(I₁-I₂)
    if ϵ < tol * max(abs(I₂), 1.0)
        return I₂
    else
        return adaptive_quadrature(N,f,a,xm;tol=tol,depth=depth+1,max_depth=max_depth) + adaptive_quadrature(N,f,xm,b;tol=tol,depth=depth+1,max_depth=max_depth)
    end
end
