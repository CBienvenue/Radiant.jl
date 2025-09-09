"""
    legendre_polynomials(L::Int64,x::Float64)

Calculate the Legendre polynomials Pl(x) values for l=0,L.

# Input Argument(s)
- `L::Int64`: truncation order.
- `x::Float64`: evaluation points.

# Output Argument(s)
- `Pl::Vector{Float64}`: Legendre polynomials.

# Reference(s)
- Weisstein (2023), Legendre Polynomial [MathWorld - A Wolfram Web Resource].

"""
function legendre_polynomials(L::Int64,x::Float64)

    # Verification of input paramters
    if L < 0 error("Legendre order is greater or equal to zero.") end
    if abs(x) > 1 error("Invalid evaluation point.") end

    # Calculations of Legendre polynomials using Bonnet's recursion formula
    if l == 0
        return 1.0
    elseif l == 1
        return x
    else
        p0 = 1.0
        p1 = x
        for n in 2:L
            p2 = ((2n-1)*x*p1 - (n-1)*p0)/n
            p0,p1 = p1,p2
        end
        return p1
    end
end

"""
    legendre_polynomials_up_to_L(L::Int64,x::Float64)

Calculate the Legendre polynomials up to order L.

# Input Argument(s)
- `L::Int64`: maximum Legendre order.
- `x::Float64`: direction cosine.

# Output Argument(s)
- `Pl::Vector{Float64}`: Legendre polynomials up to L, evaluated at μ and ϕ.

"""
function legendre_polynomials_up_to_L(L::Int64,x::Float64)
    if L < 0 error("Legendre order is greater or equal to zero.") end
    if ~(-1 ≤ x ≤ 1) error("Invalid x value (should be between -1 and 1).") end
    Pl = zeros(L+1)

    # Endpoint cases
    if x == 1
        for l in 0:L
            Pl[l+1] = 1
        end
    elseif x == -1
        for l in 0:L
            Pl[l+1] = (-1)^l
        end

    # General formula
    else
        Pl[1] = 1
        if L ≥ 1 Pl[2] = x end
        for l in 2:L
            Pl[l+1] = ((2*l-1)*x*Pl[l] - (l-1)*Pl[l-1])/l
        end
    end
    return Pl
end

"""
    half_range_legendre_polynomials_up_to_L(L::Int64,x::Float64)

Calculate the half-range Legendre polynomials up to order L.

# Input Argument(s)
- `L::Int64`: maximum Legendre order.
- `x::Float64`: direction cosine.

# Output Argument(s)
- `Pl::Vector{Float64}`: half-range Legendre polynomials up to L, evaluated at μ and ϕ.

"""
function half_range_legendre_polynomials_up_to_L(L::Int64,x::Float64)
    return sqrt.(2*(0:L).+1) .* legendre_polynomials_up_to_L(L,2*x-1)
end


