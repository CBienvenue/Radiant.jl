"""
    natural_cubic_spline(x,y)

Produce the natural cubic spline function for (x,y) dataset.

# Input Argument(s)
- 'x::Vector{Float64}': x.
- 'y::Vector{Float64}': y.

# Output Argument(s)
- 'f::function': monotone cubic hermite spline function

# Author(s)
Charles Bienvenue

# Reference(s)

"""
function natural_cubic_spline(x,y)
    
    # Initialize
    n = length(x)-1
    A = zeros(4*n,4*n)
    b = zeros(4*n)
    
    # 2n equations for function intersecting with (xᵢ,yᵢ) data points for i ∈ [1,n+1]
    for i in range(1,2n)
        index_spaces = div(i+1,2)
        index_points = div(i+2,2)
        for j in range(0,3)
            A[i,4*index_spaces-3+j] = x[index_points]^(3-j)
        end
        b[i] = y[index_points]
    end

    # (n-1) equations for 1st derivative continuity at (xᵢ,yᵢ) data points for i ∈ [2,n]
    for i in range(1,n-1)
        for j in range(0,2)
            A[2*n+i,4*i-3+j] = (3-j)*x[i+1]^(2-j)
            A[2*n+i,4*(i+1)-3+j] = -(3-j)*x[i+1]^(2-j)
        end
    end

    # (n-1) equations for 2nd derivative continuity at (xᵢ,yᵢ) data points for i ∈ [2,n]
    for i in range(1,n-1)
        for j in range(0,1)
            A[2*n+(n-1)+i,4*i-3+j] = (3-j)*(2-j)*x[i+1]^(1-j)
            A[2*n+(n-1)+i,4*(i+1)-3+j] = -(3-j)*(2-j)*x[i+1]^(1-j)
        end
    end

    # 2 equations for boundary conditions (natural spline) at (xᵢ,yᵢ) data points for i = 1 and i = n+1
    for j in range(0,1)
        A[2*n+2*(n-1)+1,1+j] = (3-j)*(2-j)*x[1]^(1-j)
        A[2*n+2*(n-1)+2,4*n-3+j] = (3-j)*(2-j)*x[n+1]^(1-j)
    end

    # Solve to find coefficients
    coeff = A\b

    # Return the resulting cubic spline interpolation function
    return function natural_cubic_spline(xi)

        i = searchsortedfirst(x,xi)-1
        if (i == 0 && x[1] ≈ xi) i=1; xi=x[1] end
        if (i == length(x) && x[end] ≈ xi) i=length(x)-1; xi=x[end] end
        if (i == 0 || i == length(x)) error("Interpolation value is outside the interpolation vector.") end
        p = coeff[4*(i-1)+1:4*(i-1)+4]
        return p[1]*xi^3 + p[2]*xi^2 + p[3]*xi + p[4]

    end 
end

"""
    cubic_hermite_spline(x,y)

Produce the monotone cubic hermite spline function for (x,y) dataset.

# Input Argument(s)
- 'x::Vector{Float64}': x.
- 'y::Vector{Float64}': y.

# Output Argument(s)
- 'f::function': monotone cubic hermite spline function

# Author(s)
Charles Bienvenue

# Reference(s)
- Fritsch (1980) : Monotone piecewise cubic interpolation.

"""
function cubic_hermite_spline(x,y) 

    # Initialize
    n = length(x)
    if x[1] > x[2] x = reverse(x); y = reverse(y) end

    # Compute parameters
    m = zeros(n)
    δ = zeros(n-1)
    for i in range(1,n-1) δ[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) end
    m[1] = δ[1]; m[n] = δ[n-1]
    for i in range(2,n-1) m[i] = (δ[i-1]+δ[i])/2 end
    for i in range(1,n-1)
        if δ[i] == 0
            m[i] = 0
            m[i+1] = 0
        else
            α = m[i]/δ[i]
            β = m[i+1]/δ[i]
            h = sqrt(α^2+β^2)
            if h > 3
                t = 3/h
                m[i] = t*α*δ[i]
                m[i+1] = t*β*δ[i]
            end
        end
    end
    
    # Return the resulting cubic spline interpolation function
    return function cubic_hermite_spline(xi)

        i = searchsortedfirst(x,xi)-1
        if (i == 0 && x[1] ≈ xi) i=1; xi=x[1] end
        if (i == length(x) && x[end] ≈ xi) i=length(x)-1; xi=x[end] end
        if (i == 0 || i == length(x)) error("Interpolation value is outside the interpolation vector.") end

        h = (x[i+1]-x[i])
        t = (xi-x[i])/h
        h00 = (1+2*t)*(1-t)^2
        h10 = t*(1-t)^2
        h01 = t^2*(3-2*t)
        h11 = t^2*(t-1)
        return h00*y[i] + h01*y[i+1] + h10*h*m[i] + h11*h*m[i+1]

    end 
end