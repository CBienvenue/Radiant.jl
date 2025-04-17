"""
    bessel_I₀(x::Float64)

Compute the modified Bessel function of the first kind (order 0).

# Input Argument(s)
- `x::Float64` : x

# Output Argument(s)
- `I₀::Float64` : Modified Bessel function of the first kind (order 0).

# Reference(s)
- Press et al. (2007), Numerical Recipes 3rd Edition: The Art of Scientific Computing.

"""
function bessel_I₀(x::Float64)
    if abs(x) < 3.75
        y = (x/3.75)^2
        return 1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))))
    else
        y = 3.75/abs(x)
        return (exp(abs(x))/sqrt(abs(x)))*(0.39894228+y*(0.1328592e-1+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1+y*0.392377e-2))))))))
    end
end

"""
    bessel_I₁(x::Float64)

Compute the modified Bessel function of the first kind (order 1).

# Input Argument(s)
- `x::Float64` : x

# Output Argument(s)
- `I₁::Float64` : Modified Bessel function of the first kind (order 1).

# Reference(s)
- Press et al. (2007), Numerical Recipes 3rd Edition: The Art of Scientific Computing.

"""
function bessel_I₁(x::Float64)
    if abs(x) < 3.75
        y = (x/3.75)^2
        return sign(x) * abs(x)*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))))
    else
        y = 3.75/abs(x)
        return sign(x) * (exp(abs(x))/sqrt(abs(x)))*(0.39894228+y*(-0.03988024+y*(-0.00362018+y*(0.00163801+y*(-0.01031555+y*(0.02282967+y*(-0.02895312+y*(0.01787654+y*-0.00420059))))))))
    end
end

"""
    bessel_K₀(x::Float64)

Compute the modified Bessel function of the second kind (order 0).

# Input Argument(s)
- `x::Float64` : x

# Output Argument(s)
- `K₀::Float64` : Modified Bessel function of the second kind (order 0).

# Reference(s)
- Press et al. (2007), Numerical Recipes 3rd Edition: The Art of Scientific Computing.

"""
function bessel_K₀(x::Float64)
    if x ≤ 2
        y = x^2/4
        return (-log(x/2.0)*bessel_I₀(x))+(-0.57721566+y*(0.42278420+y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2+y*(0.10750e-3+y*0.74e-5))))))
    else
        y = 2/x
        return (exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1+y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2+y*(-0.251540e-2+y*0.53208e-3))))))
    end
end

"""
    bessel_K₁(x::Float64)

Compute the modified Bessel function of the second kind (order 1).

# Input Argument(s)
- `x::Float64` : x

# Output Argument(s)
- `K₁::Float64` : Modified Bessel function of the second kind (order 1).

# Reference(s)
- Press et al. (2007), Numerical Recipes 3rd Edition: The Art of Scientific Computing.

"""
function bessel_K₁(x::Float64)
    if x ≤ 2
        y = x^2/4
        return (log(x/2.0)*bessel_I₁(x))+(1.0/x)*(1.0+y*(0.15443144+y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1+y*(-0.110404e-2+y*(-0.4686e-4)))))))
    else
        y = 2/x
        return (exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619+y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2+y*(0.325614e-2+y*(-0.68245e-3)))))))
    end
end