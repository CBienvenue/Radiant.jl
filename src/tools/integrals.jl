"""
    𝒢₁(x::Float64,N::Int64,a::Float64,b::Float64)

Compute the integral I = ∫(xⁿ/R²)dx, where R = a + bx, evaluated at x.

# Input Argument(s)
- 'x::Float64': evaluation point.
- 'n::Int64': exponent.
- 'a::Float64': coefficient.
- 'b::Float64': coefficient.

# Output Argument(s)
- 'v::Vector{Float64}': integral evaluated at x for 0 ≤ n ≤ N.

# Reference(s)
- Gradshteyn (2014) : Table of integrals, series, and products.

"""
function 𝒢₁(x::Real,N::Int64,a::Real,b::Real)

    # Initialize
    v = zeros(N+1)
    if x == 0 return v end
    R = a+b*x

    # Eq. 2.111 (4)
    @inbounds for n in range(0,N)
        v[n+1] += (-1)^(n-1)*a^n/(b^(n+1)*R) + (-1)^(n+1)*n*a^(n-1)/b^(n+1)*log(R)
        for g in range(1,n-1)
            v[n+1] += (-1)^(g-1) * (g*a^(g-1)*x^(n-g))/((n-g)*b^(g+1))
        end
    end
    
    return v
end

"""
    𝒢₂(x::Float64,N::Int64,a::Float64,b::Float64)

Compute the integral I = ∫(xⁿ/R²)dx, where R = a + bx^2, evaluated at x.

# Input Argument(s)
- 'x::Float64': evaluation point.
- 'n::Int64': exponent.
- 'a::Float64': coefficient.
- 'b::Float64': coefficient.

# Output Argument(s)
- 'v::Vector{Float64}': integral evaluated at x for 2 ≤ n ≤ N and n even only.

# Author(s)
Charles Bienvenue

# Reference(s)
- Gradshteyn (2014) : Table of integrals, series, and products.

"""
function 𝒢₂(x::Real,N::Int64,a::Real,b::Real)

    # Initialize
    v = zeros(N)
    if x == 0 return v end
    R = a+b*x^2

    # Eq. 2.172 and Eq. 2.173 (1)
    v₀ = x/(2*a*R) + 1/(2*a*sqrt(a*b))*atan(b*x/sqrt(a*b))

    # Eq. 2.174 (1)
    v[2] = -x/(b*R) + a/b*v₀
    @inbounds for n in range(4,N,step=2)
        v[n] = -x^(n-1)/((3-n)*b*R) + (n-1)*a/((3-n)*b)*v[n-2]
    end

    return v
end

"""
    𝒢₃(n::Int64,m::Int64,a::Real,b::Real,α::Real,β::Real,x::Real)

Compute the integral I = ∫tⁿzᵐdx, where z = a + bx and t = α + βx, evaluated at x.

# Input Argument(s)
- 'x::Float64': evaluation point.
- 'n::Int64': exponent.
- 'm::Int64': exponent.
- 'a::Float64': coefficient.
- 'b::Float64': coefficient.
- 'α::Float64': coefficient.
- 'β::Float64': coefficient.

# Output Argument(s)
- 'v::Vector{Float64}': integral evaluated at x for (n,m).

# Author(s)
Charles Bienvenue

# Reference(s)
- Gradshteyn (2014) : Table of integrals, series, and products.

"""
function 𝒢₃(n::Int64,m::Int64,a::Real,b::Real,α::Real,β::Real,x::Real)

    Δ = a*β-b*α
    z = a+b*x
    t = α+β*x
    if β == 0 || b == 0 error() end

    if Δ == 0
        return (a/α)^m*𝒢₃(n+m,0,a,b,α,β,x)
    elseif n == 0 && m == 0 
        return x
    elseif n == 0
        if m != -1 return z^(m+1)/(b*(m+1)) else return log(abs(z))/b end
    elseif m == 0
        if n != -1 return t^(n+1)/(β*(n+1)) else return log(abs(t))/β end
    elseif m < 0 && n < 0
        if  m == -1 && n == -1 # Gradshteyn - Sect. 2.154
            return log(t/z)/Δ
        elseif n != -1 && m < 0 # Gradshteyn - Sect. 2.155
            return 1/(Δ*(n+1))*( t^(n+1)*z^(m+1) - b*(m+n+2)*𝒢₃(n+1,m,a,b,α,β,x) )
        elseif m != -1 && n < 0 # Gradshteyn - Sect. 2.155
            return -1/(Δ*(m+1))*( t^(n+1)*z^(m+1) - β*(m+n+2)*𝒢₃(n,m+1,a,b,α,β,x) )
        else
            error()
        end
    elseif m > 0 && n > 0 # Gradshteyn - Sect. 2.151 & Sect. 2.153
        return 1/(b*(m+n+1))*( z^(m+1)*t^n - n*Δ*𝒢₃(n-1,m,a,b,α,β,x) )
    elseif n > 0 && m < 0 # Gradshteyn - Sect. 2.153
        if m+n != -1
            return 1/(b*(m+n+1))*( z^(m+1)*t^n - n*Δ*𝒢₃(n-1,m,a,b,α,β,x) )
        else
            return  -1/((m+1)*Δ) * ( z^(m+1)*t^(n+1) - (m+n+2)*β*𝒢₃(n,m+1,a,b,α,β,x) )
        end
    elseif m > 0 && n < 0 # Gradshteyn - Sect. 2.153
        if m+n != -1
            return 1/(β*(m+n+1))*( t^(n+1)*z^m + m*Δ*𝒢₃(n,m-1,a,b,α,β,x) )
        else
            return  1/((n+1)*Δ) * ( z^(m+1)*t^(n+1) - (m+n+2)*b*𝒢₃(n+1,m,a,b,α,β,x) )
        end
    else
        error()
    end

end

"""
    𝒢₄(n::Int64,m::Int64,a::Real,b::Real,α::Real,β::Real,x::Real)

Compute the integral I = ∫(1/(tⁿzᵐ√z))dx, where z = a + bx and t = α + βx, evaluated at x.

# Input Argument(s)
- 'x::Float64': evaluation point.
- 'n::Int64': exponent.
- 'm::Int64': exponent.
- 'a::Float64': coefficient.
- 'b::Float64': coefficient.
- 'α::Float64': coefficient.
- 'β::Float64': coefficient.

# Output Argument(s)
- 'v::Vector{Float64}': integral evaluated at x for (n,m).

# Author(s)
Charles Bienvenue

# Reference(s)
- Gradshteyn (2014) : Table of integrals, series, and products.
- Zwillinger (2003) : Standard mathematical tables and formulae.


"""
function 𝒢₄(n::Int64,m::Int64,a::Real,b::Real,α::Real,β::Real,x::Real)

    Δ = a*β-b*α
    z = a+b*x
    t = α+β*x
    if m < 0 error("Negative m index") end
    if n < 0 || m < 0 || β == 0 || b == 0 error("Negative n index") end

    if Δ == 0
        return (α/a)^m*𝒢₄(n+m,0,a,b,α,β,x)
    elseif n == 0
        return z^(0.5-m)/(b*(0.5-m))
    elseif m == 0 && n == 1 # Gradshteyn - Sect. 2.246
        if Δ*β > 0
            return 1/sqrt(β*Δ) * log((β*sqrt(z)-sqrt(β*Δ))/(β*sqrt(z)+sqrt(β*Δ)))
        elseif Δ*β < 0
            return 2/sqrt(-β*Δ) * atan(β*sqrt(z)/sqrt(-β*Δ))
        else
            return -2*sqrt(z)/(b*t)
        end 
    elseif m == 0 && n != 1 # Zwillinger - Sect. 5.4.10 (147)
        return -1/((n-1)*Δ) * (sqrt(z)/t^(n-1) + (n-3/2)*b*𝒢₄(n-1,0,a,b,α,β,x))
    elseif m > 0 && n == 1 # Gradshteyn - Sect. 2.249 (2)
        return 2/((2*m-1)*Δ) * 1/(z^(m-1)*sqrt(z)) + β/Δ*𝒢₄(1,m-1,a,b,α,β,x) 
    elseif m > 0 && n != 1 # Gradshteyn - Sect. 2.249 (1)
        return -1/((n-1)*Δ)*sqrt(z)/(z^m*t^(n-1)) - (2*n+2*m-3)*b/(2*(n-1)*Δ)*𝒢₄(n-1,m,a,b,α,β,x)
    else
        error()
    end
end

"""
    𝒢₅(n::Int64,a::Real,b::Real,c::Real,x::Real)

Compute the integral I = ∫(xⁿ/√R)dx, where R = a + bx + cx², evaluated at x.

# Input Argument(s)
- 'x::Real': evaluation point.
- 'n::Int64': exponent.
- 'a::Real': coefficient.
- 'b::Real': coefficient.
- 'c::Real': coefficient.

# Output Argument(s)
- 'v::Vector{Float64}': integral evaluated at x up to order n.

# Author(s)
Charles Bienvenue

# Reference(s)
- Gradshteyn (2014) : Table of integrals, series, and products.


"""
function 𝒢₅(n::Int64,a::Real,b::Real,c::Real,x::Real)
    R = a + b*x + c*x^2
    Δ = 4*a*c-b^2
    if n == 0 # Gradshteyn - Sect. 2.261
        if Δ > 0 && c > 0
            return 1/sqrt(c)*asinh((2*c*x+b)/sqrt(Δ))
        elseif Δ == 0
            return 1/sqrt(c)*log(2*c*x+b)
        elseif Δ < 0 && c < 0
            return -1/sqrt(-c)*asinh((2*c*x+b)/sqrt(-Δ))
        elseif Δ < 0 && c > 0 && 2*c*x+b > sqrt(-Δ)
            return 1/sqrt(c)*log(2*sqrt(c*R)+2*c*x+b)
        elseif Δ < 0 && c > 0 && 2*c*x+b < -sqrt(-Δ)
            return -1/sqrt(c)*log(2*sqrt(c*R)-2*c*x-b)
        else
            error("Unknown case.")
        end
    elseif n == 1 # Gradshteyn - Sect. 2.264 (2)
        return sqrt(R)/c - b/(2*c)*𝒢₅(0,a,b,c,x)
    elseif n ≥ 2 # Gradshteyn - Sect. 2.263 (1)
        return x^(n-1)/(n*c)*sqrt(R) - (2*n-1)*b/(2*n*c)*𝒢₅(n-1,a,b,c,x) - (n-1)*a/(n*c)*𝒢₅(n-2,a,b,c,x)
    else
        error("Negative n index")
    end
end

"""
    𝒢₆(j₁::Int64,j₂::Int64,j₃::Int64)

Compute the integral I = 0.5 ∫ Pⱼ₁(x)Pⱼ₂(x)Pⱼ₃(x) dx evaluated between -1 and 1.

# Input Argument(s)
- 'j₁::Int64': first Legendre polynomial index.
- 'j₂::Int64': second Legendre polynomial index.
- 'j₃::Int64': third Legendre polynomial index.

# Output Argument(s)
- '𝒲::Float64': integral evaluated between -1 and 1.

# Author(s)
Charles Bienvenue

# Reference(s)

"""
function 𝒢₆(j₁::Int64,j₂::Int64,j₃::Int64)

    Cjk = zeros(3,div(max(j₁,j₂,j₃),2)+1) 
    for k in range(0,div(j₁,2)) Cjk[1,k+1] = (-1)^k * exp( sum(log.(1:2*j₁-2*k)) - sum(log.(1:k)) - sum(log.(1:j₁-k)) - sum(log.(1:j₁-2*k)) ) end
    for k in range(0,div(j₂,2)) Cjk[2,k+1] = (-1)^k * exp( sum(log.(1:2*j₂-2*k)) - sum(log.(1:k)) - sum(log.(1:j₂-k)) - sum(log.(1:j₂-2*k)) ) end
    for k in range(0,div(j₃,2)) Cjk[3,k+1] = (-1)^k * exp( sum(log.(1:2*j₃-2*k)) - sum(log.(1:k)) - sum(log.(1:j₃-k)) - sum(log.(1:j₃-2*k)) ) end
    
    𝒲 = 0
    J = j₁ + j₂ + j₃
    for k₁ in range(0,div(j₁,2)), k₂ in range(0,div(j₂,2)), k₃ in range(0,div(j₃,2))
        K = k₁ + k₂ + k₃
        if mod(J-2*K,2) == 0
            𝒲 += Cjk[1,k₁+1] * Cjk[2,k₂+1] * Cjk[3,k₃+1] * 1/(J-2*K+1)
        end
    end 
    𝒲 /= 2^J
    
    return 𝒲
end