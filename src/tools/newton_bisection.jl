"""
    newton_bisection(f::Function,dfdx::Function,x⁻::Real,x⁺::Real)
⁻
Mixed Newton-bissection iterative root-finding method.

# Input Argument(s)
- `f::Function`: function f(x).
- `dfdx::Function`: derivative of f(x).
- `x₁::Real`: x-value such as f(x⁻) < 0.
- `x₂::Real`: x-value such as f(x⁺) > 0.
- `ϵ₀::Float64`: convergence criterion.
- `N_max::Int64`: maximum number of iterations.

# Output Argument(s)
- `x::Real`: root of f(x).

# Reference(s)

"""
function newton_bisection(f::Function, dfdx::Function, x₁::Real, x₂::Real, ϵ₀::Float64 = eps(), N_max::Int = 100)
    
    # Evaluate function at interval endpoints
    f₁ = f(x₁)
    f₂ = f(x₂)

    # Ensure that the root is bracketed
    if f₁ * f₂ > 0.0
        error("Root must be bracketed: f(x₁) and f(x₂) must have opposite signs.")
    elseif f₁ == 0.0
        return x₁
    elseif f₂ == 0.0
        return x₂
    end

    # Choose interval so that f(xl) < 0 and f(xh) > 0
    xl, xh = f₁ < 0 ? (x₁, x₂) : (x₂, x₁)
    
    # Initial guess: midpoint of interval
    rts = 0.5 * (x₁ + x₂)

    # Track step sizes for convergence control
    dxold = abs(x₂ - x₁)
    dx = dxold

    # Evaluate function and derivative at initial guess
    f_rts = f(rts)
    df_rts = dfdx(rts)

    for j in 1:N_max

        # Determine if Newton step is acceptable
        # If step would go outside [xl, xh] or Newton step is too large, fall back to bisection
        if ((rts - xh) * df_rts - f_rts) * ((rts - xl) * df_rts - f_rts) > 0.0 || abs(2.0 * f_rts) > abs(dxold * df_rts)

            # Bisection step
            dxold = dx
            dx = 0.5 * (xh - xl)
            rts = xl + dx
        else
            # Newton step
            dxold = dx
            dx = f_rts / df_rts
            rts -= dx
        end

        # Convergence check: step size small enough
        if abs(dx) < ϵ₀
            return rts
        end

        # Re-evaluate function and derivative at new guess
        f_rts = f(rts)
        df_rts = dfdx(rts)

        # Update bracketing interval
        if f_rts < 0.0
            xl = rts
        else
            xh = rts
        end
    end

     # If loop exits, convergence failed
    error("Maximum number of iterations ($N_max) exceeded in newton_bisection.")
end

"""
    find_zero_bisection(f::Function, a::Float64, b::Float64; tol=1e-12, max_iter=100)

Finds a scalar root of `f` in the interval `[a, b]` using the bisection method.

# Input Argument(s)
- `f::Function` : scalar function to solve.
- `a::Float64` : lower endpoint of the bracketing interval.
- `b::Float64` : upper endpoint of the bracketing interval.
- `tol::Float64` : absolute tolerance on the residual or interval width.
- `max_iter::Int` : maximum number of bisection iterations.

# Output Argument(s)
- `x::Float64` : estimated root of `f`.
"""
function find_zero_bisection(f::Function, a::Float64, b::Float64;
                             tol::Float64=1.0e-12, max_iter::Int=100)
    fa = f(a)
    fb = f(b)
    if fa == 0.0
        return a
    elseif fb == 0.0
        return b
    elseif signbit(fa) == signbit(fb)
        error("Bisection root finder requires opposite signs at the interval endpoints.")
    end

    lo = a
    hi = b
    for _ in 1:max_iter
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if abs(fm) <= tol || (hi - lo) <= tol
            return mid
        end
        if signbit(fa) == signbit(fm)
            lo = mid
            fa = fm
        else
            hi = mid
            fb = fm
        end
    end
    return 0.5 * (lo + hi)
end