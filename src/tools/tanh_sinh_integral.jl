"""
    tanh_sinh_integral(f::Function, a::Float64, b::Float64;
                       h::Float64=0.1, n::Int=60, tol::Float64=1.0e-10,
                       max_iter::Int=6)

Tanh-sinh (double-exponential) quadrature on [a,b].
"""
function tanh_sinh_integral(f::Function, a::Float64, b::Float64;
                            h::Float64=0.1, n::Int=60, tol::Float64=1.0e-10,
                            max_iter::Int=6)
    mid = 0.5 * (a + b)
    half = 0.5 * (b - a)
    function integrate(n_local::Int)
        acc = 0.0
        for k in -n_local:n_local
            t = k * h
            sh = sinh(t)
            ch = cosh(t)
            u = (pi / 2.0) * sh
            x = tanh(u)
            dxdt = (pi / 2.0) * ch / (cosh(u)^2)
            mu = mid + half * x
            acc += f(mu) * dxdt
        end
        return acc * half * h
    end

    prev = integrate(n)
    for _ in 1:max_iter
        n = Int(round(n * 1.5))
        curr = integrate(n)
        if abs(curr - prev) <= tol * max(1.0, abs(prev))
            return curr
        end
        prev = curr
    end
    return prev
end
