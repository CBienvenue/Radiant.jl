function double_factorial(n::Integer)
    n < -1 && throw(DomainError(n, "n must be ≥ -1"))
    n ≤ 0 && return one(n)
    r = one(n)
    for k in n:-2:1
        r *= k
    end
    return r
end