"""
    factorial_factor(num_terms::Vector{Int64},denom_terms::Vector{Int64},
    exp_terms::Vector{<:Tuple{<:Real,<:Real}}=Tuple{Real,Real}[])

Compute a factor involving factorials and exponentials in a numerically stable way.

# Input Argument(s)
- `num_terms::Vector{Int64}`: list of integers in the numerator of the factorial term.
- `denom_terms::Vector{Int64}`: list of integers in the denominator of the factorial term.
- `exp_terms::Vector{Tuple{<:Real,<:Real}}`: list of tuples (x,m) where x is the base and m 
  the exponent of the exponential term.

# Output Argument(s)
- `factor::Float64`: value of the factorial term.

"""
function factorial_factor(num_terms::Vector{Int64},denom_terms::Vector{Int64},exp_terms::Vector{<:Tuple{<:Real,<:Real}}=Tuple{Real,Real}[])
    factor = 0
    for n in num_terms, i in 1:n
        factor += log(i)
    end
    for n in denom_terms, i in 1:n
        factor -= log(i)
    end
    for (x,m) in exp_terms
        if m ≈ 0 continue end
        if x ≈ 0
            if m > 0 return 0 else return Inf end
        end
        factor += m * log(x)  
    end
    return exp(factor)
end