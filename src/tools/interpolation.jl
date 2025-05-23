"""
    linear_interpolation(xi::Float64,x::Vector{Float64},f::Vector{Float64})

Linear interpolation.

# Input Argument(s)
- `xi::Float64`: x-value at which to estimate fi.
- `x::Vector{Float64}`: vector of x-values.
- `f::Vector{Float64}`: vector of f-values.

# Output Argument(s)
- `fi::Float64`: interpolated f-value.

"""
function linear_interpolation(xi::Float64,x::Vector{Float64},f::Vector{Float64})
    
    # Find the index
    i = searchsortedfirst(x,xi)
    if (i == 1 && x[1] ≈ xi) i=2; xi=x[1] end
    if (i == length(x)+1 && x[end] ≈ xi) i=length(x); xi=x[end] end
    if (i == 1 || i > length(x)) error("Interpolation value is outside the interpolation vector.") end

    # Compute the value at xi
    fi = (f[i-1]*(x[i]-xi) + f[i]*(xi-x[i-1]))/(x[i]-x[i-1])

    return fi

end

"""
    linear_interpolation(xi::Float64,x::Vector{Float64},f::Vector{Float64})

Bilinear interpolation.

# Input Argument(s)
- `xi::Float64`: x-value at which to estimate fi.
- `xi::Float64`: y-value at which to estimate fi.
- `x::Vector{Float64}`: vector of x-values.
- `y::Vector{Float64}`: vector of y-values.
- `f::Array{Float64}`: array of f-values.

# Output Argument(s)
- `fi::Float64`: interpolated f-value.

"""
function linear_interpolation(xi::Float64,yi::Float64,x::Vector{Float64},y::Vector{Float64},f::Array{Float64})
    
    # Find the index
    i = searchsortedfirst(x,xi)
    j = searchsortedfirst(y,yi)
    if (i == 1 && x[1] ≈ xi) i=2; xi=x[1] end
    if (i == length(x)+1 && x[end] ≈ xi) j=length(x); xi=x[end] end
    if (j == 1 && y[1] ≈ yi) i=2; yi=y[1] end
    if (j == length(y)+1 && y[end] ≈ yi) j=length(y); yi=y[end] end
    if (i == 1 || j == 1 || i > length(x) || j > length(y)) error("Interpolation value is outside the interpolation vector.") end

    # Compute the value at xi
    fi = ( f[i-1,j-1]*(x[i]-xi)*(y[j]-yi) + f[i,j-1]*(xi-x[i-1])*(y[j]-yi) + f[i-1,j]*(x[i]-xi)*(yi-y[j-1]) + f[i,j]*(xi-x[i-1])*(yi-y[j-1]) ) / ((x[i]-x[i-1])*(y[j]-y[j-1]))

    return fi

end

"""
    log_interpolation(xi::Float64,x::Vector{Float64},f::Vector{Float64},
    is_negative::Bool=true)

Logarithmic interpolation.

# Input Argument(s)
- `xi::Float64`: x-value at which to estimate fi.
- `x::Vector{Float64}`: vector of x-values.
- `f::Vector{Float64}`: vector of f-values.
- `is_negative::Bool` : enable negative input values, which are dealt with linear
  interpolation.

# Output Argument(s)
- `fi::Float64`: interpolated f-value.

"""
function log_interpolation(xi::Float64,x::Vector{Float64},f::Vector{Float64},is_negative::Bool=true)

    # Validation of input parameters
    if xi ≤ 0 || x[1] ≤ 0 || x[2] ≤ 0 || f[1] ≤ 0 || f[2] ≤ 0
        if is_negative # Change to linear interpolation if negative values
            return linear_interpolation(xi,x,f)
        else # Error if negative values
            error("Unexpected negative values.")
        end
    end
    
    # Find the index
    i = searchsortedfirst(x,xi)
    if (i == 1 && x[1] ≈ xi) i=2; xi=x[1] end
    if (i == length(x)+1 && x[end] ≈ xi) i=length(x); xi=x[end] end
    if (i == 1 || i > length(x)) error("Interpolation value is outside the interpolation vector.") end
    
    # Compute interpolation factor in log space
    t = (log(xi) - log(x[i-1])) / (log(x[i]) - log(x[i-1]))
    
    # Linear interpolation in log-space
    fi = exp( log(f[i-1]) + t * (log(f[i]) - log(f[i-1])) )

    return fi
end
