
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

