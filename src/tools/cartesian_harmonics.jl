function cartesian_harmonics_number_basis(L)
    return div((L+1)*(L+2)*(L+3),6)
end

function cartesian_harmonics_indices(L)
    Np = cartesian_harmonics_number_basis(L)
    pl = zeros(Int64,Np)
    pa = zeros(Int64,Np)
    pb = zeros(Int64,Np)
    pc = zeros(Int64,Np)
    p = 1
    for l in 0:L
        for a in 0:l
            for b in 0:(l - a)
                c = l - a - b
                pl[p] = l
                pa[p] = a
                pb[p] = b
                pc[p] = c
                p += 1
            end
        end
    end
    return pl, pa, pb, pc
end

function cartesian_harmonics_index(l,a,b)
    return div(l*(l+1)*(l+2),6)+(l+1)*a-div(a*(a-1),2)+b+1
end

function cartesian_harmonics_normalization(l::Int64,a::Int64,b::Int64,c::Int64,is_half_range::Bool=true)
    if l < 0 || a < 0 || b < 0 || c < 0 error("l, a, b, and c should be non-negative.") end
    if l != a + b + c error("l must be equal to a + b + c.") end
    H² = 0
    for i in 0:div(a,2), j in 0:div(b,2), k in 0:div(c,2), i2 in 0:div(a,2), j2 in 0:div(b,2), k2 in 0:div(c,2)
        if is_half_range
            H² += cartesian_harmonics_β(l,i,j,k,a,b,c) * cartesian_harmonics_β(l,i2,j2,k2,a,b,c) * cartesian_harmonics_I(2*(a-i-i2),0,2*(b-j-j2),2*(c-k-k2))
        else
            H² += cartesian_harmonics_β(l,i,j,k,a,b,c) * cartesian_harmonics_β(l,i2,j2,k2,a,b,c) * cartesian_harmonics_J(2*(b-j-j2),2*(c-k-k2)) * cartesian_harmonics_K(2*(a-i-i2),2*(b-j-j2+c-k-k2))
        end
    end
    return 1 / sqrt(H²)
end

function cartesian_harmonics_β(l::Int64,i::Int64,j::Int64,k::Int64,a::Int64,b::Int64,c::Int64)
    if l < 0 || a < 0 || b < 0 || c < 0 || i < 0 || j < 0 || k < 0 error("l, a, b, c, i, j, and k should be non-negative.") end
    if l != a + b + c error("l must be equal to a + b + c.") end
    return (-1)^(i+j+k) * double_factorial(2l-2*(i+j+k)-1)/double_factorial(2l-1)/2^(i+j+k) * factorial_factor([a,b,c],[i,j,k,a-2i,b-2j,c-2k])
end

function cartesian_harmonics_I(a::Int64,b::Int64,c::Int64,d::Int64,is_half_range::Bool=true)
    if a < 0 || b < 0 || c < 0 || d < 0 error("a, b, c and d should be non-negative.") end
    I = 0
    for n in 0:a
        I += 2^n * binomial(a,n) * (-1)^(a-n) * cartesian_harmonics_K(n+b,c+d)
    end
    I *= cartesian_harmonics_J(c,d)
    return I
end

function cartesian_harmonics_J(n::Int64,m::Int64)
    if n < 0 || m < 0 error("n and m should be non-negative.") end
    if iseven(n) && iseven(m)
        return 2*π*double_factorial(n-1)*double_factorial(m-1)/double_factorial(n+m)
    else
        return 0
    end
end

function cartesian_harmonics_K(n::Int64,m::Int64,is_half_range::Bool=true)
    if n < 0 || m < 0 error("n and m should be non-negative.") end
    if is_half_range
        return gamma((n+1)/2)*gamma((m+2)/2)/gamma((n+m+3)/2)/2
    elseif iseven(n)
        return gamma((n+1)/2)*gamma((m+2)/2)/gamma((n+m+3)/2)
    else
        return 0
    end
end

# À déménager !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function spherical_harmonics_number_basis(L)
    return (L+1)^2
end

function spherical_harmonics_indices(L)
    Np = spherical_harmonics_number_basis(L)
    pl = zeros(Int64,Np)
    pm = zeros(Int64,Np)
    p = 1
    for l in range(0,L), m in range(-l,l)
        pl[p] = l
        pm[p] = m
        p += 1
    end
    return pl, pm
end