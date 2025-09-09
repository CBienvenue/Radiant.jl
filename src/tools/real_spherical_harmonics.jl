"""
    real_spherical_harmonics(l::Int64,m::Int64,μ::Float64,ϕ::Float64)

Calculate the real spherical harmonics components. 

# Input Argument(s)
- `l::Int64`: order of the spherical harmonic.
- `m::Int64`: degree of the spherical harmonic.
- `μ::Float64`: direction cosine.
- `ϕ::Float64`: azimuthal angle.

# Output Argument(s)
- `Rlm::Float64`: real spherical harmonics of order l and degree m evaluated at μ and ϕ.

# Reference(s)
- Hébert (2016), Applied Reactor Physics.

"""
function real_spherical_harmonics(l::Int64,m::Int64,μ::Float64,ϕ::Float64)

    # Verification of input paramters
    if l < 0 error("Legendre order is greater or equal to zero.") end
    if abs(m) > l error("Legendre order smaller than m-order (|m| > l).") end
    if ~(-1 ≤ μ ≤ 1) error("Invalid direction cosine (should be between -1 and 1).") end

    # Compute the real spherical harmonics
    if (m ≥ 0) 𝓣m = cos(m*ϕ) else 𝓣m = sin(abs(m)*ϕ) end
    Plm = ferrer_associated_legendre(l,abs(m),μ)
    Clm = sqrt((2-(m == 0)) * exp( sum(log.(1:l-abs(m))) - sum(log.(1:l+abs(m))) ))
    Rlm = Clm * Plm * 𝓣m

    return Rlm
end

"""
    real_spherical_harmonics_up_to_L(L::Int64,μ::Float64,ϕ::Float64)

Calculate the real spherical harmonics components up to order L.

# Input Argument(s)
- `L::Int64`: maximum Legendre order.
- `μ::Float64`: direction cosine.
- `ϕ::Float64`: azimuthal angle.

# Output Argument(s)
- `Rlm::Vector{Float64}`: real spherical harmonics up to L, evaluated at μ and ϕ.

# Reference(s)
- Hébert (2016), Applied Reactor Physics.

"""
function real_spherical_harmonics_up_to_L(L::Int64,μ::Float64,ϕ::Float64)

    # Verification of input parameters
    if L < 0 error("Legendre order is greater or equal to zero.") end
    if ~(-1 ≤ μ ≤ 1) error("Invalid direction cosine (should be between -1 and 1).") end

    # Compute the associated Legendre polynomials
    Plm = [zeros(l+1) for l in 0:L]
    if μ == -1 || μ == 1
        Pl = jacobi_polynomials_up_to_L(L,0,0,μ)
        for l in 0:L
            Plm[l+1][1] = Pl[l+1]
        end
    else
        for m in 0:L
            Pl = jacobi_polynomials_up_to_L(L-m,m,m,μ)
            for l in m:L
                Plm[l+1][m+1] = factorial_factor([l+m],[l],[(2,-m),(1-μ^2,m/2)]) * Pl[l-m+1]
            end
        end
    end

    # Validation
    for l in 0:L, m in 0:l
        if isnan(Plm[l+1][m+1]) error("NaN for Plm (l = $l and m = $m)") end
        if isinf(Plm[l+1][m+1]) error("Inf for Plm (l = $l and m = $m)") end
    end

    # Compute real spherical harmonics
    Rlm = [zeros(2*l+1) for l in 0:L]
    for l in 0:L, m in -l:l
        if (m ≥ 0) 𝓣m = cos(m*ϕ) else 𝓣m = sin(abs(m)*ϕ) end
        Clm = sqrt((2-(m == 0)) * factorial_factor([l-abs(m)],[l+abs(m)]))
        Rlm[l+1][l+m+1] = Clm * Plm[l+1][abs(m)+1] * 𝓣m
    end
    return Rlm
end

"""
    real_half_range_spherical_harmonics_up_to_L(L::Int64,μ::Float64,ϕ::Float64)

Calculate the half-range real spherical harmonics up to order L.

# Input Argument(s)
- `L::Int64`: maximum Legendre order.
- `μ::Float64`: direction cosine.
- `ϕ::Float64`: azimuthal angle.

# Output Argument(s)
- `Pl::Vector{Float64}`: half-range real spherical harmonics up to L, evaluated at μ and ϕ.

"""
function real_half_range_spherical_harmonics_up_to_L(L::Int64,μ::Float64,ϕ::Float64)

    # Verification of input paramters
    if L < 0 error("Legendre order is greater or equal to zero.") end
    if ~(0 ≤ μ ≤ 1) error("Invalid direction cosine (should be between 0 and 1).") end
    
    # Compute the associated Jacobi polynomials
    Plm_01 = [zeros(l+1) for l in 0:L]
    if μ == 0 || μ == 1
        Pl = jacobi_polynomials_up_to_L(L,0,1,2*μ-1)
        for l in 0:L
            Plm_01[l+1][1] = Pl[l+1]
        end
    else
        for m in 0:L
            Pl = jacobi_polynomials_up_to_L(L-m,m,m+1,2*μ-1)
            for l in m:L
                Plm_01[l+1][m+1] = factorial_factor([l+m+1],[l+1],[(2,-m),(1-(2*μ-1)^2,m/2)]) * Pl[l-m+1]
            end
        end
    end

    # Validation
    for l in 0:L, m in 0:l
        if isnan(Plm_01[l+1][m+1]) error("NaN for Plm_01 (l = $l and m = $m)") end
        if isinf(Plm_01[l+1][m+1]) error("Inf for Plm_01 (l = $l and m = $m)") end
    end
    
    # Compute half-range spherical harmonics
    ψlm = [zeros(2*l+1) for l in 0:L]
    for l in range(0,L), m in range(-l,l)
        if (m ≥ 0) 𝓣m = cos(m*ϕ) else 𝓣m = sin(abs(m)*ϕ) end
        Clm = sqrt((2-(m == 0))/π * (l+1)^2 * factorial_factor([l-abs(m)],[l+abs(m)+1]))
        ψlm[l+1][l+m+1] = Clm * Plm_01[l+1][abs(m)+1] * 𝓣m
    end
    return ψlm
end