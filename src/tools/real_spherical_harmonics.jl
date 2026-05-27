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
- `ψlm::Vector{Float64}`: half-range real spherical harmonics up to L, evaluated at μ and ϕ.

"""
function real_half_range_spherical_harmonics_up_to_L(L::Int64,μ::Float64,ϕ::Float64)

    # Verification of input paramters
    if L < 0 error("Legendre order is greater or equal to zero.") end
    if ~(0 ≤ μ ≤ 1) error("Invalid direction cosine (should be between 0 and 1).") end

    Plm = [zeros(l+1) for l in 0:L]
    if μ == -1 || μ == 1
        Pl = jacobi_polynomials_up_to_L(L,0,0,2*μ-1)
        for l in 0:L
            Plm[l+1][1] = Pl[l+1]
        end
    else
        for m in 0:L
            Pl = jacobi_polynomials_up_to_L(L-m,m,m,2*μ-1)
            for l in m:L
                Plm[l+1][m+1] = factorial_factor([l+m],[l],[(2,-m),(1-(2*μ-1)^2,m/2)]) * Pl[l-m+1]
            end
        end
    end

    # Validation
    for l in 0:L, m in 0:l
        if isnan(Plm[l+1][m+1]) error("NaN for Plm (l = $l and m = $m)") end
        if isinf(Plm[l+1][m+1]) error("Inf for Plm (l = $l and m = $m)") end
    end
    
    # Compute half-range spherical harmonics
    ψlm = [zeros(2*l+1) for l in 0:L]
    for l in range(0,L), m in range(-l,l)
        if (m ≥ 0) 𝓣m = cos(m*ϕ) else 𝓣m = sin(abs(m)*ϕ) end
        Clm = sqrt((2-(m == 0))/2π * (2*l+1) * factorial_factor([l-abs(m)],[l+abs(m)]))
        ψlm[l+1][l+m+1] = Clm * Plm[l+1][abs(m)+1] * 𝓣m
    end
    return ψlm
end

"""
    real_half_range_spherical_harmonics_up_to_L!(Yhalf::AbstractVector{Float64},
        Plm::AbstractMatrix{Float64}, Pl::AbstractVector{Float64},
        L::Int64, μ::Float64, ϕ::Float64)

In-place variant of `real_half_range_spherical_harmonics_up_to_L` that writes the
flat (length (L+1)²) result into `Yhalf`, ordered with `(l,m)` running outer-to-inner
over `l=0..L, m=-l..l`. Uses `Plm` and `Pl` as workspaces; no internal allocations.

Layout: `Yhalf[(l)² + l + m + 1]` holds the (l,m) component.
"""
function real_half_range_spherical_harmonics_up_to_L!(Yhalf::AbstractVector{Float64},Plm::AbstractMatrix{Float64},Pl::AbstractVector{Float64},L::Int64,μ::Float64,ϕ::Float64)

    # Associated Legendre functions Plm(l, m, 2μ-1), 0 ≤ m ≤ l ≤ L
    # Stored as Plm[l+1, m+1]; entries with m > l are unused.
    fill!(Plm, 0.0)
    twoμm1 = 2*μ - 1
    if μ == -1 || μ == 1
        jacobi_polynomials_up_to_L!(Pl, L, 0, 0, twoμm1)
        for l in 0:L
            Plm[l+1, 1] = Pl[l+1]
        end
    else
        sq = 1 - twoμm1*twoμm1
        for m in 0:L
            jacobi_polynomials_up_to_L!(Pl, L-m, m, m, twoμm1)
            # factor = (l+m)!/l! * 2^(-m) * (1-(2μ-1)²)^(m/2)
            # We compute (l+m)!/l! incrementally as `pochhammer`.
            two_pow_m = 2.0^m
            sq_pow_m_half = m == 0 ? 1.0 : sq^(m/2)
            for l in m:L
                pochhammer = 1.0
                @inbounds for i in (l+1):(l+m)
                    pochhammer *= i
                end
                Plm[l+1, m+1] = pochhammer / two_pow_m * sq_pow_m_half * Pl[l-m+1]
            end
        end
    end

    # Half-range spherical harmonics, flat ordering
    p = 0
    inv_two_π = 1/(2π)
    @inbounds for l in 0:L
        twol1 = 2*l + 1
        for m in -l:l
            p += 1
            am = abs(m)
            𝓣m = m ≥ 0 ? cos(m*ϕ) : sin(am*ϕ)
            # factor = (l-am)!/(l+am)! = 1 / ∏_{i=l-am+1..l+am} i  (for am ≥ 1; 1 for am=0)
            inv_pochhammer = 1.0
            for i in (l-am+1):(l+am)
                inv_pochhammer /= i
            end
            Clm = sqrt((2 - (m == 0)) * inv_two_π * twol1 * inv_pochhammer)
            Yhalf[p] = Clm * Plm[l+1, am+1] * 𝓣m
        end
    end
    return nothing
end

"""
    jacobi_polynomials_up_to_L!(Pl::AbstractVector{Float64}, L::Int64,
        α::Int64, β::Int64, x::Real)

In-place variant: writes Jacobi polynomials P_l^{(α,β)}(x), l = 0..L, into `Pl`.
"""
function jacobi_polynomials_up_to_L!(Pl::AbstractVector{Float64},L::Int64,α::Int64,β::Int64,x::Real)
    if x == 1
        for l in 0:L
            Pl[l+1] = factorial_factor([l+α],[l,α])
        end
        return nothing
    elseif x == -1
        for l in 0:L
            Pl[l+1] = (-1)^l * factorial_factor([l+β],[l,β])
        end
        return nothing
    end
    Pl[1] = 1
    if L ≥ 1
        Pl[2] = (α+1) + 0.5*(α+β+2)*(x-1)
    end
    @inbounds for l in 2:L
        if l < 125
            Pl[l+1] = ((2*l+α+β-1)*((2*l+α+β)*(2*l+α+β-2)*x+α^2-β^2)*Pl[l] - (2*(l+α-1)*(l+β-1)*(2*l+α+β))*Pl[l-1])/(2*l*(l+α+β)*(2*l+α+β-2))
        else
            θ = acos(x)
            kθ = 1/(sqrt(π)*sin(0.5*θ)^(α+0.5)*cos(0.5*θ)^(β+0.5))
            N = l + 0.5 * (α + β + 1)
            γ = -0.5*π * (α+0.5)
            Pl[l+1] = 1/sqrt(l) * kθ * cos(N*θ+γ)
        end
    end
    return nothing
end

"""
    real_octant_range_spherical_harmonics_up_to_L(L::Int64,μ::Float64,ϕ::Float64)

Calculate the octant-range real spherical harmonics up to order L.

# Input Argument(s)
- `L::Int64`: maximum Legendre order.
- `μ::Float64`: direction cosine.
- `ϕ::Float64`: azimuthal angle.

# Output Argument(s)
- `ψlm::Vector{Float64}`: octant-range real spherical harmonics up to L, evaluated at μ and ϕ.

"""
function real_octant_range_spherical_harmonics_up_to_L(L::Int64,μ::Float64,ϕ::Float64)

    # Verification of input paramters
    if L < 0 error("Legendre order is greater or equal to zero.") end
    if ~(0 ≤ μ ≤ 1) error("Invalid direction cosine (should be between 0 and 1).") end
    if ~(0 ≤ ϕ ≤ π/2) error("Invalid azimuthal angle (should be between 0 and π/2).") end

    Plm = [zeros(l+1) for l in 0:L]
    if μ == -1 || μ == 1
        Pl = jacobi_polynomials_up_to_L(L,0,0,2*μ-1)
        for l in 0:L
            Plm[l+1][1] = Pl[l+1]
        end
    else
        for m in 0:L
            Pl = jacobi_polynomials_up_to_L(L-m,m,m,2*μ-1)
            for l in m:L
                Plm[l+1][m+1] = factorial_factor([l+m],[l],[(2,-m),(1-(2*μ-1)^2,m/2)]) * Pl[l-m+1]
            end
        end
    end

    # Validation
    for l in 0:L, m in 0:l
        if isnan(Plm[l+1][m+1]) error("NaN for Plm (l = $l and m = $m)") end
        if isinf(Plm[l+1][m+1]) error("Inf for Plm (l = $l and m = $m)") end
    end
    
    # Compute half-range spherical harmonics
    ψlm = [zeros(2*l+1) for l in 0:L]
    for l in range(0,L), m in range(-l,l)
        if (m ≥ 0) 𝓣m = cos(4*m*ϕ) else 𝓣m = sin(4*abs(m)*ϕ) end
        Clm = sqrt(2*(2-(m == 0))/π * (2*l+1) * factorial_factor([l-abs(m)],[l+abs(m)]))
        ψlm[l+1][l+m+1] = Clm * Plm[l+1][abs(m)+1] * 𝓣m
    end
    return ψlm
end

"""
    real_patch_range_spherical_harmonics_up_to_L(L::Int64,μ::Float64,ϕ::Float64,u::Int64,v::Int64,w::Int64)

Calculate the real spherical harmonics over a patch (u,v,w) up to order L.

# Input Argument(s)
- `L::Int64`: maximum Legendre order.
- `μ::Float64`: direction cosine.
- `ϕ::Float64`: azimuthal angle.

# Output Argument(s)
- `ψlm::Vector{Float64}`: patch-range real spherical harmonics up to L, evaluated at μ and ϕ.

"""
function real_patch_range_spherical_harmonics_up_to_L(L::Int64,μ::Float64,ϕ::Float64,Nv::Int64,u::Int64,v::Int64,w::Int64)

    # Verification of input paramters
    if L < 0 error("Legendre order is greater or equal to zero.") end
    if ~(1 ≤ u ≤ 8) error("Invalid patch index u (should be between 1 and 8).") end
    if v <= 0 error("Invalid patch index v (should be greater than zero).") end

    # Direction cosines signs
    if u ∈ [1,2,3,4] sx = 1 else sx = -1 end
    if u ∈ [1,2,5,6] sy = 1 else sy = -1 end
    if u ∈ [1,3,5,7] sz = 1 else sz = -1 end

    # Compute patch-range spherical harmonics
    Nw = -sx*v + (sx+1)/2*(Nv+1)
    if ~(1 ≤ w ≤ Nw) error("Invalid patch index w (should be between 1 and Nw).") end
    μ_v(v) = sx*(1-(1-v+(sx+1)/2*Nv)*(-v+(sx+1)/2*(Nv+2))/(Nv*(Nv+1)))
    ϕ_w(w) = (π/2)*(w-1)/Nw + π/2 * (2 + (sy+1)/2 - (sz+1)/2 - (sy+1)*(sz+1)/2)
    Δμ = μ_v(v+1) - μ_v(v)
    Δϕ = ϕ_w(w+1) - ϕ_w(w)
    ψlm = sqrt(π/(2*Δμ*Δϕ)) * real_octant_range_spherical_harmonics_up_to_L(L,(μ-μ_v(v))/Δμ,π/2 * (ϕ-ϕ_w(w))/Δϕ)
    return ψlm
end

"""
    real_patch_range_spherical_harmonics_up_to_L_symmetric(L::Int64,μ::Float64,ϕ::Float64,Nv::Int64,u::Int64,i::Int64,j::Int64)

Calculate the real spherical harmonics over a sub-triangle patch from the
barycentric subdivision of octant `u`, up to order `L`.

The sub-triangle `K` (with vertices `A, B, C`) is mapped to the reference
octant triangle `T₀ = ((1,0,0), (0,1,0), (0,0,1))` via the barycentric
correspondence; the octant-range spherical harmonics are evaluated at the
image point `(μ', ϕ')` and rescaled by the patch area.

# Input Argument(s)
- `L::Int64`: maximum Legendre order.
- `μ::Float64`, `ϕ::Float64`: spherical coordinates inside patch `K`.
- `Nv::Int64`: barycentric subdivision order.
- `u::Int64`: octant index.
- `i::Int64`, `j::Int64`: row and position indices of the sub-triangle inside the octant.

# Output Argument(s)
- `ψlm::Vector{Vector{Float64}}`: patch-range real spherical harmonics up to L.
"""
function real_patch_range_spherical_harmonics_up_to_L_symmetric(L::Int64,μ::Float64,ϕ::Float64,Nv::Int64,u::Int64,i::Int64,j::Int64)
    if L < 0 error("Legendre order must be ≥ 0.") end
    A, B, C = octant_subtriangle_vertices(u, i, j, Nv)
    μp, ϕp, _, _ = muphi_to_muprime_phiprime(μ, ϕ, A, B, C)
    μp = clamp(μp, 0.0, 1.0)
    ϕp = clamp(ϕp, 0.0, π/2)
    area_K = spherical_triangle_area(A, B, C)
    ψlm = real_octant_range_spherical_harmonics_up_to_L(L, μp, ϕp)
    norm_factor = sqrt(π / (2 * area_K))
    for l in 0:L
        for m_idx in eachindex(ψlm[l+1])
            ψlm[l+1][m_idx] *= norm_factor
        end
    end
    return ψlm
end

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