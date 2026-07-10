"""
    cpm_exponential_integral(n::Int64,x::T) where {T<:AbstractFloat}

Generalized exponential integral of order `n`,
``E_n(x) = \\int_0^1 \\mu^{n-2} e^{-x/\\mu}\\,d\\mu = \\int_1^\\infty t^{-n} e^{-xt}\\,dt``,
used to close the one-dimensional collision probabilities (see TR-02). It coincides with
`SpecialFunctions.expint(n,x)` for `x > 0`; the finite value `E_n(0) = 1/(n-1)` (``n > 1``)
is returned at the origin. The argument type is preserved (`Float64` on the fast path,
`BigFloat` when the surface reductions need extended precision).

# Input Argument(s)
- `n::Int64` : order of the exponential integral (`n ≥ 1`).
- `x::T` : argument (`x ≥ 0`).

# Output Argument(s)
- `En::T` : value of `E_n(x)`.

"""
function cpm_exponential_integral(n::Int64,x::T) where {T<:AbstractFloat}
    if x ≤ 0
        return n > 1 ? one(T)/(n-1) : T(Inf)
    end
    return expint(n,x)
end

"""
    cpm_optical_thickness(Σ::Vector{Float64},Δx::Vector{Float64},i::Int64,j::Int64)

Accumulated optical thickness ``F_{i,j} = \\sum_{k=i}^{j} \\Delta x_k \\Sigma_k`` between
regions `i` and `j` (`0` when `j < i`).

# Input Argument(s)
- `Σ::Vector{Float64}` : total cross-section per region.
- `Δx::Vector{Float64}` : width per region.
- `i::Int64`, `j::Int64` : region indices.

# Output Argument(s)
- `F::Float64` : accumulated optical thickness.

"""
function cpm_optical_thickness(Σ::Vector{Float64},Δx::Vector{Float64},i::Int64,j::Int64)
    if j < i return 0.0 end
    F = 0.0
    for k in range(i,j) F += Δx[k]*Σ[k] end
    return F
end

"""
    cpm_collision_coefficient(ℓp::Int64,ℓq::Int64,k₁::Int64,k₂::Int64)

Reduction coefficient ``C^{0,k_1,k_2}_{\\ell_p,\\ell_q,0}`` of the azimuthally-symmetric
(`m = 0`) volume--volume collision probability (TR-02, Eq. for ``C^{k_0,k_1,k_2}`` with
`m = 0`, `k_0 = 0`). The associated streaming exponent is
``\\kappa = \\ell_p + \\ell_q - 2k_1 - 2k_2``.

# Input Argument(s)
- `ℓp::Int64`, `ℓq::Int64` : Legendre degrees of the two flux moments.
- `k₁::Int64`, `k₂::Int64` : summation indices.

# Output Argument(s)
- `C::Float64` : reduction coefficient.

"""
function cpm_collision_coefficient(ℓp::Int64,ℓq::Int64,k₁::Int64,k₂::Int64)
    return (-1)^(k₁+k₂)/2.0^(ℓp+ℓq+1) * binomial(ℓp,k₁)*binomial(ℓq,k₂)*binomial(2ℓp-2k₁,ℓp)*binomial(2ℓq-2k₂,ℓq)
end

"""
    cpm_reduced_collision_probability(ℓp::Int64,ℓq::Int64,i::Int64,j::Int64,
    Σ::Vector{Float64},Δx::Vector{Float64})

Geometry-reduced volume--volume collision probability ``\\tilde{p}^{(p,q)}_{i,j}`` for the
azimuthally-symmetric (`m = 0`) case, in one-dimensional Cartesian geometry (TR-02,
Section on the finite-domain method). All cases (adjacent/cross-cell, diagonal, void
regions and reciprocity for `i > j`) are handled in terms of the exponential integrals
``E_n(F_{i,j})``.

# Input Argument(s)
- `ℓp::Int64`, `ℓq::Int64` : Legendre degrees of the flux (`p`) and source (`q`) moments.
- `i::Int64`, `j::Int64` : region indices of the flux (`i`) and source (`j`) regions.
- `Σ::Vector{Float64}` : total cross-section per region.
- `Δx::Vector{Float64}` : width per region.

# Output Argument(s)
- `p̃::Float64` : reduced collision probability ``\\tilde{p}^{(p,q)}_{i,j}``.

"""
function cpm_reduced_collision_probability(ℓp::Int64,ℓq::Int64,i::Int64,j::Int64,Σ::Vector{Float64},Δx::Vector{Float64})

    F(a,b) = cpm_optical_thickness(Σ,Δx,a,b)
    s = 0.0

    if i != j
        # Cross-cell collision. The magnitude is evaluated with the μ > 0 convention on the
        # ordered pair (lo < hi); the physical directional sign is applied afterwards. The ray
        # from the source cell j to the field cell i points toward -x when i < j, so
        # R_p R_q = P_p(μ)P_q(μ) → (-1)^{ℓp+ℓq}, which is +1 for even ℓp+ℓq and -1 otherwise.
        lo = min(i,j); hi = max(i,j)
        if Σ[lo] != 0.0 && Σ[hi] != 0.0
            for k₁ in range(0,ℓp÷2), k₂ in range(0,ℓq÷2)
                κ = ℓp+ℓq-2k₁-2k₂
                s += cpm_collision_coefficient(ℓp,ℓq,k₁,k₂)*(cpm_exponential_integral(κ+3,F(lo+1,hi-1))
                     - cpm_exponential_integral(κ+3,F(lo+1,hi)) - cpm_exponential_integral(κ+3,F(lo,hi-1))
                     + cpm_exponential_integral(κ+3,F(lo,hi)))
            end
            s /= (Σ[lo]*Σ[hi])
        elseif Σ[lo] != 0.0 && Σ[hi] == 0.0
            for k₁ in range(0,ℓp÷2), k₂ in range(0,ℓq÷2)
                κ = ℓp+ℓq-2k₁-2k₂
                s += cpm_collision_coefficient(ℓp,ℓq,k₁,k₂)*(cpm_exponential_integral(κ+2,F(lo+1,hi-1)) - cpm_exponential_integral(κ+2,F(lo,hi-1)))
            end
            s *= Δx[hi]/Σ[lo]
        elseif Σ[lo] == 0.0 && Σ[hi] != 0.0
            for k₁ in range(0,ℓp÷2), k₂ in range(0,ℓq÷2)
                κ = ℓp+ℓq-2k₁-2k₂
                s += cpm_collision_coefficient(ℓp,ℓq,k₁,k₂)*(cpm_exponential_integral(κ+2,F(lo+1,hi-1)) - cpm_exponential_integral(κ+2,F(lo+1,hi)))
            end
            s *= Δx[lo]/Σ[hi]
        else
            for k₁ in range(0,ℓp÷2), k₂ in range(0,ℓq÷2)
                κ = ℓp+ℓq-2k₁-2k₂
                s += cpm_collision_coefficient(ℓp,ℓq,k₁,k₂)*cpm_exponential_integral(κ+1,F(lo+1,hi-1))
            end
            s *= Δx[lo]*Δx[hi]
        end
        return (i < j) ? (-1)^(ℓp+ℓq)*s : s
    else
        # Diagonal collision (i = j) : both streaming directions add
        if Σ[i] != 0.0
            for k₁ in range(0,ℓp÷2), k₂ in range(0,ℓq÷2)
                κ = ℓp+ℓq-2k₁-2k₂
                s += (1+(-1)^κ)*cpm_collision_coefficient(ℓp,ℓq,k₁,k₂)*(
                        (1/(κ+1))*(Δx[i]/Σ[i]) - (1/(κ+2))/Σ[i]^2 + cpm_exponential_integral(κ+3,Δx[i]*Σ[i])/Σ[i]^2)
            end
            return s
        else
            # Void diagonal term multiplies a vanishing source density (set to zero).
            return 0.0
        end
    end
end

"""
    cpm_collision_matrix(Σ::Vector{Float64},Δx::Vector{Float64},Lp::Int64)

Assemble the (un-normalized) volume--volume collision-probability matrix ``P_{vv}`` coupling
the flux moments to the source moments for the azimuthally-symmetric (`m = 0`) 1D CPM. The
matrix has size ``(N_x P) \\times (N_x P)`` with `P = Lp+1` Legendre moments per region and is
indexed so that the flux relation reads ``\\vec{\\phi} = P_{vv}\\,\\vec{Q}`` (before the
boundary and scattering closures). The un-normalized probability is recovered from the
reduced one through ``p^{(p,q)}_{i,j} = \\frac{2\\ell_q+1}{\\Delta x_i}\\tilde{p}^{(p,q)}_{i,j}``
(the `(2\\ell_q+1)` factor is that of the source moment order `q`, from the source expansion).

# Input Argument(s)
- `Σ::Vector{Float64}` : total cross-section per region.
- `Δx::Vector{Float64}` : width per region.
- `Lp::Int64` : maximum Legendre degree of the volume flux expansion.

# Output Argument(s)
- `Pvv::Matrix{Float64}` : collision-probability matrix, size `(Nx*(Lp+1), Nx*(Lp+1))`.

"""
function cpm_collision_matrix(Σ::Vector{Float64},Δx::Vector{Float64},Lp::Int64)
    Nx = length(Σ)
    P = Lp+1
    Pvv = zeros(Nx*P,Nx*P)
    for i in range(1,Nx), p in range(1,P), j in range(1,Nx), q in range(1,P)
        ℓp = p-1; ℓq = q-1
        row = (i-1)*P + p
        col = (j-1)*P + q
        # The source expansion Q(Ω) = Σ_q (2ℓ_q+1)/4π R_q(Ω) Q^{(q)} carries the (2ℓ_q+1) factor
        # of the *source* moment order q (not the flux moment order p); this is what makes the
        # anisotropic (p ≠ q) couplings and the particle balance correct.
        Pvv[row,col] = (2ℓq+1)/Δx[i] * cpm_reduced_collision_probability(ℓp,ℓq,i,j,Σ,Δx)
    end
    return Pvv
end
