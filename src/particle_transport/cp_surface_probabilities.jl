"""
    cp_legendre_polynomial_coefficients(ℓ::Int64,::Type{T}=Float64) where {T<:AbstractFloat}

Monomial coefficients of the Legendre polynomial ``P_\\ell(\\mu) = \\sum_r c_r \\mu^r`` (the
returned vector has `c_r` at index `r+1`), in the requested floating type `T`. Binomials are
formed exactly as `BigInt` before conversion, so no overflow occurs at high order.
"""
function cp_legendre_polynomial_coefficients(ℓ::Int64,::Type{T}=Float64) where {T<:AbstractFloat}
    c = zeros(T,ℓ+1)
    for k in range(0,ℓ÷2)
        c[ℓ-2k+1] += (-1)^k * T(binomial(big(ℓ),big(k))) * T(binomial(big(2ℓ-2k),big(ℓ))) / T(2)^ℓ
    end
    return c
end

"""
    cp_shifted_legendre_polynomial_coefficients(ℓ::Int64,::Type{T}=Float64) where {T<:AbstractFloat}

Monomial coefficients of the shifted Legendre polynomial ``P_\\ell(2\\mu-1) = \\sum_r c_r \\mu^r``,
in the requested floating type `T`. The half-range surface basis of Radiant is
``\\bar{R}_\\ell(\\hat\\mu) = \\sqrt{2\\ell+1}\\,P_\\ell(2\\hat\\mu-1)``, so its coefficients are
`sqrt(2ℓ+1)` times the ones returned here.
"""
function cp_shifted_legendre_polynomial_coefficients(ℓ::Int64,::Type{T}=Float64) where {T<:AbstractFloat}
    c = zeros(T,ℓ+1)
    for k in range(0,ℓ)
        c[k+1] = (-1)^(ℓ+k) * T(binomial(big(ℓ),big(k))) * T(binomial(big(ℓ+k),big(k)))
    end
    return c
end

"""
    cp_polynomial_product(a::Vector{T},b::Vector{T}) where {T<:AbstractFloat}

Monomial coefficients of the product of two polynomials given by their monomial coefficients
(discrete convolution).
"""
function cp_polynomial_product(a::Vector{T},b::Vector{T}) where {T<:AbstractFloat}
    c = zeros(T,length(a)+length(b)-1)
    for i in eachindex(a), j in eachindex(b)
        c[i+j-1] += a[i]*b[j]
    end
    return c
end

"""
    cp_ray_reduction(c::Vector{T},Σi::T,Δxi::T,Fout::T,w::Int64) where {T<:AbstractFloat}

Closed form of the surface ray integral
``\\int_0^1 \\mu^{w}\\,[\\sum_r c_r \\mu^r]\\, e^{-F_{\\text{out}}/\\mu}\\, \\varepsilon_i(\\mu)\\,d\\mu``,
where ``\\varepsilon_i(\\mu) = \\int_0^{\\Delta x_i} e^{-\\xi\\Sigma_i/\\mu}d\\xi`` is the escape
across region `i`. Each monomial reduces to the generalized exponential integral through
``\\int_0^1 \\mu^n e^{-F/\\mu}d\\mu = E_{n+2}(F)``, giving

- ``\\Sigma_i > 0`` : ``\\Sigma_i^{-1} \\sum_r c_r \\left[E_{r+w+2}(F_{\\text{out}}) - E_{r+w+2}(F_{\\text{out}}+\\Sigma_i\\Delta x_i)\\right]``,
- ``\\Sigma_i = 0`` (void) : ``\\Delta x_i \\sum_r c_r\\, E_{r+w+1}(F_{\\text{out}})``.

`w = 1` is used for the volume-to-surface leakage (the escape carries an extra `μ`) and `w = 0`
for the surface-to-volume leakage (the `1/μ` streaming factor cancels one power). The computation
is carried out in the type `T` of the coefficients (`BigFloat` when extended precision is needed
to resolve the cancellation between the high-order monomials).
"""
function cp_ray_reduction(c::Vector{T},Σi::T,Δxi::T,Fout::T,w::Int64) where {T<:AbstractFloat}
    s = zero(T)
    if Σi > 0
        Fin = Fout + Σi*Δxi
        for r in range(0,length(c)-1)
            if c[r+1] != 0
                s += c[r+1]*(cp_exponential_integral(r+w+2,Fout) - cp_exponential_integral(r+w+2,Fin))
            end
        end
        return s/Σi
    else
        for r in range(0,length(c)-1)
            if c[r+1] != 0
                s += c[r+1]*cp_exponential_integral(r+w+1,Fout)
            end
        end
        return Δxi*s
    end
end

"""
    cp_surface_matrices(Σ::Vector{Float64},Δx::Vector{Float64},Lp::Int64,Nν::Int64)

Assemble the boundary-coupling probability matrices of the 1D CP in Radiant's half-range
surface basis, the normalized shifted Legendre polynomials
``\\bar{R}_\\ell(\\hat\\mu) = \\sqrt{2\\ell+1}\\,P_\\ell(2\\hat\\mu-1)`` on ``\\hat\\mu = |\\vec\\Omega\\cdot\\vec N| \\in [0,1]``
(the same basis used by [`half_range_legendre_polynomials_up_to_L`](@ref) for surface sources).
The two faces (`f = 1` left, `f = 2` right) each carry `Nν+1` surface moments; the surface
degrees of freedom are stacked as `[face 1 : ℓ = 0…Nν, face 2 : ℓ = 0…Nν]`.

- `Pvs` : incoming boundary flux moments → volume flux moments, size `(Nx·(Lp+1), 2(Nν+1))`.
- `Psv` : volume source moments → outgoing boundary flux moments, size `(2(Nν+1), Nx·(Lp+1))`.
- `Pss` : incoming boundary flux moments → outgoing boundary flux moments (transmission across
  the slab; the same-face block vanishes for the uncollided flux), size `(2(Nν+1), 2(Nν+1))`.

The volume moments are indexed as `(i-1)(Lp+1)+p`, matching [`cp_collision_matrix`](@ref). Each
angular integral over ``\\mu \\in [0,1]`` is evaluated in closed form as a finite combination of
generalized exponential integrals ``E_n`` (see [`cp_ray_reduction`](@ref)).

The monomial expansion of the shifted Legendre polynomials has coefficients that grow like
``4^{\\ell}``; the closed-form sum is therefore well-conditioned in `Float64` up to `Nν ≈ 8`
(all that a smooth interface flux requires), and is carried out in scaled-precision `BigFloat`
for larger orders to defeat the catastrophic cancellation. The matrices are returned in `Float64`.

# Input Argument(s)
- `Σ::Vector{Float64}` : total cross-section per region.
- `Δx::Vector{Float64}` : width per region.
- `Lp::Int64` : Legendre order of the volume flux.
- `Nν::Int64` : order of the half-range surface flux expansion.

# Output Argument(s)
- `Pvs::Matrix{Float64}`, `Psv::Matrix{Float64}`, `Pss::Matrix{Float64}`.

"""
function cp_surface_matrices(Σ::Vector{Float64},Δx::Vector{Float64},Lp::Int64,Nν::Int64)
    if Nν ≤ 8
        return cp_surface_matrices_impl(Σ,Δx,Lp,Nν,Float64)
    else
        # Extra bits ≈ log2(4^{2Nν}) = 4Nν for the worst (transmission) product, plus a margin.
        return setprecision(BigFloat, 64 + 8*Nν) do
            cp_surface_matrices_impl(Σ,Δx,Lp,Nν,BigFloat)
        end
    end
end

function cp_surface_matrices_impl(Σ::Vector{Float64},Δx::Vector{Float64},Lp::Int64,Nν::Int64,::Type{T}) where {T<:AbstractFloat}
    Nx = length(Σ)
    P = Lp+1
    Sv = Nν+1
    NS = 2*Sv
    Ftot = T(cp_optical_thickness(Σ,Δx,1,Nx))

    # Monomial coefficients of the volume (Legendre) and surface (shifted Legendre) bases in T.
    Pcoef = [cp_legendre_polynomial_coefficients(p-1,T) for p in range(1,P)]                          # P_{p-1}(μ)
    Rcoef = [sqrt(T(2*l+1)) .* cp_shifted_legendre_polynomial_coefficients(l,T) for l in range(0,Nν)] # R̄_ℓ(μ)

    vidx(i,p) = (i-1)*P + p
    sidx(f,l) = (f-1)*Sv + (l+1)   # l = 0…Nν

    Pvs = zeros(T,Nx*P,NS)
    Psv = zeros(T,NS,Nx*P)
    Pss = zeros(T,NS,NS)

    # Pvs : incoming boundary moment → volume flux moment (escape carries an extra μ, w = 1).
    for i in range(1,Nx)
        FL = T(cp_optical_thickness(Σ,Δx,1,i-1))
        FR = T(cp_optical_thickness(Σ,Δx,i+1,Nx))
        Σi = T(Σ[i]); Δxi = T(Δx[i])
        for p in range(1,P), l in range(0,Nν)
            base = cp_polynomial_product(Pcoef[p],Rcoef[l+1])       # P_p(μ) R̄_ℓ(μ)
            Pvs[vidx(i,p),sidx(1,l)] = (2π/Δx[i])*cp_ray_reduction(base,Σi,Δxi,FL,1)
            # Right face : the incoming ray travels toward -x, P_p(-μ) = (-1)^{ℓp} P_p(μ).
            Pvs[vidx(i,p),sidx(2,l)] = (-1)^(p-1)*(2π/Δx[i])*cp_ray_reduction(base,Σi,Δxi,FR,1)
        end
    end

    # Psv : volume source moment → outgoing boundary moment (1/μ streaming cancels one power, w = 0).
    for j in range(1,Nx)
        FL = T(cp_optical_thickness(Σ,Δx,1,j-1))
        FR = T(cp_optical_thickness(Σ,Δx,j+1,Nx))
        Σj = T(Σ[j]); Δxj = T(Δx[j])
        for q in range(1,P), l in range(0,Nν)
            base = cp_polynomial_product(Rcoef[l+1],Pcoef[q])       # R̄_ℓ(μ) P_q(μ)
            # Left face : outgoing toward -x uses P_q(-μ) = (-1)^{ℓq} P_q(μ) ; right : P_q(μ).
            Psv[sidx(1,l),vidx(j,q)] = (2(q-1)+1)/(4π)*(-1)^(q-1)*cp_ray_reduction(base,Σj,Δxj,FL,0)
            Psv[sidx(2,l),vidx(j,q)] = (2(q-1)+1)/(4π)*cp_ray_reduction(base,Σj,Δxj,FR,0)
        end
    end

    # Pss : transmission across the whole slab (same-face block is zero for the uncollided flux).
    for l in range(0,Nν), lp in range(0,Nν)
        g = cp_polynomial_product(Rcoef[l+1],Rcoef[lp+1])
        t = zero(T)
        for r in range(0,length(g)-1)
            if g[r+1] != 0 t += g[r+1]*cp_exponential_integral(r+2,Ftot) end
        end
        Pss[sidx(2,l),sidx(1,lp)] = t   # left incoming → right outgoing
        Pss[sidx(1,l),sidx(2,lp)] = t   # right incoming → left outgoing
    end

    return Matrix{Float64}(Pvs), Matrix{Float64}(Psv), Matrix{Float64}(Pss)
end
