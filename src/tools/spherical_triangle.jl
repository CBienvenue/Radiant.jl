"""
    barycentric_inverse(P, P1, P2, P3)

Solve `(1 - ξ - η)·P1 + ξ·P2 + η·P3 = s·P` for the barycentric coordinates `(ξ, η)`
on the planar triangle through `(P1, P2, P3)` and the scaling factor `s`.

# Input Argument(s)
- `P::AbstractVector{<:Real}`: target direction (3-vector).
- `P1, P2, P3::AbstractVector{<:Real}`: vertices of the planar triangle (3-vectors).

# Output Argument(s)
- `(ξ, η, s)::NTuple{3,Float64}`.
"""
function barycentric_inverse(P::AbstractVector{<:Real}, P1::AbstractVector{<:Real}, P2::AbstractVector{<:Real}, P3::AbstractVector{<:Real})
    a1 = P2[1] - P1[1]; a2 = P2[2] - P1[2]; a3 = P2[3] - P1[3]
    b1 = P3[1] - P1[1]; b2 = P3[2] - P1[2]; b3 = P3[3] - P1[3]
    M = @inbounds [a1 b1 -P[1]; a2 b2 -P[2]; a3 b3 -P[3]]
    rhs = @inbounds [-P1[1], -P1[2], -P1[3]]
    sol = M \ rhs
    return sol[1], sol[2], sol[3]
end

"""
    jacobian_from_affine(P1, P2, P3, ξ, η)

Surface-area Jacobian of the map `(ξ, η) → V(ξ,η)/‖V(ξ,η)‖` where
`V(ξ,η) = (1-ξ-η)·P1 + ξ·P2 + η·P3`.

Returns `|∂Q/∂ξ × ∂Q/∂η|` with `Q = V/‖V‖`. Used to convert the planar
barycentric parameterization into a spherical surface measure.

# Input Argument(s)
- `P1, P2, P3::AbstractVector{<:Real}`: vertices on the unit sphere.
- `ξ, η::Real`: barycentric coordinates.

# Output Argument(s)
- `J::Float64`: spherical surface Jacobian at `(ξ, η)`.
"""
function jacobian_from_affine(P1::AbstractVector{<:Real}, P2::AbstractVector{<:Real}, P3::AbstractVector{<:Real}, ξ::Real, η::Real)
    a = (P2[1] - P1[1], P2[2] - P1[2], P2[3] - P1[3])
    b = (P3[1] - P1[1], P3[2] - P1[2], P3[3] - P1[3])
    V = ((1 - ξ - η)*P1[1] + ξ*P2[1] + η*P3[1],
         (1 - ξ - η)*P1[2] + ξ*P2[2] + η*P3[2],
         (1 - ξ - η)*P1[3] + ξ*P2[3] + η*P3[3])
    r2 = V[1]*V[1] + V[2]*V[2] + V[3]*V[3]
    if r2 == 0.0 error("Affine combination is zero (P1+P2+P3 = 0).") end
    axb = (a[2]*b[3] - a[3]*b[2], a[3]*b[1] - a[1]*b[3], a[1]*b[2] - a[2]*b[1])
    Vxa = (V[2]*a[3] - V[3]*a[2], V[3]*a[1] - V[1]*a[3], V[1]*a[2] - V[2]*a[1])
    Vxb = (V[2]*b[3] - V[3]*b[2], V[3]*b[1] - V[1]*b[3], V[1]*b[2] - V[2]*b[1])
    Va  = V[1]*a[1] + V[2]*a[2] + V[3]*a[3]
    Vb  = V[1]*b[1] + V[2]*b[2] + V[3]*b[3]
    n1 = r2*r2*axb[1] + r2*(Vb*Vxa[1] - Va*Vxb[1])
    n2 = r2*r2*axb[2] + r2*(Vb*Vxa[2] - Va*Vxb[2])
    n3 = r2*r2*axb[3] + r2*(Vb*Vxa[3] - Va*Vxb[3])
    return sqrt(n1*n1 + n2*n2 + n3*n3) / (r2*r2*r2)
end

"""
    spherical_triangle_area(A, B, C)

Area of the spherical triangle with unit-vector vertices `(A, B, C)` on the unit
sphere, using l'Huilier's formula. Equivalent to the spherical excess.

# Input Argument(s)
- `A, B, C::AbstractVector{<:Real}`: unit-vector vertices.

# Output Argument(s)
- `area::Float64`: area in steradians.
"""
function spherical_triangle_area(A::AbstractVector{<:Real}, B::AbstractVector{<:Real}, C::AbstractVector{<:Real})
    dot_BC = clamp(B[1]*C[1] + B[2]*C[2] + B[3]*C[3], -1.0, 1.0)
    dot_AC = clamp(A[1]*C[1] + A[2]*C[2] + A[3]*C[3], -1.0, 1.0)
    dot_AB = clamp(A[1]*B[1] + A[2]*B[2] + A[3]*B[3], -1.0, 1.0)
    a = acos(dot_BC)
    b = acos(dot_AC)
    c = acos(dot_AB)
    s = (a + b + c) / 2
    t = tan(s/2) * tan((s-a)/2) * tan((s-b)/2) * tan((s-c)/2)
    t = max(t, 0.0)
    return 4 * atan(sqrt(t))
end

"""
    muphi_to_muprime_phiprime(μ, ϕ, A, B, C; tol=1e-10)

Convert spherical coordinates `(μ, ϕ)` of a point `P` lying inside the spherical
triangle `K = (A, B, C)` to the coordinates `(μ', ϕ')` of the corresponding point
on the reference octant triangle `T₀ = ((1,0,0), (0,1,0), (0,0,1))` sharing the
same barycentric coordinates `(ξ, η)`.

# Input Argument(s)
- `μ::Real`, `ϕ::Real`: spherical coordinates on `K`.
- `A, B, C::AbstractVector{<:Real}`: unit-vector vertices of `K`.

# Output Argument(s)
- `(μ', ϕ', ξ, η)::NTuple{4,Float64}`.
"""
function muphi_to_muprime_phiprime(μ::Real, ϕ::Real, A::AbstractVector{<:Real}, B::AbstractVector{<:Real}, C::AbstractVector{<:Real}; tol::Float64=1e-10)
    r_perp = sqrt(max(0.0, 1 - μ^2))
    P = (μ, r_perp * cos(ϕ), r_perp * sin(ϕ))
    ξ, η, _ = barycentric_inverse([P[1], P[2], P[3]], A, B, C)
    if ξ < -tol || η < -tol || ξ + η > 1 + tol
        error("(μ,ϕ) is outside spherical triangle K (ξ=$ξ, η=$η).")
    end
    ξc = clamp(ξ, 0.0, 1.0)
    ηc = clamp(η, 0.0, 1.0 - ξc)
    Vt = (1 - ξc - ηc, ξc, ηc)
    r = sqrt(Vt[1]^2 + Vt[2]^2 + Vt[3]^2)
    Qx = Vt[1] / r; Qy = Vt[2] / r; Qz = Vt[3] / r
    μp = Qx
    ϕp = atan(Qz, Qy)
    if ϕp < 0 ϕp += 2π end
    return μp, ϕp, ξc, ηc
end

"""
    muprime_phiprime_to_mu_phi(μp, ϕp, A, B, C; tol=1e-10)

Inverse of `muphi_to_muprime_phiprime`: given `(μ', ϕ')` on the reference
octant triangle `T₀`, return the `(μ, ϕ)` of the corresponding point on the
spherical triangle `K = (A, B, C)`.

# Output Argument(s)
- `(μ, ϕ, ξ, η)::NTuple{4,Float64}`.
"""
function muprime_phiprime_to_mu_phi(μp::Real, ϕp::Real, A::AbstractVector{<:Real}, B::AbstractVector{<:Real}, C::AbstractVector{<:Real}; tol::Float64=1e-10)
    r_perp = sqrt(max(0.0, 1 - μp^2))
    Q = [μp, r_perp * cos(ϕp), r_perp * sin(ϕp)]
    e1 = [1.0, 0.0, 0.0]; e2 = [0.0, 1.0, 0.0]; e3 = [0.0, 0.0, 1.0]
    ξ, η, _ = barycentric_inverse(Q, e1, e2, e3)
    if ξ < -tol || η < -tol || ξ + η > 1 + tol
        error("(μ',ϕ') is outside reference triangle T₀ (ξ=$ξ, η=$η).")
    end
    ξc = clamp(ξ, 0.0, 1.0)
    ηc = clamp(η, 0.0, 1.0 - ξc)
    Px = (1 - ξc - ηc)*A[1] + ξc*B[1] + ηc*C[1]
    Py = (1 - ξc - ηc)*A[2] + ξc*B[2] + ηc*C[2]
    Pz = (1 - ξc - ηc)*A[3] + ξc*B[3] + ηc*C[3]
    r = sqrt(Px*Px + Py*Py + Pz*Pz)
    Px /= r; Py /= r; Pz /= r
    μ = Px
    ϕ = atan(Pz, Py)
    if ϕ < 0 ϕ += 2π end
    return μ, ϕ, ξc, ηc
end

"""
    octant_subtriangle_vertices(u, i, j, Nv)

Vertices `(A, B, C)` on the unit sphere of the sub-triangle indexed by `(i, j)`
inside octant `u` of the barycentric subdivision of order `Nv`.

The reference octant triangle is `(e1, e2, e3) = ((1,0,0), (0,1,0), (0,0,1))`
in octant `u = 1`. For other octants the signs are applied component-wise.
Row `i ∈ 1:Nv` contains `2i-1` sub-triangles: odd `j` → upward triangle,
even `j` → downward triangle. Total: `Nv²` sub-triangles per octant.

# Input Argument(s)
- `u::Int64`: octant index in `1:8`.
- `i::Int64`: row index in `1:Nv`.
- `j::Int64`: position index in `1:(2i-1)`.
- `Nv::Int64`: barycentric subdivision order.

# Output Argument(s)
- `(A, B, C)::NTuple{3,Vector{Float64}}`: unit-vector vertices of the sub-triangle.
"""
function octant_subtriangle_vertices(u::Int64, i::Int64, j::Int64, Nv::Int64)
    if ~(1 ≤ u ≤ 8) error("Octant index u must be in 1:8.") end
    if ~(1 ≤ i ≤ Nv) error("Row index i must be in 1:Nv.") end
    if ~(1 ≤ j ≤ 2i - 1) error("Position index j must be in 1:(2i-1).") end
    sx = (u ∈ (1,2,3,4)) ? 1 : -1
    sy = (u ∈ (1,2,5,6)) ? 1 : -1
    sz = (u ∈ (1,3,5,7)) ? 1 : -1
    if isodd(j)
        k = (j + 1) ÷ 2
        a1, b1, c1 = Nv - i + 1, i - k,     k - 1
        a2, b2, c2 = Nv - i,     i - k + 1, k - 1
        a3, b3, c3 = Nv - i,     i - k,     k
    else
        k = j ÷ 2
        a1, b1, c1 = Nv - i + 1, i - k,     k - 1
        a2, b2, c2 = Nv - i + 1, i - k - 1, k
        a3, b3, c3 = Nv - i,     i - k,     k
    end
    A = _normalize_signed(sx, sy, sz, a1, b1, c1)
    B = _normalize_signed(sx, sy, sz, a2, b2, c2)
    C = _normalize_signed(sx, sy, sz, a3, b3, c3)
    return A, B, C
end

function _normalize_signed(sx::Int, sy::Int, sz::Int, a::Int, b::Int, c::Int)
    r = sqrt(float(a*a + b*b + c*c))
    return [sx * a / r, sy * b / r, sz * c / r]
end

"""
    symmetric_patch_orthonormal_basis(Lq, Nv, u, i, j, x_gl, weight_gl)

For sub-triangle `(u, i, j)` of the symmetric tiling, build an *exactly orthonormal*
basis (via Gram-Schmidt / Cholesky) on the curved spherical triangle. Returns the
basis values at the 32-point Duffy-chart Gauss-Legendre quadrature nodes, together
with the direction cosines `(Px, Py, Pz)` and the spherical-area quadrature weight
`JW = J_K · w_Duffy` at each node.

The orthonormal basis `Ψ_orth` satisfies `Σ_k Ψ_orth[p,k] · Ψ_orth[q,k] · JW[k] = δ_{pq}`,
which restores the assumption `mass matrix = I` used by the Galerkin sweep routines.

# Input Argument(s)
- `Lq::Int64`: patch-local Legendre order.
- `Nv::Int64`: barycentric subdivision order.
- `u, i, j::Int64`: patch indices.
- `x_gl, weight_gl::Vector{Float64}`: 1-D Gauss-Legendre nodes and weights on `[-1,1]`.

# Output Argument(s)
- `Ψ_orth::Matrix{Float64}`: orthonormal basis of size `(Nq, length(x_gl)^2)`.
- `Px, Py, Pz::Vector{Float64}`: direction cosines at quadrature nodes.
- `JW::Vector{Float64}`: spherical-area quadrature weight at each node.
"""
function symmetric_patch_orthonormal_basis(Lq::Int64, Nv::Int64, u::Int64, i::Int64, j::Int64, x_gl::Vector{Float64}, weight_gl::Vector{Float64})
    A, B, C = octant_subtriangle_vertices(u, i, j, Nv)
    Nq = (Lq + 1)^2
    N = length(x_gl)
    Nquad = N * N

    Ψ_old = Matrix{Float64}(undef, Nq, Nquad)
    Px = Vector{Float64}(undef, Nquad)
    Py = Vector{Float64}(undef, Nquad)
    Pz = Vector{Float64}(undef, Nquad)
    JW = Vector{Float64}(undef, Nquad)

    let k = 0
        for m in 1:N, n in 1:N
            k += 1
            α = x_gl[n]; β = x_gl[m]
            ξ = (α + 1.0)/4.0 * (1.0 - β)
            η = (β + 1.0)/2.0
            w_k = weight_gl[n] * weight_gl[m] * (1.0 - β)/8.0

            Vx = (1 - ξ - η)*A[1] + ξ*B[1] + η*C[1]
            Vy = (1 - ξ - η)*A[2] + ξ*B[2] + η*C[2]
            Vz = (1 - ξ - η)*A[3] + ξ*B[3] + η*C[3]
            r = sqrt(Vx*Vx + Vy*Vy + Vz*Vz)
            Px[k] = Vx/r; Py[k] = Vy/r; Pz[k] = Vz/r

            # Full-range real spherical harmonics evaluated at the global (μ, ϕ)
            μ_g = Px[k]
            ϕ_g = atan(Pz[k], Py[k])
            if ϕ_g < 0.0 ϕ_g += 2π end
            Yf = real_spherical_harmonics_up_to_L(Lq, μ_g, ϕ_g)
            q = 0
            for lq in 0:Lq, mq in -lq:lq
                q += 1
                Ψ_old[q, k] = Yf[lq+1][lq+mq+1]
            end

            J_K = jacobian_from_affine(A, B, C, ξ, η)
            JW[k] = J_K * w_k
        end
    end

    # Mass matrix: M[p,q] = Σ_k Ψ_old[p,k] · Ψ_old[q,k] · JW[k]
    M_mass = zeros(Nq, Nq)
    for k in 1:Nquad
        wk = JW[k]
        for q in 1:Nq
            Ψqk = Ψ_old[q, k]
            for p in 1:Nq
                M_mass[p, q] += Ψ_old[p, k] * Ψqk * wk
            end
        end
    end
    # Symmetrize for robustness against floating-point roundoff
    @inbounds for q in 1:Nq, p in 1:q
        avg = 0.5 * (M_mass[p, q] + M_mass[q, p])
        M_mass[p, q] = avg
        M_mass[q, p] = avg
    end

    # Löwdin (symmetric) orthogonalization: M^{-1/2} = V · Λ^{-1/2} · V^T.
    # Symmetric under any orthogonal transformation U: if M' = U M U^T then
    # M'^{-1/2} = U M^{-1/2} U^T. Combined with the full-range basis (which
    # transforms by U under reflections), this restores exact x↔y symmetry.
    E = LinearAlgebra.eigen(LinearAlgebra.Symmetric(M_mass))
    sqrt_inv_vals = 1.0 ./ sqrt.(max.(E.values, eps(Float64)))
    M_inv_sqrt = E.vectors * LinearAlgebra.Diagonal(sqrt_inv_vals) * LinearAlgebra.transpose(E.vectors)
    Ψ_orth = M_inv_sqrt * Ψ_old

    return Ψ_orth, Px, Py, Pz, JW
end
