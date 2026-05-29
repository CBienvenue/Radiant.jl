"""
    anderson!(z,fixedpoint!;depth,maxit,tol,β)

Depth-`m` Anderson acceleration of the in-group fixed-point iteration `z = g(z)`, where the map
`g` (one full source-iteration pass `T`, fixed sources included) is supplied through the in-place
call `fixedpoint!(out,z)` (`out ← g(z)`) and all vectors are Krylov states (`KState`, see
[krylov_state.jl](krylov_state.jl)).

Anderson mixing replaces the plain update `z ← g(z)` with a least-squares combination of the last
`m` iterates that minimizes the linearized residual `f(z) = g(z) − z`; the depth-1 case recovers
the spirit of the two-point [livolant.jl](livolant.jl) extrapolation, while larger depths are
markedly more robust on scattering-dominated, optically-thick problems. Because it accelerates
the exact fixed-point sequence, it applies unchanged to every boundary condition and physics mode
(BTE / BFP / CSD / Fokker-Planck / electromagnetic).

The small `mₖ×mₖ` least-squares system is solved through its (Tikhonov-regularized) normal
equations with an in-place Gaussian elimination ([`anderson_solve`](@ref)), avoiding any heap
allocation from a factorization object.

# Input Argument(s)
- `z::KState` : initial guess, overwritten with the accelerated solution.
- `fixedpoint!` : function `fixedpoint!(out::KState,z::KState)` computing `out ← g(z)`.
- `depth::Int64` : Anderson memory `m` (number of stored residual differences).
- `maxit::Int64` : maximum number of `fixedpoint!` applications.
- `tol::Float64` : relative residual tolerance on `‖g(z) − z‖`.
- `β::Float64` : mixing/damping parameter (`1.0` = no damping).

# Output Argument(s)
- `iter::Int64` : number of `fixedpoint!` applications performed.
- `resid::Float64` : final relative residual.
- `converged::Bool` : whether `resid < tol` was reached.
- `ρ::Float64` : estimated spectral radius (see [`spectral_radius_estimate`](@ref)).

# Reference(s)
- Anderson (1965), Iterative Procedures for Nonlinear Integral Equations.
- Walker and Ni (2011), Anderson Acceleration for Fixed-Point Iterations.
"""
function anderson!(z::KState,fixedpoint!;depth::Int64,maxit::Int64,tol::Float64,β::Float64=1.0)

    m = max(1,depth)
    g     = state_similar(z)
    f     = state_similar(z)
    f_old = state_similar(z)
    g_old = state_similar(z)
    ΔF = KState[]
    ΔG = KState[]

    iter = 0
    resid = 1.0
    resid0 = NaN
    converged = false

    for it in range(1,maxit)

        fixedpoint!(g,z); iter += 1

        # Residual f = g(z) - z
        state_copy!(f,g); state_axpy!(-1.0,z,f)
        fnorm = state_norm(f)
        scale = max(state_norm(z),state_norm(g),1e-16)
        resid = fnorm/scale
        if isnan(resid0) resid0 = resid end
        if resid < tol converged = true; break end

        if it == 1
            # Plain damped fixed-point update: z ← z + β f = (1-β) z + β g(z)
            state_axpy!(β,f,z)
        else
            # Append newest residual / iterate differences to the history
            dF = state_clone(f); state_axpy!(-1.0,f_old,dF)
            dG = state_clone(g); state_axpy!(-1.0,g_old,dG)
            push!(ΔF,dF); push!(ΔG,dG)
            if length(ΔF) > m popfirst!(ΔF); popfirst!(ΔG) end
            mk = length(ΔF)

            # Normal equations (ΔFᵀΔF) γ = ΔFᵀ f, with Tikhonov regularization
            Nmat = zeros(mk,mk)
            bb = zeros(mk)
            for i in range(1,mk)
                bb[i] = state_dot(ΔF[i],f)
                for j in range(1,mk)
                    Nmat[i,j] = state_dot(ΔF[i],ΔF[j])
                end
            end
            tr = 0.0; for i in range(1,mk) tr += Nmat[i,i] end
            λ = 1e-12 * (tr/mk + 1e-300)
            for i in range(1,mk) Nmat[i,i] += λ end
            γ = anderson_solve(Nmat,bb)

            # z ← z + β f - Σ γ_i (ΔG_i + (β-1) ΔF_i), the damped Anderson update written in
            # terms of the stored g-differences ΔG (= Δx + Δf): for β = 1 this is z + f - ΔG γ.
            state_axpy!(β,f,z)
            for i in range(1,mk)
                state_axpy!(-γ[i],ΔG[i],z)
                state_axpy!(-(β-1.0)*γ[i],ΔF[i],z)
            end
        end

        state_copy!(f_old,f); state_copy!(g_old,g)
    end

    return iter, resid, converged, spectral_radius_estimate(resid0,resid,iter)
end

"""
    anderson_solve(A::Matrix{Float64},b::Vector{Float64})

Solve the small dense system `A x = b` by Gaussian elimination with partial pivoting, returning
`x`. Used for the Anderson mixing coefficients; kept allocation-light (no factorization object)
to suit the tiny systems involved.
"""
function anderson_solve(A::Matrix{Float64},b::Vector{Float64})
    n = length(b)
    M = copy(A)
    x = copy(b)
    for k in range(1,n)
        # Partial pivoting
        p = k; mx = abs(M[k,k])
        for i in range(k+1,n)
            if abs(M[i,k]) > mx mx = abs(M[i,k]); p = i end
        end
        if p != k
            for j in range(k,n) M[k,j],M[p,j] = M[p,j],M[k,j] end
            x[k],x[p] = x[p],x[k]
        end
        piv = M[k,k]
        if piv == 0.0 continue end
        for i in range(k+1,n)
            factor = M[i,k]/piv
            for j in range(k,n) M[i,j] -= factor*M[k,j] end
            x[i] -= factor*x[k]
        end
    end
    for k in n:-1:1
        s = x[k]
        for j in range(k+1,n) s -= M[k,j]*x[j] end
        x[k] = M[k,k] != 0.0 ? s/M[k,k] : 0.0
    end
    return x
end
