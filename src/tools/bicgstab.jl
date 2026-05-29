"""
    bicgstab!(z,matvec!,c;maxit,tol)

Matrix-free BiCGStab solving the linear system `M z = c`, where `M` is given only through the
in-place product `matvec!(out,v)` (`out ← M·v`) and all vectors are Krylov states (`KState`, see
[krylov_state.jl](krylov_state.jl)). A low-memory alternative to [gmres.jl](gmres.jl) for the
in-group discrete-ordinates iteration (`M = I − A`): short recurrences, no restart, a fixed set
of work vectors, two `matvec!` (i.e. two transport sweeps) per iteration.

Convergence is measured relative to `‖c‖`. Breakdowns (`ρ ≈ 0` or `ω ≈ 0`) terminate the solve
and are reported through `converged = false`.

# Input Argument(s)
- `z::KState` : initial guess, overwritten with the solution.
- `matvec!` : function `matvec!(out::KState,v::KState)` computing `out ← M·v`.
- `c::KState` : right-hand side.
- `maxit::Int64` : maximum number of BiCGStab iterations.
- `tol::Float64` : relative residual tolerance.

# Output Argument(s)
- `iter::Int64` : number of `matvec!` applications performed.
- `resid::Float64` : final relative residual `‖c − M z‖ / ‖c‖`.
- `converged::Bool` : whether `resid < tol` was reached.
- `ρ::Float64` : estimated spectral radius (see [`spectral_radius_estimate`](@ref)).

# Reference(s)
- van der Vorst (1992), Bi-CGSTAB: A Fast and Smoothly Converging Variant of Bi-CG for the
  Solution of Nonsymmetric Linear Systems.
"""
function bicgstab!(z::KState,matvec!,c::KState;maxit::Int64,tol::Float64)

    bnorm = state_norm(c)
    if bnorm == 0.0
        state_zero!(z)
        return 0, 0.0, true, NaN
    end

    # Work vectors
    r  = state_similar(z)
    r̂0 = state_similar(z)
    p  = state_similar(z)
    v  = state_similar(z)
    s  = state_similar(z)
    t  = state_similar(z)

    # Initial residual r = c - M z
    matvec!(r,z)
    state_scale!(-1.0,r); state_axpy!(1.0,c,r)
    iter = 1
    resid = state_norm(r) / bnorm
    resid0 = resid
    if resid < tol return iter, resid, true, NaN end

    state_copy!(r̂0,r)
    ρ_old = 1.0; α = 1.0; ω = 1.0
    state_zero!(p); state_zero!(v)
    converged = false
    tiny = 1e-300

    while iter < maxit
        ρ_new = state_dot(r̂0,r)
        if abs(ρ_new) < tiny break end                       # breakdown

        β = (ρ_new/ρ_old) * (α/ω)
        # p = r + β (p - ω v)
        state_axpy!(-ω,v,p)        # p ← p - ω v
        state_scale!(β,p)          # p ← β (p - ω v)
        state_axpy!(1.0,r,p)       # p ← r + β (p - ω v)

        matvec!(v,p); iter += 1
        denom = state_dot(r̂0,v)
        if abs(denom) < tiny break end
        α = ρ_new / denom

        # s = r - α v
        state_copy!(s,r); state_axpy!(-α,v,s)
        if state_norm(s)/bnorm < tol
            state_axpy!(α,p,z)     # z ← z + α p
            resid = state_norm(s)/bnorm
            converged = true
            break
        end

        matvec!(t,s); iter += 1
        tt = state_dot(t,t)
        if tt < tiny break end
        ω = state_dot(t,s) / tt

        # z = z + α p + ω s
        state_axpy!(α,p,z); state_axpy!(ω,s,z)
        # r = s - ω t
        state_copy!(r,s); state_axpy!(-ω,t,r)

        resid = state_norm(r)/bnorm
        if resid < tol converged = true; break end
        if abs(ω) < tiny break end                           # breakdown
        ρ_old = ρ_new
    end

    return iter, resid, converged, spectral_radius_estimate(resid0,resid,iter)
end
