"""
    gmres!(z,matvec!,c;restart,maxit,tol)

Restarted, matrix-free GMRES(`restart`) solving the linear system `M z = c`, where `M` is given
only through the in-place product `matvec!(out,v)` (`out ← M·v`) and all vectors are Krylov
states (`KState`, see [krylov_state.jl](krylov_state.jl)). Used to accelerate the in-group
discrete-ordinates iteration, where `M = I − A` (`A` = scattering + transport sweep) and one
`matvec!` costs one transport sweep over all ordinates.

Standard Arnoldi process with modified Gram-Schmidt orthogonalization and incremental Givens
rotations on the Hessenberg matrix; the running residual norm is read from the rotated
right-hand side, so no extra sweep is needed to monitor convergence. Convergence is measured
relative to `‖c‖`.

# Input Argument(s)
- `z::KState` : initial guess, overwritten with the solution.
- `matvec!` : function `matvec!(out::KState,v::KState)` computing `out ← M·v`.
- `c::KState` : right-hand side.
- `restart::Int64` : Krylov subspace size before restart.
- `maxit::Int64` : maximum number of `matvec!` applications.
- `tol::Float64` : relative residual tolerance.

# Output Argument(s)
- `iter::Int64` : number of `matvec!` applications performed.
- `resid::Float64` : final relative residual `‖c − M z‖ / ‖c‖`.
- `converged::Bool` : whether `resid < tol` was reached.
- `ρ::Float64` : estimated spectral radius (see [`spectral_radius_estimate`](@ref)).

# Reference(s)
- Saad and Schultz (1986), GMRES: A Generalized Minimal Residual Algorithm for Solving
  Nonsymmetric Linear Systems.
- Warsa, Wareing and Morel (2004), Krylov Iterative Methods and the Degraded Effectiveness 
  of Diffusion Synthetic Acceleration for Multidimensional SN Calculations in Problems with Material Discontinuities.
"""
function gmres!(z::KState,matvec!,c::KState;restart::Int64,maxit::Int64,tol::Float64)

    m = max(1,restart)
    bnorm = state_norm(c)
    if bnorm == 0.0
        state_zero!(z)
        return 0, 0.0, true, NaN
    end

    # Preallocated Krylov basis and work buffers
    V = [state_similar(z) for _ in range(1,m+1)]
    w = state_similar(z)
    r = state_similar(z)
    H = zeros(m+1,m)
    cs = zeros(m)
    sn = zeros(m)
    g  = zeros(m+1)

    iter = 0
    resid = 1.0
    resid0 = NaN
    converged = false

    while iter < maxit && ~converged

        # Initial residual r = c - M z
        matvec!(r,z); iter += 1
        state_scale!(-1.0,r); state_axpy!(1.0,c,r)
        β = state_norm(r)
        resid = β / bnorm
        if isnan(resid0) resid0 = resid end
        if resid < tol converged = true; break end

        # First Arnoldi vector
        state_copy!(V[1],r); state_scale!(1.0/β,V[1])
        fill!(g,0.0); g[1] = β

        k = 0
        for j in range(1,m)
            if iter ≥ maxit break end
            matvec!(w,V[j]); iter += 1

            # Modified Gram-Schmidt orthogonalization
            for i in range(1,j)
                H[i,j] = state_dot(w,V[i])
                state_axpy!(-H[i,j],V[i],w)
            end
            hjp1 = state_norm(w)            # subdiagonal entry H[j+1,j] before rotation
            H[j+1,j] = hjp1
            if hjp1 > 0.0
                state_copy!(V[j+1],w); state_scale!(1.0/hjp1,V[j+1])
            end

            # Apply previous Givens rotations to the new column
            for i in range(1,j-1)
                temp     =  cs[i]*H[i,j] + sn[i]*H[i+1,j]
                H[i+1,j] = -sn[i]*H[i,j] + cs[i]*H[i+1,j]
                H[i,j]   =  temp
            end

            # New Givens rotation eliminating H[j+1,j]
            denom = hypot(H[j,j],H[j+1,j])
            if denom == 0.0
                cs[j] = 1.0; sn[j] = 0.0
            else
                cs[j] = H[j,j]/denom
                sn[j] = H[j+1,j]/denom
            end
            H[j,j]   = cs[j]*H[j,j] + sn[j]*H[j+1,j]
            H[j+1,j] = 0.0
            g[j+1]   = -sn[j]*g[j]
            g[j]     =  cs[j]*g[j]

            k = j
            resid = abs(g[j+1]) / bnorm
            if resid < tol converged = true; break end
            if hjp1 == 0.0 break end        # happy breakdown: Krylov space exhausted
        end

        # Solve the k×k upper-triangular least-squares system H y = g
        y = zeros(k)
        for i in k:-1:1
            s = g[i]
            for l in range(i+1,k)
                s -= H[i,l]*y[l]
            end
            y[i] = H[i,i] != 0.0 ? s/H[i,i] : 0.0
        end

        # Update the solution z += Σ y[i] V[i]
        for i in range(1,k)
            state_axpy!(y[i],V[i],z)
        end

        if converged break end
    end

    return iter, resid, converged, spectral_radius_estimate(resid0,resid,iter)
end
