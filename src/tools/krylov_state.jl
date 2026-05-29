"""
    krylov_state.jl

Vector-space algebra over a *state* — a `Vector{Array{Float64}}` gathering the unknowns of the
in-group source iteration into a single abstract vector. The first entry is always the Legendre
flux moments `𝚽l` (shape `(Np,Nm[5],Ns[1],Ns[2],Ns[3])`); the remaining entries are the incoming
boundary angular fluxes (`𝚽x12_in`, and in 2D/3D `𝚽y12_in`, `𝚽z12_in`), which are part of the
fixed-point state whenever a face is reflective or periodic.

These helpers let the hand-rolled Krylov ([gmres.jl](gmres.jl), [bicgstab.jl](bicgstab.jl)) and
Anderson ([anderson.jl](anderson.jl)) solvers operate directly on the heterogeneously-shaped
arrays — without flattening to a single `Vector{Float64}` — by iterating over each array with
`eachindex`. All mutating operations are in-place to keep allocations out of the iteration loops,
matching the rest of the transport solver.
"""

# Type alias for a Krylov state vector.
const KState = Vector{Array{Float64}}

"""
    state_dot(a::KState,b::KState)

Euclidean inner product `⟨a,b⟩` summed over every component array.
"""
function state_dot(a::KState,b::KState)
    s = 0.0
    @inbounds for k in eachindex(a,b)
        ak = a[k]; bk = b[k]
        for i in eachindex(ak,bk)
            s += ak[i] * bk[i]
        end
    end
    return s
end

"""
    state_norm(a::KState)

Euclidean norm `√⟨a,a⟩`.
"""
state_norm(a::KState) = sqrt(state_dot(a,a))

"""
    state_axpy!(α::Float64,x::KState,y::KState)

In-place `y ← y + α·x`. Returns `y`.
"""
function state_axpy!(α::Float64,x::KState,y::KState)
    @inbounds for k in eachindex(x,y)
        xk = x[k]; yk = y[k]
        for i in eachindex(xk,yk)
            yk[i] += α * xk[i]
        end
    end
    return y
end

"""
    state_scale!(α::Float64,x::KState)

In-place `x ← α·x`. Returns `x`.
"""
function state_scale!(α::Float64,x::KState)
    @inbounds for k in eachindex(x)
        xk = x[k]
        for i in eachindex(xk)
            xk[i] *= α
        end
    end
    return x
end

"""
    state_copy!(dst::KState,src::KState)

In-place copy `dst ← src`. Returns `dst`.
"""
function state_copy!(dst::KState,src::KState)
    @inbounds for k in eachindex(dst,src)
        copyto!(dst[k],src[k])
    end
    return dst
end

"""
    state_zero!(x::KState)

In-place `x ← 0`. Returns `x`.
"""
function state_zero!(x::KState)
    @inbounds for k in eachindex(x)
        fill!(x[k],0.0)
    end
    return x
end

"""
    state_similar(x::KState)

Allocate a new zero-filled state with the same component shapes as `x`.
"""
state_similar(x::KState) = Array{Float64}[zeros(Float64,size(xk)) for xk in x]

"""
    state_clone(x::KState)

Allocate a new state that is a copy of `x`.
"""
state_clone(x::KState) = Array{Float64}[copy(xk) for xk in x]

"""
    spectral_radius_estimate(resid0::Float64,resid::Float64,iter::Int64)

Estimate the iteration's spectral radius from the geometric-average relative residual reduction
per pass, `ρ ≈ (resid/resid0)^{1/(iter-1)}`, where `resid0`/`resid` are the first/last residuals
over `iter` applications. Returns `NaN` when too few iterations or degenerate residuals make the
estimate meaningless. Shared by the iterative solvers so every method reports a comparable `ρ`
(≈ the scattering ratio for plain source iteration, much smaller for the accelerated methods).
"""
function spectral_radius_estimate(resid0::Float64,resid::Float64,iter::Int64)
    if iter > 1 && resid0 > 0.0 && resid > 0.0
        return (resid/resid0)^(1.0/(iter-1))
    end
    return NaN
end
