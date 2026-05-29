"""
    livolant(𝚽₂::Array{Float64},𝚽₁::Array{Float64},𝚽₀::Array{Float64})

Estimate the solution using estimate from Livolant acceleration method. 

# Input Argument(s)
- `𝚽₂::Array{Float64}`: flux at iteration (i)
- `𝚽₁::Array{Float64}`: flux at iteration (i-1)
- `𝚽₀::Array{Float64}`: flux at iteration (i-2)

# Output Argument(s)
- `𝚽::Array{Float64}`: flux estimated by Livolant acceleration method for iteration (i+1).

# Reference(s)
- Hébert (2016),  Applied Reactor Physics (Sect. C.1.3 - Iterative approach).

"""
function livolant(𝚽₂::Array{Float64},𝚽₁::Array{Float64},𝚽₀::Array{Float64})

# Compute acceleration factor
sum_eΔe = 0.0
sum_Δe2 = 0.0
for i in eachindex(𝚽₂, 𝚽₁, 𝚽₀)
    e₀ = 𝚽₁[i] - 𝚽₀[i]
    e₁ = 𝚽₂[i] - 𝚽₁[i]
    Δe = e₁ - e₀
    sum_eΔe += e₀ * Δe
    sum_Δe2 += Δe * Δe
end
μj = -sum_eΔe / sum_Δe2
if μj ≤ 0 μj = 1 end

# Compute the next iteration
𝚽 = similar(𝚽₂)
for i in eachindex(𝚽₂, 𝚽₁)
    𝚽[i] = μj * 𝚽₂[i] + (1-μj) * 𝚽₁[i]
end

return 𝚽

end

"""
    livolant_factor(z₂::KState,z₁::KState,z₀::KState)

Livolant extrapolation factor `μ` computed over a full Krylov state (`KState`, see
[krylov_state.jl](krylov_state.jl)) rather than a single array: the block-vector generalization of
the `μ` formula used in [`livolant`](@ref). With `e₀ = z₁ − z₀`, `e₁ = z₂ − z₁` and `Δe = e₁ − e₀`,

`μ = − ⟨e₀,Δe⟩ / ⟨Δe,Δe⟩` ,

guarded so that a non-positive or undefined value falls back to `μ = 1` (no extrapolation).
"""
function livolant_factor(z₂::KState,z₁::KState,z₀::KState)
    sum_eΔe = 0.0
    sum_Δe2 = 0.0
    @inbounds for k in eachindex(z₂,z₁,z₀)
        a = z₂[k]; b = z₁[k]; c = z₀[k]
        for i in eachindex(a,b,c)
            e₀ = b[i] - c[i]
            e₁ = a[i] - b[i]
            Δe = e₁ - e₀
            sum_eΔe += e₀ * Δe
            sum_Δe2 += Δe * Δe
        end
    end
    if sum_Δe2 == 0.0 return 1.0 end
    μ = -sum_eΔe / sum_Δe2
    if μ ≤ 0 μ = 1.0 end
    return μ
end

"""
    livolant!(z,fixedpoint!;maxit,tol,period)

Source iteration of the in-group fixed point `z = g(z)` with periodic Livolant two-point
extrapolation, written in the same closure-driven, allocation-light form as
[anderson.jl](anderson.jl): the map `g` (one full source-iteration pass `T`, fixed sources
included) is supplied through the in-place call `fixedpoint!(out,z)` (`out ← g(z)`) and all vectors
are Krylov states (`KState`, see [krylov_state.jl](krylov_state.jl)).

Every `period` iterations (once two prior iterates are available) the last three iterates are
recombined through the Livolant factor [`livolant_factor`](@ref). Setting `period = typemax(Int64)`
disables the extrapolation, recovering plain (unaccelerated) source iteration — this is the
`"none"` acceleration option. Because it accelerates the exact fixed-point sequence, it applies
unchanged to every boundary condition and physics mode (BTE / BFP / CSD / Fokker-Planck /
electromagnetic).

# Input Argument(s)
- `z::KState` : initial guess, overwritten with the (accelerated) solution.
- `fixedpoint!` : function `fixedpoint!(out::KState,z::KState)` computing `out ← g(z)`.
- `maxit::Int64` : maximum number of `fixedpoint!` applications.
- `tol::Float64` : relative residual tolerance on `‖g(z) − z‖`.
- `period::Int64` : number of iterations between extrapolations (`typemax(Int64)` ⇒ none).

# Output Argument(s)
- `iter::Int64` : number of `fixedpoint!` applications performed.
- `resid::Float64` : final relative residual.
- `converged::Bool` : whether `resid < tol` was reached.
- `ρ::Float64` : estimated spectral radius (see [`spectral_radius_estimate`](@ref)).

# Reference(s)
- Hébert (2016), Applied Reactor Physics (Sect. C.1.3 - Iterative approach).
"""
function livolant!(z::KState,fixedpoint!;maxit::Int64,tol::Float64,period::Int64=3)

    g  = state_similar(z)
    f  = state_similar(z)
    z1 = state_similar(z)   # post-pass iterate (i-1)
    z0 = state_similar(z)   # post-pass iterate (i-2)
    nhist = 0

    iter = 0
    resid = 1.0
    resid0 = NaN
    converged = false

    for it in range(1,maxit)

        fixedpoint!(g,z); iter += 1

        # Residual f = g(z) - z on the full state
        state_copy!(f,g); state_axpy!(-1.0,z,f)
        fnorm = state_norm(f)
        scale = max(state_norm(z),state_norm(g),1e-16)
        resid = fnorm/scale
        if isnan(resid0) resid0 = resid end
        if resid < tol state_copy!(z,g); converged = true; break end

        if period != typemax(Int64) && it % period == 0 && nhist ≥ 2
            # Livolant extrapolation z ← μ g + (1-μ) z1, using the last three iterates
            μ = livolant_factor(g,z1,z0)
            state_copy!(z,g); state_scale!(μ,z); state_axpy!(1.0-μ,z1,z)
        else
            state_copy!(z,g)
        end

        # History stores the un-extrapolated post-pass iterates
        state_copy!(z0,z1); state_copy!(z1,g); nhist += 1
    end

    return iter, resid, converged, spectral_radius_estimate(resid0,resid,iter)
end