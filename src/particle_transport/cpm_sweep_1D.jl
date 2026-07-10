"""
    cpm_cell_response(Σi::Float64,Δxi::Float64,Lp::Int64,Nν::Int64)

Response blocks of a single voxel for the sweeping CPM (TR-02, sweeping method). The unit-cell
collision, leakage and transmission probabilities are exactly the single-region (`Nx = 1`)
outputs of [`cpm_collision_matrix`](@ref) and [`cpm_surface_matrices`](@ref); this routine
extracts, from those, the per-cell maps used by the interface-current sweep. With the surface
degrees of freedom stacked as `[left ℓ=0…Nν, right ℓ=0…Nν]` (`S = Nν+1` per face), the cell
relates its in-cell flux moments `φ` and outgoing face currents `J⁺` to the in-cell source `Q`
and incoming face currents `J⁻` through `φ = Cvv Q + Cvs J⁻` and `J⁺ = Csv Q + Css J⁻`.

# Output Argument(s)
- `Cvv::Matrix{Float64}` : source → volume flux, size `(P, P)` with `P = Lp+1`.
- `Csv_L`, `Csv_R::Matrix{Float64}` : source → outgoing left / right current, size `(S, P)`.
- `Css_LR`, `Css_RL::Matrix{Float64}` : transmission right-incoming → left-outgoing and
  left-incoming → right-outgoing, size `(S, S)` (the same-face blocks vanish).
- `Cvs_L`, `Cvs_R::Matrix{Float64}` : incoming left / right current → volume flux, size `(P, S)`.
"""
function cpm_cell_response(Σi::Float64,Δxi::Float64,Lp::Int64,Nν::Int64)
    P = Lp+1
    S = Nν+1
    Cvv = cpm_collision_matrix([Σi],[Δxi],Lp)
    Cvs,Csv,Css = cpm_surface_matrices([Σi],[Δxi],Lp,Nν)
    L = 1:S ; R = S+1:2S
    return Cvv, Csv[R,:], Csv[L,:], Css[L,R], Css[R,L], Cvs[:,L], Cvs[:,R]
end

"""
    cpm_sweep_pass!(work::KState,Qext::Matrix{Float64},cells,Σs_in::Matrix{Float64},
    βL::Float64,βR::Float64,P::Int64,S::Int64,Nx::Int64;homogeneous::Bool=false)

One source-iteration pass of the 1D sweeping CPM. The state `work = [φ (P×Nx), Jbc (S×2)]` holds
the in-cell flux moments and the incoming boundary currents (left, right); the pass rebuilds the
in-cell source `Q = Q_{ext} + Σ_s φ` (only the scattering part when `homogeneous`), propagates the
interface currents by a forward (rightward) and a backward (leftward) sweep, updates the flux, and
sets the new incoming boundary currents from the reflected outgoing ones (`βL`, `βR` albedos).
"""
function cpm_sweep_pass!(work::KState,Qext::Matrix{Float64},cells,Σs_in::Matrix{Float64},βL::Float64,βR::Float64,S::Int64,Nx::Int64;homogeneous::Bool=false)

    φ = work[1]
    Jbc = work[2]

    # In-cell source (scattering + external, external dropped for the homogeneous operator).
    Q = Σs_in .* φ
    if ~homogeneous Q .+= Qext end

    # Per-cell blocks are stored as (Cvv, Csv_R, Csv_L, Css_LR, Css_RL, Cvs_L, Cvs_R).
    # Forward sweep : rightward interface currents Jp[:,k+1] at interface k (k = 0…Nx).
    Jp = zeros(S,Nx+1)
    Jp[:,1] = Jbc[:,1]
    for i in range(1,Nx)
        Jp[:,i+1] = cells[i][2]*Q[:,i] .+ cells[i][5]*Jp[:,i]   # Csv_R Q + Css_RL Jp
    end

    # Backward sweep : leftward interface currents Jm[:,k+1] at interface k.
    Jm = zeros(S,Nx+1)
    Jm[:,Nx+1] = Jbc[:,2]
    for i in range(Nx,1,step=-1)
        Jm[:,i] = cells[i][3]*Q[:,i] .+ cells[i][4]*Jm[:,i+1]   # Csv_L Q + Css_LR Jm
    end

    # In-cell flux : collision + incoming currents on both faces.
    for i in range(1,Nx)
        φ[:,i] = cells[i][1]*Q[:,i] .+ cells[i][6]*Jp[:,i] .+ cells[i][7]*Jm[:,i+1]   # Cvv Q + Cvs_L Jp + Cvs_R Jm
    end

    # New incoming boundary currents = albedo × outgoing boundary currents.
    Jbc[:,1] = βL .* Jm[:,1]        # left  : incoming = βL × outgoing-left  (Jm at interface 0)
    Jbc[:,2] = βR .* Jp[:,Nx+1]     # right : incoming = βR × outgoing-right (Jp at interface Nx)

    return work
end

"""
    cpm_sweep_1D(Qext::Matrix{Float64},cells,Σs_in::Matrix{Float64},βL::Float64,βR::Float64,
    P::Int64,S::Int64,Nx::Int64,𝒜::String,I_max::Int64,ϵ_max::Float64,ig::Int64,
    gmres_restart::Int64,anderson_depth::Int64)

Solve one energy group of the sweeping CPM by source iteration over the interface-current sweep
([`cpm_sweep_pass!`](@ref)), accelerated with the shared in-group solvers (`livolant!`, `anderson!`,
`gmres!`, `bicgstab!`) exactly as [`sn_one_speed`](@ref). Returns the in-cell flux moments
`φ (P×Nx)` and the estimated spectral radius.
"""
function cpm_sweep_1D(Qext::Matrix{Float64},cells,Σs_in::Matrix{Float64},βL::Float64,βR::Float64,P::Int64,S::Int64,Nx::Int64,𝒜::String,I_max::Int64,ϵ_max::Float64,ig::Int64,gmres_restart::Int64,anderson_depth::Int64)

    # State : in-cell flux moments and incoming boundary currents (left, right).
    work = Array{Float64}[zeros(P,Nx),zeros(S,2)]

    # Trivial zero solution when there is no source and no reflected inflow is possible.
    if ~any(x->x!=0,Qext) && βL == 0.0 && βR == 0.0
        println(">>>Group ",ig," has converged ( ϵ = ",@sprintf("%.4E",0.0)," , N = ",1," , ρ = ",@sprintf("%.2f",NaN)," )")
        return work[1], NaN
    end

    pass!(homogeneous) = cpm_sweep_pass!(work,Qext,cells,Σs_in,βL,βR,S,Nx;homogeneous=homogeneous)

    function load_and_pass!(out::KState,zin::KState,homogeneous::Bool)
        state_copy!(work,zin)
        pass!(homogeneous)
        state_copy!(out,work)
    end
    fixedpoint!(out::KState,zin::KState) = load_and_pass!(out,zin,false)
    matvec!(out::KState,zin::KState) = (load_and_pass!(out,zin,true); state_scale!(-1.0,out); state_axpy!(1.0,zin,out))

    local niter, resid, conv, ρ_in
    if 𝒜 == "none"
        z = state_clone(work)
        niter,resid,conv,ρ_in = livolant!(z,fixedpoint!;maxit=I_max,tol=ϵ_max,period=typemax(Int64))
    elseif 𝒜 == "livolant"
        z = state_clone(work)
        niter,resid,conv,ρ_in = livolant!(z,fixedpoint!;maxit=I_max,tol=ϵ_max,period=3)
    elseif 𝒜 == "anderson"
        z = state_clone(work)
        niter,resid,conv,ρ_in = anderson!(z,fixedpoint!;depth=anderson_depth,maxit=I_max,tol=ϵ_max,β=1.0)
    elseif 𝒜 == "gmres"
        state_zero!(work); pass!(false)
        c = state_clone(work)
        z = state_similar(work)
        niter,resid,conv,ρ_in = gmres!(z,matvec!,c;restart=gmres_restart,maxit=I_max,tol=ϵ_max)
    elseif 𝒜 == "bicgstab"
        state_zero!(work); pass!(false)
        c = state_clone(work)
        z = state_similar(work)
        niter,resid,conv,ρ_in = bicgstab!(z,matvec!,c;maxit=I_max,tol=ϵ_max)
    else
        error("Unknown acceleration method: $𝒜.")
    end

    # Reconstruction pass on the converged state.
    state_copy!(work,z); pass!(false)

    if conv
        println(">>>Group $ig has converged ( ϵ = ",@sprintf("%.4E",resid)," , N = ",niter," , ρ = ",@sprintf("%.4f",ρ_in)," )")
    else
        println(">>>Group $ig has not converged ( ϵ = ",@sprintf("%.4E",resid)," , N = ",niter," , ρ = ",@sprintf("%.4f",ρ_in)," )")
    end

    return work[1], ρ_in
end
