"""
    gn_one_speed(𝚽l::Array{Float64},Qlout::Array{Float64},Σt::Vector{Float64},
    Σs::Array{Float64},mat::Array{Int64,3},Ndims::Int64,ig::Int64,Ns::Vector{Int64},
    Δs::Vector{Vector{Float64}},Np::Int64,pl::Vector{Int64},pm::Vector{Int64},
    Np_surf::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},isFC::Bool,C::Vector{Float64},
    ω::Vector{Array{Float64}},I_max::Int64,ϵ_max::Float64,
    sources::Array{Union{Array{Float64},Float64}},isCSD::Bool,solver::Int64,
    𝚽E12::Array{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},
    T::Vector{Float64},ℳ::Array{Float64},𝒜::String,Ntot::Int64,𝒲::Array{Float64},
    Mll::Array{Float64},is_SPH::Bool,𝒩::Vector{Matrix{Float64}},
    boundary_conditions::Vector{Int64},Np_source::Int64)

Solve the one-speed transport equation for a given particle.  

# Input Argument(s)
- `𝚽l::Array{Float64}`: Legendre components of the in-cell flux.
- `Qlout::Array{Float64}`: Legendre components of the out-of-group in-cell source.
- `Σt::Vector{Float64}`: total cross-sections.
- `Σs::Array{Float64}`: Legendre moments of the scattering differential cross-sections.
- `mat::Array{Int64,3}`: material identifier per voxel.
- `Ndims::Int64`: dimension of the geometry.
- `ig::Int64`: energy group index.
- `Ns::Vector{Int64}`: number of voxels per axis.
- `Δs::Vector{Vector{Float64}}`: size of each voxels per axis.
- `Np::Int64`: number of angular interpolation basis.
- `pl::Vector{Int64}`: legendre order associated with each interpolation basis. 
- `pm::Vector{Int64}`: associate legendre order associated with each interpolation basis. 
- `Np_surf::Int64`: number of angular interpolation basis for each geometry surface.
- `𝒪::Vector{Int64}`: spatial and/or energy closure relation order.
- `Nm::Vector{Int64}`: number of spatial and/or energy moments.
- `isFC::Bool`: boolean indicating if the high-order incoming moments are fully coupled.
- `C::Vector{Float64}`: constants related to the spatial and energy normalized
   Legendre expansion.
- `ω::Vector{Array{Float64}}`: weighting factors of the closure relations.
- `I_max::Int64`: maximum number of iterations of inner iterations.
- `ϵ_max::Float64`: convergence criterion on the flux solution.
- `sources::Array{Union{Array{Float64},Float64}}`: surface sources intensities.
- `isCSD::Bool`: boolean to indicate if continuous slowing-down term is treated in 
   calculations.
- `solver::Int64`: indicate the type of solver to execute.
- `𝚽E12::Array{Float64}`: incoming flux along the energy axis.
- `S⁻::Vector{Float64}`: stopping power at higher energy group boundary.
- `S⁺::Vector{Float64}`: stopping power at lower energy group boundary.
- `S::Array{Float64}` : stopping powers.
- `T::Vector{Float64}`: momentum transfer.
- `ℳ::Array{Float64}`: Fokker-Planck scattering matrix.
- `𝒜::String` : acceleration method for in-group iterations ("none", "livolant", "anderson",
   "gmres" or "bicgstab").
- `Ntot::Int64` : accumulator for the total number of in-group iterations.
- `𝒲::Array{Float64}` : weighting constants.
- `Mll::Array{Float64}` : transformation matrix from restrict-angle to full-range fluxes.
- `is_SPH::Bool`: boolean indicating if spherical harmonics or Legendre basis is used.
- `𝒩::Vector{Matrix{Float64}}`: weights matrices for Legendre/spherical harmonics basis.
- `boundary_conditions::Vector{Int64}`: boundary conditions types for each axis.

# Output Argument(s)
- `𝚽l::Array{Float64}`: Legendre components of the in-cell flux.
- `𝚽E12::Array{Float64}`: outgoing flux along the energy axis.
- `ρ_in::Float64`: estimated spectral radius.
- `Ntot::Int64` : accumulator for the total number of in-group iterations.

# Reference(s)

"""
function gn_one_speed(𝚽l::Array{Float64},Qlout::Array{Float64},Σt::Vector{Float64},Σs::Array{Float64},mat::Array{Int64,3},Ndims::Int64,ig::Int64,Ns::Vector{Int64},Δs::Vector{Vector{Float64}},Np::Int64,Nq::Int64,pl::Vector{Int64},pm::Vector{Int64},Np_surf::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},isFC::Bool,C::Vector{Float64},ω::Vector{Vector{Float64}},I_max::Int64,ϵ_max::Float64,sources::Array{Union{Array{Float64},Float64}},isCSD::Bool,solver::Int64,𝚽E12::Array{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},T::Vector{Float64},ℳ::Array{Float64},𝒜::String,Ntot::Int64,𝒲::Array{Float64},Mll::Array{Float64},is_SPH::Bool,𝒩::Array{Float64},boundary_conditions::Vector{Int64},Np_source::Int64,Nv::Int64,Mll_surf::Array{Float64},Rpq::Array{Float64},tiling::String="polar-anchored",gmres_restart::Int64=30,anderson_depth::Int64=3,fold::Bool=false)

    # Flux Initialization
    𝚽E12_temp = Array{Float64}(undef)
    if isCSD
        𝚽E12_temp = zeros(Np,Nm[4],Ns[1],Ns[2],Ns[3])
    end
    sx = [1,1,1,1,-1,-1,-1,-1]
    if (Ndims > 1) sy = [1,1,-1,-1,1,1,-1,-1] end
    if (Ndims > 2) sz = [1,-1,1,-1,1,-1,1,-1] end

    # Patch indexing helper: number of patches along the "w" axis for octant u,
    # row v, subdivision Nv, and the chosen angular tiling. In 1D (azimuthally
    # symmetric, both Legendre and spherical-harmonics bases) the angular domain
    # collapses to μ-bands: a single azimuthal slot and only octants u ∈ {1, 5}
    # carry patches.
    # In 1D the angular domain collapses to two half-spheres (u ∈ {1,5}, full
    # azimuth) for the Legendre basis and for the folded spherical-harmonics basis;
    # the unfolded spherical-harmonics basis keeps the full octant tiling. In 2D the
    # z-symmetry fold skips the even octants (four quadrants {1,3,5,7}).
    azim_collapsed = (Ndims == 1) && (!is_SPH || fold)
    z_fold_2D = (Ndims == 2) && (Nv == 1) && fold
    Nw_max = azim_collapsed ? 1 : ((tiling == "symmetric") ? (2*Nv - 1) : Nv)
    Nw_of = azim_collapsed ? ((u, v) -> (u == 1 || u == 5) ? 1 : 0) :
            ((u, v) -> (z_fold_2D && iseven(u)) ? 0 :
                       ((tiling == "symmetric") ? (2*v - 1) : ((sx[u] == 1) ? (Nv + 1 - v) : v)))

    # Fixed boundary sources
    if Ndims == 1
        sources_q = zeros(Nq,2*Ndims,8,Nv,Nw_max)
    else
        sources_q = Array{Union{Float64,Array{Float64}}}(undef,Nq,2*Ndims,8,Nv,Nw_max)
        for q in range(1,Nq), u in range(1,8), v in range(1,Nv)
            Nw = Nw_of(u, v)
            for w in range(1,Nw), ib in range(1,2)
                if Ndims == 2
                    sources_q[q,ib,u,v,w] = zeros(Ns[2])
                    sources_q[q,ib+2,u,v,w] = zeros(Ns[1])
                elseif Ndims == 3
                    sources_q[q,ib,u,v,w] = zeros(Ns[2],Ns[3])
                    sources_q[q,ib+2,u,v,w] = zeros(Ns[1],Ns[3])
                    sources_q[q,ib+4,u,v,w] = zeros(Ns[1],Ns[2])
                else
                    error("Invalid number of dimensions.")
                end
            end
        end
    end
    for p in range(1,Np_source), q in range(1,Nq), u in range(1,8), v in range(1,Nv)
        Nw = Nw_of(u, v)
        for w in range(1,Nw)
            for ib in range(1,2)
                if Ndims == 1
                    sources_q[q,ib,u,v,w] += sources[p,ib] * Mll_surf[p,q,u,v,w,ib,1]
                elseif Ndims == 2
                    for iy in range(1,Ns[2])
                        sources_q[q,ib,u,v,w][iy] += sources[p,ib][iy] * Mll_surf[p,q,u,v,w,ib,1]
                    end
                    for ix in range(1,Ns[1])
                        sources_q[q,ib+2,u,v,w][ix] += sources[p,ib+2][ix] * Mll_surf[p,q,u,v,w,ib+2,1]
                    end
                elseif Ndims == 3
                    for iy in range(1,Ns[2]), iz in range(1,Ns[3])
                        sources_q[q,ib,u,v,w][iy,iz] += sources[p,ib][iy,iz] * Mll_surf[p,q,u,v,w,ib,1]
                    end
                    for ix in range(1,Ns[1]), iz in range(1,Ns[3])
                        sources_q[q,ib+2,u,v,w][ix,iz] += sources[p,ib+2][ix,iz] * Mll_surf[p,q,u,v,w,ib+2,1]
                    end
                    for ix in range(1,Ns[1]), iy in range(1,Ns[2])
                        sources_q[q,ib+4,u,v,w][ix,iy] += sources[p,ib+4][ix,iy] * Mll_surf[p,q,u,v,w,ib+4,1]
                    end
                else
                    error("Invalid number of dimensions.")
                end
            end
        end
    end

    # Zeroed clone of the surface source, same shape as sources_q, used by the homogeneous
    # operator A·z (the GN surface source enters the sweeps through sources_q, not Np_source).
    if Ndims == 1
        sources_q_zero = zeros(Nq,2*Ndims,8,Nv,Nw_max)
    else
        sources_q_zero = Array{Union{Float64,Array{Float64}}}(undef,Nq,2*Ndims,8,Nv,Nw_max)
        for q in range(1,Nq), u in range(1,8), v in range(1,Nv)
            Nw = Nw_of(u, v)
            for w in range(1,Nw), ib in range(1,2)
                if Ndims == 2
                    sources_q_zero[q,ib,u,v,w] = zeros(Ns[2])
                    sources_q_zero[q,ib+2,u,v,w] = zeros(Ns[1])
                elseif Ndims == 3
                    sources_q_zero[q,ib,u,v,w] = zeros(Ns[2],Ns[3])
                    sources_q_zero[q,ib+2,u,v,w] = zeros(Ns[1],Ns[3])
                    sources_q_zero[q,ib+4,u,v,w] = zeros(Ns[1],Ns[2])
                else
                    error("Invalid number of dimensions.")
                end
            end
        end
    end

    # Boundary fluxes initialization (unused axes get empty placeholders, never accessed)
    𝚽y12⁻ = zeros(0); 𝚽y12⁺ = zeros(0); 𝚽z12⁻ = zeros(0); 𝚽z12⁺ = zeros(0)
    if Ndims == 1
        𝚽x12⁻ = zeros(Np_surf,Nm[1],2)
        𝚽x12⁺ = zeros(Np_surf,Nm[1],2)
    elseif Ndims == 2
        𝚽x12⁻ = zeros(Np_surf,Nm[1],Ns[2],2)
        𝚽x12⁺ = zeros(Np_surf,Nm[1],Ns[2],2)
        𝚽y12⁻ = zeros(Np_surf,Nm[2],Ns[1],2)
        𝚽y12⁺ = zeros(Np_surf,Nm[2],Ns[1],2)
    else
        𝚽x12⁻ = zeros(Np_surf,Nm[1],Ns[2],Ns[3],2)
        𝚽x12⁺ = zeros(Np_surf,Nm[1],Ns[2],Ns[3],2)
        𝚽y12⁻ = zeros(Np_surf,Nm[2],Ns[1],Ns[3],2)
        𝚽y12⁺ = zeros(Np_surf,Nm[2],Ns[1],Ns[3],2)
        𝚽z12⁻ = zeros(Np_surf,Nm[3],Ns[1],Ns[2],2)
        𝚽z12⁺ = zeros(Np_surf,Nm[3],Ns[1],Ns[2],2)
    end

    # Pre-allocate restricted-angle buffers and precompute scaled Mll outside the
    # source iteration loop to avoid repeated allocations; BLAS mul! replaces loops.
    # (Unused axes get empty placeholders, never accessed.)
    𝚽y12_q = zeros(0); 𝚽z12_q = zeros(0)
    inv_4π = 1/(4*π)
    if Ndims == 3
        𝚽_q = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3],8,Nv,Nw_max)
        Q_q = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3],8,Nv,Nw_max)
        𝚽E12_q = zeros(Nq,Nm[4],Ns[1],Ns[2],Ns[3],8,Nv,Nw_max)
        𝚽x12_q = zeros(Nq,Nm[1],Ns[2],Ns[3],2,8,Nv,Nw_max)
        𝚽y12_q = zeros(Nq,Nm[2],Ns[1],Ns[3],2,8,Nv,Nw_max)
        𝚽z12_q = zeros(Nq,Nm[3],Ns[1],Ns[2],2,8,Nv,Nw_max)
        Mll_factored = similar(Mll)
        for u in 1:8, v in 1:Nv, w in 1:Nw_max, q in 1:Nq, p in 1:Np
            Mll_factored[p,q,u,v,w] = (2*pl[p]+1) * inv_4π * Mll[p,q,u,v,w]
        end
    elseif Ndims == 2
        𝚽_q = zeros(Nq,Nm[5],Ns[1],Ns[2],8,Nv,Nw_max)
        Q_q = zeros(Nq,Nm[5],Ns[1],Ns[2],8,Nv,Nw_max)
        𝚽E12_q = zeros(Nq,Nm[4],Ns[1],Ns[2],8,Nv,Nw_max)
        𝚽x12_q = zeros(Nq,Nm[1],Ns[2],2,8,Nv,Nw_max)
        𝚽y12_q = zeros(Nq,Nm[2],Ns[1],2,8,Nv,Nw_max)
        Mll_factored = similar(Mll)
        for u in 1:8, v in 1:Nv, w in 1:Nw_max, q in 1:Nq, p in 1:Np
            Mll_factored[p,q,u,v,w] = (2*pl[p]+1) * inv_4π * Mll[p,q,u,v,w]
        end
    elseif Ndims == 1
        𝚽_q = zeros(Nq,Nm[5],Ns[1],8,Nv,Nw_max)
        Q_q = zeros(Nq,Nm[5],Ns[1],8,Nv,Nw_max)
        𝚽E12_q = zeros(Nq,Nm[4],Ns[1],8,Nv,Nw_max)
        𝚽x12_q = zeros(Nq,Nm[1],2,8,Nv,Nw_max)
        Mll_factored = similar(Mll)
        for u in 1:8, v in 1:Nv, w in 1:Nw_max, q in 1:Nq, p in 1:Np
            fac = is_SPH ? (2*pl[p]+1)*inv_4π : (2*pl[p]+1)/2
            Mll_factored[p,q,u,v,w] = fac * Mll[p,q,u,v,w]
        end
    end

    # Pre-allocate per-voxel workspaces (𝒮, Q, 𝚽) and moving-boundary buffers
    # to avoid allocations inside the deep sweep loops. The cell-system
    # dimension matches gn_3D_BTE!/gn_3D_BFP! exactly: Nm[5] in-cell moments × Nq angular.
    Nm_solve = Nm[5] * Nq
    if Ndims == 3
        𝚽x12_buf = zeros(Nq, Nm[1], Ns[2], Ns[3])
        𝚽y12_buf = zeros(Nq, Nm[2], Ns[3])
        𝚽z12_buf = zeros(Nq, Nm[3])
    elseif Ndims == 2
        𝚽x12_buf = zeros(Nq, Nm[1], Ns[2])
        𝚽y12_buf = zeros(Nq, Nm[2])
        𝚽z12_buf = zeros(0,0)
    else
        𝚽x12_buf = zeros(Nq, Nm[1])
        𝚽y12_buf = zeros(0,0)
        𝚽z12_buf = zeros(0,0)
    end
    𝒮_ws = zeros(Nm_solve, Nm_solve)
    Q_ws = zeros(Nm_solve)
    𝚽_ws = zeros(Nm_solve)

    # Source scratch
    Ql = similar(Qlout)
    ρ_in = NaN

    # If there is no source anywhere, the in-group solution is trivially zero
    if ~any(x->x!=0,sources) && ~any(x->x!=0,Qlout) && (~isCSD || ~any(x->x!=0,𝚽E12))
        𝚽l .= 0
        println(">>>Group ",ig," has converged ( ϵ = ",@sprintf("%.4E",0.0)," , N = ",1," , ρ = ",@sprintf("%.2f",ρ_in)," )")
        return 𝚽l,𝚽E12_temp,ρ_in,Ntot
    end

    # Shorthand wrapper around one source-iteration pass
    pass!(homogeneous) = gn_inner_pass!(𝚽l,Qlout,Σt,Σs,mat,Ndims,Ns,Δs,Np,Nq,pl,Np_surf,𝒪,Nm,isFC,C,ω,isCSD,solver,𝚽E12,S⁻,S⁺,S,T,ℳ,𝒲,𝒩,boundary_conditions,Np_source,Nv,Mll,Mll_surf,Rpq,Mll_factored,tiling,is_SPH,fold,Ql,𝚽E12_temp,sources_q,sources_q_zero,𝚽x12⁻,𝚽x12⁺,𝚽y12⁻,𝚽y12⁺,𝚽z12⁻,𝚽z12⁺,Q_q,𝚽_q,𝚽E12_q,𝚽x12_q,𝚽y12_q,𝚽z12_q,𝒮_ws,Q_ws,𝚽_ws,𝚽x12_buf,𝚽y12_buf,𝚽z12_buf;homogeneous=homogeneous)

    # State vector z = (𝚽l, incoming boundary angular fluxes on each active axis)
    if Ndims == 1
        work = Array{Float64}[𝚽l,𝚽x12⁻]
    elseif Ndims == 2
        work = Array{Float64}[𝚽l,𝚽x12⁻,𝚽y12⁻]
    elseif Ndims == 3
        work = Array{Float64}[𝚽l,𝚽x12⁻,𝚽y12⁻,𝚽z12⁻]
    end
    Nref = Ref(Ntot)

    # One application of the affine map T (homogeneous=false) or the linear operator A (=true):
    # load the state into the working arrays, run one pass, read the result back.
    function load_and_pass!(out::KState,zin::KState,homogeneous::Bool)
        state_copy!(work,zin)
        pass!(homogeneous)
        state_copy!(out,work)
        Nref[] += 1
    end
    # Closure that applies the fixed-point map T (used by the source-iteration solvers), and the
    # matrix-vector product zin ↦ (I − A)·zin = zin − A·zin (used by the Krylov solvers).
    fixedpoint!(out::KState,zin::KState) = load_and_pass!(out,zin,false)
    matvec!(out::KState,zin::KState) = (load_and_pass!(out,zin,true); state_scale!(-1.0,out); state_axpy!(1.0,zin,out))

    # One branch per acceleration method, all built on the single pass above.
    local niter, resid, conv
    if 𝒜 == "none"
        # Plain source iteration: livolant! with the extrapolation disabled.
        z = state_clone(work)
        niter,resid,conv,ρ_in = livolant!(z,fixedpoint!;maxit=I_max,tol=ϵ_max,period=typemax(Int64))

    elseif 𝒜 == "livolant"
        # Source iteration with the periodic two-point Livolant extrapolation (every 3 passes).
        z = state_clone(work)
        niter,resid,conv,ρ_in = livolant!(z,fixedpoint!;maxit=I_max,tol=ϵ_max,period=3)

    elseif 𝒜 == "anderson"
        # Depth-m Anderson acceleration of the fixed-point iteration.
        z = state_clone(work)
        niter,resid,conv,ρ_in = anderson!(z,fixedpoint!;depth=anderson_depth,maxit=I_max,tol=ϵ_max,β=1.0)

    elseif 𝒜 == "gmres"
        # Restarted GMRES on (I − A) z = c, with the right-hand side c = T(0) from a zero state.
        state_zero!(work); pass!(false); Nref[] += 1
        c = state_clone(work)
        z = state_similar(work)
        niter,resid,conv,ρ_in = gmres!(z,matvec!,c;restart=gmres_restart,maxit=I_max,tol=ϵ_max)

    elseif 𝒜 == "bicgstab"
        # BiCGStab on (I − A) z = c, with the right-hand side c = T(0) from a zero state.
        state_zero!(work); pass!(false); Nref[] += 1
        c = state_clone(work)
        z = state_similar(work)
        niter,resid,conv,ρ_in = bicgstab!(z,matvec!,c;maxit=I_max,tol=ϵ_max)

    else
        error("Unknown acceleration method: $𝒜.")
    end

    # Reconstruction pass: a final non-homogeneous pass on the converged state refills the physical
    # 𝚽l, the outgoing energy flux 𝚽E12_temp and the boundary fluxes (z is a fixed point).
    state_copy!(work,z); pass!(false); Nref[] += 1
    Ntot = Nref[]

    if conv
        println(">>>Group $ig has converged ( ϵ = ",@sprintf("%.4E",resid)," , N = ",niter," , ρ = ",@sprintf("%.4f",ρ_in)," )")
    else
        println(">>>Group $ig has not converged ( ϵ = ",@sprintf("%.4E",resid)," , N = ",niter," , ρ = ",@sprintf("%.4f",ρ_in)," )")
    end
    return 𝚽l,𝚽E12_temp,ρ_in,Ntot
end