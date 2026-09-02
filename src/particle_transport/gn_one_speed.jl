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

"""
    gn_patch_list(Ndims,Nv,tiling,is_SPH,fold)

Flat list of the `(u,v,w)` triples carrying an angular patch, in the order the reference
`gn_inner_pass!` visits them.

The reference dimensions its per-patch arrays on the full `(8,Nv,Nw_max)` box even though
only `4·Nv·(Nv+1)` of those slots hold a patch (48 of 72 at `Nv = 3`, with the
polar-anchored tiling). Enumerating them once here lets the fast chain allocate exactly the
patches it sweeps, and gives each patch the linear id used to index the context and the
per-patch workspaces.
"""
function gn_patch_list(Ndims::Int64,Nv::Int64,tiling::String,is_SPH::Bool,fold::Bool)
    sx = [1,1,1,1,-1,-1,-1,-1]
    azim_collapsed = (Ndims == 1) && (!is_SPH || fold)
    z_fold_2D = (Ndims == 2) && (Nv == 1) && fold
    Nw_of(u,v) = azim_collapsed ? ((u == 1 || u == 5) ? 1 : 0) :
                 (z_fold_2D && iseven(u)) ? 0 :
                 ((tiling == "symmetric") ? (2*v - 1) : ((sx[u] == 1) ? (Nv + 1 - v) : v))
    list = NTuple{3,Int64}[]
    for u in 1:8, v in 1:Nv, w in 1:Nw_of(u,v)
        push!(list,(u,v,w))
    end
    return list
end

"""
    gn_one_speed_fast(...)

Optimized counterpart of `gn_one_speed`, restricted to 3D. Solves the one-speed transport
equation for one energy group, with the same arguments, the same acceleration schemes and
the same return values as the reference.

Three things change, all of them organizational:

- the cell systems are assembled and factorized once per (material, mesh-width triple,
  angular patch) in `gn_fast_context`, instead of once per voxel per patch per pass;
- the per-patch arrays are dimensioned on the flat patch list (`gn_patch_list`) rather than
  the `(8,Nv,Nw_max)` box, and each patch gets its own `GNFastWorkspace`;
- the boundary sources are folded onto the patches once, into concrete arrays.

The workspaces the reference reallocates for every energy group are the caller's
responsibility here — `compute_flux` hoists them out of the group loop and passes them in.

See `set_fast_path`.
"""
function gn_one_speed_fast(𝚽l::Array{Float64},Qlout::Array{Float64},Σt::Vector{Float64},Σs::Array{Float64},mat::Array{Int64,3},Ndims::Int64,ig::Int64,Ns::Vector{Int64},Δs::Vector{Vector{Float64}},Np::Int64,Nq::Int64,pl::Vector{Int64},Np_surf::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},isFC::Bool,C::Vector{Float64},ω::Vector{Vector{Float64}},I_max::Int64,ϵ_max::Float64,sources::Array{Union{Array{Float64},Float64}},isCSD::Bool,solver::Int64,𝚽E12::Array{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},T::Vector{Float64},ℳ::Array{Float64},𝒜::String,Ntot::Int64,𝒲::Array{Float64},Mll::Array{Float64},is_SPH::Bool,𝒩::Array{Float64},boundary_conditions::Vector{Int64},Np_source::Int64,Nv::Int64,Mll_surf::Array{Float64},Rpq::Array{Float64},need_boundary_flux::Bool,has_surface::Bool,tiling::String="polar-anchored",gmres_restart::Int64=30,anderson_depth::Int64=3,fold::Bool=false)

    Nmat  = length(Σt)
    plist = gn_patch_list(Ndims,Nv,tiling,is_SPH,fold)
    Npatch = length(plist)

    # Flux Initialization
    𝚽E12_temp = Array{Float64}(undef)
    if isCSD
        𝚽E12_temp = zeros(Np,Nm[4],Ns[1],Ns[2],Ns[3])
    end
    sxs = [1,1,1,1,-1,-1,-1,-1]
    sys = [1,1,-1,-1,1,1,-1,-1]
    szs = [1,-1,1,-1,1,-1,1,-1]

    # Precomputed cell systems and index tables for this group
    ctx = gn_fast_context(Ndims,Δs,Nmat,Σt,S⁻,S⁺,S,𝒪,isFC,isCSD,C,ω,𝒲,Nq,𝒩,plist)

    # Fixed boundary sources, folded onto the patches. Only the incoming face of each axis
    # is ever read by the sweep, so one array per axis is enough (the reference carries all
    # six and selects at every voxel column).
    srcx = zeros(Nq,Ns[2],Ns[3],Npatch)
    srcy = zeros(Nq,Ns[1],Ns[3],Npatch)
    srcz = zeros(Nq,Ns[1],Ns[2],Npatch)
    for (k,(u,v,w)) in enumerate(plist)
        has_surface || break
        fx = sxs[u] > 0 ? 1 : 2
        fy = sys[u] > 0 ? 3 : 4
        fz = szs[u] > 0 ? 5 : 6
        for p in range(1,Np_source), q in range(1,Nq)
            # A geometry has 2·Ndims faces, so the y and z slots do not exist below 3D.
            ax = Mll_surf[p,q,u,v,w,fx,1]
            Sx = sources[p,fx]
            if Ndims == 3
                ay = Mll_surf[p,q,u,v,w,fy,1]; Sy = sources[p,fy]
                az = Mll_surf[p,q,u,v,w,fz,1]; Sz = sources[p,fz]
                for iz in range(1,Ns[3]), iy in range(1,Ns[2])
                    srcx[q,iy,iz,k] += Sx[iy,iz] * ax
                end
                for iz in range(1,Ns[3]), ix in range(1,Ns[1])
                    srcy[q,ix,iz,k] += Sy[ix,iz] * ay
                end
                for iy in range(1,Ns[2]), ix in range(1,Ns[1])
                    srcz[q,ix,iy,k] += Sz[ix,iy] * az
                end
            elseif Ndims == 2
                # A 2D geometry has four faces, each carrying a 1-D array.
                ay = Mll_surf[p,q,u,v,w,fy,1]; Sy = sources[p,fy]
                for iy in range(1,Ns[2]); srcx[q,iy,1,k] += Sx[iy] * ax end
                for ix in range(1,Ns[1]); srcy[q,ix,1,k] += Sy[ix] * ay end
            elseif Ndims == 1
                # A 1D geometry has two faces, each carrying a single scalar.
                srcx[q,1,1,k] += Sx * ax
            else
                error("gn_one_speed_fast: unsupported dimension $(Ndims).")
            end
        end
    end
    # Zeroed clones used by the homogeneous operator A·z (the GN surface source enters the
    # sweeps through these arrays, not Np_source).
    srcx0 = zeros(Nq,Ns[2],Ns[3],Npatch)
    srcy0 = zeros(Nq,Ns[1],Ns[3],Npatch)
    srcz0 = zeros(Nq,Ns[1],Ns[2],Npatch)

    # Boundary fluxes initialization. With every face void nothing can ever make these
    # non-zero (only reflective and periodic conditions feed them), so they collapse to
    # empty arrays: the transforms that move them are skipped, and their contribution to the
    # Krylov/Livolant inner products — exactly zero — is unchanged.
    nb = need_boundary_flux ? 1 : 0
    𝚽x12⁻ = zeros(Np_surf,Nm[1],Ns[2]*nb,Ns[3]*nb,2)
    𝚽x12⁺ = zeros(Np_surf,Nm[1],Ns[2]*nb,Ns[3]*nb,2)
    𝚽y12⁻ = zeros(Np_surf,Nm[2],Ns[1]*nb,Ns[3]*nb,2)
    𝚽y12⁺ = zeros(Np_surf,Nm[2],Ns[1]*nb,Ns[3]*nb,2)
    𝚽z12⁻ = zeros(Np_surf,Nm[3],Ns[1]*nb,Ns[2]*nb,2)
    𝚽z12⁺ = zeros(Np_surf,Nm[3],Ns[1]*nb,Ns[2]*nb,2)

    # Restricted-angle buffers, dimensioned on the flat patch list, and the scaled Mll.
    NmE_face = isCSD ? Nm[4] : 0
    𝚽_q     = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3],Npatch)
    Q_q     = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3],Npatch)
    𝚽E12_q  = zeros(Nq,NmE_face,Ns[1],Ns[2],Ns[3],Npatch)
    𝚽x12_q  = zeros(Nq,Nm[1],Ns[2],Ns[3],2,Npatch)
    𝚽y12_q  = zeros(Nq,Nm[2],Ns[1],Ns[3],2,Npatch)
    𝚽z12_q  = zeros(Nq,Nm[3],Ns[1],Ns[2],2,Npatch)
    inv_4π  = 1/(4*π)
    Mll_factored = similar(Mll)
    for u in 1:8, v in 1:Nv, w in 1:size(Mll,5), q in 1:Nq, p in 1:Np
        Mll_factored[p,q,u,v,w] = (2*pl[p]+1) * inv_4π * Mll[p,q,u,v,w]
    end
    # Patch-major copies of the transform matrices. With Nq = 1 they turn the whole patch
    # loop of each volume transform into a single gemm (see gn_inner_pass_fast!).
    Mll_all   = zeros(Np,Nq*Npatch)
    Mfact_all = zeros(Np,Nq*Npatch)
    for (k,(u,v,w)) in enumerate(plist), q in 1:Nq, p in 1:Np
        Mll_all[p,(k-1)*Nq+q]   = Mll[p,q,u,v,w]
        Mfact_all[p,(k-1)*Nq+q] = Mll_factored[p,q,u,v,w]
    end

    # One workspace per patch, so a sweep never has to clear state left by the previous one.
    Nmf = [Nm[1],Nm[2],Nm[3]]
    wss = [GNFastWorkspace(Ndims,Nm[5]*Nq,Nq,Nmf,Ns) for _ in 1:Npatch]

    # Staging buffers for the batched angular transforms (Nq > 1 only; with Nq == 1 the
    # transform already lands in the sweep's layout and needs none). Blocked so the buffer
    # stays cache-sized whatever the mesh — see GN_FAST_BLOCK.
    tbufs = (Nq == 1) ? [Matrix{Float64}(undef,0,0)] :
                        [zeros(Nq*Npatch,GN_FAST_BLOCK)]

    # Source scratch
    Ql = similar(Qlout)
    ρ_in = NaN

    # If there is no source anywhere, the in-group solution is trivially zero
    if ~any(x->x!=0,sources) && ~any(x->x!=0,Qlout) && (~isCSD || ~any(x->x!=0,𝚽E12))
        𝚽l .= 0
        println(">>>Group ",ig," has converged ( ϵ = ",@sprintf("%.4E",0.0)," , N = ",1," , ρ = ",@sprintf("%.2f",ρ_in)," )")
        return 𝚽l,𝚽E12_temp,ρ_in,Ntot
    end

    # Tracks which incoming energy flux 𝚽E12_q currently holds (0 nothing, 1 the affine map's,
    # 2 the homogeneous operator's zeros), so the pass only rebuilds it when the caller
    # switches between the two.
    E12_state = Ref(0)

    # Shorthand wrapper around one source-iteration pass
    pass!(homogeneous;reconstruct::Bool=false) = gn_inner_pass_fast!(𝚽l,Qlout,Σs,mat,Ndims,Ns,Np,Nq,pl,Np_surf,Nm,isCSD,solver,𝚽E12,T,ℳ,boundary_conditions,Mll,Mll_surf,Rpq,Mll_factored,Mll_all,Mfact_all,ctx,plist,wss,Ql,𝚽E12_temp,homogeneous ? srcx0 : srcx,homogeneous ? srcy0 : srcy,homogeneous ? srcz0 : srcz,𝚽x12⁻,𝚽x12⁺,𝚽y12⁻,𝚽y12⁺,𝚽z12⁻,𝚽z12⁺,Q_q,𝚽_q,𝚽E12_q,𝚽x12_q,𝚽y12_q,𝚽z12_q,E12_state,need_boundary_flux,tbufs;homogeneous=homogeneous,reconstruct=reconstruct)

    # State vector z = (𝚽l, incoming boundary angular fluxes on each active axis)
    work = Ndims == 3 ? Array{Float64}[𝚽l,𝚽x12⁻,𝚽y12⁻,𝚽z12⁻] :
           Ndims == 2 ? Array{Float64}[𝚽l,𝚽x12⁻,𝚽y12⁻] :
                        Array{Float64}[𝚽l,𝚽x12⁻]
    Nref = Ref(Ntot)

    function load_and_pass!(out::KState,zin::KState,homogeneous::Bool)
        state_copy!(work,zin)
        pass!(homogeneous)
        state_copy!(out,work)
        Nref[] += 1
    end
    fixedpoint!(out::KState,zin::KState) = load_and_pass!(out,zin,false)
    matvec!(out::KState,zin::KState) = (load_and_pass!(out,zin,true); state_scale!(-1.0,out); state_axpy!(1.0,zin,out))

    local niter, resid, conv
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
        state_zero!(work); pass!(false); Nref[] += 1
        c = state_clone(work)
        z = state_similar(work)
        niter,resid,conv,ρ_in = gmres!(z,matvec!,c;restart=gmres_restart,maxit=I_max,tol=ϵ_max)
    elseif 𝒜 == "bicgstab"
        state_zero!(work); pass!(false); Nref[] += 1
        c = state_clone(work)
        z = state_similar(work)
        niter,resid,conv,ρ_in = bicgstab!(z,matvec!,c;maxit=I_max,tol=ϵ_max)
    else
        error("Unknown acceleration method: $𝒜.")
    end

    # Reconstruction pass on the converged state. This is the only pass that rebuilds the
    # outgoing energy flux 𝚽E12_temp, which is read after this call and never during the
    # iteration.
    state_copy!(work,z); pass!(false;reconstruct=true); Nref[] += 1
    Ntot = Nref[]

    if conv
        println(">>>Group $ig has converged ( ϵ = ",@sprintf("%.4E",resid)," , N = ",niter," , ρ = ",@sprintf("%.4f",ρ_in)," )")
    else
        println(">>>Group $ig has not converged ( ϵ = ",@sprintf("%.4E",resid)," , N = ",niter," , ρ = ",@sprintf("%.4f",ρ_in)," )")
    end
    return 𝚽l,𝚽E12_temp,ρ_in,Ntot
end