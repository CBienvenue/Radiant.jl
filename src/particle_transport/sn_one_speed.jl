"""
    compute_one_speed(𝚽l::Array{Float64},Qlout::Array{Float64},Σt::Vector{Float64},
    Σs::Array{Float64},mat::Array{Int64,3},ndims::Int64,Nd::Int64,ig::Int64,
    Ns::Vector{Int64},Δs::Vector{Vector{Float64}},Ω::Vector{Vector{Float64}},
    Mn::Array{Float64,2},Dn::Array{Float64,2},Np::Int64,pl::Vector{Int64},
    Mn_surf::Vector{Array{Float64}},Dn_surf::Vector{Array{Float64}},Np_surf::Int64,
    n_to_n⁺::Vector{Vector{Int64}},𝒪::Vector{Int64},Nm::Vector{Int64},isFC::Bool,
    C::Vector{Float64},ω::Vector{Array{Float64}},I_max::Int64,ϵ_max::Float64,
    sources::Array{Union{Array{Float64},Float64}},isAdapt::Bool,isCSD::Bool,solver::Int64,
    E::Float64,ΔE::Float64,𝚽E12::Array{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64},
    S::Array{Float64},T::Vector{Float64},ℳ::Array{Float64},𝒜::String,Ntot::Int64,
    is_EM::Bool,ℳ_EM::Array{Float64},𝒲::Array{Float64})

Solve the one-speed transport equation for a given particle.  

# Input Argument(s)
- `𝚽l::Array{Float64}`: Legendre components of the in-cell flux.
- `Qlout::Array{Float64}`: Legendre components of the out-of-group in-cell source.
- `Σt::Vector{Float64}`: total cross-sections.
- `Σs::Array{Float64}`: Legendre moments of the scattering differential cross-sections.
- `mat::Array{Int64,3}`: material identifier per voxel.
- `ndims::Int64`: dimension of the geometry.
- `Nd::Int64`: number of discrete ordinates.
- `ig::Int64`: energy group index.
- `Ns::Vector{Int64}`: number of voxels per axis.
- `Δs::Vector{Vector{Float64}}`: size of each voxels per axis.
- `Ω::Union{Vector{Vector{Float64}},Vector{Float64}}`: director cosines.
- `Mn::Array{Float64,2}`: moment-to-discrete matrix.
- `Dn::Array{Float64,2}`: discrete-to-moment matrix.
- `Np::Int64`: number of angular interpolation basis.
- `Mn_surf::Vector{Array{Float64}}`: moment-to-discrete matrix for each geometry surface.
- `Dn_surf::Vector{Array{Float64}}`: discrete-to-moment matrix for each geometry surface.
- `Np_surf::Int64`: number of angular interpolation basis for each geometry surface.
- `n_to_n⁺::Vector{Vector{Int64}}`: mapping from full-range to half-range indices.
- `pl::Vector{Int64}`: legendre order associated with each interpolation basis. 
- `𝒪::Vector{Int64}`: spatial and/or energy closure relation order.
- `Nm::Vector{Int64}`: number of spatial and/or energy moments.
- `isFC::Bool`: boolean indicating if the high-order incoming moments are fully coupled.
- `C::Vector{Float64}`: constants related to the spatial and energy normalized
   Legendre expansion.
- `ω::Vector{Array{Float64}}`: weighting factors of the closure relations.
- `I_max::Int64`: maximum number of iterations of inner iterations.
- `ϵ_max::Float64`: convergence criterion on the flux solution.
- `sources::Array{Union{Array{Float64},Float64}}`: surface sources intensities.
- `isAdapt::Bool`: boolean for adaptive calculations.
- `isCSD::Bool`: boolean to indicate if continuous slowing-down term is treated in 
   calculations.
- `solver::Int64`: indicate the type of solver to execute.
- `E::Float64`: group midpoint energy.
- `ΔE::Float64`: energy group width.
- `𝚽E12::Array{Float64}`: incoming flux along the energy axis.
- `S⁻::Vector{Float64}`: stopping power at higher energy group boundary.
- `S⁺::Vector{Float64}`: stopping power at lower energy group boundary.
- `S::Array{Float64}` : stopping powers.
- `T::Vector{Float64}`: momentum transfer.
- `ℳ::Array{Float64}`: Fokker-Planck scattering matrix.
- `𝒜::String` : acceleration method for in-group iterations ("none", "livolant", "anderson",
   "gmres" or "bicgstab").
- `Ntot::Int64` : accumulator for the total number of in-group iterations.
- `is_EM::Bool` : boolean for electromagnetic fields.
- `ℳ_EM::Array{Float64}` : electromagnetic scattering matrix.
- `𝒲::Array{Float64}` : weighting constants.
- `boundary_conditions::Vector{Int64}` : boundary condition per geometry surface.
- `Np_source::Int64` : number of angular basis of the surface sources.
- `gmres_restart::Int64` : Krylov subspace size before restart (GMRES only).
- `anderson_depth::Int64` : Anderson memory depth.

# Output Argument(s)
- `𝚽l::Array{Float64}`: Legendre components of the in-cell flux.
- `𝚽E12::Array{Float64}`: outgoing flux along the energy axis.
- `ρ_in::Float64`: estimated spectral radius (NaN for the Krylov/Anderson methods, where it is
   replaced in the printout by the achieved relative residual).
- `Ntot::Int64` : accumulator for the total number of in-group iterations.

# Reference(s)
- Larsen and Morel (2010) : Advances in Discrete-Ordinates Methodology.
- Warsa, Wareing and Morel (2004) : Krylov Iterative Methods Applied to Multidimensional Sn
  Calculations.

"""
function sn_one_speed(𝚽l::Array{Float64},Qlout::Array{Float64},Σt::Vector{Float64},Σs::Array{Float64},mat::Array{Int64,3},ndims::Int64,Nd::Int64,ig::Int64,Ns::Vector{Int64},Δs::Vector{Vector{Float64}},Ω::Vector{Vector{Float64}},Mn::Array{Float64,2},Dn::Array{Float64,2},Np::Int64,pl::Vector{Int64},Mn_surf::Vector{Array{Float64}},Dn_surf::Vector{Array{Float64}},Np_surf::Int64,n_to_n⁺::Vector{Vector{Int64}},𝒪::Vector{Int64},Nm::Vector{Int64},isFC::Bool,C::Vector{Float64},ω::Vector{Array{Float64}},I_max::Int64,ϵ_max::Float64,sources::Array{Union{Array{Float64},Float64}},isAdapt::Bool,isCSD::Bool,solver::Int64,ΔE::Float64,𝚽E12::Array{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},T::Vector{Float64},ℳ::Array{Float64},𝒜::String,Ntot::Int64,is_EM::Bool,ℳ_EM::Array{Float64},𝒲::Array{Float64},boundary_conditions::Vector{Int64},Np_source::Int64,gmres_restart::Int64=30,anderson_depth::Int64=3)

    # Flux Initialization
    𝚽E12_temp = Array{Float64}(undef)
    if isCSD
        𝚽E12_temp = zeros(Nd,Nm[4],Ns[1],Ns[2],Ns[3])
    end

    # Boundary conditions initialization (unused axes get empty placeholders, never accessed)
    𝚽y12_in = zeros(0); 𝚽z12_in = zeros(0)
    𝚽y12_temp = zeros(0); 𝚽z12_temp = zeros(0)
    if ndims == 1
        𝚽x12_in = zeros(Np_surf,Nm[1],2)
        𝚽x12_temp = copy(𝚽x12_in)
    elseif ndims == 2
        𝚽x12_in = zeros(Np_surf,Nm[1],2,Ns[2])
        𝚽y12_in = zeros(Np_surf,Nm[2],2,Ns[1])
        𝚽x12_temp = copy(𝚽x12_in)
        𝚽y12_temp = copy(𝚽y12_in)
    elseif ndims == 3
        𝚽x12_in = zeros(Np_surf,Nm[1],2,Ns[2],Ns[3])
        𝚽y12_in = zeros(Np_surf,Nm[2],2,Ns[1],Ns[3])
        𝚽z12_in = zeros(Np_surf,Nm[3],2,Ns[1],Ns[2])
        𝚽x12_temp = copy(𝚽x12_in)
        𝚽y12_temp = copy(𝚽y12_in)
        𝚽z12_temp = copy(𝚽z12_in)
    else
        error("Dimension is not 1, 2 or 3.")
    end

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
    pass!(homogeneous) = sn_inner_pass!(𝚽l,𝚽x12_in,𝚽y12_in,𝚽z12_in,𝚽x12_temp,𝚽y12_temp,𝚽z12_temp,𝚽E12_temp,Ql,Qlout,𝚽E12,Σt,Σs,mat,ndims,Nd,Ns,Δs,Ω,Mn,Dn,Np,pl,Mn_surf,Dn_surf,Np_surf,n_to_n⁺,𝒪,Nm,isFC,C,ω,sources,isAdapt,isCSD,solver,ΔE,S⁻,S⁺,S,T,ℳ,is_EM,ℳ_EM,𝒲,boundary_conditions,Np_source;homogeneous=homogeneous)
    
    # State vector z = (𝚽l, incoming boundary angular fluxes on each active axis)
    if ndims == 1
        work = Array{Float64}[𝚽l,𝚽x12_in]
    elseif ndims == 2
        work = Array{Float64}[𝚽l,𝚽x12_in,𝚽y12_in]
    elseif ndims == 3
        work = Array{Float64}[𝚽l,𝚽x12_in,𝚽y12_in,𝚽z12_in]
    else
        error("Dimension is not 1, 2 or 3.")
    end
    Nref = Ref(Ntot)

    # One application of the affine map T (homogeneous=false) or the linear operator A (=true): load the state into the working arrays, run one pass, read the result back.
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

    # Reconstruction pass: a final non-homogeneous pass on the converged state fills the physical
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