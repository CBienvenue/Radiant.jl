"""
    pn_one_speed(𝚽l::Array{Float64},Qlout::Array{Float64},Σt::Vector{Float64},
    Σs::Array{Float64},mat::Array{Int64,3},ndims::Int64,ig::Int64,Ns::Vector{Int64},
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
- `ndims::Int64`: dimension of the geometry.
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
- `𝒜::String` : acceleration method for in-group iterations.
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
function pn_one_speed(𝚽l::Array{Float64},Qlout::Array{Float64},Σt::Vector{Float64},Σs::Array{Float64},mat::Array{Int64,3},ndims::Int64,ig::Int64,Ns::Vector{Int64},Δs::Vector{Vector{Float64}},Np::Int64,pl::Vector{Int64},pm::Vector{Int64},Np_surf::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},isFC::Bool,C::Vector{Float64},ω::Vector{Array{Float64}},I_max::Int64,ϵ_max::Float64,sources::Array{Union{Array{Float64},Float64}},isCSD::Bool,solver::Int64,𝚽E12::Array{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},T::Vector{Float64},ℳ::Array{Float64},𝒜::String,Ntot::Int64,𝒲::Array{Float64},Mll::Array{Float64},is_SPH::Bool,𝒩::Vector{Matrix{Float64}},boundary_conditions::Vector{Int64},Np_source::Int64)

    # Flux Initialization
    𝚽E12_temp = Array{Float64}(undef)
    if isCSD
        𝚽E12_temp = zeros(Np,Nm[4],Ns[1],Ns[2],Ns[3])
    end
    N⁻ = 2
    𝚽l⁻ = zeros(N⁻,Np,Nm[5],Ns[1],Ns[2],Ns[3])

    # Boundary conditions initialization
    if ndims == 1
        𝚽x12⁺ = zeros(Np,Nm[1],2)
        𝚽x12⁻ = zeros(Np,Nm[1],2)
    elseif ndims == 2
        𝚽x12⁺⁺ = zeros(Np,Nm[1],2,Ns[2])
        𝚽x12⁺⁻ = zeros(Np,Nm[1],2,Ns[2])
        𝚽x12⁻⁺ = zeros(Np,Nm[1],2,Ns[2])
        𝚽x12⁻⁻ = zeros(Np,Nm[1],2,Ns[2])
        𝚽y12⁺⁺ = zeros(Np,Nm[2],2,Ns[1])
        𝚽y12⁺⁻ = zeros(Np,Nm[2],2,Ns[1])
        𝚽y12⁻⁺ = zeros(Np,Nm[2],2,Ns[1])
        𝚽y12⁻⁻ = zeros(Np,Nm[2],2,Ns[1])
    elseif ndims == 3
        𝚽x12⁺⁺⁺ = zeros(Np,Nm[1],2,Ns[2],Ns[3])
        𝚽x12⁺⁺⁻ = zeros(Np,Nm[1],2,Ns[2],Ns[3])
        𝚽x12⁺⁻⁺ = zeros(Np,Nm[1],2,Ns[2],Ns[3])
        𝚽x12⁺⁻⁻ = zeros(Np,Nm[1],2,Ns[2],Ns[3])
        𝚽x12⁻⁺⁺ = zeros(Np,Nm[1],2,Ns[2],Ns[3])
        𝚽x12⁻⁺⁻ = zeros(Np,Nm[1],2,Ns[2],Ns[3])
        𝚽x12⁻⁻⁺ = zeros(Np,Nm[1],2,Ns[2],Ns[3])
        𝚽x12⁻⁻⁻ = zeros(Np,Nm[1],2,Ns[2],Ns[3])
        𝚽y12⁺⁺⁺ = zeros(Np,Nm[2],2,Ns[1],Ns[3])
        𝚽y12⁺⁺⁻ = zeros(Np,Nm[2],2,Ns[1],Ns[3])
        𝚽y12⁺⁻⁺ = zeros(Np,Nm[2],2,Ns[1],Ns[3])
        𝚽y12⁺⁻⁻ = zeros(Np,Nm[2],2,Ns[1],Ns[3])
        𝚽y12⁻⁺⁺ = zeros(Np,Nm[2],2,Ns[1],Ns[3])
        𝚽y12⁻⁺⁻ = zeros(Np,Nm[2],2,Ns[1],Ns[3])
        𝚽y12⁻⁻⁺ = zeros(Np,Nm[2],2,Ns[1],Ns[3])
        𝚽y12⁻⁻⁻ = zeros(Np,Nm[2],2,Ns[1],Ns[3])
        𝚽z12⁺⁺⁺ = zeros(Np,Nm[3],2,Ns[1],Ns[2])
        𝚽z12⁺⁺⁻ = zeros(Np,Nm[3],2,Ns[1],Ns[2])
        𝚽z12⁺⁻⁺ = zeros(Np,Nm[3],2,Ns[1],Ns[2])
        𝚽z12⁺⁻⁻ = zeros(Np,Nm[3],2,Ns[1],Ns[2])
        𝚽z12⁻⁺⁺ = zeros(Np,Nm[3],2,Ns[1],Ns[2])
        𝚽z12⁻⁺⁻ = zeros(Np,Nm[3],2,Ns[1],Ns[2])
        𝚽z12⁻⁻⁺ = zeros(Np,Nm[3],2,Ns[1],Ns[2])
        𝚽z12⁻⁻⁻ = zeros(Np,Nm[3],2,Ns[1],Ns[2])
    else
        error("Dimension is not 1, 2 or 3.")
    end

    # Source iteration loop until convergence
    i_in = 1
    ϵ_in = 0.0
    ρ_in = NaN
    isInnerConv=false
    while ~(isInnerConv)

        # Calculation of the Legendre components of the source (in-scattering)
        Ql = copy(Qlout)
        if solver ∉ [4,5,6] Ql = scattering_source(Ql,𝚽l,Σs,mat,Np,pl,Nm[5],Ns) end

        # Finite element treatment of the angular Fokker-Planck term
        if solver ∈ [2,4] Ql = fokker_planck_source(Np,Nm[5],T,𝚽l,Ql,Ns,mat,ℳ) end

        # If there is no source
        if ~any(x->x!=0,sources) && ~any(x->x!=0,Ql) && (~isCSD || (isCSD && ~any(x->x!=0,𝚽E12)))
            𝚽l = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            ϵ_in = 0.0; i_in = 1
            println(">>>Group ",ig," has converged ( ϵ = ",@sprintf("%.4E",ϵ_in)," , Nd = ",i_in," , ρ = ",@sprintf("%.2f",ρ_in)," )")
            break
        end

        #----
        # Loop over all discrete ordinates
        #----
        𝚽l .= 0
        𝚽E12_temp .= 0
        if ndims == 1
            𝚽⁺ = zeros(Np,Nm[5],Ns[1])
            𝚽⁻ = zeros(Np,Nm[5],Ns[1])
            Q⁺ = zeros(Np,Nm[5],Ns[1])
            Q⁻ = zeros(Np,Nm[5],Ns[1])
            𝚽E12⁺ = zeros(Np,Nm[4],Ns[1])
            𝚽E12⁻ = zeros(Np,Nm[4],Ns[1])
            for p in range(1,Np), q in range(1,Np)
                if is_SPH factor = (2*pl[p]+1)/(4*π) else factor = (2*pl[p]+1)/2 end
                for is in range(1,Nm[5]), ix in range(1,Ns[1])
                    Q⁺[q,is,ix] += factor * Ql[p,is,ix,1,1] * Mll[p,q]
                    Q⁻[q,is,ix] += factor * (-1)^pl[p] * Ql[p,is,ix,1,1] * Mll[p,q]
                end
                if isCSD
                    for is in range(1,Nm[4]), ix in range(1,Ns[1])
                        𝚽E12⁺[q,is,ix] += factor * 𝚽E12[p,is,ix,1,1] * Mll[p,q]
                        𝚽E12⁻[q,is,ix] += factor * (-1)^pl[p] * 𝚽E12[p,is,ix,1,1] * Mll[p,q]
                    end
                end
            end
            𝚽⁺,𝚽E12⁺,𝚽x12⁺ = pn_sweep_1D(1,𝚽⁺,Q⁺,Σt,mat[:,1,1],Ns[1],Δs[1],Np,Np_source,Np_surf,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁺,𝒲,isFC,isCSD,𝒩[1],𝚽x12⁺,𝚽x12⁻,boundary_conditions)
            𝚽⁻,𝚽E12⁻,𝚽x12⁻ = pn_sweep_1D(-1,𝚽⁻,Q⁻,Σt,mat[:,1,1],Ns[1],Δs[1],Np,Np_source,Np_surf,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁻,𝒲,isFC,isCSD,𝒩[1],𝚽x12⁺,𝚽x12⁻,boundary_conditions)
            for p in range(1,Np), q in range(1,Np)
                for is in range(1,Nm[5]), ix in range(1,Ns[1])
                    𝚽l[p,is,ix,1,1] += Mll[p,q] * (𝚽⁺[q,is,ix] + (-1)^pl[p]*𝚽⁻[q,is,ix])
                end
                if isCSD
                    for is in range(1,Nm[4]), ix in range(1,Ns[1])
                        𝚽E12_temp[p,is,ix,1,1] += Mll[p,q] * (𝚽E12⁺[q,is,ix] + (-1)^pl[p]*𝚽E12⁻[q,is,ix])
                    end
                end
            end
        elseif ndims == 2
            𝚽⁺⁺ = zeros(Np,Nm[5],Ns[1],Ns[2])
            𝚽⁺⁻ = zeros(Np,Nm[5],Ns[1],Ns[2])
            𝚽⁻⁺ = zeros(Np,Nm[5],Ns[1],Ns[2])
            𝚽⁻⁻ = zeros(Np,Nm[5],Ns[1],Ns[2])
            Q⁺⁺ = zeros(Np,Nm[5],Ns[1],Ns[2])
            Q⁺⁻ = zeros(Np,Nm[5],Ns[1],Ns[2])
            Q⁻⁺ = zeros(Np,Nm[5],Ns[1],Ns[2])
            Q⁻⁻ = zeros(Np,Nm[5],Ns[1],Ns[2])
            𝚽E12⁺⁺ = zeros(Np,Nm[5],Ns[1],Ns[2])
            𝚽E12⁺⁻ = zeros(Np,Nm[5],Ns[1],Ns[2])
            𝚽E12⁻⁺ = zeros(Np,Nm[5],Ns[1],Ns[2])
            𝚽E12⁻⁻ = zeros(Np,Nm[5],Ns[1],Ns[2])
            for ix in range(1,Ns[1]), iy in range(1,Ns[2]), p in range(1,Np), q in range(1,Np)
                factor = (2*pl[p]+1)/(4*π)
                for is in range(1,Nm[5])
                    Q⁺⁺[q,is,ix,iy] += factor * Ql[p,is,ix,iy,1] * Mll[p,q]
                    Q⁺⁻[q,is,ix,iy] += factor * (-1)^pm[p] * Ql[p,is,ix,iy,1] * Mll[p,q]
                    Q⁻⁺[q,is,ix,iy] += factor * (-1)^pl[p] * Ql[p,is,ix,iy,1] * Mll[p,q]
                    Q⁻⁻[q,is,ix,iy] += factor * (-1)^(pl[p]+pm[p]) * Ql[p,is,ix,iy,1] * Mll[p,q]
                end
                if isCSD
                    for is in range(1,Nm[4])
                        𝚽E12⁺⁺[q,is,ix,iy] += factor * 𝚽E12[p,is,ix,iy,1] * Mll[p,q]
                        𝚽E12⁺⁻[q,is,ix,iy] += factor * (-1)^pm[p] * 𝚽E12[p,is,ix,iy,1] * Mll[p,q]
                        𝚽E12⁻⁺[q,is,ix,iy] += factor * (-1)^pl[p] * 𝚽E12[p,is,ix,iy,1] * Mll[p,q]
                        𝚽E12⁻⁻[q,is,ix,iy] += factor * (-1)^(pl[p]+pm[p]) * 𝚽E12[p,is,ix,iy,1] * Mll[p,q]
                    end
                end
            end
            𝚽⁺⁺,𝚽E12⁺⁺,𝚽x12⁺⁺,𝚽y12⁺⁺ = pn_sweep_2D(1,1,𝚽⁺⁺,Q⁺⁺,Σt,mat[:,:,1],Ns[1],Ns[2],Δs[1],Δs[2],Np,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁺⁺,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝚽x12⁺⁺,𝚽x12⁻⁺,𝚽y12⁺⁺,𝚽y12⁺⁻,boundary_conditions)
            𝚽⁺⁻,𝚽E12⁺⁻,𝚽x12⁺⁻,𝚽y12⁺⁻ = pn_sweep_2D(1,-1,𝚽⁺⁻,Q⁺⁻,Σt,mat[:,:,1],Ns[1],Ns[2],Δs[1],Δs[2],Np,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁺⁻,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝚽x12⁺⁻,𝚽x12⁻⁻,𝚽y12⁺⁻,𝚽y12⁺⁺,boundary_conditions)
            𝚽⁻⁺,𝚽E12⁻⁺,𝚽x12⁻⁺,𝚽y12⁻⁺ = pn_sweep_2D(-1,1,𝚽⁻⁺,Q⁻⁺,Σt,mat[:,:,1],Ns[1],Ns[2],Δs[1],Δs[2],Np,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁻⁺,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝚽x12⁻⁺,𝚽x12⁺⁺,𝚽y12⁻⁺,𝚽y12⁻⁻,boundary_conditions)
            𝚽⁻⁻,𝚽E12⁻⁻,𝚽x12⁻⁻,𝚽y12⁻⁻ = pn_sweep_2D(-1,-1,𝚽⁻⁻,Q⁻⁻,Σt,mat[:,:,1],Ns[1],Ns[2],Δs[1],Δs[2],Np,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁻⁻,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝚽x12⁻⁻,𝚽x12⁺⁻,𝚽y12⁻⁻,𝚽y12⁻⁺,boundary_conditions)
            for ix in range(1,Ns[1]), iy in range(1,Ns[2]), p in range(1,Np), q in range(1,Np)
                for is in range(1,Nm[5])
                    𝚽l[p,is,ix,iy,1] += Mll[p,q] * (𝚽⁺⁺[q,is,ix,iy] + (-1)^pm[p] * 𝚽⁺⁻[q,is,ix,iy] + (-1)^pl[p] * 𝚽⁻⁺[q,is,ix,iy] + (-1)^(pl[p]+pm[p]) * 𝚽⁻⁻[q,is,ix,iy])
                end
                if isCSD
                    for is in range(1,Nm[4])
                        𝚽E12_temp[p,is,ix,iy,1] += Mll[p,q] * (𝚽E12⁺⁺[q,is,ix,iy] + (-1)^pm[p] * 𝚽E12⁺⁻[q,is,ix,iy] + (-1)^pl[p] * 𝚽E12⁻⁺[q,is,ix,iy] + (-1)^(pl[p]+pm[p]) * 𝚽E12⁻⁻[q,is,ix,iy])
                    end
                end
            end
        elseif ndims == 3
            𝚽⁺⁺⁺ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            𝚽⁺⁺⁻ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            𝚽⁺⁻⁺ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            𝚽⁺⁻⁻ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            𝚽⁻⁺⁺ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            𝚽⁻⁺⁻ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            𝚽⁻⁻⁺ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            𝚽⁻⁻⁻ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            Q⁺⁺⁺ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            Q⁺⁺⁻ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            Q⁺⁻⁺ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            Q⁺⁻⁻ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            Q⁻⁺⁺ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            Q⁻⁺⁻ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            Q⁻⁻⁺ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            Q⁻⁻⁻ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            𝚽E12⁺⁺⁺ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            𝚽E12⁺⁺⁻ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            𝚽E12⁺⁻⁺ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            𝚽E12⁺⁻⁻ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            𝚽E12⁻⁺⁺ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            𝚽E12⁻⁺⁻ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            𝚽E12⁻⁻⁺ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            𝚽E12⁻⁻⁻ = zeros(Np,Nm[5],Ns[1],Ns[2],Ns[3])
            for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3]), p in range(1,Np), q in range(1,Np)
                factor = (2*pl[p]+1)/(4*π)
                for is in range(1,Nm[5])
                    Q⁺⁺⁺[q,is,ix,iy,iz] += factor * Ql[p,is,ix,iy,iz] * Mll[p,q,1]
                    Q⁺⁺⁻[q,is,ix,iy,iz] += factor * (-1)^pm[p] * Ql[p,is,ix,iy,iz] * Mll[p,q,2]
                    Q⁺⁻⁺[q,is,ix,iy,iz] += factor * Ql[p,is,ix,iy,iz] * Mll[p,q,2]
                    Q⁺⁻⁻[q,is,ix,iy,iz] += factor * (-1)^pm[p] * Ql[p,is,ix,iy,iz] * Mll[p,q,1]
                    Q⁻⁺⁺[q,is,ix,iy,iz] += factor * (-1)^pl[p] * Ql[p,is,ix,iy,iz] * Mll[p,q,1]
                    Q⁻⁺⁻[q,is,ix,iy,iz] += factor * (-1)^(pl[p]+pm[p]) * Ql[p,is,ix,iy,iz] * Mll[p,q,2]
                    Q⁻⁻⁺[q,is,ix,iy,iz] += factor * (-1)^pl[p] * Ql[p,is,ix,iy,iz] * Mll[p,q,2]
                    Q⁻⁻⁻[q,is,ix,iy,iz] += factor * (-1)^(pl[p]+pm[p]) * Ql[p,is,ix,iy,iz] * Mll[p,q,1]
                end
                if isCSD
                    for is in range(1,Nm[4])
                        𝚽E12⁺⁺⁺[q,is,ix,iy,iz] += factor * 𝚽E12[p,is,ix,iy,iz] * Mll[p,q,1]
                        𝚽E12⁺⁺⁻[q,is,ix,iy,iz] += factor * (-1)^pm[p] * 𝚽E12[p,is,ix,iy,iz] * Mll[p,q,2]
                        𝚽E12⁺⁻⁺[q,is,ix,iy,iz] += factor * 𝚽E12[p,is,ix,iy,iz] * Mll[p,q,2]
                        𝚽E12⁺⁻⁻[q,is,ix,iy,iz] += factor * (-1)^pm[p] * 𝚽E12[p,is,ix,iy,iz] * Mll[p,q,1]
                        𝚽E12⁻⁺⁺[q,is,ix,iy,iz] += factor * (-1)^pl[p] * 𝚽E12[p,is,ix,iy,iz] * Mll[p,q,1]
                        𝚽E12⁻⁺⁻[q,is,ix,iy,iz] += factor * (-1)^(pl[p]+pm[p]) * 𝚽E12[p,is,ix,iy,iz] * Mll[p,q,2]
                        𝚽E12⁻⁻⁺[q,is,ix,iy,iz] += factor * (-1)^pl[p] * 𝚽E12[p,is,ix,iy,iz] * Mll[p,q,2]
                        𝚽E12⁻⁻⁻[q,is,ix,iy,iz] += factor * (-1)^(pl[p]+pm[p]) * 𝚽E12[p,is,ix,iy,iz] * Mll[p,q,1]
                    end
                end
            end
            𝚽⁺⁺⁺,𝚽E12⁺⁺⁺,𝚽x12⁺⁺⁺,𝚽y12⁺⁺⁺,𝚽z12⁺⁺⁺ = pn_sweep_3D(1,1,1,𝚽⁺⁺⁺,Q⁺⁺⁺,Σt,mat,Ns[1],Ns[2],Ns[3],Δs[1],Δs[2],Δs[3],Np,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁺⁺⁺,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝒩[3],𝚽x12⁺⁺⁺,𝚽x12⁻⁺⁺,𝚽y12⁺⁺⁺,𝚽y12⁺⁻⁺,𝚽z12⁺⁺⁺,𝚽z12⁺⁺⁻,boundary_conditions)
            𝚽⁺⁺⁻,𝚽E12⁺⁺⁻,𝚽x12⁺⁺⁻,𝚽y12⁺⁺⁻,𝚽z12⁺⁺⁻ = pn_sweep_3D(1,1,-1,𝚽⁺⁺⁻,Q⁺⁺⁻,Σt,mat,Ns[1],Ns[2],Ns[3],Δs[1],Δs[2],Δs[3],Np,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁺⁺⁻,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝒩[3],𝚽x12⁺⁺⁻,𝚽x12⁻⁺⁻,𝚽y12⁺⁺⁻,𝚽y12⁺⁻⁻,𝚽z12⁺⁺⁻,𝚽z12⁺⁺⁺,boundary_conditions)
            𝚽⁺⁻⁺,𝚽E12⁺⁻⁺,𝚽x12⁺⁻⁺,𝚽y12⁺⁻⁺,𝚽z12⁺⁻⁺ = pn_sweep_3D(1,-1,1,𝚽⁺⁻⁺,Q⁺⁻⁺,Σt,mat,Ns[1],Ns[2],Ns[3],Δs[1],Δs[2],Δs[3],Np,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁺⁻⁺,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝒩[3],𝚽x12⁺⁻⁺,𝚽x12⁻⁻⁺,𝚽y12⁺⁻⁺,𝚽y12⁺⁺⁺,𝚽z12⁺⁻⁺,𝚽z12⁺⁻⁻,boundary_conditions)
            𝚽⁺⁻⁻,𝚽E12⁺⁻⁻,𝚽x12⁺⁻⁻,𝚽y12⁺⁻⁻,𝚽z12⁺⁻⁻ = pn_sweep_3D(1,-1,-1,𝚽⁺⁻⁻,Q⁺⁻⁻,Σt,mat,Ns[1],Ns[2],Ns[3],Δs[1],Δs[2],Δs[3],Np,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁺⁻⁻,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝒩[3],𝚽x12⁺⁻⁻,𝚽x12⁻⁻⁻,𝚽y12⁺⁻⁻,𝚽y12⁺⁺⁻,𝚽z12⁺⁻⁻,𝚽z12⁺⁻⁺,boundary_conditions)
            𝚽⁻⁺⁺,𝚽E12⁻⁺⁺,𝚽x12⁻⁺⁺,𝚽y12⁻⁺⁺,𝚽z12⁻⁺⁺ = pn_sweep_3D(-1,1,1,𝚽⁻⁺⁺,Q⁻⁺⁺,Σt,mat,Ns[1],Ns[2],Ns[3],Δs[1],Δs[2],Δs[3],Np,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁻⁺⁺,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝒩[3],𝚽x12⁻⁺⁺,𝚽x12⁺⁺⁺,𝚽y12⁻⁺⁺,𝚽y12⁻⁻⁺,𝚽z12⁻⁺⁺,𝚽z12⁻⁺⁻,boundary_conditions)
            𝚽⁻⁺⁻,𝚽E12⁻⁺⁻,𝚽x12⁻⁺⁻,𝚽y12⁻⁺⁻,𝚽z12⁻⁺⁻ = pn_sweep_3D(-1,1,-1,𝚽⁻⁺⁻,Q⁻⁺⁻,Σt,mat,Ns[1],Ns[2],Ns[3],Δs[1],Δs[2],Δs[3],Np,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁻⁺⁻,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝒩[3],𝚽x12⁻⁺⁻,𝚽x12⁺⁺⁻,𝚽y12⁻⁺⁻,𝚽y12⁻⁻⁻,𝚽z12⁻⁺⁻,𝚽z12⁻⁺⁺,boundary_conditions)
            𝚽⁻⁻⁺,𝚽E12⁻⁻⁺,𝚽x12⁻⁻⁺,𝚽y12⁻⁻⁺,𝚽z12⁻⁻⁺ = pn_sweep_3D(-1,-1,1,𝚽⁻⁻⁺,Q⁻⁻⁺,Σt,mat,Ns[1],Ns[2],Ns[3],Δs[1],Δs[2],Δs[3],Np,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁻⁻⁺,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝒩[3],𝚽x12⁻⁻⁺,𝚽x12⁺⁻⁺,𝚽y12⁻⁻⁺,𝚽y12⁻⁺⁺,𝚽z12⁻⁻⁺,𝚽z12⁻⁻⁻,boundary_conditions)
            𝚽⁻⁻⁻,𝚽E12⁻⁻⁻,𝚽x12⁻⁻⁻,𝚽y12⁻⁻⁻,𝚽z12⁻⁻⁻ = pn_sweep_3D(-1,-1,-1,𝚽⁻⁻⁻,Q⁻⁻⁻,Σt,mat,Ns[1],Ns[2],Ns[3],Δs[1],Δs[2],Δs[3],Np,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁻⁻⁻,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝒩[3],𝚽x12⁻⁻⁻,𝚽x12⁺⁻⁻,𝚽y12⁻⁻⁻,𝚽y12⁻⁺⁻,𝚽z12⁻⁻⁻,𝚽z12⁻⁻⁺,boundary_conditions)
            for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3]), p in range(1,Np), q in range(1,Np)
                for is in range(1,Nm[5])
                    𝚽l[p,is,ix,iy,iz] += Mll[p,q,1] * (𝚽⁺⁺⁺[q,is,ix,iy,iz] + (-1)^pm[p] * 𝚽⁺⁻⁻[q,is,ix,iy,iz] + (-1)^pl[p] * 𝚽⁻⁺⁺[q,is,ix,iy,iz] + (-1)^(pl[p]+pm[p]) * 𝚽⁻⁻⁻[q,is,ix,iy,iz]) + Mll[p,q,2] * (𝚽⁺⁻⁺[q,is,ix,iy,iz] + (-1)^pm[p] * 𝚽⁺⁺⁻[q,is,ix,iy,iz] + (-1)^pl[p] * 𝚽⁻⁻⁺[q,is,ix,iy,iz] + (-1)^(pl[p]+pm[p]) * 𝚽⁻⁺⁻[q,is,ix,iy,iz])
                end
                if isCSD
                    for is in range(1,Nm[4])
                        𝚽E12_temp[p,is,ix,iy,iz] += Mll[p,q,1] * (𝚽E12⁺⁺⁺[q,is,ix,iy,iz] + (-1)^pm[p] * 𝚽E12⁺⁺⁻[q,is,ix,iy,iz] + (-1)^pl[p] * 𝚽E12⁺⁻⁺[q,is,ix,iy,iz] + (-1)^(pl[p]+pm[p]) * 𝚽E12⁺⁻⁻[q,is,ix,iy,iz]) + Mll[p,q,2] * (𝚽E12⁻⁺⁺[q,is,ix,iy,iz] + (-1)^pm[p] * 𝚽E12⁻⁺⁻[q,is,ix,iy,iz] + (-1)^pl[p] * 𝚽E12⁻⁻⁺[q,is,ix,iy,iz] + (-1)^(pl[p]+pm[p]) * 𝚽E12⁻⁻⁻[q,is,ix,iy,iz])
                    end
                end
            end
        else
            error("Multidimensional PN method not implemented yet.")
        end
        
        #----
        # Verification of convergence of the one-group flux
        #----  
        ϵ_in = 0.0
        if (solver ∉ [5,6]) ϵ_in = norm(𝚽l .- 𝚽l⁻[1,:,:,:,:,:]) / max(norm(𝚽l), 1e-16) end
        if (ϵ_in < ϵ_max) || i_in >= I_max

            # Convergence or maximum iterations reach
            isInnerConv = true
            Ntot += i_in
            if i_in ≥ 3 ρ_in = sqrt(sum(( vec(𝚽l[1,1,:,:,:]) .- vec(𝚽l⁻[1,1,1,:,:,:]) ).^2))/sqrt(sum(( vec(𝚽l⁻[1,1,1,:,:,:]) .- vec(𝚽l⁻[2,1,1,:,:,:]) ).^2)) end
            if ~(i_in >= I_max)
                println(">>>Group $ig has converged ( ϵ = ",@sprintf("%.4E",ϵ_in)," , N = ",i_in," , ρ = ",@sprintf("%.2f",ρ_in)," )")
            else
                println(">>>Group $ig has not converged ( ϵ = ",@sprintf("%.4E",ϵ_in)," , N = ",i_in," , ρ = ",@sprintf("%.2f",ρ_in)," )")
            end

        else

            # Livolant acceleration
            if 𝒜 == "livolant" && mod(i_in,3) == 0
                𝚽l⁺ = livolant(𝚽l,𝚽l⁻[1,:,:,:,:,:],𝚽l⁻[2,:,:,:,:,:])
                𝚽l⁻[2,:,:,:,:,:] = 𝚽l⁻[1,:,:,:,:,:]
                𝚽l⁻[1,:,:,:,:,:] = 𝚽l
                𝚽l .= 𝚽l⁺
            else
                𝚽l⁻[2,:,:,:,:,:] = 𝚽l⁻[1,:,:,:,:,:]
                𝚽l⁻[1,:,:,:,:,:] = 𝚽l
            end
            
            # Save flux solution and go to next iteration
            i_in += 1

        end
    end
    return 𝚽l,𝚽E12_temp,ρ_in,Ntot
end