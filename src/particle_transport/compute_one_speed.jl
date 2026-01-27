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
- `𝒜::String` : acceleration method for in-group iterations.
- `is_CUDA::Bool` : boolean for CUDA parallelism.
- `Ntot::Int64` : accumulator for the total number of in-group iterations.
- `is_EM::Bool` : boolean for electromagnetic fields.
- `ℳ_EM::Array{Float64}` : electromagnetic scattering matrix.
- `𝒲::Array{Float64}` : weighting constants.

# Output Argument(s)
- `𝚽l::Array{Float64}`: Legendre components of the in-cell flux.
- `𝚽E12::Array{Float64}`: outgoing flux along the energy axis.
- `ρ_in::Float64`: estimated spectral radius.
- `Ntot::Int64` : accumulator for the total number of in-group iterations.

# Reference(s)
- Larsen and Morel (2010) : Advances in Discrete-Ordinates Methodology.

"""
function compute_one_speed(𝚽l::Array{Float64},Qlout::Array{Float64},Σt::Vector{Float64},Σs::Array{Float64},mat::Array{Int64,3},ndims::Int64,Nd::Int64,ig::Int64,Ns::Vector{Int64},Δs::Vector{Vector{Float64}},Ω::Vector{Vector{Float64}},Mn::Array{Float64,2},Dn::Array{Float64,2},Np::Int64,Nq::Int64,pl::Vector{Int64},pm::Vector{Int64},Mn_surf::Vector{Array{Float64}},Dn_surf::Vector{Array{Float64}},Np_surf::Int64,n_to_n⁺::Vector{Vector{Int64}},𝒪::Vector{Int64},Nm::Vector{Int64},isFC::Bool,C::Vector{Float64},ω::Vector{Array{Float64}},I_max::Int64,ϵ_max::Float64,sources::Array{Union{Array{Float64},Float64}},isAdapt::Bool,isCSD::Bool,solver::Int64,E::Float64,ΔE::Float64,𝚽E12::Array{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},T::Vector{Float64},ℳ::Array{Float64},𝒜::String,Ntot::Int64,is_EM::Bool,ℳ_EM::Array{Float64},𝒲::Array{Float64},Mll::Array{Float64},is_SN::Bool,is_PN::Bool,is_SPH::Bool,𝒩,boundary_conditions::Vector{Int64},Np_source)

    # Flux Initialization
    𝚽E12_temp = Array{Float64}(undef)
    if isCSD
        if is_SN
            𝚽E12_temp = zeros(Nd,Nm[4],Ns[1],Ns[2],Ns[3])
        elseif is_PN
            𝚽E12_temp = zeros(Np,Nm[4],Ns[1],Ns[2],Ns[3])
        else
            error("Unknown angular discretization method.")
        end
    end
    N⁻ = 2
    𝚽l⁻ = zeros(N⁻,Np,Nm[5],Ns[1],Ns[2],Ns[3])

    # Boundary conditions initialization
    if is_SN
        if ndims == 1
            𝚽x12_in = zeros(Np_surf,Nm[1],2)
            𝚽x12_out = copy(𝚽x12_in)
            𝚽x12_temp = copy(𝚽x12_in)
        elseif ndims == 2
            𝚽x12_in = zeros(Np_surf,Nm[1],2,Ns[2])
            𝚽y12_in = zeros(Np_surf,Nm[2],2,Ns[1])
            𝚽x12_out = copy(𝚽x12_in)
            𝚽y12_out = copy(𝚽y12_in)
            𝚽x12_temp = copy(𝚽x12_in)
            𝚽y12_temp = copy(𝚽y12_in)
        elseif ndims == 3
            𝚽x12_in = zeros(Np_surf,Nm[1],2,Ns[2],Ns[3])
            𝚽y12_in = zeros(Np_surf,Nm[2],2,Ns[1],Ns[3])
            𝚽z12_in = zeros(Np_surf,Nm[3],2,Ns[1],Ns[2])
            𝚽x12_out = copy(𝚽x12_in)
            𝚽y12_out = copy(𝚽y12_in)
            𝚽z12_out = copy(𝚽z12_in)
            𝚽x12_temp = copy(𝚽x12_in)
            𝚽y12_temp = copy(𝚽y12_in)
            𝚽z12_temp = copy(𝚽z12_in)
        else
            error("Dimension is not 1, 2 or 3.")
        end
    else
        if ndims == 1
            𝚽x12⁺ = zeros(Nq,Nm[1],2)
            𝚽x12⁻ = zeros(Nq,Nm[1],2)
        elseif ndims == 2
            𝚽x12⁺⁺ = zeros(Nq,Nm[1],2,Ns[2])
            𝚽x12⁺⁻ = zeros(Nq,Nm[1],2,Ns[2])
            𝚽x12⁻⁺ = zeros(Nq,Nm[1],2,Ns[2])
            𝚽x12⁻⁻ = zeros(Nq,Nm[1],2,Ns[2])
            𝚽y12⁺⁺ = zeros(Nq,Nm[2],2,Ns[1])
            𝚽y12⁺⁻ = zeros(Nq,Nm[2],2,Ns[1])
            𝚽y12⁻⁺ = zeros(Nq,Nm[2],2,Ns[1])
            𝚽y12⁻⁻ = zeros(Nq,Nm[2],2,Ns[1])
        elseif ndims == 3
            𝚽x12⁺⁺⁺ = zeros(Nq,Nm[1],2,Ns[2],Ns[3])
            𝚽x12⁺⁺⁻ = zeros(Nq,Nm[1],2,Ns[2],Ns[3])
            𝚽x12⁺⁻⁺ = zeros(Nq,Nm[1],2,Ns[2],Ns[3])
            𝚽x12⁺⁻⁻ = zeros(Nq,Nm[1],2,Ns[2],Ns[3])
            𝚽x12⁻⁺⁺ = zeros(Nq,Nm[1],2,Ns[2],Ns[3])
            𝚽x12⁻⁺⁻ = zeros(Nq,Nm[1],2,Ns[2],Ns[3])
            𝚽x12⁻⁻⁺ = zeros(Nq,Nm[1],2,Ns[2],Ns[3])
            𝚽x12⁻⁻⁻ = zeros(Nq,Nm[1],2,Ns[2],Ns[3])
            𝚽y12⁺⁺⁺ = zeros(Nq,Nm[2],2,Ns[1],Ns[3])
            𝚽y12⁺⁺⁻ = zeros(Nq,Nm[2],2,Ns[1],Ns[3])
            𝚽y12⁺⁻⁺ = zeros(Nq,Nm[2],2,Ns[1],Ns[3])
            𝚽y12⁺⁻⁻ = zeros(Nq,Nm[2],2,Ns[1],Ns[3])
            𝚽y12⁻⁺⁺ = zeros(Nq,Nm[2],2,Ns[1],Ns[3])
            𝚽y12⁻⁺⁻ = zeros(Nq,Nm[2],2,Ns[1],Ns[3])
            𝚽y12⁻⁻⁺ = zeros(Nq,Nm[2],2,Ns[1],Ns[3])
            𝚽y12⁻⁻⁻ = zeros(Nq,Nm[2],2,Ns[1],Ns[3])
            𝚽z12⁺⁺⁺ = zeros(Nq,Nm[3],2,Ns[1],Ns[2])
            𝚽z12⁺⁺⁻ = zeros(Nq,Nm[3],2,Ns[1],Ns[2])
            𝚽z12⁺⁻⁺ = zeros(Nq,Nm[3],2,Ns[1],Ns[2])
            𝚽z12⁺⁻⁻ = zeros(Nq,Nm[3],2,Ns[1],Ns[2])
            𝚽z12⁻⁺⁺ = zeros(Nq,Nm[3],2,Ns[1],Ns[2])
            𝚽z12⁻⁺⁻ = zeros(Nq,Nm[3],2,Ns[1],Ns[2])
            𝚽z12⁻⁻⁺ = zeros(Nq,Nm[3],2,Ns[1],Ns[2])
            𝚽z12⁻⁻⁻ = zeros(Nq,Nm[3],2,Ns[1],Ns[2])
        else
            error("Dimension is not 1, 2 or 3.")
        end
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

        # Electromagnetic source
        if is_EM
            for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3]), is in range(1,Nm[5])
                for p in range(1,Np), q in range(1,Np)
                    Ql[p,is,ix,iy,iz] += ℳ_EM[p,q] * 𝚽l[q,is,ix,iy,iz]
                end
            end
        end

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
        if is_PN
            𝚽E12_temp .= 0
            if ndims == 1
                𝚽⁺ = zeros(Nq,Nm[5],Ns[1])
                𝚽⁻ = zeros(Nq,Nm[5],Ns[1])
                Q⁺ = zeros(Nq,Nm[5],Ns[1])
                Q⁻ = zeros(Nq,Nm[5],Ns[1])
                𝚽E12⁺ = zeros(Nq,Nm[4],Ns[1])
                𝚽E12⁻ = zeros(Nq,Nm[4],Ns[1])
                for p in range(1,Np), q in range(1,Nq)
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
                𝚽⁺,𝚽E12⁺,𝚽x12⁺ = pn_sweep_1D(1,𝚽⁺,Q⁺,Σt,mat[:,1,1],Ns[1],Δs[1],Nq,Np_source,Np_surf,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁺,𝒲,isFC,isCSD,𝒩[1],𝚽x12⁺,𝚽x12⁻,boundary_conditions)
                𝚽⁻,𝚽E12⁻,𝚽x12⁻ = pn_sweep_1D(-1,𝚽⁻,Q⁻,Σt,mat[:,1,1],Ns[1],Δs[1],Nq,Np_source,Np_surf,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁻,𝒲,isFC,isCSD,𝒩[1],𝚽x12⁺,𝚽x12⁻,boundary_conditions)
                for p in range(1,Np), q in range(1,Nq)
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
                𝚽⁺⁺ = zeros(Nq,Nm[5],Ns[1],Ns[2])
                𝚽⁺⁻ = zeros(Nq,Nm[5],Ns[1],Ns[2])
                𝚽⁻⁺ = zeros(Nq,Nm[5],Ns[1],Ns[2])
                𝚽⁻⁻ = zeros(Nq,Nm[5],Ns[1],Ns[2])
                Q⁺⁺ = zeros(Nq,Nm[5],Ns[1],Ns[2])
                Q⁺⁻ = zeros(Nq,Nm[5],Ns[1],Ns[2])
                Q⁻⁺ = zeros(Nq,Nm[5],Ns[1],Ns[2])
                Q⁻⁻ = zeros(Nq,Nm[5],Ns[1],Ns[2])
                𝚽E12⁺⁺ = zeros(Nq,Nm[5],Ns[1],Ns[2])
                𝚽E12⁺⁻ = zeros(Nq,Nm[5],Ns[1],Ns[2])
                𝚽E12⁻⁺ = zeros(Nq,Nm[5],Ns[1],Ns[2])
                𝚽E12⁻⁻ = zeros(Nq,Nm[5],Ns[1],Ns[2])
                for ix in range(1,Ns[1]), iy in range(1,Ns[2]), p in range(1,Np), q in range(1,Nq)
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
                𝚽⁺⁺,𝚽E12⁺⁺,𝚽x12⁺⁺,𝚽y12⁺⁺ = pn_sweep_2D(1,1,𝚽⁺⁺,Q⁺⁺,Σt,mat[:,:,1],Ns[1],Ns[2],Δs[1],Δs[2],Nq,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁺⁺,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝚽x12⁺⁺,𝚽x12⁻⁺,𝚽y12⁺⁺,𝚽y12⁺⁻,boundary_conditions)
                𝚽⁺⁻,𝚽E12⁺⁻,𝚽x12⁺⁻,𝚽y12⁺⁻ = pn_sweep_2D(1,-1,𝚽⁺⁻,Q⁺⁻,Σt,mat[:,:,1],Ns[1],Ns[2],Δs[1],Δs[2],Nq,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁺⁻,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝚽x12⁺⁻,𝚽x12⁻⁻,𝚽y12⁺⁻,𝚽y12⁺⁺,boundary_conditions)
                𝚽⁻⁺,𝚽E12⁻⁺,𝚽x12⁻⁺,𝚽y12⁻⁺ = pn_sweep_2D(-1,1,𝚽⁻⁺,Q⁻⁺,Σt,mat[:,:,1],Ns[1],Ns[2],Δs[1],Δs[2],Nq,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁻⁺,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝚽x12⁻⁺,𝚽x12⁺⁺,𝚽y12⁻⁺,𝚽y12⁻⁻,boundary_conditions)
                𝚽⁻⁻,𝚽E12⁻⁻,𝚽x12⁻⁻,𝚽y12⁻⁻ = pn_sweep_2D(-1,-1,𝚽⁻⁻,Q⁻⁻,Σt,mat[:,:,1],Ns[1],Ns[2],Δs[1],Δs[2],Nq,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁻⁻,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝚽x12⁻⁻,𝚽x12⁺⁻,𝚽y12⁻⁻,𝚽y12⁻⁺,boundary_conditions)
                for ix in range(1,Ns[1]), iy in range(1,Ns[2]), p in range(1,Np), q in range(1,Nq)
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
                𝚽⁺⁺⁺ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                𝚽⁺⁺⁻ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                𝚽⁺⁻⁺ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                𝚽⁺⁻⁻ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                𝚽⁻⁺⁺ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                𝚽⁻⁺⁻ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                𝚽⁻⁻⁺ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                𝚽⁻⁻⁻ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                Q⁺⁺⁺ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                Q⁺⁺⁻ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                Q⁺⁻⁺ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                Q⁺⁻⁻ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                Q⁻⁺⁺ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                Q⁻⁺⁻ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                Q⁻⁻⁺ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                Q⁻⁻⁻ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                𝚽E12⁺⁺⁺ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                𝚽E12⁺⁺⁻ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                𝚽E12⁺⁻⁺ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                𝚽E12⁺⁻⁻ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                𝚽E12⁻⁺⁺ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                𝚽E12⁻⁺⁻ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                𝚽E12⁻⁻⁺ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                𝚽E12⁻⁻⁻ = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3])
                for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3]), p in range(1,Np), q in range(1,Nq)
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
                𝚽⁺⁺⁺,𝚽E12⁺⁺⁺,𝚽x12⁺⁺⁺,𝚽y12⁺⁺⁺,𝚽z12⁺⁺⁺ = pn_sweep_3D(1,1,1,𝚽⁺⁺⁺,Q⁺⁺⁺,Σt,mat,Ns[1],Ns[2],Ns[3],Δs[1],Δs[2],Δs[3],Nq,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁺⁺⁺,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝒩[3],𝚽x12⁺⁺⁺,𝚽x12⁻⁺⁺,𝚽y12⁺⁺⁺,𝚽y12⁺⁻⁺,𝚽z12⁺⁺⁺,𝚽z12⁺⁺⁻,boundary_conditions)
                𝚽⁺⁺⁻,𝚽E12⁺⁺⁻,𝚽x12⁺⁺⁻,𝚽y12⁺⁺⁻,𝚽z12⁺⁺⁻ = pn_sweep_3D(1,1,-1,𝚽⁺⁺⁻,Q⁺⁺⁻,Σt,mat,Ns[1],Ns[2],Ns[3],Δs[1],Δs[2],Δs[3],Nq,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁺⁺⁻,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝒩[3],𝚽x12⁺⁺⁻,𝚽x12⁻⁺⁻,𝚽y12⁺⁺⁻,𝚽y12⁺⁻⁻,𝚽z12⁺⁺⁻,𝚽z12⁺⁺⁺,boundary_conditions)
                𝚽⁺⁻⁺,𝚽E12⁺⁻⁺,𝚽x12⁺⁻⁺,𝚽y12⁺⁻⁺,𝚽z12⁺⁻⁺ = pn_sweep_3D(1,-1,1,𝚽⁺⁻⁺,Q⁺⁻⁺,Σt,mat,Ns[1],Ns[2],Ns[3],Δs[1],Δs[2],Δs[3],Nq,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁺⁻⁺,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝒩[3],𝚽x12⁺⁻⁺,𝚽x12⁻⁻⁺,𝚽y12⁺⁻⁺,𝚽y12⁺⁺⁺,𝚽z12⁺⁻⁺,𝚽z12⁺⁻⁻,boundary_conditions)
                𝚽⁺⁻⁻,𝚽E12⁺⁻⁻,𝚽x12⁺⁻⁻,𝚽y12⁺⁻⁻,𝚽z12⁺⁻⁻ = pn_sweep_3D(1,-1,-1,𝚽⁺⁻⁻,Q⁺⁻⁻,Σt,mat,Ns[1],Ns[2],Ns[3],Δs[1],Δs[2],Δs[3],Nq,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁺⁻⁻,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝒩[3],𝚽x12⁺⁻⁻,𝚽x12⁻⁻⁻,𝚽y12⁺⁻⁻,𝚽y12⁺⁺⁻,𝚽z12⁺⁻⁻,𝚽z12⁺⁻⁺,boundary_conditions)
                𝚽⁻⁺⁺,𝚽E12⁻⁺⁺,𝚽x12⁻⁺⁺,𝚽y12⁻⁺⁺,𝚽z12⁻⁺⁺ = pn_sweep_3D(-1,1,1,𝚽⁻⁺⁺,Q⁻⁺⁺,Σt,mat,Ns[1],Ns[2],Ns[3],Δs[1],Δs[2],Δs[3],Nq,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁻⁺⁺,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝒩[3],𝚽x12⁻⁺⁺,𝚽x12⁺⁺⁺,𝚽y12⁻⁺⁺,𝚽y12⁻⁻⁺,𝚽z12⁻⁺⁺,𝚽z12⁻⁺⁻,boundary_conditions)
                𝚽⁻⁺⁻,𝚽E12⁻⁺⁻,𝚽x12⁻⁺⁻,𝚽y12⁻⁺⁻,𝚽z12⁻⁺⁻ = pn_sweep_3D(-1,1,-1,𝚽⁻⁺⁻,Q⁻⁺⁻,Σt,mat,Ns[1],Ns[2],Ns[3],Δs[1],Δs[2],Δs[3],Nq,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁻⁺⁻,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝒩[3],𝚽x12⁻⁺⁻,𝚽x12⁺⁺⁻,𝚽y12⁻⁺⁻,𝚽y12⁻⁻⁻,𝚽z12⁻⁺⁻,𝚽z12⁻⁺⁺,boundary_conditions)
                𝚽⁻⁻⁺,𝚽E12⁻⁻⁺,𝚽x12⁻⁻⁺,𝚽y12⁻⁻⁺,𝚽z12⁻⁻⁺ = pn_sweep_3D(-1,-1,1,𝚽⁻⁻⁺,Q⁻⁻⁺,Σt,mat,Ns[1],Ns[2],Ns[3],Δs[1],Δs[2],Δs[3],Nq,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁻⁻⁺,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝒩[3],𝚽x12⁻⁻⁺,𝚽x12⁺⁻⁺,𝚽y12⁻⁻⁺,𝚽y12⁻⁺⁺,𝚽z12⁻⁻⁺,𝚽z12⁻⁻⁻,boundary_conditions)
                𝚽⁻⁻⁻,𝚽E12⁻⁻⁻,𝚽x12⁻⁻⁻,𝚽y12⁻⁻⁻,𝚽z12⁻⁻⁻ = pn_sweep_3D(-1,-1,-1,𝚽⁻⁻⁻,Q⁻⁻⁻,Σt,mat,Ns[1],Ns[2],Ns[3],Δs[1],Δs[2],Δs[3],Nq,Np_source,𝒪,Nm,C,ω,sources,S⁻,S⁺,S,𝚽E12⁻⁻⁻,𝒲,isFC,isCSD,𝒩[1],𝒩[2],𝒩[3],𝚽x12⁻⁻⁻,𝚽x12⁺⁻⁻,𝚽y12⁻⁻⁻,𝚽y12⁻⁺⁻,𝚽z12⁻⁻⁻,𝚽z12⁻⁻⁺,boundary_conditions)
                for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3]), p in range(1,Np), q in range(1,Nq)
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
        elseif is_SN
            for n in range(1,Nd)
                if isCSD 𝚽E12_out = 𝚽E12[n,:,:,:,:] else 𝚽E12_out = Array{Float64}(undef) end
                if ndims == 1
                    nx⁻ = n_to_n⁺[1][n]
                    nx⁺ = n_to_n⁺[2][n]
                    if nx⁻ != 0
                        Mnx⁻ = Mn_surf[1][nx⁻,:]
                        Dnx⁻ = Dn_surf[1][:,nx⁻]
                    elseif nx⁺ != 0
                        Mnx⁻ = Mn_surf[2][nx⁺,:]
                        Dnx⁻ = Dn_surf[2][:,nx⁺]
                    else
                        Mnx⁻ = zeros(Np)
                        Dnx⁻ = zeros(Np)
                    end
                    𝚽l[:,:,:,1,1], 𝚽E12_out,𝚽x12_out = sn_sweep_1D(𝚽l[:,:,:,1,1],Ql[:,:,:,1,1],Σt,mat[:,1,1],Ns[1],Δs[1],Ω[1][n],Mn[n,:],Dn[:,n],Np,Mnx⁻,Dnx⁻,Np_surf,𝒪,Nm,C,ω,sources,isAdapt,isCSD,ΔE,𝚽E12_out,S⁻,S⁺,S,𝒲,isFC,𝚽x12_in,boundary_conditions,Np_source)
                elseif ndims == 2
                    nx⁻ = n_to_n⁺[1][n]
                    nx⁺ = n_to_n⁺[2][n]
                    ny⁻ = n_to_n⁺[3][n]
                    ny⁺ = n_to_n⁺[4][n]
                    if nx⁻ != 0
                        Mnx⁻ = Mn_surf[1][nx⁻,:]
                        Dnx⁻ = Dn_surf[1][:,nx⁻]
                    elseif nx⁺ != 0
                        Mnx⁻ = Mn_surf[2][nx⁺,:]
                        Dnx⁻ = Dn_surf[2][:,nx⁺]
                    else
                        Mnx⁻ = zeros(Np)
                        Dnx⁻ = zeros(Np)
                    end
                    if ny⁻ != 0
                        Mny⁻ = Mn_surf[3][ny⁻,:]
                        Dny⁻ = Dn_surf[3][:,ny⁻]
                    elseif ny⁺ != 0
                        Mny⁻ = Mn_surf[4][ny⁺,:]
                        Dny⁻ = Dn_surf[4][:,ny⁺]
                    else
                        Mny⁻ = zeros(Np)
                        Dny⁻ = zeros(Np)
                    end
                    𝚽l[:,:,:,:,1],𝚽E12_out,𝚽x12_out,𝚽y12_out = sn_sweep_2D(𝚽l[:,:,:,:,1],Ql[:,:,:,:,1],Σt,mat[:,:,1],Ns[1:2],Δs[1:2],[Ω[1][n],Ω[2][n]],Mn[n,:],Dn[:,n],Np,Mnx⁻,Dnx⁻,Mny⁻,Dny⁻,Np_surf,𝒪,Nm,C,ω,sources,isAdapt,isCSD,ΔE,𝚽E12_out,S⁻,S⁺,S,𝒲,isFC,𝚽x12_in,𝚽y12_in,boundary_conditions,Np_source)
                elseif ndims == 3
                    nx⁻ = n_to_n⁺[1][n]
                    nx⁺ = n_to_n⁺[2][n]
                    ny⁻ = n_to_n⁺[3][n]
                    ny⁺ = n_to_n⁺[4][n]
                    nz⁻ = n_to_n⁺[5][n]
                    nz⁺ = n_to_n⁺[6][n]
                    if nx⁻ != 0
                        Mnx⁻ = Mn_surf[1][nx⁻,:]
                        Dnx⁻ = Dn_surf[1][:,nx⁻]
                    elseif nx⁺ != 0
                        Mnx⁻ = Mn_surf[2][nx⁺,:]
                        Dnx⁻ = Dn_surf[2][:,nx⁺]
                    else
                        Mnx⁻ = zeros(Np)
                        Dnx⁻ = zeros(Np)
                    end
                    if ny⁻ != 0
                        Mny⁻ = Mn_surf[3][ny⁻,:]
                        Dny⁻ = Dn_surf[3][:,ny⁻]
                    elseif ny⁺ != 0
                        Mny⁻ = Mn_surf[4][ny⁺,:]
                        Dny⁻ = Dn_surf[4][:,ny⁺]
                    else
                        Mny⁻ = zeros(Np)
                        Dny⁻ = zeros(Np)
                    end
                    if nz⁻ != 0
                        Mnz⁻ = Mn_surf[5][nz⁻,:]
                        Dnz⁻ = Dn_surf[5][:,nz⁻]
                    elseif nz⁺ != 0
                        Mnz⁻ = Mn_surf[6][nz⁺,:]
                        Dnz⁻ = Dn_surf[6][:,nz⁺]
                    else
                        Mnz⁻ = zeros(Np)
                        Dnz⁻ = zeros(Np)
                    end
                    𝚽l,𝚽E12_out,𝚽x12_out,𝚽y12_out,𝚽z12_out = sn_sweep_3D(𝚽l,Ql,Σt,mat,Ns,Δs,[Ω[1][n],Ω[2][n],Ω[3][n]],Mn[n,:],Dn[:,n],Np,Mnx⁻,Dnx⁻,Mny⁻,Dny⁻,Mnz⁻,Dnz⁻,Np_surf,𝒪,Nm,C,ω,sources,isAdapt,isCSD,ΔE,𝚽E12_out,S⁻,S⁺,S,𝒲,isFC,𝚽x12_in,𝚽y12_in,𝚽z12_in,boundary_conditions,Np_source)
                else
                    error("Dimension is not 1, 2 or 3.")
                end
                𝚽x12_temp += 𝚽x12_out
                if ndims >= 2 𝚽y12_temp += 𝚽y12_out end
                if ndims >= 3 𝚽z12_temp += 𝚽z12_out end
                if isCSD 𝚽E12_temp[n,:,:,:,:] = 𝚽E12_out end
            end
            𝚽x12_in = copy(𝚽x12_temp); 𝚽x12_temp .= 0.0
            if ndims >= 2 𝚽y12_in = copy(𝚽y12_temp); 𝚽y12_temp .= 0.0 end
            if ndims >= 3 𝚽z12_in = copy(𝚽z12_temp); 𝚽z12_temp .= 0.0 end
        else
            error("Unknown angular discretization method.")
        end
        
        #----
        # Verification of convergence of the one-group flux
        #----  
        ϵ_in = 0.0
        if (solver ∉ [5,6]) ϵ_in = maximum(vec(abs.((𝚽l[1,1,:,:,:] .- 𝚽l⁻[1,1,1,:,:,:])./max.(abs.(𝚽l[1,1,:,:,:]),1e-16)))) end
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



