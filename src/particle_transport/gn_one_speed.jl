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
function gn_one_speed(𝚽l::Array{Float64},Qlout::Array{Float64},Σt::Vector{Float64},Σs::Array{Float64},mat::Array{Int64,3},Ndims::Int64,ig::Int64,Ns::Vector{Int64},Δs::Vector{Vector{Float64}},Np::Int64,Nq::Int64,pl::Vector{Int64},pm::Vector{Int64},Np_surf::Int64,𝒪::Vector{Int64},Nm::Vector{Int64},isFC::Bool,C::Vector{Float64},ω::Vector{Vector{Float64}},I_max::Int64,ϵ_max::Float64,sources::Array{Union{Array{Float64},Float64}},isCSD::Bool,solver::Int64,𝚽E12::Array{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},T::Vector{Float64},ℳ::Array{Float64},𝒜::String,Ntot::Int64,𝒲::Array{Float64},Mll::Array{Float64},is_SPH::Bool,𝒩::Array{Float64},boundary_conditions::Vector{Int64},Np_source::Int64,Nv::Int64,Mll_surf::Array{Float64})

    # Flux Initialization
    𝚽E12_temp = Array{Float64}(undef)
    if isCSD
        𝚽E12_temp = zeros(Np,Nm[4],Ns[1],Ns[2],Ns[3])
    end
    N⁻ = 2
    𝚽l⁻ = zeros(N⁻,Np,Nm[5],Ns[1],Ns[2],Ns[3])
    sx = [1,1,1,1,-1,-1,-1,-1]
    if (Ndims > 1) sy = [1,1,-1,-1,1,1,-1,-1] end
    if (Ndims > 2) sz = [1,-1,1,-1,1,-1,1,-1] end

    # Boundary sources
    if Ndims == 1
        sources_q = zeros(Nq,2*Ndims,8,Nv,Nv)
    else
        sources_q = Array{Union{Float64,Array{Float64}}}(undef,Nq,2*Ndims,8,Nv,Nv)
        for q in range(1,Nq), u in range(1,8), v in range(1,Nv)
            Nw = Int(-sx[u]*v + (sx[u]+1)/2*(Nv+1))
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
        Nw = Int(-sx[u]*v + (sx[u]+1)/2*(Nv+1))
        for w in range(1,Nw)
            for ib in range(1,2)
                if Ndims == 1
                    sources_q[q,ib,u,v,w] += sources[p,ib] * Mll_surf[p,q,u,v,w,ib]
                elseif Ndims == 2
                    for iy in range(1,Ns[2])
                        sources_q[q,ib,u,v,w][iy] += sources[p,ib][iy] * Mll_surf[p,q,u,v,w,ib]
                    end
                    for ix in range(1,Ns[1])
                        sources_q[q,ib+2,u,v,w][ix] += sources[p,ib+2][ix] * Mll_surf[p,q,u,v,w,ib+2]
                    end
                elseif Ndims == 3
                    for iy in range(1,Ns[2]), iz in range(1,Ns[3])
                        sources_q[q,ib,u,v,w][iy,iz] += sources[p,ib][iy,iz] * Mll_surf[p,q,u,v,w,ib]
                    end
                    for ix in range(1,Ns[1]), iz in range(1,Ns[3])
                        sources_q[q,ib+2,u,v,w][ix,iz] += sources[p,ib+2][ix,iz] * Mll_surf[p,q,u,v,w,ib+2]
                    end
                    for ix in range(1,Ns[1]), iy in range(1,Ns[2])
                        sources_q[q,ib+4,u,v,w][ix,iy] += sources[p,ib+4][ix,iy] * Mll_surf[p,q,u,v,w,ib+4]
                    end
                else
                    error("Invalid number of dimensions.")
                end
            end
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
        if Ndims == 1
            𝚽_q = zeros(Nq,Nm[5],Ns[1],8,Nv,Nv)
            Q_q = zeros(Nq,Nm[5],Ns[1],8,Nv,Nv)
            𝚽E12_q = zeros(Nq,Nm[4],Ns[1],8,Nv,Nv)
            for p in range(1,Np), q in range(1,Nq), u in range(1,8), v in range(1,Nv)
                Nw = Int(-sx[u]*v + (sx[u]+1)/2*(Nv+1))
                for w in range(1,Nw)
                    if is_SPH factor = (2*pl[p]+1)/(4*π) else factor = (2*pl[p]+1)/2 end
                    for is in range(1,Nm[5]), ix in range(1,Ns[1])
                        Q_q[q,is,ix,u,v,w] += factor * Ql[p,is,ix,1,1] * Mll[p,q,u,v,w]
                    end
                    if isCSD
                        for is in range(1,Nm[4]), ix in range(1,Ns[1])
                            𝚽E12_q[q,is,ix,u,v,w] += factor * 𝚽E12[p,is,ix,1,1] * Mll[p,q,u,v,w]
                        end
                    end
                end
            end
            for u in range(1,8), v in range(1,Nv)
                Nw = Int(-sx[u]*v + (sx[u]+1)/2*(Nv+1))
                for w in range(1,Nw)
                    𝚽_q[:,:,:,u,v,w],𝚽E12_q[:,:,:,u,v,w] = gn_sweep_1D(sx[u],𝚽_q[:,:,:,u,v,w],Q_q[:,:,:,u,v,w],Σt,mat[:,1,1],Ns[1],Δs[1],Nq,Np_source,Np_surf,𝒪,Nm,C,ω,sources_q[:,:,u,v,w],S⁻,S⁺,S,𝚽E12_q[:,:,:,u,v,w],𝒲,isFC,isCSD,𝒩[:,:,1,u,v,w])
                end
            end
            for p in range(1,Np), q in range(1,Nq), u in range(1,8), v in range(1,Nv)
                Nw = Int(-sx[u]*v + (sx[u]+1)/2*(Nv+1))
                for w in range(1,Nw)
                    for is in range(1,Nm[5]), ix in range(1,Ns[1])
                        𝚽l[p,is,ix,1,1] += Mll[p,q,u,v,w] * 𝚽_q[q,is,ix,u,v,w]
                    end
                    if isCSD
                        for is in range(1,Nm[4]), ix in range(1,Ns[1])
                            𝚽E12_temp[p,is,ix,1,1] += Mll[p,q,u,v,w] * 𝚽E12_q[q,is,ix,u,v,w]
                        end
                    end
                end
            end
        elseif Ndims == 2
            𝚽_q = zeros(Nq,Nm[5],Ns[1],Ns[2],8,Nv,Nv)
            Q_q = zeros(Nq,Nm[5],Ns[1],Ns[2],8,Nv,Nv)
            𝚽E12_q = zeros(Nq,Nm[4],Ns[1],Ns[2],8,Nv,Nv)
            for p in range(1,Np), q in range(1,Nq), u in range(1,8), v in range(1,Nv)
                Nw = Int(-sx[u]*v + (sx[u]+1)/2*(Nv+1))
                for w in range(1,Nw)
                    factor = (2*pl[p]+1)/(4*π)
                    for is in range(1,Nm[5]), ix in range(1,Ns[1]), iy in range(1,Ns[2])
                        Q_q[q,is,ix,iy,u,v,w] += factor * Ql[p,is,ix,iy,1] * Mll[p,q,u,v,w]
                    end
                    if isCSD
                        for is in range(1,Nm[4]), ix in range(1,Ns[1]), iy in range(1,Ns[2])
                            𝚽E12_q[q,is,ix,iy,u,v,w] += factor * 𝚽E12[p,is,ix,iy,1] * Mll[p,q,u,v,w]
                        end
                    end
                end
            end
            for u in range(1,8), v in range(1,Nv)
                Nw = Int(-sx[u]*v + (sx[u]+1)/2*(Nv+1))
                for w in range(1,Nw)
                    𝚽_q[:,:,:,:,u,v,w],𝚽E12_q[:,:,:,:,u,v,w] = gn_sweep_2D(sx[u],sy[u],𝚽_q[:,:,:,:,u,v,w],Q_q[:,:,:,:,u,v,w],Σt,mat[:,:,1],Ns[1],Ns[2],Δs[1],Δs[2],Nq,Np_source,𝒪,Nm,C,ω,sources_q[:,:,u,v,w],S⁻,S⁺,S,𝚽E12_q[:,:,:,:,u,v,w],𝒲,isFC,isCSD,𝒩[:,:,1,u,v,w],𝒩[:,:,2,u,v,w])
                end
            end
            for p in range(1,Np), q in range(1,Nq), u in range(1,8), v in range(1,Nv)
                Nw = Int(-sx[u]*v + (sx[u]+1)/2*(Nv+1))
                for w in range(1,Nw)
                    for is in range(1,Nm[5]), ix in range(1,Ns[1]), iy in range(1,Ns[2])
                        𝚽l[p,is,ix,iy,1] += Mll[p,q,u,v,w] * 𝚽_q[q,is,ix,iy,u,v,w]
                    end
                    if isCSD
                        for is in range(1,Nm[4]), ix in range(1,Ns[1]), iy in range(1,Ns[2])
                            𝚽E12_temp[p,is,ix,iy,1] += Mll[p,q,u,v,w] * 𝚽E12_q[q,is,ix,iy,u,v,w]
                        end
                    end
                end
            end
        elseif Ndims == 3
            𝚽_q = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3],8,Nv,Nv)
            Q_q = zeros(Nq,Nm[5],Ns[1],Ns[2],Ns[3],8,Nv,Nv)
            𝚽E12_q = zeros(Nq,Nm[4],Ns[1],Ns[2],Ns[3],8,Nv,Nv)
            for p in range(1,Np), q in range(1,Nq), u in range(1,8), v in range(1,Nv)
                Nw = Int(-sx[u]*v + (sx[u]+1)/2*(Nv+1))
                for w in range(1,Nw)
                    factor = (2*pl[p]+1)/(4*π)
                    for is in range(1,Nm[5]), ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3])
                        Q_q[q,is,ix,iy,iz,u,v,w] += factor * Ql[p,is,ix,iy,iz] * Mll[p,q,u,v,w]
                    end
                    if isCSD
                        for is in range(1,Nm[4]), ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3])
                            𝚽E12_q[q,is,ix,iy,iz,u,v,w] += factor * 𝚽E12[p,is,ix,iy,iz] * Mll[p,q,u,v,w]
                        end
                    end
                end
            end
            for u in range(1,8), v in range(1,Nv)
                Nw = Int(-sx[u]*v + (sx[u]+1)/2*(Nv+1))
                for w in range(1,Nw)
                    𝚽_q[:,:,:,:,:,u,v,w],𝚽E12_q[:,:,:,:,:,u,v,w] = gn_sweep_3D(sx[u],sy[u],sz[u],𝚽_q[:,:,:,:,:,u,v,w],Q_q[:,:,:,:,:,u,v,w],Σt,mat,Ns[1],Ns[2],Ns[3],Δs[1],Δs[2],Δs[3],Nq,Np_source,𝒪,Nm,C,ω,sources_q[:,:,u,v,w],S⁻,S⁺,S,𝚽E12_q[:,:,:,:,:,u,v,w],𝒲,isFC,isCSD,𝒩[:,:,1,u,v,w],𝒩[:,:,2,u,v,w],𝒩[:,:,3,u,v,w])
                end
            end
            for p in range(1,Np), q in range(1,Nq), u in range(1,8), v in range(1,Nv)
                Nw = Int(-sx[u]*v + (sx[u]+1)/2*(Nv+1))
                for w in range(1,Nw)
                    for is in range(1,Nm[5]), ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3])
                        𝚽l[p,is,ix,iy,iz] += Mll[p,q,u,v,w] * 𝚽_q[q,is,ix,iy,iz,u,v,w]
                    end
                    if isCSD
                        for is in range(1,Nm[4]), ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3])
                            𝚽E12_temp[p,is,ix,iy,iz] += Mll[p,q,u,v,w] * 𝚽E12_q[q,is,ix,iy,iz,u,v,w]
                        end
                    end
                end
            end
        else
            error("Geometry dimension is either 1D, 2D or 3D.")
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