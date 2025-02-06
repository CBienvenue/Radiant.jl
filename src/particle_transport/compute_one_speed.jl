"""
    compute_one_speed(𝚽ℓ::Array{Float64},Qℓout::Array{Float64},Σt::Vector{Float64},
    Σs::Array{Float64},mat::Array{Int64,3},ndims::Int64,N::Int64,ig::Int64,
    Ns::Vector{Int64},Δs::Vector{Vector{Float64}},Ω::Union{Vector{Vector{Float64}},
    Vector{Float64}},Mn::Array{Float64,2},Dn::Array{Float64,2},P::Int64,pℓ::Vector{Int64},
    𝒪::Vector{Int64},Nm::Vector{Int64},isFC::Bool,C::Vector{Vector{Float64}},
    ω::Vector{Array{Float64}},I_max::Int64,ϵ_max::Float64,
    sources::Array{Union{Array{Float64},Float64}},isAdapt::Vector{Bool},isCSD::Bool,
    solver::Int64,E::Float64,ΔE::Float64,𝚽E12::Array{Float64},S⁻::Vector{Float64},
    S⁺::Vector{Float64},α::Vector{Float64},ℳ::Array{Float64})

Solve the one-speed transport equation for a given particle.  

# Input Argument(s)
- '𝚽ℓ::Array{Float64}': Legendre components of the in-cell flux.
- 'Qℓout::Array{Float64}': Legendre components of the out-of-group in-cell source.
- 'Σt::Vector{Float64}': total cross-sections.
- 'Σs::Array{Float64}': Legendre moments of the scattering differential cross-sections.
- 'mat::Array{Int64,3}': material identifier per voxel.
- 'ndims::Int64': dimension of the geometry.
- 'N::Int64': number of discrete ordinates.
- 'ig::Int64': energy group index.
- 'Ns::Vector{Int64}': number of voxels per axis.
- 'Δs::Vector{Vector{Float64}}': size of each voxels per axis.
- 'Ω::Union{Vector{Vector{Float64}},Vector{Float64}}': director cosines.
- 'Mn::Array{Float64,2}': moment-to-discrete matrix.
- 'Dn::Array{Float64,2}': discrete-to-moment matrix.
- 'P::Int64': number of angular interpolation basis.
- 'pℓ::Vector{Int64}': legendre order associated with each interpolation basis. 
- '𝒪::Vector{Int64}': spatial and/or energy closure relation order.
- 'Nm::Vector{Int64}': number of spatial and/or energy moments.
- 'isFC::Bool': boolean to indicate if full coupling or not.
- 'C::Vector{Vector{Float64}}': constants related to the spatial and energy normalized
   Legendre expansion.
- 'ω::Vector{Array{Float64}}': weighting factors of the closure relations.
- 'I_max::Int64': maximum number of iterations of inner iterations.
- 'ϵ_max::Float64': convergence criterion on the flux solution.
- 'sources::Array{Union{Array{Float64},Float64}}': surface sources intensities.
- 'isAdapt::Vector{Bool}': boolean for adaptive calculations.
- 'isCSD::Bool': boolean to indicate if continuous slowing-down term is treated in 
   calculations.
- 'solver::Int64': indicate the type of solver to execute.
- 'E::Float64': group midpoint energy.
- 'ΔE::Float64': energy group width.
- '𝚽E12::Array{Float64}': incoming flux along the energy axis.
- 'S⁻::Vector{Float64}': restricted stopping power at higher energy group boundary.
- 'S⁺::Vector{Float64}': restricted stopping power at lower energy group boundary.
- 'α::Vector{Float64}': restricted momentum transfer.
- 'ℳ::Array{Float64}': Fokker-Planck scattering matrix.

# Output Argument(s)
- '𝚽ℓ::Array{Float64}': Legendre components of the in-cell flux.
- '𝚽E12::Array{Float64}': outgoing flux along the energy axis.
- 'ρ_in::Float64': estimated spectral radius.

# Reference(s)
- Larsen (2010) : Advances in Discrete-Ordinates Methodology.

"""
function compute_one_speed(𝚽ℓ::Array{Float64},Qℓout::Array{Float64},Σt::Vector{Float64},Σs::Array{Float64},mat::Array{Int64,3},ndims::Int64,N::Int64,ig::Int64,Ns::Vector{Int64},Δs::Vector{Vector{Float64}},Ω::Vector{Vector{Float64}},Mn::Array{Float64,2},Dn::Array{Float64,2},P::Int64,pℓ::Vector{Int64},𝒪::Vector{Int64},Nm::Vector{Int64},isFC::Bool,C::Vector{Float64},ω::Vector{Array{Float64}},I_max::Int64,ϵ_max::Float64,sources::Array{Union{Array{Float64},Float64}},isAdapt::Bool,isCSD::Bool,solver::Int64,E::Float64,ΔE::Float64,𝚽E12::Array{Float64},S⁻::Vector{Float64},S⁺::Vector{Float64},S::Array{Float64},α::Vector{Float64},ℳ::Array{Float64},Mn_FP::Array{Float64},Dn_FP::Array{Float64},N_FP::Int64,𝒜::String,is_CUDA::Bool,Ntot::Int64,is_EM::Bool,ℳ_EM::Array{Float64},𝒲::Array{Float64})

# Flux Initialization
𝚽E12_temp = Array{Float64}(undef)
if isCSD
    𝚽E12_temp = zeros(N,Nm[4],Ns[1],Ns[2],Ns[3])
end
N⁻ = 2
𝚽ℓ⁻ = zeros(N⁻,P,Nm[5],Ns[1],Ns[2],Ns[3])

# Source iteration loop until convergence
i_in = 1
ϵ_in = 0.0
ρ_in = NaN
isInnerConv=false
@inbounds while ~(isInnerConv)

    # Calculation of the Legendre components of the source (in-scattering)
    Qℓ = copy(Qℓout)
    if solver ∉ [4,5,6] Qℓ = scattering_source(Qℓ,𝚽ℓ,Σs,mat,P,pℓ,Nm[5],Ns) end

    # Finite element treatment of the angular Fokker-Planck term
    if solver ∈ [2,4] Qℓ = fokker_planck_source(N_FP,P,Nm[5],α,𝚽ℓ,Qℓ,Ns,mat,ℳ,Mn_FP,Dn_FP) end

    # Electromagnetic source
    if is_EM
        for ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3]), is in range(1,Nm[5])
            for n in range(1,P), m in range(1,P)
                Qℓ[n,is,ix,iy,iz] += ℳ_EM[n,m] * 𝚽ℓ[m,is,ix,iy,iz]
            end
        end
    end

    # If there is no source
    if ~any(x->x!=0,sources) && ~any(x->x!=0,Qℓ) && (~isCSD || (isCSD && ~any(x->x!=0,𝚽E12)))
        𝚽ℓ = zeros(P,Nm[5],Ns[1],Ns[2],Ns[3])
        ϵ_in = 0.0; i_in = 1
        println(">>>Group ",ig," has converged ( ϵ = ",@sprintf("%.4E",ϵ_in)," , N = ",i_in," , ρ = ",@sprintf("%.2f",ρ_in)," )")
        break
    end

    #----
    # Loop over all discrete ordinates
    #----
    #println(string(i_in," ",ϵ_in))
    𝚽ℓ .= 0
    @inbounds for n in range(1,N)
        if isCSD 𝚽E12ⁿ = 𝚽E12[n,:,:,:,:] else 𝚽E12ⁿ = Array{Float64}(undef) end
        if ndims == 1
            𝚽ℓ[:,:,:,1,1], 𝚽E12ⁿ = compute_sweep_1D(𝚽ℓ[:,:,:,1,1],Qℓ[:,:,:,1,1],Σt,mat[:,1,1],Ns[1],Δs[1],Ω[1][n],Mn[n,:],Dn[:,n],P,𝒪,Nm,isFC,C,ω,sources[n,:],isAdapt,isCSD,ΔE,𝚽E12ⁿ,S⁻,S⁺,S,𝒲)
        elseif ndims == 2
            if is_CUDA
                error()
            else
                𝚽ℓ[:,:,:,:,1],𝚽E12ⁿ = compute_sweep_2D(𝚽ℓ[:,:,:,:,1],Qℓ[:,:,:,:,1],Σt,mat[:,:,1],Ns[1:2],Δs[1:2],[Ω[1][n],Ω[2][n]],Mn[n,:],Dn[:,n],P,𝒪,Nm,C,ω,sources[n,:],isAdapt,isCSD,ΔE,𝚽E12ⁿ,S⁻,S⁺,S,𝒲)
            end
        elseif ndims == 3
            if is_CUDA
                error()
            else
                𝚽ℓ,𝚽E12ⁿ = compute_sweep_3D(𝚽ℓ,Qℓ,Σt,mat,Ns,Δs,[Ω[1][n],Ω[2][n],Ω[3][n]],Mn[n,:],Dn[:,n],P,𝒪,Nm,C,ω,sources[n,:],isAdapt,isCSD,ΔE,𝚽E12ⁿ,S⁻,S⁺,S,𝒲)
            end
        else
            error("Error in computeOneSpeed.jl: Dimension is not 1, 2 or 3.")
        end
        if isCSD 𝚽E12_temp[n,:,:,:,:] = 𝚽E12ⁿ end
    end
    
    #----
    # Verification of convergence of the one-group flux
    #----  
    ϵ_in = 0.0
    if (solver ∉ [5,6]) ϵ_in = maximum(vec(abs.((𝚽ℓ[1,1,:,:,:] .- 𝚽ℓ⁻[1,1,1,:,:,:])./max.(abs.(𝚽ℓ[1,1,:,:,:]),1e-16)))) end
    if (ϵ_in < ϵ_max) || i_in >= I_max

        # Convergence or maximum iterations reach
        isInnerConv = true
        Ntot += i_in
        if i_in ≥ 3 ρ_in = sqrt(sum(( vec(𝚽ℓ[1,1,:,:,:]) .- vec(𝚽ℓ⁻[1,1,1,:,:,:]) ).^2))/sqrt(sum(( vec(𝚽ℓ⁻[1,1,1,:,:,:]) .- vec(𝚽ℓ⁻[2,1,1,:,:,:]) ).^2)) end
        if ~(i_in >= I_max)
            println(">>>Group $ig has converged ( ϵ = ",@sprintf("%.4E",ϵ_in)," , N = ",i_in," , ρ = ",@sprintf("%.2f",ρ_in)," )")
        else
            println(">>>Group $ig has not converged ( ϵ = ",@sprintf("%.4E",ϵ_in)," , N = ",i_in," , ρ = ",@sprintf("%.2f",ρ_in)," )")
        end

    else

        # Livolant acceleration
        if 𝒜 == "livolant" && mod(i_in,3) == 0
            𝚽ℓ⁺ = livolant(𝚽ℓ,𝚽ℓ⁻[1,:,:,:,:,:],𝚽ℓ⁻[2,:,:,:,:,:])
            𝚽ℓ⁻[2,:,:,:,:,:] = 𝚽ℓ⁻[1,:,:,:,:,:]
            𝚽ℓ⁻[1,:,:,:,:,:] = 𝚽ℓ
            𝚽ℓ .= 𝚽ℓ⁺
        else
            𝚽ℓ⁻[2,:,:,:,:,:] = 𝚽ℓ⁻[1,:,:,:,:,:]
            𝚽ℓ⁻[1,:,:,:,:,:] = 𝚽ℓ
        end
        
        # Save flux solution and go to next iteration
        i_in += 1

    end

end

return 𝚽ℓ,𝚽E12_temp,ρ_in,Ntot

end



