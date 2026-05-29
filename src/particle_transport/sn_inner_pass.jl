"""
    sn_inner_pass!(...)

Apply one in-group source-iteration pass for the discrete-ordinates solver: build the source
from the current flux `𝚽l`, sweep all `Nd` ordinates, and update `𝚽l` and the incoming
boundary fluxes (`𝚽x12_in`, plus `𝚽y12_in`/`𝚽z12_in` in 2D/3D). This is the building block 
every in-group acceleration scheme is built on (see `sn_one_speed`).

The pass is the affine map `T(z) = A·z + c` over the state `z = (𝚽l, boundary fluxes)`:
- `homogeneous = false` gives `T(z)` (all fixed sources active);
- `homogeneous = true` gives the linear part `A·z`, by dropping the fixed sources 
   (`Qlout = 0`,`Np_source = 0`, zeroed incoming energy flux `𝚽E12`).

`𝚽E12_temp` receives the outgoing energy flux (when `isCSD`); `Ql` and `𝚽*12_temp` are
scratch. Remaining arguments mirror `sn_one_speed`.
"""
function sn_inner_pass!(𝚽l,𝚽x12_in,𝚽y12_in,𝚽z12_in,𝚽x12_temp,𝚽y12_temp,𝚽z12_temp,𝚽E12_temp,Ql,Qlout,𝚽E12_in,Σt,Σs,mat,ndims,Nd,Ns,Δs,Ω,Mn,Dn,Np,pl,Mn_surf,Dn_surf,Np_surf,n_to_n⁺,𝒪,Nm,isFC,C,ω,sources,isAdapt,isCSD,solver,ΔE,S⁻,S⁺,S,T,ℳ,is_EM,ℳ_EM,𝒲,boundary_conditions,Np_source;homogeneous::Bool)

    # Calculation of the Legendre components of the source (in-scattering)
    if homogeneous Ql .= 0.0 else Ql .= Qlout end
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

    # Fixed-source switches: drop the external surface source and incoming energy flux for the
    # homogeneous operator (out-of-group source already dropped above via Ql .= 0).
    Np_source_eff = homogeneous ? 0 : Np_source
    𝚽E12_eff = (isCSD && homogeneous) ? zero(𝚽E12_in) : 𝚽E12_in

    #----
    # Loop over all discrete ordinates
    #----
    𝚽l .= 0
    for n in range(1,Nd)
        if isCSD 𝚽E12_out = 𝚽E12_eff[n,:,:,:,:] else 𝚽E12_out = Array{Float64}(undef) end
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
            𝚽l[:,:,:,1,1], 𝚽E12_out,𝚽x12_out = sn_sweep_1D(𝚽l[:,:,:,1,1],Ql[:,:,:,1,1],Σt,mat[:,1,1],Ns[1],Δs[1],Ω[1][n],Mn[n,:],Dn[:,n],Np,Mnx⁻,Dnx⁻,Np_surf,𝒪,Nm,C,ω,sources,isAdapt,isCSD,ΔE,𝚽E12_out,S⁻,S⁺,S,𝒲,isFC,𝚽x12_in,boundary_conditions,Np_source_eff)
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
            𝚽l[:,:,:,:,1],𝚽E12_out,𝚽x12_out,𝚽y12_out = sn_sweep_2D(𝚽l[:,:,:,:,1],Ql[:,:,:,:,1],Σt,mat[:,:,1],Ns[1:2],Δs[1:2],[Ω[1][n],Ω[2][n]],Mn[n,:],Dn[:,n],Np,Mnx⁻,Dnx⁻,Mny⁻,Dny⁻,Np_surf,𝒪,Nm,C,ω,sources,isAdapt,isCSD,ΔE,𝚽E12_out,S⁻,S⁺,S,𝒲,isFC,𝚽x12_in,𝚽y12_in,boundary_conditions,Np_source_eff)
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
            𝚽l,𝚽E12_out,𝚽x12_out,𝚽y12_out,𝚽z12_out = sn_sweep_3D(𝚽l,Ql,Σt,mat,Ns,Δs,[Ω[1][n],Ω[2][n],Ω[3][n]],Mn[n,:],Dn[:,n],Np,Mnx⁻,Dnx⁻,Mny⁻,Dny⁻,Mnz⁻,Dnz⁻,Np_surf,𝒪,Nm,C,ω,sources,isAdapt,isCSD,ΔE,𝚽E12_out,S⁻,S⁺,S,𝒲,isFC,𝚽x12_in,𝚽y12_in,𝚽z12_in,boundary_conditions,Np_source_eff)
        else
            error("Dimension is not 1, 2 or 3.")
        end
        𝚽x12_temp .+= 𝚽x12_out
        if ndims >= 2 𝚽y12_temp .+= 𝚽y12_out end
        if ndims >= 3 𝚽z12_temp .+= 𝚽z12_out end
        if isCSD 𝚽E12_temp[n,:,:,:,:] = 𝚽E12_out end
    end
    𝚽x12_in .= 𝚽x12_temp; 𝚽x12_temp .= 0.0
    if ndims >= 2 𝚽y12_in .= 𝚽y12_temp; 𝚽y12_temp .= 0.0 end
    if ndims >= 3 𝚽z12_in .= 𝚽z12_temp; 𝚽z12_temp .= 0.0 end

    return 𝚽l
end
