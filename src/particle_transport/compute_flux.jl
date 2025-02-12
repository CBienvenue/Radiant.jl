"""
    compute_flux(cross_sections::Cross_Sections,geometry::Geometry,discrete_ordinates::Discrete_Ordinates,
    source::Source)

Solve the transport equation for a given particle.  

# Input Argument(s)
- 'cross_sections::Cross_Sections' : cross section informations.
- 'geometry::Geometry' : geometry informations.
- 'discrete_ordinates::Discrete_Ordinates' : discrete_ordinates informations.
- 'source::Source' : source informations.

# Output Argument(s)
- 'flux::Flux_Per_Particle': flux informations.

# Reference(s)
N/A

"""
function compute_flux(cross_sections::Cross_Sections,geometry::Geometry,discrete_ordinates::Discrete_Ordinates,source::Source)

#----
# Geometry data
#----
Ndims = geometry.get_dimension()
geo_type = geometry.get_type()
if geo_type != "cartesian" error("Transport of particles in",geo_type," is unavailable.") end
Ns = geometry.get_number_of_voxels()
Δs = geometry.get_voxels_width()
mat = geometry.get_material_per_voxel()

#----
# Preparation of angular discretisation
#----

L = discrete_ordinates.get_legendre_order()
N = discrete_ordinates.get_quadrature_order()
quadrature_type = discrete_ordinates.get_quadrature_type()
SN_type = discrete_ordinates.get_angular_boltzmann()
Qdims = discrete_ordinates.get_quadrature_dimension(Ndims)

# Compute quadrature weights and abscissae
Ω,w = quadrature(N,quadrature_type,Ndims,Qdims)
if typeof(Ω) == Vector{Float64} Ω = [Ω,0*Ω,0*Ω] end
P,Mn,Dn,pℓ,pm = angular_polynomial_basis(Ndims,Ω,w,L,N,SN_type,Qdims)
Nd = length(w)

#----
# Preparation of cross sections
#----

part = discrete_ordinates.get_particle()
solver,isCSD = discrete_ordinates.get_solver_type()
Nmat = cross_sections.get_number_of_materials()
Ng = cross_sections.get_number_of_groups(part)
if isCSD
    ΔE = cross_sections.get_energy_width(part)
    E = cross_sections.get_energies(part)
    Eb = cross_sections.get_energy_boundaries(part)
end

is_full_coupling = true
schemes,𝒪,Nm = discrete_ordinates.get_schemes(geometry,is_full_coupling)
ω,𝒞,is_adaptive,𝒲 = scheme_weights(𝒪,schemes,Ndims,isCSD)

println(">>>Particle: $(get_type(part)) <<<")

# Total cross sections
Σtot = zeros(Ng,Nmat)
if solver ∈ [4,5] 
    Σtot = cross_sections.get_absorption(part)
else
    Σtot = cross_sections.get_total(part)
end

# Scattering cross sections
Σs = zeros(Nmat,Ng,Ng,L+1)
if solver ∉ [4,5]
    Σs = cross_sections.get_scattering(part,part,L)
end

# Stopping powers
if isCSD
    S⁻ = zeros(Ng,Nmat); S⁺ = zeros(Ng,Nmat)
    Sb = cross_sections.get_stopping_powers(part)
    for n in range(1,Nmat)
        S⁻[:,n] = Sb[1:Ng,n] ; S⁺[:,n] = Sb[2:Ng+1,n]
    end
    S = zeros(Ng,Nmat,𝒪[4])
    for n in range(1,Nmat), ig in range(1,Ng)
        S[ig,n,1] = (S⁻[ig,n]+S⁺[ig,n])/2
        if (𝒪[4] > 1) S[ig,n,2] = (S⁻[ig,n]-S⁺[ig,n])/(2*sqrt(3)) end
    end
end

# Momentum transfer
if solver ∈ [2,4]
    α = zeros(Ng,Nmat)
    α = cross_sections.get_momentum_transfer(part)
    fokker_planck_type = discrete_ordinates.get_angular_fokker_planck()
    ℳ,λ₀ = fokker_planck_scattering_matrix(N,Nd,quadrature_type,Ndims,fokker_planck_type,Mn,Dn,pℓ,P,Qdims)
    Σtot .+= α .* λ₀/2
end

# Elastic-free approximation
if solver == 6
    for n in range(1,Nmat), ig in range(1,Ng)
        Σtot[ig,n] -= Σs[n,ig,ig,1]
    end
end

# External electro-magnetic fields
𝓔 = [0.0,0.0,0.0]; 𝓑 = [0.0,0.0,0.0]
is_EM = false
if is_EM
    q = part.get_charge()
    ℳ_EM = electromagnetic_scattering_matrix(𝓔,𝓑,q,Ω,w,Ndims,Mn,Dn,pℓ,pm,P,Ng,Eb,ΔE,Qdims)
else
    ℳ_EM = zeros(Ng,P,P);
end

#----
# Acceleration discrete_ordinates
#----

𝒜 = discrete_ordinates.get_acceleration()

#----
# Fixed sources
#----

surface_sources = source.get_surface_sources()
volume_sources = source.get_volume_sources()

#----
# Flux calculations
#----

ϵ_max = discrete_ordinates.get_convergence_criterion()
I_max = discrete_ordinates.get_maximum_iteration()

# Initialization flux
𝚽ℓ = zeros(Ng,P,Nm[5],Ns[1],Ns[2],Ns[3])
if isCSD 𝚽cutoff = zeros(P,Nm[5],Ns[1],Ns[2],Ns[3]) end

# All-group iteration
i_out = 1
is_outer_convergence = false
ϵ_out = Inf
is_outer_iteration = false
Ntot = 0
if is_outer_iteration 𝚽ℓ⁻ = zeros(Ng,Ns[1],Ns[2],Ns[3]) end

@inbounds while ~(is_outer_convergence)

    ρ_in = -ones(Ng) # In-group spectral radius
    isCSD ? 𝚽E12 = zeros(Nd,Nm[4],Ns[1],Ns[2],Ns[3]) : 𝚽E12 = Array{Float64}(undef)

    # Loop over energy group
    @inbounds for ig in range(1,Ng)

        # Calculation of the Legendre components of the source (out-scattering)
        Qℓout = zeros(P,Nm[5],Ns[1],Ns[2],Ns[3])
        if solver ∉ [4,5] Qℓout = scattering_source(Qℓout,𝚽ℓ,Σs[:,:,ig,:],mat,P,pℓ,Nm[5],Ns,Ng,ig) end

        # Fixed volumic sources
        Qℓout .+= volume_sources[ig,:,:,:,:,:]

        # Calculation of the group flux
        if isCSD
            if (ig != 1) 𝚽E12 = 𝚽E12 .* ΔE[ig]/ΔE[ig-1] end
            Eg = E[ig]
            ΔEg = ΔE[ig]
            Sg⁻ = S⁻[ig,:]/ΔEg
            Sg⁺ = S⁺[ig,:]/ΔEg
            Sg = S[ig,:,:]/ΔEg
            if solver ∈ [2,4]
                αg = α[ig,:]
            else
                αg = Vector{Float64}()
                ℳ = Array{Float64}(undef)
            end
        else
            Eg = 0.0
            ΔEg = 0.0
            Sg⁻ = Vector{Float64}()
            Sg⁺ = Vector{Float64}()
            Sg = Vector{Float64}()
            αg = Vector{Float64}()
            ℳ = Array{Float64}(undef)
        end
        𝚽ℓ[ig,:,:,:,:,:],𝚽E12,ρ_in[ig],Ntot = compute_one_speed(𝚽ℓ[ig,:,:,:,:,:],Qℓout,Σtot[ig,:],Σs[:,ig,ig,:],mat,Ndims,Nd,ig,Ns,Δs,Ω,Mn,Dn,P,pℓ,𝒪,Nm,is_full_coupling,𝒞,ω,I_max,ϵ_max,surface_sources[ig,:,:],is_adaptive,isCSD,solver,Eg,ΔEg,𝚽E12,Sg⁻,Sg⁺,Sg,αg,ℳ,𝒜,Ntot,is_EM,ℳ_EM[ig,:,:],𝒲)
        
    end

    # Verification of convergence in all energy groups
    if is_outer_iteration
        ϵ_out = maximum(vec(abs.(𝚽ℓ[:,1,1,:,:,:] .- 𝚽ℓ⁻)))/maximum(vec(abs.(𝚽ℓ[:,1,1,:,:,:])))
        𝚽ℓ⁻ = 𝚽ℓ[:,1,1,:,:,:]
    end
    if (ϵ_out < ϵ_max || i_out >= I_max) || ~is_outer_iteration
        is_outer_convergence = true
        # Calculate the flux at the cutoff energy
        if isCSD
            @inbounds for n in range(1,Nd), ix in range(1,Ns[1]), iy in range(1,Ns[2]), iz in range(1,Ns[3]), is in range(1,Nm[4]), p in range(1,P)
                𝚽cutoff[p,is,ix,iy,iz] += Dn[p,n] * 𝚽E12[n,is,ix,iy,iz]
            end
        end
    else
        i_out += 1
    end
    
end

# Save flux
flux = Flux_Per_Particle(part)
flux.add_flux(𝚽ℓ)
if isCSD flux.add_flux_cutoff(𝚽cutoff) end

return flux

end