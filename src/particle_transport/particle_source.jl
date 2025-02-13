"""
    particle_source(flux::Flux_Per_Particle,cross_sections::Cross_Sections,
    geometry::Geometry,discrete_ordinates_in::Discrete_Ordinates,
    discrete_ordinates_out::Discrete_Ordinates)

Compute the source of particle produced by interaction of another type of particle with
matter.

# Input Argument(s)
- `flux::Flux_Per_Particle`: flux informations of the incoming particle.
- `cross_sections::Cross_Sections`: cross section informations.
- `geometry::Geometry`: geometry informations.
- `discrete_ordinates_in::Discrete_Ordinates`: method informations of the incoming particle.
- `discrete_ordinates_out::Discrete_Ordinates`: method informations of the outgoing particle.

# Output Argument(s)
- `ps::Source`: source information of the outgoing particle.

# Reference(s)
N/A

"""
function particle_source(flux::Flux_Per_Particle,cross_sections::Cross_Sections,geometry::Geometry,discrete_ordinates_in::Discrete_Ordinates,discrete_ordinates_out::Discrete_Ordinates)

# Geometry data
Ndims = geometry.get_dimension()
Ns = geometry.get_number_of_voxels()
mat = geometry.get_material_per_voxel()

# Discrete ordinates data
_,isCSD = discrete_ordinates_in.get_solver_type()
particle_in = discrete_ordinates_in.get_particle()
particle_out = discrete_ordinates_out.get_particle()
L_in = discrete_ordinates_in.get_legendre_order()
L_out = discrete_ordinates_out.get_legendre_order()
_,𝒪_in,Nm_in = discrete_ordinates_in.get_schemes(geometry,true)
_,𝒪_out,Nm_out = discrete_ordinates_out.get_schemes(geometry,true)
Nm_in = Nm_in[5]; Nm_out = Nm_out[5]
Qdims_in = discrete_ordinates_in.get_quadrature_dimension(Ndims)
Qdims_out = discrete_ordinates_out.get_quadrature_dimension(Ndims)
Ω_in,w_in = quadrature(discrete_ordinates_in.get_quadrature_order(),discrete_ordinates_in.get_quadrature_type(),Ndims,Qdims_in)
Ω_out,w_out = quadrature(discrete_ordinates_out.get_quadrature_order(),discrete_ordinates_out.get_quadrature_type(),Ndims,Qdims_out)
if typeof(Ω_in) == Vector{Float64} Ω_in = [Ω_in,0*Ω_in,0*Ω_in] end
if typeof(Ω_out) == Vector{Float64} Ω_out = [Ω_out,0*Ω_out,0*Ω_out] end

# Compute transfer matrix
P_in,_,_,pℓ_in,pm_in = angular_polynomial_basis(Ndims,Ω_in,w_in,L_in,discrete_ordinates_in.get_quadrature_order(),discrete_ordinates_in.get_angular_boltzmann(),Qdims_in)
P_out,_,Dn_out,_,_ = angular_polynomial_basis(Ndims,Ω_out,w_out,L_out,discrete_ordinates_out.get_quadrature_order(),discrete_ordinates_out.get_angular_boltzmann(),Qdims_out)
if discrete_ordinates_in.get_angular_boltzmann() == discrete_ordinates_out.get_angular_boltzmann() && length(w_out) == length(w_in) && w_out == w_in
    type_scat = discrete_ordinates_in.get_angular_boltzmann()
else
    type_scat = "standard"
end
if (Qdims_in == 1 && Qdims_out ∈ [2,3]) || (Qdims_in ∈ [2,3] && Qdims_out == 1)
    Nd = length(w_out)
    μ = Ω_out[1]
    Pℓ = zeros(Nd,maximum(pℓ_in)+1,1)
    @inbounds for n in range(1,Nd)
        Pℓ[n,:] = legendre_polynomials(maximum(pℓ_in),μ[n])
    end
    P = length(pℓ_in)
    Mn_tr = zeros(Nd,P)
    for p in range(1,P)
        for n in range(1,Nd)
            if pm_in[p] == 0 || Qdims_in == 1
                Mn_tr[n,pℓ_in[p]+1] = (2*pℓ_in[p]+1)/2 * Pℓ[n,pℓ_in[p]+1]
            end
        end
    end
elseif Qdims_in == Qdims_out
    _,Mn_tr,_,_,_ = angular_polynomial_basis(Ndims,Ω_out,w_out,L_in,discrete_ordinates_out.get_quadrature_order(),type_scat,Qdims_out)
else
    error("Unknown particle transfer.")
end

# Cross-sections data
Nmat = cross_sections.get_number_of_materials()
Ng_in = cross_sections.get_number_of_groups(particle_in)
Ng_out = cross_sections.get_number_of_groups(particle_out)
Σs = zeros(Nmat,Ng_in+1,Ng_out,L_in+1)
Σs = cross_sections.get_scattering(particle_in,particle_out,L_in)

# Flux data
𝚽ℓ = flux.get_flux()
if isCSD 𝚽cutoff = flux.get_flux_cutoff() else 𝚽cutoff = zeros(P_in,Nm_in,Ns[1],Ns[2],Ns[3]) end

# Compute the scattered particle source
Qℓ_in = zeros(Ng_in,P_in,Nm_in,Ns[1],Ns[2],Ns[3])
Qℓ_out = zeros(Ng_out,P_out,Nm_out,Ns[1],Ns[2],Ns[3])
particle_source(Qℓ_in,𝚽ℓ,Σs,mat,P_in,pℓ_in,Nm_in,Ns,Ng_in,Ng_out)

# Adapt the source to the new particle flux expansions
map = map_moments(𝒪_in,𝒪_out)
T = Dn_out*Mn_tr
for i in range(1,length(map))
    m = map[i]
    if m != 0
        for p_in in range(1,P_in), p_out in range(1,P_out)
            Qℓ_out[:,p_out,i,:,:,:] .+= T[p_out,p_in] .* Qℓ_in[:,p_in,m,:,:,:]
        end
    end
end

# Save source informations
ps = Source(particle_out,cross_sections,geometry,discrete_ordinates_out)
ps.add_volume_source(Qℓ_out)
return ps

end