"""
    particle_source(flux::Flux_Per_Particle,cross_sections::Cross_Sections,
    geometry::Geometry,solver_in::Solver,
    solver_out::Solver)

Compute the source of particle produced by interaction of another type of particle with
matter.

# Input Argument(s)
- `flux::Flux_Per_Particle`: flux informations of the incoming particle.
- `cross_sections::Cross_Sections`: cross section informations.
- `geometry::Geometry`: geometry informations.
- `solver_in::Solver`: method informations of the incoming particle.
- `solver_out::Solver`: method informations of the outgoing particle.

# Output Argument(s)
- `ps::Source`: source information of the outgoing particle.

# Reference(s)
N/A

"""
function particle_source(flux::Flux_Per_Particle,cross_sections::Cross_Sections,geometry::Geometry,solver_in::Solver,solver_out::Solver)

# Geometry data
Ndims = geometry.get_dimension()
Ns = geometry.get_number_of_voxels()
mat = geometry.get_material_per_voxel()

# Discrete ordinates data
_,isCSD = solver_in.get_solver_type()
particle_in = solver_in.get_particle()
particle_out = solver_out.get_particle()
L_in = solver_in.get_legendre_order()
L_out = solver_out.get_legendre_order()
isFC_in = solver_in.get_is_full_coupling()
isFC_out = solver_out.get_is_full_coupling()
_,𝒪_in,Nm_in = solver_in.get_schemes(geometry,isFC_in)
_,𝒪_out,Nm_out = solver_out.get_schemes(geometry,isFC_out)
Nm_in = Nm_in[5]
Nm_out = Nm_out[5]
if solver_in isa Discrete_Ordinates
    Qdims_in = solver_in.get_quadrature_dimension(Ndims)
    Ω_in,w_in = quadrature(solver_in.get_quadrature_order(),solver_in.get_quadrature_type(),Ndims,Qdims_in)
    if typeof(Ω_in) == Vector{Float64} Ω_in = [Ω_in,0*Ω_in,0*Ω_in] end
    P_in,_,_,pl_in,pm_in = angular_polynomial_basis(Ω_in,w_in,L_in,solver_in.get_angular_boltzmann(),Qdims_in)
elseif solver_in isa Spherical_Harmonics
    P_in = L_in+1
    if Ndims == 1
        pl_in = collect(0:L_in)
        pm_in = zeros(Int64,L_in+1)
    else
        error("Unsupported spatial dimension.")
    end
else
    error("Unknown angular discretization method for the incoming particle.")
end

if solver_out isa Discrete_Ordinates
    Qdims_out = solver_out.get_quadrature_dimension(Ndims)
    Ω_out,w_out = quadrature(solver_out.get_quadrature_order(),solver_out.get_quadrature_type(),Ndims,Qdims_out)
    if typeof(Ω_out) == Vector{Float64} Ω_out = [Ω_out,0*Ω_out,0*Ω_out] end
    P_out,_,Dn_out,_,_ = angular_polynomial_basis(Ω_out,w_out,L_out,solver_out.get_angular_boltzmann(),Qdims_out)
elseif solver_out isa Spherical_Harmonics
    P_out = L_out+1
else
    error("Unknown angular discretization method for the outgoing particle.")
end


#_,isCSD = solver_in.get_solver_type()
#particle_in = solver_in.get_particle()
#particle_out = solver_out.get_particle()
#L_in = solver_in.get_legendre_order()
#L_out = solver_out.get_legendre_order()
#isFC_in = solver_in.get_is_full_coupling()
#isFC_out = solver_out.get_is_full_coupling()
#_,𝒪_in,Nm_in = solver_in.get_schemes(geometry,isFC_in)
#_,𝒪_out,Nm_out = solver_out.get_schemes(geometry,isFC_out)
#Nm_in = Nm_in[5]; Nm_out = Nm_out[5]
#Qdims_in = solver_in.get_quadrature_dimension(Ndims)
#Qdims_out = solver_out.get_quadrature_dimension(Ndims)
#Ω_in,w_in = quadrature(solver_in.get_quadrature_order(),solver_in.get_quadrature_type(),Ndims,Qdims_in)
#Ω_out,w_out = quadrature(solver_out.get_quadrature_order(),solver_out.get_quadrature_type(),Ndims,Qdims_out)
# if typeof(Ω_in) == Vector{Float64} Ω_in = [Ω_in,0*Ω_in,0*Ω_in] end
# if typeof(Ω_out) == Vector{Float64} Ω_out = [Ω_out,0*Ω_out,0*Ω_out] end

# Compute transfer matrix
#P_in,_,_,pl_in,pm_in = angular_polynomial_basis(Ω_in,w_in,L_in,solver_in.get_angular_boltzmann(),Qdims_in)
#P_out,_,Dn_out,_,_ = angular_polynomial_basis(Ω_out,w_out,L_out,solver_out.get_angular_boltzmann(),Qdims_out)

if solver_in isa Discrete_Ordinates && solver_out isa Discrete_Ordinates

    if solver_in.get_angular_boltzmann() == solver_out.get_angular_boltzmann() && length(w_out) == length(w_in) && w_out == w_in
        type_scat = solver_in.get_angular_boltzmann()
    else
        type_scat = "standard"
    end

    if (Qdims_in == 1 && Qdims_out ∈ [2,3]) || (Qdims_in ∈ [2,3] && Qdims_out == 1)
        Nd = length(w_out)
        μ = Ω_out[1]
        Pl = zeros(Nd,maximum(pl_in)+1,1)
        for n in range(1,Nd)
            Pl[n,:] = legendre_polynomials_up_to_L(maximum(pl_in),μ[n])
        end
        P = length(pl_in)
        Mn_tr = zeros(Nd,P)
        for p in range(1,P)
            for n in range(1,Nd)
                if pm_in[p] == 0 || Qdims_in == 1
                    Mn_tr[n,pl_in[p]+1] = (2*pl_in[p]+1)/2 * Pl[n,pl_in[p]+1]
                end
            end
        end
    elseif Qdims_in == Qdims_out
        _,Mn_tr,_,_,_ = angular_polynomial_basis(Ω_out,w_out,L_in,type_scat,Qdims_out)
    else
        error("Unknown particle transfer.")
    end
    T = Dn_out*Mn_tr

elseif solver_in isa Spherical_Harmonics && solver_out isa Spherical_Harmonics
    T = ones(P_out,P_in)
else
    error("Unknown angular discretization method for particle transfer.")
end

# Cross-sections data
Nmat = cross_sections.get_number_of_materials()
Ng_in = cross_sections.get_number_of_groups(particle_in)
Ng_out = cross_sections.get_number_of_groups(particle_out)
Σs = zeros(Nmat,Ng_in+1,Ng_out,L_in+1)
Σs = cross_sections.get_scattering(particle_in,particle_out,L_in)

# Flux data
𝚽l = flux.get_flux()
if isCSD 𝚽cutoff = flux.get_flux_cutoff() else 𝚽cutoff = zeros(P_in,Nm_in,Ns[1],Ns[2],Ns[3]) end

# Compute the scattered particle source
Ql_in = zeros(Ng_out,P_in,Nm_in,Ns[1],Ns[2],Ns[3])
Ql_out = zeros(Ng_out,P_out,Nm_out,Ns[1],Ns[2],Ns[3])
particle_sources(Ql_in,𝚽l,Σs,mat,P_in,pl_in,Nm_in,Ns,Ng_in,Ng_out)

# Adapt the source to the new particle flux expansions
map = map_moments(𝒪_in,𝒪_out,isFC_in,isFC_out)
for i in range(1,length(map))
    m = map[i]
    if m != 0
        for p_in in range(1,P_in), p_out in range(1,P_out)
            Ql_out[:,p_out,i,:,:,:] .+= T[p_out,p_in] .* Ql_in[:,p_in,m,:,:,:]
        end
    end
end

# Save source informations
ps = Source(particle_out,cross_sections,geometry,solver_out)
ps.add_volume_source(Ql_out)
return ps

end