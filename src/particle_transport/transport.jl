"""
    transport(cross_sections::Cross_Sections,geometry::Geometry,coupled_transport_solvers::Coupled_Transport_Solvers,
    sources::Fixed_Sources)

Solve transport calculations.

See also [`compute_flux`](@ref).

# Input Argument(s)
- 'cross_sections::Cross_Sections': cross section informations.
- 'geometry::Geometry': geometry informations.
- 'coupled_transport_solvers::Coupled_Transport_Solvers': coupled_transport_solvers informations.
- 'sources::Fixed_Sources': sources informations.

# Output Argument(s)
N/A

# Author(s)
Charles Bienvenue

# Reference(s)
N/A

"""
function transport(cross_sections::Cross_Sections,geometry::Geometry,coupled_transport_solvers::Coupled_Transport_Solvers,sources::Fixed_Sources,is_CUDA::Bool)

#----
# Initialization
#----
Npart = coupled_transport_solvers.get_number_of_particles()
flux = Flux()

#---
# One particles transport
#---
if Npart == 1

    # Initialization
    particle = coupled_transport_solvers.get_particles()[1]

    # Transport
    method = coupled_transport_solvers.get_method(particle)
    fixed_source = sources.get_source(particle)

    particle_flux = compute_flux(cross_sections,geometry,method,fixed_source,is_CUDA)
    flux.add_flux(particle_flux)

#----
# N-particles coupled transport
#----
else

    # Initialization
    particles = coupled_transport_solvers.get_particles()
    fixed_source = Vector{Source}(undef,Npart)
    particle_sources = Vector{Source}(undef,Npart)
    method = Vector{Method}(undef,Npart) 
    for i in range(1,Npart)
        method[i] = coupled_transport_solvers.get_method(particles[i])
        fixed_source[i] = sources.get_source(particles[i])
        particle_sources[i] = Source(particles[i],cross_sections,geometry,method[i])
    end

    # Ordering the particles
    particle_index = zero(Npart)
    for i in range(1,Npart)
        if particles[i] ∈ sources.get_particles()
            i0 = i
            if i0 == 1
                particle_index = collect(1:Npart)
            else
                particle_index = append!(collect(i0:Npart),collect(1:i0-1))
            end
            break
        end
        if i == Npart error("No fixed sources are defined.") end
    end

    # Coupled transport
    Ngen = coupled_transport_solvers.get_number_of_generations()
    for n in range(1,Ngen), p in range(1,Npart)

        i = particle_index[p]

        # Combinaison of fixed external and scattered particles sources 
        source = particle_sources[i]
        if n == 1 source += fixed_source[i] end

        # Transport
        particle_flux = compute_flux(cross_sections,geometry,method[i],source,is_CUDA)
        flux.add_flux(particle_flux)

        # Compute scattered particles sources
        for q in range(1,Npart)
            
            j = particle_index[q]
            if (n == Ngen) && (q ≤ p) continue end

            if i != j
                particle_sources[j] += particle_source(particle_flux,cross_sections,geometry,method[i],method[j])
            else
                particle_sources[i] = Source(particles[i],cross_sections,geometry,method[i])
            end
        end
    end
end

return flux
end