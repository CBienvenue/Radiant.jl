
mutable struct Flux

    # Variable(s)
    number_of_particles             ::Int64
    particles                       ::Vector{String}
    flux_per_particle               ::Vector{Flux_Per_Particle}

    # Constructor(s)
    function Flux()

        this = new()
        this.number_of_particles = 0
        this.particles = Vector{String}()
        this.flux_per_particle = Vector{Flux_Per_Particle}()

        return this
    end
end

# Method(s)
Base.propertynames(::Flux) = 
(
    fieldnames(Flux)...,
    :add_flux,
    :get_particles,
    :get_flux,
    :get_flux_cutoff
)

function add_flux(this::Flux,flux_per_particle::Flux_Per_Particle)
    particle = flux_per_particle.particle
    if particle ∈ this.particles 
        index = findfirst(x -> x == particle,this.particles)
        this.flux_per_particle[index].add_flux(flux_per_particle)
    else
        this.number_of_particles += 1
        push!(this.particles,particle)
        push!(this.flux_per_particle,flux_per_particle)
    end
end

function get_particles(this::Flux)
    return this.particles
end

function get_flux(this::Flux,particle::String)
    if particle ∉ this.particles error("No data for the specified particle.") end
    index = findfirst(x -> x == particle,this.particles)
    return this.flux_per_particle[index].get_flux()
end

function get_flux_cutoff(this::Flux,particle::String)
    if particle ∉ this.particles error("No data for the specified particle.") end
    index = findfirst(x -> x == particle,this.particles)
    return this.flux_per_particle[index].get_flux_cutoff()
end