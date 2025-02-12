"""
    Flux

Structure used to contain flux information for all particles.

"""
mutable struct Flux

    # Variable(s)
    number_of_particles             ::Int64
    particles                       ::Vector{Particle}
    flux_per_particle               ::Vector{Flux_Per_Particle}

    # Constructor(s)
    function Flux()

        this = new()
        this.number_of_particles = 0
        this.particles = Vector{Particle}()
        this.flux_per_particle = Vector{Flux_Per_Particle}()

        return this
    end
end

# Method(s)
"""
    add_flux(this::Flux_Per_Particle,flux::Array{Float64,6})

Add flux solution for a particle.

# Input Argument(s)
- `this::Flux` : structure to contain flux solutions.
- `flux_per_particle::Flux_Per_Particle` : flux solution for a given particle.

# Output Argument(s)
N/A

"""
function add_flux(this::Flux,flux_per_particle::Flux_Per_Particle)
    particle = flux_per_particle.particle
    if get_id(particle) ∈ get_id.(this.particles) 
        index = findfirst(x -> get_id(x) == get_id(particle),this.particles)
        this.flux_per_particle[index].add_flux(flux_per_particle)
    else
        this.number_of_particles += 1
        push!(this.particles,particle)
        push!(this.flux_per_particle,flux_per_particle)
    end
end

"""
    get_particles(this::Flux)

Get particle list.

# Input Argument(s)
- `this::Flux` : structure to contain flux solutions.

# Output Argument(s)
- `particles::Vector{Particle}` : particle list.

"""
function get_particles(this::Flux)
    return this.particles
end

"""
    get_flux(this::Flux,particle::Particle)

Get the total flux solution for a given particle.

# Input Argument(s)
- `this::Flux` : structure to contain flux solutions.
- `particle::Particle` : particle

# Output Argument(s)
- `flux_per_particle::Array{Float64}` : flux solution.

"""
function get_flux(this::Flux,particle::Particle)
    if get_id(particle) ∉ get_id.(this.particles) error("No data for the specified particle.") end
    index = findfirst(x -> get_id(x) == get_id(particle),this.particles)
    return this.flux_per_particle[index].get_flux()
end

"""
    get_flux_cutoff(this::Flux,particle::Particle)

Get the total flux solution at cutoff for a given particle.

# Input Argument(s)
- `this::Flux` : structure to contain flux solutions.
- `particle::Particle` : particle

# Output Argument(s)
- `flux_per_particle::Array{Float64}` : flux solution at cutoff.

"""
function get_flux_cutoff(this::Flux,particle::Particle)
    if get_id(particle) ∉ get_id.(this.particles) error("No data for the specified particle.") end
    index = findfirst(x -> get_id(x) == get_id(particle),this.particles)
    return this.flux_per_particle[index].get_flux_cutoff()
end