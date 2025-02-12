"""
    Source

Structure used to describe sources for every particles.

"""
mutable struct Sources

    # Variable(s)
    name                       ::Union{Missing,String}
    number_of_particles        ::Int64
    particles                  ::Vector{Particle}
    sources_names              ::Vector{String}
    sources_list               ::Vector{Source}

    # Constructor(s)
    function Sources()

        this = new()

        this.name = missing
        this.number_of_particles = 0
        this.particles = Vector{Particle}()
        this.sources_names = Vector{String}()
        this.sources_list = Vector{Source}()

        return this
    end
end

# Method(s)
"""
    add_source(this::Sources,source::Source)

Add sources associated to a particle to the Sources structure.

# Input Argument(s)
- `this::Sources` : sources structure.
- `source::Source` : sources for a given particle.

# Output Argument(s)
N/A

"""
function add_source(this::Sources,source::Source)
    this.number_of_particles += 1
    push!(this.particles,source.particle)
    push!(this.sources_list,source)
end

"""
    get_source(this::Sources,particle::String,cross_sections::Cross_Sections,
    geometry::Geometry,discrete_ordinates::Discrete_Ordinates)

Get sources for a given particle.

# Input Argument(s)
- `this::Sources` : sources structure.
- `particle::String` : particle.
- `cross_sections::Cross_Sections` : cross-sections library.
- `geometry::Geometry` : geometry.
- `discrete_ordinates::Discrete_Ordinates` : discrete ordinates solver.

# Output Argument(s)
- `source::Source` : sources for the given particle.

"""
function get_source(this::Sources,particle::String,cross_sections::Cross_Sections,geometry::Geometry,discrete_ordinates::Discrete_Ordinates)
    index = findfirst(x -> get_id(x) == get_id(particle),this.particles)
    if isnothing(index)
        source = Source(particle,cross_sections,geometry,discrete_ordinates)
        return source
    else
        return this.sources_list[index]
    end
end

"""
    get_particles(this::Sources)

Get the list of particles.

# Input Argument(s)
- `this::Sources` : sources structure.

# Output Argument(s)
- `particles::Vector{Particles}` : particle list.

"""
function get_particles(this::Sources)
    return this.particles
end