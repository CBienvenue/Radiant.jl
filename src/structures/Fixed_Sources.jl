"""
    Fixed_Sources

Structure used to define a collection of fixed sources and their properties.

# Mandatory field(s)
- `sources_list::Vector{Source}` : list of Volume_Source or Surface_Source sources.
- `cross_sections` : Cross_Sections structure.
- `geometry` : Geometry structure.
- `solvers` : Solvers structure.

# Optional field(s) - with default values
- N/A

"""
mutable struct Fixed_Sources

    # Variable(s)
    number_of_particles        ::Int64
    particles                  ::Vector{Particle}
    normalization_factor       ::Float64
    sources_names              ::Vector{String}
    sources_list               ::Vector{Source}
    cross_sections             ::Cross_Sections
    geometry                   ::Geometry
    solvers                    ::Solvers

    # Constructor(s)
    function Fixed_Sources(cross_sections,geometry,solvers)

        this = new()

        this.number_of_particles = 0
        this.normalization_factor = 0
        this.particles = Vector{Particle}()
        this.sources_names = Vector{String}()
        this.sources_list = Vector{Source}()
        this.cross_sections = cross_sections
        this.geometry = geometry
        this.solvers = solvers

        return this
    end
end

# Method(s)
"""
    add_source(this::Fixed_Sources,fixed_source::Union{Surface_Source,Volume_Source})

To add a surface or volume source to the Fixed_Sources structure.

# Input Argument(s)
- `this::Fixed_Sources` : collection of fixed sources.
- `fixed_source::Union{Surface_Source,Volume_Source}` : volume or surface source.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> vs = Volume_Source()
julia> ... # Define the volume source properties
julia> fs = Fixed_Sources()
julia> fs.add_source(vs)
```
"""
function add_source(this::Fixed_Sources,fixed_source::Union{Surface_Source,Volume_Source})
    particle = fixed_source.get_particle()
    method = this.solvers.get_method(particle)
    if get_id(particle) ∈ get_id.(this.particles)
        source = Source(particle,this.cross_sections,this.geometry,method)
        source.add_source(fixed_source)
        index = findfirst(x -> get_id(x) == get_id(particle),this.particles)
        this.sources_list[index] += source
    else
        this.number_of_particles += 1
        push!(this.particles,particle)
        source = Source(particle,this.cross_sections,this.geometry,method)
        source.add_source(fixed_source)
        push!(this.sources_list,source)
    end
    this.normalization_factor += fixed_source.get_normalization_factor()
end

"""
    get_source(this::Fixed_Sources,particle::Particle)

Get the sources for a given particle.

# Input Argument(s)
- `this::Fixed_Sources` : collection of fixed sources.
- `particle::Particle` : particle.

# Output Argument(s)
- `source::Source` : sources for the given particle.

"""
function get_source(this::Fixed_Sources,particle::Particle)
    index = findfirst(x -> get_id(x) == get_id(particle),this.particles)
    method = this.solvers.get_method(particle)
    if isnothing(index)
        source = Source(particle,this.cross_sections,this.geometry,method)
        return source
    else
        return this.sources_list[index]
    end
end

"""
    get_particles(this::Fixed_Sources)

Get the particle list.

# Input Argument(s)
- `this::Fixed_Sources` : collection of fixed sources.

# Output Argument(s)
- `particles::Vector{Particle}` : particle list.

"""
function get_particles(this::Fixed_Sources)
    return this.particles
end

"""
    get_normalization_factor(this::Fixed_Sources)

Get the normalization factor.

# Input Argument(s)
- `this::Fixed_Sources` : collection of fixed sources.

# Output Argument(s)
- `normalization_factor::Float64` : normalization factor.

"""
function get_normalization_factor(this::Fixed_Sources)
    return this.normalization_factor
end