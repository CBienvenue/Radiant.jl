
mutable struct Fixed_Sources

    # Variable(s)
    name                       ::Union{Missing,String}
    number_of_particles        ::Int64
    particles                  ::Vector{String}
    normalization_factor       ::Float64
    sources_names              ::Vector{String}
    sources_list               ::Vector{Source}
    cross_sections             ::Cross_Sections
    geometry                   ::Geometry
    methods                    ::Methods

    # Constructor(s)
    function Fixed_Sources(cross_sections,geometry,methods)

        this = new()

        this.name = missing
        this.number_of_particles = 0
        this.normalization_factor = 0
        this.particles = Vector{String}()
        this.sources_names = Vector{String}()
        this.sources_list = Vector{Source}()
        this.cross_sections = cross_sections
        this.geometry = geometry
        this.methods = methods

        return this
    end
end

# Method(s)
Base.propertynames(::Fixed_Sources) = 
(
    fieldnames(Fixed_Sources)...,
    :add_source,
    :get_source,
    :get_particles,
    :get_normalization_factor
)

function add_source(this::Fixed_Sources,fixed_source::Union{Surface_Source,Volume_Source})
    particle = fixed_source.get_particle()
    method = this.methods.get_method(particle)
    if particle ∈ this.particles
        source = Source(particle,this.cross_sections,this.geometry,method)
        source.add_source(fixed_source)
        index = findfirst(x -> x == particle,this.particles)
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

function get_source(this::Fixed_Sources,particle::String)
    index = findfirst(x -> x == particle,this.particles)
    method = this.methods.get_method(particle)
    if isnothing(index)
        source = Source(particle,this.cross_sections,this.geometry,method)
        return source
    else
        return this.sources_list[index]
    end
end

function get_particles(this::Fixed_Sources)
    return this.particles
end

function get_normalization_factor(this::Fixed_Sources)
    return this.normalization_factor
end