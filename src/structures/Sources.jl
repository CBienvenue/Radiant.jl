
mutable struct Sources

    # Variable(s)
    name                       ::Union{Missing,String}
    number_of_particles        ::Int64
    particles                  ::Vector{String}
    sources_names              ::Vector{String}
    sources_list               ::Vector{Source}

    # Function(s)
    add_source                 ::Function
    get_source                 ::Function
    get_particles              ::Function

    # Constructor(s)
    function Sources()

        this = new()

        this.name = missing
        this.number_of_particles = 0
        this.particles = Vector{String}()
        this.sources_names = Vector{String}()
        this.sources_list = Vector{Source}()

        this.add_source = function (source) add_source!(this,source) end
        this.get_source = function (particle,cross_sections,geometry,discrete_ordinates) get_source!(this,particle,cross_sections,geometry,discrete_ordinates) end
        this.get_particles = function () get_particles(this) end

        return this
    end
end

# Method(s)
function println(this::Sources)
    entries = ["Name","Number of particles","Particles","Sources"]
    values = [this.name,this.number_of_particles,this.particles,this.sources_names]
    N = length(entries); L = length.(entries); Lmax = maximum(L)
    println("Sources")
    for n in range(1,N)
        println(string("   ",entries[n]," "^(Lmax-L[n])),"  :  ",values[n])
    end
end

function add_source!(this::Sources,source::Source)
    this.number_of_particles += 1
    push!(this.particles,source.particle)
    #push!(this.sources_names,source.name)
    push!(this.sources_list,source)
end

function get_source!(this::Sources,particle::String,cross_sections::Cross_Sections,geometry::Geometry,discrete_ordinates::Discrete_Ordinates)
    index = findfirst(x -> x == particle,this.particles)
    if isnothing(index)
        source = Source(particle,cross_sections,geometry,discrete_ordinates)
        return source
    else
        return this.sources_list[index]
    end
end

function get_particles(this::Sources)
    return this.particles
end