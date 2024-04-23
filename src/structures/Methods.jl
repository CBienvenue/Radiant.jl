
mutable struct Methods

    # Variable(s)
    name                       ::Union{Missing,String}
    number_of_particles        ::Int64
    particles                  ::Vector{String}
    methods_names              ::Vector{String}
    methods_list               ::Vector{Method}
    number_of_generations      ::Int64
    
    # Constructor(s)
    function Methods()

        this = new()

        this.name = missing
        this.number_of_particles = 0
        this.particles = Vector{String}()
        this.methods_names = Vector{String}()
        this.methods_list = Vector{Method}()
        this.number_of_generations = 1

        return this
    end
end

# Method(s)
Base.propertynames(::Methods) = 
(
    fieldnames(Methods)...,
    :add_method,
    :set_number_of_generations,
    :get_method,
    :get_number_of_generations,
    :get_particles,
    :get_number_of_particles
)

function println(this::Methods)
    entries = ["Name","Number of particles","Particles","Methods"]
    values = [this.name,this.number_of_particles,this.particles,this.methods_names]
    N = length(entries); L = length.(entries); Lmax = maximum(L)
    println("Source")
    for n in range(1,N)
        println(string("   ",entries[n]," "^(Lmax-L[n])),"  :  ",values[n])
    end
end

function add_method(this::Methods,method::Method)
    this.number_of_particles += 1
    push!(this.particles,method.particle)
    #push!(this.methods_names,method.name)
    push!(this.methods_list,method)
end

function set_number_of_generations(this::Methods,number_of_generations::Int64)
    if number_of_generations < 1 error("Number of particle generations should be at least 1.") end
    this.number_of_generations = number_of_generations
end

function get_method(this::Methods,particle::String)
    index = findfirst(x -> x == particle,this.particles)
    return this.methods_list[index]
end

function get_number_of_generations(this::Methods)
    return this.number_of_generations
end

function get_particles(this::Methods)
    return this.particles
end

function get_number_of_particles(this::Methods)
    return this.number_of_particles
end