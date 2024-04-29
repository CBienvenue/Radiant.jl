"""
    Methods

Structure used to define the collection of discretization methods for transport calculations associated with each of the particle and additionnal coupled transport informations.

# User-defined field(s)

- ## Mandatory field(s)
    - `methods_list::Vector{Method}`: list of the particle methods
    - `number_of_generations::Int64`: number of particle generations to transport.

- ## Optional field(s) - with default values
    N/A

# System-defined field(s)
- `number_of_particles::Int64`: number of particle for which methods are described
- `particles::Vector{String}`: vector of the particle for which methods are described
- `methods_names::Vector{String}`: vector of the names of the particle methods

"""
mutable struct Methods

    # Variable(s)
    number_of_particles        ::Int64
    particles                  ::Vector{String}
    methods_names              ::Vector{String}
    methods_list               ::Vector{Method}
    number_of_generations      ::Int64
    
    # Constructor(s)
    function Methods()

        this = new()

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

"""
    add_method(this::Methods,method::Method)

To add a particle and is associated methods to the Methods structure.

# Input Argument(s)
- `this::Methods`: collection of discretization method.
- `method::Method`: discretization method.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = Method()
julia> ... # Define the methods properties
julia> ms = Methods()
julia> ms.add_method(m)
```
"""
function add_method(this::Methods,method::Method)
    this.number_of_particles += 1
    push!(this.particles,method.particle)
    push!(this.methods_list,method)
end

"""
    set_number_of_generations(this::Methods,number_of_generations::Int64)

To set the number of particle generation to transport during calculations

# Input Argument(s)
- `this::Methods`: collection of discretization method.
- `number_of_generations::Int64`: number of particle generation to transport during calculations.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ms = Methods()
julia> ms.set_number_of_generations(2)
```
"""
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