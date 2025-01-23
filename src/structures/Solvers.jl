"""
    Solvers

Structure used to define the collection of discretization methods for transport calculations associated with each of the particle and additionnal coupled transport informations.

# Mandatory field(s)
- `methods_list::Vector{Discrete_Ordinates}` : list of the particle methods
- `number_of_generations::Int64` : number of particle generations to transport.

# Optional field(s) - with default values
- N/A

"""
mutable struct Solvers

    # Variable(s)
    number_of_particles        ::Int64
    particles                  ::Vector{Particle}
    methods_names              ::Vector{String}
    methods_list               ::Vector{Discrete_Ordinates}
    number_of_generations      ::Int64
    
    # Constructor(s)
    function Solvers()

        this = new()

        this.number_of_particles = 0
        this.particles = Vector{Particle}()
        this.methods_names = Vector{String}()
        this.methods_list = Vector{Discrete_Ordinates}()
        this.number_of_generations = 1

        return this
    end
end

# Method(s)
"""
    add_solver(this::Solvers,method::Discrete_Ordinates)

To add a particle and is associated methods to the Solvers structure.

# Input Argument(s)
- `this::Solvers` : collection of discretization method.
- `method::Discrete_Ordinates` : discretization method.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = Discrete_Ordinates()
julia> ... # Define the methods properties
julia> ms = Solvers()
julia> ms.add_solver(m)
```
"""
function add_solver(this::Solvers,method::Discrete_Ordinates)
    this.number_of_particles += 1
    push!(this.particles,method.particle)
    push!(this.methods_list,method)
end

"""
    set_number_of_generations(this::Solvers,number_of_generations::Int64)

To set the number of particle generation to transport during calculations

# Input Argument(s)
- `this::Solvers` : collection of discretization method.
- `number_of_generations::Int64` : number of particle generation to transport during calculations.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ms = Solvers()
julia> ms.set_number_of_generations(2)
```
"""
function set_number_of_generations(this::Solvers,number_of_generations::Int64)
    if number_of_generations < 1 error("Number of particle generations should be at least 1.") end
    this.number_of_generations = number_of_generations
end

function get_method(this::Solvers,particle::Particle)
    index = findfirst(x -> get_id(x) == get_id(particle),this.particles)
    if isnothing(index) error("Solvers don't contain data for the given particle.") end
    return this.methods_list[index]
end

function get_number_of_generations(this::Solvers)
    return this.number_of_generations
end

function get_particles(this::Solvers)
    return this.particles
end

function get_number_of_particles(this::Solvers)
    return this.number_of_particles
end