"""
    Coupled_Transport_Solvers

Structure used to define the collection of discretization methods for transport calculations associated with each of the particle and additionnal coupled transport informations.

# User-defined field(s)

- **Mandatory field(s)**
    - `methods_list::Vector{Discrete_Ordinates}`: list of the particle methods
    - `number_of_generations::Int64`: number of particle generations to transport.

- **Optional field(s) - with default values**
    N/A

# System-defined field(s)
- `number_of_particles::Int64`: number of particle for which methods are described
- `particles::Vector{String}`: vector of the particle for which methods are described
- `methods_names::Vector{String}`: vector of the names of the particle methods

"""
mutable struct Coupled_Transport_Solvers

    # Variable(s)
    number_of_particles        ::Int64
    particles                  ::Vector{String}
    methods_names              ::Vector{String}
    methods_list               ::Vector{Discrete_Ordinates}
    number_of_generations      ::Int64
    
    # Constructor(s)
    function Coupled_Transport_Solvers()

        this = new()

        this.number_of_particles = 0
        this.particles = Vector{String}()
        this.methods_names = Vector{String}()
        this.methods_list = Vector{Discrete_Ordinates}()
        this.number_of_generations = 1

        return this
    end
end

# Discrete_Ordinates(s)
Base.propertynames(::Coupled_Transport_Solvers) = 
(
    fieldnames(Coupled_Transport_Solvers)...,
    :add_method,
    :set_number_of_generations,
    :get_method,
    :get_number_of_generations,
    :get_particles,
    :get_number_of_particles
)

"""
    add_solver(this::Coupled_Transport_Solvers,method::Discrete_Ordinates)

To add a particle and is associated methods to the Coupled_Transport_Solvers structure.

# Input Argument(s)
- `this::Coupled_Transport_Solvers`: collection of discretization method.
- `method::Discrete_Ordinates`: discretization method.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = Discrete_Ordinates()
julia> ... # Define the methods properties
julia> ms = Coupled_Transport_Solvers()
julia> ms.add_solver(m)
```
"""
function add_solver(this::Coupled_Transport_Solvers,method::Discrete_Ordinates)
    this.number_of_particles += 1
    push!(this.particles,method.particle)
    push!(this.methods_list,method)
end

"""
    set_number_of_generations(this::Coupled_Transport_Solvers,number_of_generations::Int64)

To set the number of particle generation to transport during calculations

# Input Argument(s)
- `this::Coupled_Transport_Solvers`: collection of discretization method.
- `number_of_generations::Int64`: number of particle generation to transport during calculations.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ms = Coupled_Transport_Solvers()
julia> ms.set_number_of_generations(2)
```
"""
function set_number_of_generations(this::Coupled_Transport_Solvers,number_of_generations::Int64)
    if number_of_generations < 1 error("Number of particle generations should be at least 1.") end
    this.number_of_generations = number_of_generations
end

function get_method(this::Coupled_Transport_Solvers,particle::String)
    index = findfirst(x -> x == particle,this.particles)
    return this.methods_list[index]
end

function get_number_of_generations(this::Coupled_Transport_Solvers)
    return this.number_of_generations
end

function get_particles(this::Coupled_Transport_Solvers)
    return this.particles
end

function get_number_of_particles(this::Coupled_Transport_Solvers)
    return this.number_of_particles
end