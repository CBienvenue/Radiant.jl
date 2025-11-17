const Solver = Union{Discrete_Ordinates, Spherical_Harmonics}
"""
    Solvers

Structure used to define the collection of discretization methods for transport calculations associated with each of the particle and additionnal coupled transport informations.

# Mandatory field(s)
- `methods_list::Vector{Discrete_Ordinates}` : list of the particle methods
- `maximum_number_of_generations::Int64` : number of particle generations to transport.

# Optional field(s) - with default values
- N/A

"""
mutable struct Solvers

    # Variable(s)
    number_of_particles                ::Int64
    particles                          ::Vector{Particle}
    methods_names                      ::Vector{String}
    methods_list                       ::Vector{Solver}
    maximum_number_of_generations      ::Int64
    convergence_criterion              ::Real
    convergence_type                   ::String
    
    # Constructor(s)
    function Solvers()

        this = new()

        this.number_of_particles = 0
        this.particles = Vector{Particle}()
        this.methods_names = Vector{String}()
        this.methods_list = Vector{Solver}()
        this.maximum_number_of_generations = 10
        this.convergence_criterion = 1e-7
        this.convergence_type = "flux"

        return this
    end
end

# Method(s)
"""
    add_solver(this::Solvers,method::Solver)

To add a particle and is associated methods to the Solvers structure.

# Input Argument(s)
- `this::Solvers` : collection of discretization method.
- `method::Solver` : discretization method.

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
function add_solver(this::Solvers,method::Solver)
    this.number_of_particles += 1
    push!(this.particles,method.particle)
    push!(this.methods_list,method)
end

"""
    set_maximum_number_of_generations(this::Solvers,maximum_number_of_generations::Int64)

To set the maximum number of particle generation to transport during calculations.

# Input Argument(s)
- `this::Solvers` : collection of discretization method.
- `maximum_number_of_generations::Int64` : maximum number of particle generation to
   transport during calculations.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ms = Solvers()
julia> ms.set_maximum_number_of_generations(2)
```
"""
function set_maximum_number_of_generations(this::Solvers,maximum_number_of_generations::Int64)
    if maximum_number_of_generations < 1 error("Number of particle generations should be at least 1.") end
    this.maximum_number_of_generations = maximum_number_of_generations
end

"""
    set_convergence_criterion(this::Solvers,convergence_criterion::Real)

To set the convergence criterion for coupled particle transport.

# Input Argument(s)
- `this::Solvers` : collection of discretization method.
- `convergence_criterion::Real` : convergence criterion.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ms = Solvers()
julia> ms.set_convergence_criterion(1e-7)
```
"""
function set_convergence_criterion(this::Solvers,convergence_criterion::Real)
    if convergence_criterion ≤ 0 error("Convergence criterion has to be greater than 0.") end
    this.convergence_criterion = convergence_criterion
end

"""
    set_convergence_type(this::Solvers,convergence_type::String)

To set the convergence type.

# Input Argument(s)
- `this::Solvers` : collection of discretization method.
- `convergence_type::String` : convergence type, which can be:
    - `convergence_type = "flux"` : convergence of the fluxes per particle.
    - `convergence_type = "energy-deposition"` : convergence of the total energy deposition.
    - `convergence_type = "charge-deposition"` : convergence of the total charge deposition.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ms = Solvers()
julia> ms.set_convergence_type("energy-deposition")
```
"""
function set_convergence_type(this::Solvers,convergence_type::String)
    if convergence_type ∉ ["flux","energy-deposition","charge-deposition"] error("Unknown convergence type.") end
    this.convergence_type = convergence_type
end

"""
    get_method(this::Solvers,particle::Particle)

Get the method associated with a particle.

# Input Argument(s)
- `this::Solvers` : collection of discretization method.
- `particle::Particle` : particle.

# Output Argument(s)
- `method::Solver` : method associated with the particle.

"""
function get_method(this::Solvers,particle::Particle)
    index = findfirst(x -> get_id(x) == get_id(particle),this.particles)
    if isnothing(index) error("Solvers does not contain data for the given particle.") end
    return this.methods_list[index]
end

"""
    get_maximum_number_of_generations(this::Solvers)

Get the maximum number of generations.

# Input Argument(s)
- `this::Solvers` : collection of discretization method.

# Output Argument(s)
- `maximum_number_of_generations::Int64` :  maximum number of generations.

"""
function get_maximum_number_of_generations(this::Solvers)
    return this.maximum_number_of_generations
end

"""
    get_convergence_criterion(this::Solvers)

Get the convergence criterion.

# Input Argument(s)
- `this::Solvers` : collection of discretization method.

# Output Argument(s)
- `convergence_criterion::Float64` : convergence criterion.

"""
function get_convergence_criterion(this::Solvers)
    return this.convergence_criterion
end

"""
    get_convergence_type(this::Solvers)

Get the convergence type.

# Input Argument(s)
- `this::Solvers` : collection of discretization method.

# Output Argument(s)
- `convergence_type::String` : convergence type.

"""
function get_convergence_type(this::Solvers)
    return this.convergence_type
end

"""
    get_particles(this::Solvers)

Get the particle list.

# Input Argument(s)
- `this::Solvers` : collection of discretization method.

# Output Argument(s)
- `particle::Vector{Particle}` : particle list.

"""
function get_particles(this::Solvers)
    return this.particles
end

"""
    get_number_of_particles(this::Solvers)

Get the number of particles.

# Input Argument(s)
- `this::Solvers` : collection of discretization method.

# Output Argument(s)
- `number_of_particles::Int64` : number of particles.

"""
function get_number_of_particles(this::Solvers)
    return this.number_of_particles
end