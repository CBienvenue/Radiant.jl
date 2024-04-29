"""
    Material

Structure used to define a material and its properties.

# User-defined field(s)

- ## Mandatory field(s)
    - `name::String`: name (or identifier) of the Material structure.
    - `density::Float64`: density [in g/cm³].
    - `elements::Vector{String}`: vector of the element in the composition of the material.
    - `weight_fractions::Vector{Float64}`: vector of the weight fraction for each element in the composition of the material.

- ## Optional field(s) - with default values
    - `state_of_matter::String = "solid"`: state of the matter.

# System-defined field(s)
- `number_of_elements::Int64`: number of elements in the composition of the material.
- `atomic_numbers::Vector{Int64}`: vector of the atomic number corresponding to the element in the composition of the material.

"""
mutable struct Material

    # Variable(s)
    name                ::Union{Missing,String}
    density             ::Union{Missing,Float64}
    state_of_matter     ::String
    number_of_elements  ::Int64
    elements            ::Vector{String}
    atomic_numbers      ::Vector{Int64}
    weight_fractions    ::Vector{Float64}

    # Constructor(s)
    function Material()
        
        this = new()

        this.name = missing
        this.density = missing
        this.state_of_matter = "solid"
        this.number_of_elements = 0
        this.elements = Vector{String}()
        this.atomic_numbers = Vector{Int64}()
        this.weight_fractions = Vector{Float64}()

        return this
    end
end

# Method(s)
Base.propertynames(::Material) = 
(
    fieldnames(Material)...,
    :set_name,
    :set_density,
    :set_state_of_matter,
    :add_element,
    :get_density,
    :get_atomic_numbers,
    :weight_fractions
    :get_name
)

"""
    set_name(this::Material,name::String)

To set the identifier (name) of the Material structure.

# Input Argument(s)
- `this::Material`: material.
- `name::String`: identifier (name).

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> mat = Material()
julia> mat.set_name("water")
```
"""
function set_name(this::Material,name::String)
    this.name = name
end

"""
    set_density(this::Material,density::Real)

To set the density of the material.

# Input Argument(s)
- `this::Material`: material.
- `density::Real`: density [in g/cm³].

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> mat = Material()
julia> mat.set_density(19.3)
```
"""
function set_density(this::Material,density::Real)
    if density < 0 error("The density should be positive.") end
    this.density = density
end

"""
    set_state_of_matter(this::Material,state_of_matter::String)

To set the state of the matter (solid or liquid or gaz).

# Input Argument(s)
- `this::Material`: material.
- `state_of_matter::String`: state of matter, which value is given by "solid", "liquid" or "gaz".

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> mat = Material()
julia> mat.set_state_of_matter("liquid")
```
"""
function set_state_of_matter(this::Material,state_of_matter::String)
    if state_of_matter ∉ ["solid","liquid","gaz"] error("The state of matter is either solid, liquid or gaz.") end
    this.state_of_matter = state_of_matter
end

"""
    add_element(this::Material,symbol::String,weight_fraction::Real)

To add an element which is part of the composition of the material.

# Input Argument(s)
- `this::Material`: material.
- `symbol::String`: element symbol, with value such as "H", "He", etc.
- `weight_fraction::Real`: weight fraction of the element in the material (between 0 and 1).

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> water = Material()
julia> water.add_element("H",0.1111)
julia> water.add_element("O",0.8889)
```
"""
function add_element(this::Material,symbol::String,weight_fraction::Real)
    if ~(0 ≤ weight_fraction ≤ 1) error("Weight fraction should have values between 0 and 1.") end
    push!(this.elements,symbol)
    push!(this.atomic_numbers,atomic_number(lowercase(symbol)))
    push!(this.weight_fractions,weight_fraction)
    this.number_of_elements += 1
end

function get_density(this::Material)
    if ismissing(density) error("Unknown density.") end
    return this.density
end

function get_number_of_elements(this::Material)
    return this.number_of_elements
end

function get_atomic_numbers(this::Material)
    return this.atomic_numbers
end

function get_weight_fractions(this::Material)
    return this.weight_fractions
end

function get_name(this::Material)
    return this.name
end

function get_state_of_matter(this::Material)
    return this.state_of_matter
end