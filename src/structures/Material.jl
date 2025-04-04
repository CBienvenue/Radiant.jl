"""
    Material

Structure used to define a material and its properties.

# Mandatory field(s)
- `name::String` : name (or identifier) of the Material structure.
- `density::Float64` : density [in g/cm³].
- `elements::Vector{String}` : vector of the element in the composition of the material.
- `weight_fractions::Vector{Float64}` : vector of the weight fraction for each element in the composition of the material.

# Optional field(s) - with default values
- `state_of_matter::String = "solid"` : state of the matter.

"""
mutable struct Material

    # Variable(s)
    id                  ::Int64
    density             ::Union{Missing,Float64}
    state_of_matter     ::String
    number_of_elements  ::Int64
    elements            ::Vector{String}
    atomic_numbers      ::Vector{Int64}
    weight_fractions    ::Vector{Float64}

    # Constructor(s)
    function Material()
        this = new()
        this.id                 = generate_unique_id()
        this.density            = missing
        this.state_of_matter    = "solid"
        this.number_of_elements = 0
        this.elements           = Vector{String}()
        this.atomic_numbers     = Vector{Int64}()
        this.weight_fractions   = Vector{Float64}()
        return this
    end
end

# Method(s)
"""
    println(this::Material)

To print the material properties.

# Input Argument(s)
- `this::Material` : material.

# Output Argument(s)
N/A

"""
function println(this::Material)
    println("Material:")
    println("   Density (g/cm³):              $(this.density)")
    println("   Elements in the compound:     $(this.elements)")
    println("   Weight fractions:             $(this.weight_fractions)")
    println("   State of matter:              $(this.state_of_matter)")
end

"""
    set_density(this::Material,density::Real)

To set the density of the material.

# Input Argument(s)
- `this::Material` : material.
- `density::Real` : density [in g/cm³].

# Output Argument(s)
N/A

# Examples
```jldoctest
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
- `this::Material` : material.
- `state_of_matter::String` : state of matter, which value is given by "solid", "liquid" or "gaz".

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> mat.set_state_of_matter("liquid")
```
"""
function set_state_of_matter(this::Material,state_of_matter::String)
    if state_of_matter ∉ ["solid","liquid","gaz"] error("The state of matter is either solid, liquid or gaz.") end
    this.state_of_matter = state_of_matter
end

"""
    add_element(this::Material,symbol::String,weight_fraction::Real=1)

To add an element which is part of the composition of the material.

# Input Argument(s)
- `this::Material` : material.
- `symbol::String` : element symbol, with value such as "H", "He", etc.
- `weight_fraction::Real` : weight fraction of the element in the material (between 0 and 1).

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> water.add_element("H",0.1111)
julia> water.add_element("O",0.8889)
```
"""
function add_element(this::Material,symbol::String,weight_fraction::Real=1)
    if ~(0 ≤ weight_fraction ≤ 1) error("Weight fraction should have values between 0 and 1.") end
    push!(this.elements,lowercase(symbol))
    push!(this.atomic_numbers,atomic_number(lowercase(symbol)))
    push!(this.weight_fractions,weight_fraction)
    this.number_of_elements += 1
    if sum(this.weight_fractions) > 1 error("Weight fraction exceed 1.") end
    if (weight_fraction == 1) this.set_density(density(atomic_number(lowercase(symbol)))) end
end

"""
    get_density(this::Material)

To get the density of the material.

# Input Argument(s)
- `this::Material` : material.

# Output Argument(s)
- `density::Real` : density [in g/cm³].

# Examples
```jldoctest
julia> density = mat.get_density()
```
"""
function get_density(this::Material)
    if ismissing(density) error("Density is missing.") end
    return this.density
end

"""
    get_number_of_elements(this::Material)

To get the number of elements in the material composition.

# Input Argument(s)
- `this::Material` : material.

# Output Argument(s)
- `number_of_elements::Int64` : number of elements.

# Examples
```jldoctest
julia> N = mat.get_number_of_elements()
```
"""
function get_number_of_elements(this::Material)
    return this.number_of_elements
end

"""
    get_atomic_numbers(this::Material)

To get the atomic numbers of the elements in the material composition.

# Input Argument(s)
- `this::Material` : material.

# Output Argument(s)
- `atomic_numbers::Vector{Int64}` : atomic numbers.

# Examples
```jldoctest
julia> atomic_number = mat.get_atomic_numbers()
```
"""
function get_atomic_numbers(this::Material)
    return this.atomic_numbers
end

"""
    get_weight_fractions(this::Material)

To get the weight fractions associated with the elements in the material composition.

# Input Argument(s)
- `this::Material` : material.

# Output Argument(s)
- `weight_fractions::Vector{Float64}` : weight fractions.

# Examples
```jldoctest
julia> weight_fractions = mat.weight_fractions()
```
"""
function get_weight_fractions(this::Material)
    return this.weight_fractions
end

"""
    get_id(this::Material)

To get the unique ID of the material.

# Input Argument(s)
- `this::Material` : material.

# Output Argument(s)
- `ID::Int64` : ID of the material.

# Examples
```jldoctest
julia> id = mat.get_id()
```
"""
function get_id(this::Material)
    return this.id
end

"""
    get_state_of_matter(this::Material)

To get the state of matter of the material.

# Input Argument(s)
- `this::Material` : material.

# Output Argument(s)
- `state_of_matter::String` : state of matter.

# Examples
```jldoctest
julia> state_of_matter = mat.get_state_of_matter()
```
"""
function get_state_of_matter(this::Material)
    return this.state_of_matter
end