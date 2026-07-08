"""
    Electromagnetic_Field

Structure used to define an external electromagnetic field, constant over the geometry, that
acts on charged particles through the Lorentz force during transport. It is built by the user
and attached to a `Computation_Unit` (like `Cross_Sections`, `Geometry`, `Solvers` and
`Fixed_Sources`).

Only **magnetic** fields are currently supported (electric-field transport, which redistributes
particles in energy, is not yet implemented). Field components are given along the geometry axes
(x, y, z); the magnetic field is expressed in **tesla**.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `electric_field::Vector{Float64} = [0.0,0.0,0.0]` : electric field along x-, y- and z-axis.
- `magnetic_field::Vector{Float64} = [0.0,0.0,0.0]` : magnetic field [T] along x-, y- and z-axis.

"""
mutable struct Electromagnetic_Field

    # Variable(s)
    electric_field ::Vector{Float64}
    magnetic_field ::Vector{Float64}

    # Constructor(s)
    function Electromagnetic_Field()
        this = new()
        this.electric_field = [0.0,0.0,0.0]
        this.magnetic_field = [0.0,0.0,0.0]
        return this
    end
end

# Method(s)
"""
    set_magnetic_field(this::Electromagnetic_Field,magnetic_field::Vector{<:Real})

To set the external magnetic field vector [in tesla], expressed along the x-, y- and z-axis of
the geometry.

# Input Argument(s)
- `this::Electromagnetic_Field` : electromagnetic field.
- `magnetic_field::Vector{<:Real}` : magnetic field [T] along x-, y- and z-axis.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> emf = Electromagnetic_Field()
julia> emf.set_magnetic_field([0.0,0.0,1.5])
```
"""
function set_magnetic_field(this::Electromagnetic_Field,magnetic_field::Vector{<:Real})
    if length(magnetic_field) != 3 error("Magnetic field should be a 3-vector [Bx,By,Bz].") end
    this.magnetic_field = Float64.(magnetic_field)
end

"""
    set_electric_field(this::Electromagnetic_Field,electric_field::Vector{<:Real})

To set the external electric field vector, expressed along the x-, y- and z-axis of the geometry.

!!! warning
    Electric-field transport is not yet implemented (its energy-redistribution operator couples
    the energy groups). Setting any non-zero component raises an error.

# Input Argument(s)
- `this::Electromagnetic_Field` : electromagnetic field.
- `electric_field::Vector{<:Real}` : electric field along x-, y- and z-axis.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> emf = Electromagnetic_Field()
julia> emf.set_electric_field([0.0,0.0,0.0])
```
"""
function set_electric_field(this::Electromagnetic_Field,electric_field::Vector{<:Real})
    if length(electric_field) != 3 error("Electric field should be a 3-vector [Ex,Ey,Ez].") end
    if any(electric_field .!= 0) error("Electric-field transport not yet implemented; only magnetic fields are supported.") end
    this.electric_field = Float64.(electric_field)
end

"""
    get_magnetic_field(this::Electromagnetic_Field)

Get the external magnetic field vector [in tesla].

# Input Argument(s)
- `this::Electromagnetic_Field` : electromagnetic field.

# Output Argument(s)
- `magnetic_field::Vector{Float64}` : magnetic field [T] along x-, y- and z-axis.

"""
function get_magnetic_field(this::Electromagnetic_Field)
    return this.magnetic_field
end

"""
    get_electric_field(this::Electromagnetic_Field)

Get the external electric field vector.

# Input Argument(s)
- `this::Electromagnetic_Field` : electromagnetic field.

# Output Argument(s)
- `electric_field::Vector{Float64}` : electric field along x-, y- and z-axis.

"""
function get_electric_field(this::Electromagnetic_Field)
    return this.electric_field
end

"""
    get_electromagnetic_field(this::Electromagnetic_Field)

Get whether the field is active together with the electric and magnetic field vectors.

# Input Argument(s)
- `this::Electromagnetic_Field` : electromagnetic field.

# Output Argument(s)
- `is_EM::Bool` : indicate if a non-zero field is present.
- `electric_field::Vector{Float64}` : electric field along x-, y- and z-axis.
- `magnetic_field::Vector{Float64}` : magnetic field [T] along x-, y- and z-axis.

"""
function get_electromagnetic_field(this::Electromagnetic_Field)
    is_EM = any(this.magnetic_field .!= 0) || any(this.electric_field .!= 0)
    return is_EM, this.electric_field, this.magnetic_field
end
