"""
    Surface_Source

Structure used to define a directionnal boundary source and its properties.

# User-defined field(s)

- ## Mandatory field(s)
    - `name::String`: name (or identifier) of the Surface_Source structure.
    - `particle::String`: type of particle emitted.
    - `energy_group::Int64`: energy group index in which the particle are emitted.
    - `direction::Vector{Float64}`: direction cosine.
    - `location::String`: boundary at which the source is located.
    - `boundaries::Vector{Float64}`: boundaries of the source along each axis [in cm].

- ## Optional field(s) - with default values
    - `intensity::Float64=1.0`: intensity [# particles/cm⁽ᴺ⁻¹⁾, where N is the geometry dimension].

# System-defined field(s)
- `is_build::Bool`: boolean value defining if the Multigroup_Cross_Sections was build or not.
- `surface_sources::Array{Union{Array{Float64},Float64}}`: formatted object for transport calculations.
- `normalization_factor::Float64`: normalization over the number of particles.

"""
mutable struct Surface_Source

    # Variable(s)
    name                       ::Union{Missing,String}
    particle                   ::Union{Missing,String}
    intensity                  ::Union{Missing,Float64}
    energy_group               ::Union{Missing,Int64}
    direction                  ::Union{Missing,Vector{Float64}}
    location                   ::Union{Missing,String}
    boundaries                 ::Dict{String,Vector{Float64}}
    is_build                   ::Bool
    surface_sources            ::Union{Missing,Array{Union{Array{Float64},Float64}}}
    normalization_factor       ::Float64

    # Constructor(s)
    function Surface_Source()

        this = new()

        this.name = missing
        this.particle = missing
        this.intensity = 1.0
        this.energy_group = missing
        this.direction = missing
        this.location = missing
        this.boundaries = Dict{String,Vector{Float64}}()
        this.is_build = false
        this.normalization_factor = 0.0
        this.surface_sources = missing

        return this
    end
end

# Method(s)
Base.propertynames(::Surface_Source) = 
(
    fieldnames(Surface_Source)...,
    :set_particle,
    :set_intensity,
    :set_energy_group,
    :set_direction,
    :set_location,
    :set_boundaries,
    :get_particle,
    :get_normalization_factor,
    :build
)

"""
    set_particle(this::Surface_Source,particle::String)

To define the source particle.

# Input Argument(s)
- `this::Surface_Source`: surface source.
- `particle::String`: type of particle, which can takes the following values:
    - `particle = "photons"`: photons.
    - `particle = "electrons"`: electrons.
    - `particle = "positrons"`: positrons.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ss = Surface_Source()
julia> ss.set_particle("electrons")
```
"""
function set_particle(this::Surface_Source,particle::String)
    if particle ∉ ["electrons","positrons","photons"] end
    this.particle = particle
end

"""
    set_intensity(this::Surface_Source,intensity::Real)

To define the intensity of the source.

# Input Argument(s)
- `this::Surface_Source`: surface source.
- `intensity::Float64`: intensity [# particles/cmᴺ, where N is the geometry dimension]

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ss = Surface_Source()
julia> ss.set_intensity(100)
```
"""
function set_intensity(this::Surface_Source,intensity::Real)
    if intensity ≤ 0 error("The intensity should be greater than 0.") end
    this.intensity = intensity
end

"""
    set_energy_group(this::Surface_Source,energy_group::Int64)

To define the energy of the source by setting the energy group in which they are produced.

# Input Argument(s)
- `this::Surface_Source`: surface source.
- `energy_group::Int64`: energy group index.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ss = Surface_Source()
julia> ss.set_energy_group(1)
```
"""
function set_energy_group(this::Surface_Source,energy_group::Int64)
    if energy_group ≤ 0 error("The energy group number should be greater than 0.") end
    this.energy_group = energy_group
end

"""
    set_direction(this::Surface_Source,direction::Vector{Float64})

To set a direction of the source.

# Input Argument(s)
- `this::Surface_Source`: surface source.
- `direction::Vector{Float64}`: director cosines [μ,η,ξ]

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ss = Surface_Source()
julia> ss.set_direction([1.0,0.0,0.0])
```
"""
function set_direction(this::Surface_Source,direction::Vector{Float64})
    if length(direction) != 3 error("Three director cosines has to be provided.") end
    if abs(sum(direction).^2-1) > 1e-3 error("The sum of the squared three director cosines is not equal to 1.") end
    this.direction = direction
end

"""
    set_location(this::Surface_Source,location::String)

To set the location of the surface source.

# Input Argument(s)
- `this::Surface_Source`: surface source.
- `location::String`: boundary on which the surface source is, which can takes the following value:
    - `boundary = "x-"`: the lower bound along x-axis
    - `boundary = "x+"`: the upper bound along x-axis
    - `boundary = "y-"`: the lower bound along y-axis
    - `boundary = "y+"`: the upper bound along y-axis
    - `boundary = "z-"`: the lower bound along z-axis
    - `boundary = "z+"`: the upper bound along z-axis

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ss = Surface_Source()
julia> ss.set_location("x-")
```
"""
function set_location(this::Surface_Source,location::String)
    if uppercase(location) ∉ ["X-","X+","Y-","Y+","Z-","Z+"] error("Unknown location.") end
    this.location = location
end

"""
    set_boundaries(this::Surface_Source,axis::String,boundaries::Vector{Float64})

To define the boundaries of the source along the specified axis.

# Input Argument(s)
- `this::Surface_Source`: surface source.
- `axis::String`: axis along which the boundaries are defined.
- `boundaries::Vector{Float64}`: boundaries of the source in accending order [in cm]

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ss = Surface_Source()
julia> ss.set_boundaries("x",[1.0,3.0])
```
"""
function set_boundaries(this::Surface_Source,axis::String,boundaries::Vector{Float64})
    if lowercase(axis) ∉ ["x","y","z"] error("Unknown axis.") end
    if length(boundaries) != 2 error("The two boundary side should be provided.") end
    this.boundaries[axis] = boundaries
end

function get_particle(this::Surface_Source)
    return this.particle
end

function get_normalization_factor(this::Surface_Source)
    return this.normalization_factor
end