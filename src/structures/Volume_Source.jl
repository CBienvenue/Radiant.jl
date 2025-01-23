"""
    Volume_Source

Structure used to define an isotropic volume source and its properties.

# Mandatory field(s)
- `name::String` : name (or identifier) of the Volume_Source structure.
- `particle::Particle` : type of particle emitted.
- `energy_group::Int64` : energy group index in which the particle are emitted.
- `boundaries::Vector{Float64}` : boundaries of the source along each axis [in cm].

# Optional field(s) - with default values
- `intensity::Float64=1.0` : intensity [# particles/cmᴺ, where N is the geometry dimension].

"""
mutable struct Volume_Source

    # Variable(s)
    name                       ::Union{Missing,String}
    particle                   ::Union{Missing,Particle}
    intensity                  ::Union{Missing,Float64}
    energy_group               ::Union{Missing,Int64}
    boundaries                 ::Dict{String,Vector{Float64}}
    is_build                   ::Bool
    volume_sources             ::Union{Missing,Array{Float64}}
    normalization_factor       ::Float64

    # Constructor(s)
    function Volume_Source()

        this = new()

        this.name = missing
        this.particle = missing
        this.intensity = 1.0
        this.energy_group = missing
        this.boundaries = Dict{String,Vector{Float64}}()
        this.is_build = false
        this.normalization_factor = 0.0

        return this
    end
end

# Method(s)
"""
    set_particle(this::Volume_Source,particle::String)

To define the source particle.

# Input Argument(s)
- `this::Volume_Source` : volume source.
- `particle::Particle` : particle.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> vs = Volume_Source()
julia> vs.set_particle(electron)
```
"""
function set_particle(this::Volume_Source,particle::Particle)
    this.particle = particle
end

"""
    set_intensity(this::Volume_Source,intensity::Real)

To define the intensity of the source.

# Input Argument(s)
- `this::Volume_Source` : volume source.
- `intensity::Float64` : intensity [# particles/cmᴺ, where N is the geometry dimension]

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> vs = Volume_Source()
julia> vs.set_intensity(100)
```
"""
function set_intensity(this::Volume_Source,intensity::Real)
    if intensity ≤ 0 error("The intensity should be greater than 0.") end
    this.intensity = intensity
end

"""
    set_energy_group(this::Volume_Source,energy_group::Int64)

To define the energy of the source by setting the energy group in which they are produced.

# Input Argument(s)
- `this::Volume_Source` : volume source.
- `energy_group::Int64` : energy group index.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> vs = Volume_Source()
julia> vs.set_energy_group(1)
```
"""
function set_energy_group(this::Volume_Source,energy_group::Int64)
    if energy_group ≤ 0 error("The energy group number should be greater than 0.") end
    this.energy_group = energy_group
end

"""
    set_boundaries(this::Volume_Source,axis::String,boundaries::Vector{Float64})

To define the boundaries of the source along the specified axis.

# Input Argument(s)
- `this::Volume_Source` : volume source.
- `axis::String` : axis along which the boundaries are defined.
- `boundaries::Vector{Float64}` : boundaries of the source in accending order [in cm]

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> vs = Volume_Source()
julia> vs.set_boundaries("x",[1.0,3.0])
```
"""
function set_boundaries(this::Volume_Source,axis::String,boundaries::Vector{Float64})
    if lowercase(axis) ∉ ["x","y","z"] error("Unknown axis.") end
    if length(boundaries) != 2 error("The two boundary side should be provided.") end
    this.boundaries[axis] = boundaries
end

function get_particle(this::Volume_Source)
    return this.particle
end

function get_normalization_factor(this::Volume_Source)
    return this.normalization_factor
end