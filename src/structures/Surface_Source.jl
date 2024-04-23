
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
        this.intensity = missing
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

function println(this::Surface_Source)
    entries = ["Name","Intensity","Energy group","Location","Direction","Boundaries"]
    values = [this.name,this.intensity,this.energy_group,this.location,this.direction,this.boundaries]
    N = length(entries); L = length.(entries); Lmax = maximum(L)
    println("Surface_Source")
    for n in range(1,N)
        println(string("   ",entries[n]," "^(Lmax-L[n])),"  :  ",values[n])
    end
end

function set_particle(this::Surface_Source,particle::String)
    if particle ∉ ["electrons","positrons","photons"] end
    this.particle = particle
end

function set_intensity(this::Surface_Source,intensity::Float64)
    if intensity ≤ 0 error("The intensity should be greater than 0.") end
    this.intensity = intensity
end

function set_energy_group(this::Surface_Source,energy_group::Int64)
    if energy_group ≤ 0 error("The energy group number should be greater than 0.") end
    this.energy_group = energy_group
end

function set_direction(this::Surface_Source,direction::Vector{Float64})
    if length(direction) != 3 error("Three director cosines has to be provided.") end
    if abs(sum(direction).^2-1) > 1e-3 error("The sum of the squared three director cosines is not equal to 1.") end
    this.direction = direction
end

function set_location(this::Surface_Source,location::String)
    if uppercase(location) ∉ ["X-","X+","Y-","Y+","Z-","Z+"] error("Unknown location.") end
    this.location = location
end

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