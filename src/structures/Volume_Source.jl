
mutable struct Volume_Source

    # Variable(s)
    name                       ::Union{Missing,String}
    particle                   ::Union{Missing,String}
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
        this.intensity = missing
        this.energy_group = missing
        this.boundaries = Dict{String,Vector{Float64}}()
        this.is_build = false
        this.normalization_factor = 0.0

        return this
    end
end

# Method(s)
Base.propertynames(::Volume_Source) = 
(
    fieldnames(Volume_Source)...,
    :set_particle,
    :set_intensity,
    :set_energy_group,
    :set_boundaries,
    :get_particle,
    :get_normalization_factor,
    :build
)

function println(this::Volume_Source)
    entries = ["Name","Intensity","Energy group","Boundaries"]
    values = [this.name,this.intensity,this.energy_group,this.boundaries]
    N = length(entries); L = length.(entries); Lmax = maximum(L)
    println("Volume_Source")
    for n in range(1,N)
        println(string("   ",entries[n]," "^(Lmax-L[n])),"  :  ",values[n])
    end
end

function set_particle(this::Volume_Source,particle::String)
    if particle ∉ ["electrons","positrons","photons"] end
    this.particle = particle
end

function set_intensity(this::Volume_Source,intensity::Float64)
    if intensity ≤ 0 error("The intensity should be greater than 0.") end
    this.intensity = intensity
end

function set_energy_group(this::Volume_Source,energy_group::Int64)
    if energy_group ≤ 0 error("The energy group number should be greater than 0.") end
    this.energy_group = energy_group
end

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