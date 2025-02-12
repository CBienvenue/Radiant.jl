"""
    Source

Structure used to describe sources for a given particle.

"""
mutable struct Source

    # Variable(s)
    name                       ::Union{Missing,String}
    particle                   ::Union{Missing,Particle}
    volume_sources             ::Array{Float64}
    surface_sources            ::Array{Union{Array{Float64},Float64}}
    normalization_factor       ::Float64
    cross_sections             ::Cross_Sections
    geometry                   ::Geometry
    discrete_ordinates         ::Discrete_Ordinates

    # Constructor(s)
    function Source(particle::Particle,cross_sections::Cross_Sections,geometry::Geometry,discrete_ordinates::Discrete_Ordinates)
        this = new()
        this.name = missing
        this.particle = particle
        this.normalization_factor = 0
        this.cross_sections = cross_sections
        this.geometry = geometry
        this.discrete_ordinates = discrete_ordinates
        initalize_sources(this,cross_sections,geometry,discrete_ordinates)
        return this
    end
end

# Method(s)
"""
    initalize_sources(this::Source,cross_sections::Cross_Sections,geometry::Geometry,discrete_ordinates::Discrete_Ordinates)

Initialize volume and boundary sources.

# Input Argument(s)
- `this::Source` : source structure.
- `cross_sections::Cross_Sections` : cross-sections library.
- `geometry::Geometry` : geometry.
- `discrete_ordinates::Discrete_Ordinates` : discrete ordinates solver.

# Output Argument(s)
N/A

"""
function initalize_sources(this::Source,cross_sections::Cross_Sections,geometry::Geometry,discrete_ordinates::Discrete_Ordinates)
    particle = this.particle
    if get_id(particle) ∉ get_id.(cross_sections.particles) error(string("No cross sections available for ",particle," particle.")) end
    index = findfirst(x -> get_id(x) == get_id(particle),cross_sections.particles)
    Ng = cross_sections.number_of_groups[index]
    Nx = geometry.number_of_voxels["x"]
    Qdims = discrete_ordinates.get_quadrature_dimension(geometry.dimension)
    if geometry.dimension ≥ 2 Ny = geometry.number_of_voxels["y"] else Ny = 1 end
    if geometry.dimension ≥ 3 Nz = geometry.number_of_voxels["z"] else Nz = 1 end
    Ω,w = quadrature(discrete_ordinates.quadrature_order,discrete_ordinates.quadrature_type,Qdims)
    if typeof(Ω) == Vector{Float64} Ω = [Ω,0*Ω,0*Ω] end
    number_of_directions = length(w)
    P,_,_,_ = angular_polynomial_basis(geometry.dimension,Ω,w,discrete_ordinates.get_legendre_order(),discrete_ordinates.quadrature_order,discrete_ordinates.get_angular_boltzmann(),Qdims)
    _,_,Nm = discrete_ordinates.get_schemes(geometry,true)

    this.volume_sources = zeros(Ng,P,Nm[5],Nx,Ny,Nz)
    this.surface_sources = Array{Union{Array{Float64},Float64}}(undef,Ng,number_of_directions,2*geometry.dimension)
    this.surface_sources .= 0.0
end

"""
    set_particle(this::Source,particle::String)

Set the source particle.

# Input Argument(s)
- `this::Source` : source structure.
- `particle::String` : particle.

# Output Argument(s)
N/A

"""
function set_particle(this::Source,particle::String)
    if lowercase(particle) ∉ ["photons","electrons","positrons"] error("Unknown particle type") end
    this.particle = particle
end

"""
    add_source(this::Source,source::Volume_Source)

Add a volume source to the source structure.

# Input Argument(s)
- `this::Source` : source structure.
- `source::Volume_Source` : volume source.

# Output Argument(s)
N/A

"""
function add_source(this::Source,source::Volume_Source)
    particle = source.particle
    Qv,norm = volume_source(particle,source,this.cross_sections,this.geometry)
    this.volume_sources[:,1,1,:,:,:] += Qv[:,1,1,:,:,:]
    source.normalization_factor += norm
end

"""
    add_source(this::Source,surface_sources::Surface_Source)

Add a surface source to the source structure.

# Input Argument(s)
- `this::Source` : source structure.
- `source::Surface_Source` : surface source.

# Output Argument(s)
N/A

"""
function add_source(this::Source,surface_sources::Surface_Source)
    particle = surface_sources.particle
    if get_id(particle) ∉ get_id.(this.cross_sections.particles) error(string("No cross sections available for ",get_type(particle)," particle.")) end
    index = findfirst(x -> get_id(x) == get_id(particle),this.cross_sections.particles)
    Ng = this.cross_sections.number_of_groups[index]
    if this.geometry.dimension ≥ 2 Ny = this.geometry.number_of_voxels["y"] else Ny = 1 end
    if this.geometry.dimension ≥ 3 Nz = this.geometry.number_of_voxels["z"] else Nz = 1 end
    if get_id(particle) != get_id(this.discrete_ordinates.particle) error(string("No methods available for ",get_type(particle)," particle.")) end
    Qdims = this.discrete_ordinates.get_quadrature_dimension(this.geometry.dimension)
    _,w = quadrature(this.discrete_ordinates.quadrature_order,this.discrete_ordinates.quadrature_type,Qdims)
    number_of_directions = length(w)
    Q = Array{Union{Array{Float64},Float64}}(undef,Ng,number_of_directions,2*this.geometry.dimension)
    Q .= 0.0
    Q,norm = surface_source(Q,particle,surface_sources,this.cross_sections,this.geometry,this.discrete_ordinates)
    d = size(this.surface_sources)
    for i in range(1,d[1]), j in range(1,d[2]), k in range(1,d[3])
        if length(size(this.surface_sources[i,j,k])) == length(size(Q[i,j,k])) 
            this.surface_sources[i,j,k] += Q[i,j,k]
        elseif length(size(this.surface_sources[i,j,k])) < length(size(Q[i,j,k])) 
            this.surface_sources[i,j,k] = Q[i,j,k]
        end
    end
    surface_sources.normalization_factor += norm
end

"""
    add_volume_source(this::Source,source::Array{Float64})

Set a volume source.

# Input Argument(s)
- `this::Source` : source structure.
- `source::Array{Float64}` : volume source.

# Output Argument(s)
N/A

"""
function add_volume_source(this::Source,source::Array{Float64})
    this.volume_sources = source
end

"""
    get_surface_sources(this::Source)

Get the surface sources.

# Input Argument(s)
- `this::Source` : source structure.

# Output Argument(s)
- `surface_sources::Array{Array{Float64}}` : surfaces sources.

"""
function get_surface_sources(this::Source)
    return this.surface_sources
end

"""
    get_volume_sources(this::Source)

Get the volume sources.

# Input Argument(s)
- `this::Source` : source structure.

# Output Argument(s)
- `volume_sources::Array{Float64}` : volume sources.

"""
function get_volume_sources(this::Source)
    return this.volume_sources
end

"""
    get_normalization_factor(this::Source)

Get the source normalization factor.

# Input Argument(s)
- `this::Source` : source structure.

# Output Argument(s)
- `normalization_factor::Float64` : normalization factor.

"""
function get_normalization_factor(this::Source)
    return this.normalization_factor
end

"""
    get_particle(this::Source)

Get the source particle.

# Input Argument(s)
- `this::Source` : source structure.

# Output Argument(s)
- `particle::Particle` : particle.

"""
function get_particle(this::Source)
    return this.particle
end 

"""
    Base.:+(source1::Source,source2::Source)

Combination of two sources.

# Input Argument(s)
- `source1::Source` : a source structure.
- `source2::Source` : a source structure.

# Output Argument(s)
- `source1::Source` : a combination of the two sources.

"""
function Base.:+(source1::Source,source2::Source)
    if get_id(source1.get_particle()) != get_id(source2.get_particle()) error("Forbitten addition of different particle sources.") end
    source1.volume_sources += source2.volume_sources
    d = size(source1.surface_sources)
    for i in range(1,d[1]), j in range(1,d[2]), k in range(1,d[3])
        if length(size(source1.surface_sources[i,j,k])) == length(size(source2.surface_sources[i,j,k])) 
            source1.surface_sources[i,j,k] += source2.surface_sources[i,j,k]
        elseif length(size(source1.surface_sources[i,j,k])) < length(size(source2.surface_sources[i,j,k])) 
            source1.surface_sources[i,j,k] = source2.surface_sources[i,j,k]
        end
    end
    return source1
end