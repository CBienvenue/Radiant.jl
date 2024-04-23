
mutable struct Source

    # Variable(s)
    name                       ::Union{Missing,String}
    particle                   ::Union{Missing,String}
    volume_sources             ::Array{Float64}
    surface_sources            ::Array{Union{Array{Float64},Float64}}
    normalization_factor       ::Float64
    cross_sections             ::Cross_Sections
    geometry                   ::Geometry
    method                     ::Method

    # Function(s)
    set_particle               ::Function
    add_source                 ::Function
    add_volume_source          ::Function
    get_surface_sources        ::Function
    get_volume_sources         ::Function
    get_normalization_factor   ::Function
    get_particle               ::Function

    # Constructor(s)
    function Source(particle::String,cross_sections::Cross_Sections,geometry::Geometry,method::Method)

        this = new()

        this.name = missing
        this.particle = particle
        this.normalization_factor = 0
        this.cross_sections = cross_sections
        this.geometry = geometry
        this.method = method
        initalize_sources!(this,cross_sections,geometry,method)

        this.set_particle = function (particle) set_particle!(this,particle) end
        this.add_source = function (source) add_source!(this,source) end
        this.add_volume_source = function (source) add_volume_source!(this,source) end
        this.get_surface_sources = function () get_surface_sources(this) end
        this.get_volume_sources = function () get_volume_sources(this) end
        this.get_normalization_factor = function () get_normalization_factor(this) end
        this.get_particle = function () get_particle(this) end

        return this
    end
end

# Method(s)
function println(this::Source)
    entries = ["Name","Particle"]
    values = [this.name,this.particle]
    N = length(entries); L = length.(entries); Lmax = maximum(L)
    println("Source")
    for n in range(1,N)
        println(string("   ",entries[n]," "^(Lmax-L[n])),"  :  ",values[n])
    end
end

function initalize_sources!(this::Source,cross_sections::Cross_Sections,geometry::Geometry,method::Method)
    particle = this.particle
    if particle ∉ cross_sections.particles error(string("No cross sections available for ",particle," particle.")) end
    index = findfirst(x -> x == particle,cross_sections.particles)
    Ng = cross_sections.number_of_groups[index]
    Nx = geometry.number_of_voxels["x"]
    if geometry.dimension ≥ 2 Ny = geometry.number_of_voxels["y"] else Ny = 1 end
    if geometry.dimension ≥ 3 Nz = geometry.number_of_voxels["z"] else Nz = 1 end
    Ω,w = quadrature(method.quadrature_order,method.quadrature_type,geometry.dimension)
    number_of_directions = length(w)
    P,_,_,_ = angular_polynomial_basis(geometry.dimension,Ω,w,method.get_legendre_order(),method.quadrature_order,method.get_angular_boltzmann())
    _,_,Nm = method.get_schemes(geometry,true)

    this.volume_sources = zeros(Ng,P,Nm[5],Nx,Ny,Nz)
    this.surface_sources = Array{Union{Array{Float64},Float64}}(undef,Ng,number_of_directions,2*geometry.dimension)
    this.surface_sources .= 0.0
end

function set_particle!(this::Source,particle::String)
    if lowercase(particle) ∉ ["photons","electrons","positrons"] error("Unknown particle type") end
    this.particle = particle
end

function add_source!(this::Source,source::Volume_Source)
    
    particle = source.particle
    Qv,norm = volume_source(particle,source,this.cross_sections,this.geometry)
    this.volume_sources[:,1,1,:,:,:] += Qv[:,1,1,:,:,:]
    source.normalization_factor += norm

end

function add_source!(this::Source,surface_sources::Surface_Source)

    particle = surface_sources.particle
    if particle ∉ this.cross_sections.particles error(string("No cross sections available for ",particle," particle.")) end
    index = findfirst(x -> x == particle,this.cross_sections.particles)
    Ng = this.cross_sections.number_of_groups[index]
    if this.geometry.dimension ≥ 2 Ny = this.geometry.number_of_voxels["y"] else Ny = 1 end
    if this.geometry.dimension ≥ 3 Nz = this.geometry.number_of_voxels["z"] else Nz = 1 end
    if particle != this.method.particle error(string("No methods available for ",particle," particle.")) end
    _,w = quadrature(this.method.quadrature_order,this.method.quadrature_type,this.geometry.dimension)
    number_of_directions = length(w)
    Q = Array{Union{Array{Float64},Float64}}(undef,Ng,number_of_directions,2*this.geometry.dimension)
    Q .= 0.0
    Q,norm = surface_source(Q,particle,surface_sources,this.cross_sections,this.geometry,this.method)

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

function add_volume_source!(this::Source,source::Array{Float64})
    this.volume_sources = source
end

function get_surface_sources(this::Source)
    return this.surface_sources
end

function get_volume_sources(this::Source)
    return this.volume_sources
end

function get_normalization_factor(this::Source)
    return this.normalization_factor
end

function get_particle(this::Source)
    return this.particle
end 

function Base.:+(source1::Source,source2::Source)
    if source1.get_particle() != source2.get_particle() error("Forbitten addition of different particle sources.") end
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