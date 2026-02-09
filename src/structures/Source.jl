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
    solver         ::Solver

    # Constructor(s)
    function Source(particle::Particle,cross_sections::Cross_Sections,geometry::Geometry,solver::Solver)
        this = new()
        this.name = missing
        this.particle = particle
        this.normalization_factor = 0
        this.cross_sections = cross_sections
        this.geometry = geometry
        this.solver = solver
        initalize_sources(this,cross_sections,geometry,solver)
        return this
    end
end

# Method(s)
"""
    initalize_sources(this::Source,cross_sections::Cross_Sections,geometry::Geometry,solver::Solver)

Initialize volume and boundary sources.

# Input Argument(s)
- `this::Source` : source structure.
- `cross_sections::Cross_Sections` : cross-sections library.
- `geometry::Geometry` : geometry.
- `solver::Solver` : discrete ordinates solver.

# Output Argument(s)
N/A

"""
function initalize_sources(this::Source,cross_sections::Cross_Sections,geometry::Geometry,solver::Solver)

    # Data extraction and validation
    particle = this.particle
    if get_id(particle) ∉ get_id.(cross_sections.particles) error(string("No cross sections available for ",particle," particle.")) end
    index = findfirst(x -> get_id(x) == get_id(particle),cross_sections.particles)
    Ng = cross_sections.number_of_groups[index]
    Nx = geometry.number_of_voxels["x"]
    Ndims = geometry.dimension
    if Ndims ≥ 2 Ny = geometry.number_of_voxels["y"] else Ny = 1 end
    if Ndims ≥ 3 Nz = geometry.number_of_voxels["z"] else Nz = 1 end
    _,_,Nm = solver.get_schemes(geometry,solver.get_is_full_coupling())

    # Angular discretization
    if solver isa Discrete_Ordinates
        Qdims = solver.get_quadrature_dimension(Ndims)
        Ω,w = quadrature(solver.quadrature_order,solver.quadrature_type,Qdims)
        if typeof(Ω) == Vector{Float64} Ω = [Ω,0*Ω,0*Ω] end
        P,_,_,_ = angular_polynomial_basis(Ω,w,solver.get_legendre_order(),solver.get_angular_boltzmann(),Qdims)
    elseif solver isa Spherical_Harmonics
        polynomial_basis = solver.get_polynomial_basis(Ndims)
        L = solver.get_legendre_order()
        if polynomial_basis == "legendre"
            P = L+1
        else
            P = (L+1)^2
        end
    elseif solver isa Galerkin
        polynomial_basis = solver.get_polynomial_basis(Ndims)
        Lp,Lq = solver.get_legendre_order()
        if polynomial_basis == "legendre"
            error("Not available yet for Galerkin with Legendre basis.")
        else
            P = (Lp+1)^2
        end
    else
        error("No methods available for $(get_type(particle)) particle.")
    end

    # Initialize sources
    L = 0 
    this.volume_sources = zeros(Ng,P,Nm[5],Nx,Ny,Nz)
    this.surface_sources = Array{Union{Array{Float64},Float64}}(undef,Ng,L+1,2*Ndims)
    for ig in range(1,Ng), l in range(0,L)
        if Ndims == 1
            for i in range(1,2) this.surface_sources[ig,l+1,i] = 0.0 end
        elseif Ndims == 2
            for i in range(1,2) this.surface_sources[ig,l+1,i] = zeros(Ny) end
            for i in range(3,4) this.surface_sources[ig,l+1,i] = zeros(Nx) end
        elseif Ndims == 3
            for i in range(1,2) this.surface_sources[ig,l+1,i] = zeros(Ny,Nz) end
            for i in range(3,4) this.surface_sources[ig,l+1,i] = zeros(Nx,Nz) end
            for i in range(5,6) this.surface_sources[ig,l+1,i] = zeros(Nx,Ny) end
        end
    end
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

    # Data extraction and validation
    particle = surface_sources.particle
    if get_id(particle) ∉ get_id.(this.cross_sections.particles) error(string("No cross sections available for ",get_type(particle)," particle.")) end
    Ndims = this.geometry.dimension
    Nx = this.geometry.number_of_voxels["x"]
    if Ndims ≥ 2 Ny = this.geometry.number_of_voxels["y"] else Ny = 1 end
    if Ndims ≥ 3 Nz = this.geometry.number_of_voxels["z"] else Nz = 1 end
    if get_id(particle) != get_id(this.solver.particle) error(string("No methods available for ",get_type(particle)," particle.")) end

    # Compute and format the surface source for transport solver
    Q_old = this.surface_sources
    Q_new,norm = surface_source(particle,surface_sources,this.cross_sections,this.geometry,this.solver)

    # Add it to the source object
    dims_new = size(Q_new)
    dims_old = size(Q_old)
    if (dims_new[1] != dims_old[1]) || (dims_new[3] != dims_old[3]) error("Number of groups and/or dimension of sources are not coherent.") end
    Ng = dims_old[1]
    L_new = dims_new[2] - 1
    L_old = dims_old[2] - 1
    if L_new > L_old
        Q_block = Array{Union{Array{Float64},Float64}}(undef,Ng,L_new-L_old,2*Ndims)
        for ig in range(1,Ng), l in range(0,L_new-L_old-1)
            if Ndims == 1
                for i in range(1,2) Q_block[ig,l+1,i] = 0.0 end
            elseif Ndims == 2
                for i in range(1,2) Q_block[ig,l+1,i] = zeros(Ny) end
                for i in range(3,4) Q_block[ig,l+1,i] = zeros(Nx) end
            elseif Ndims == 3
                for i in range(1,2) Q_block[ig,l+1,i] = zeros(Ny,Nz) end
                for i in range(3,4) Q_block[ig,l+1,i] = zeros(Nx,Nz) end
                for i in range(5,6) Q_block[ig,l+1,i] = zeros(Nx,Ny) end
            end
        end
        Q_old = cat(Q_old,Q_block; dims=2)
    end
    for ig in range(1,Ng), l in range(0,L_new), i in range(1,2*Ndims)
        Q_old[ig,l+1,i] += Q_new[ig,l+1,i]
    end
    this.surface_sources = Q_old
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
    Ndims = source1.geometry.dimension
    Nx = source1.geometry.number_of_voxels["x"]
    if Ndims ≥ 2 Ny = source1.geometry.number_of_voxels["y"] else Ny = 1 end
    if Ndims ≥ 3 Nz = source1.geometry.number_of_voxels["z"] else Nz = 1 end

    # Add it to the source object
    Q_1 = source1.surface_sources
    Q_2 = source2.surface_sources
    dims_1 = size(Q_1)
    dims_2 = size(Q_2)
    if (dims_1[1] != dims_2[1]) || (dims_1[3] != dims_2[3]) error("Number of groups and/or dimension of sources are not coherent.") end
    Ng = dims_1[1]
    Ndims = div(dims_1[3],2)
    L_1 = dims_1[2] - 1
    L_2 = dims_2[2] - 1
    if L_2 > L_1
        Q_block = Array{Union{Array{Float64},Float64}}(undef,Ng,L_2-L_1,2*Ndims)
        for ig in range(1,Ng), l in range(0,L_2-L_1-1)
            if Ndims == 1
                for i in range(1,2) Q_block[ig,l+1,i] = 0.0 end
            elseif Ndims == 2
                for i in range(1,2) Q_block[ig,l+1,i] = zeros(Ny) end
                for i in range(3,4) Q_block[ig,l+1,i] = zeros(Nx) end
            elseif Ndims == 3
                for i in range(1,2) Q_block[ig,l+1,i] = zeros(Ny,Nz) end
                for i in range(3,4) Q_block[ig,l+1,i] = zeros(Nx,Nz) end
                for i in range(5,6) Q_block[ig,l+1,i] = zeros(Nx,Ny) end
            end
        end
        Q_1 = cat(Q_1,Q_block; dims=2)
    end
    for ig in range(1,Ng), l in range(0,L_2), i in range(1,2*Ndims)
        Q_1[ig,l+1,i] += Q_2[ig,l+1,i]
    end
    source1.surface_sources = Q_1
    return source1
end