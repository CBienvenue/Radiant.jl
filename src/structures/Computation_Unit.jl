
mutable struct Computation_Unit

    # Variable(s)
    name                     ::Union{Missing,String}
    cross_sections           ::Union{Missing,Cross_Sections}
    geometry                 ::Union{Missing,Geometry}
    methods                  ::Union{Missing,Methods}
    sources                  ::Union{Missing,Fixed_Sources}
    flux                     ::Union{Missing,Flux}

    # Constructor(s)
    function Computation_Unit()

        this = new()

        this.name = missing
        this.cross_sections = missing
        this.geometry = missing
        this.methods = missing
        this.sources = missing
        this.flux = missing

        return this
    end
end

# Method(s)
Base.propertynames(::Computation_Unit) = 
(
    fieldnames(Computation_Unit)...,
    :set_cross_sections,
    :set_geometry,
    :set_methods,
    :set_sources,
    :run,
    :get_flux,
    :get_energy_deposition,
    :get_charge_deposition,
    :get_voxels_position,
    :get_energies
)

function println(this::Computation_Unit)
    entries = ["Name","Cross-Sections","Geometry","Methods","Sources"]
    values = [this.name,this.cross_sections.name,this.geometry.name,this.methods.name,this.sources.name]
    N = length(entries); L = length.(entries); Lmax = maximum(L)
    println("Computation_Unit")
    for n in range(1,N)
        println(string("   ",entries[n]," "^(Lmax-L[n])),"  :  ",values[n])
    end
end

function set_cross_sections(this::Computation_Unit,cross_sections::Cross_Sections)
    this.cross_sections = cross_sections
end

function set_geometry(this::Computation_Unit,geometry::Geometry)
    this.geometry = geometry
end

function set_methods(this::Computation_Unit,methods::Methods)
    this.methods = methods
end

function set_sources(this::Computation_Unit,sources::Fixed_Sources)
    this.sources = sources
end

function run(this::Computation_Unit,is_CUDA::Bool=false)
    this.flux = transport(this.cross_sections,this.geometry,this.methods,this.sources,is_CUDA)
end

function get_energy_deposition(this::Computation_Unit,type::String,which_generations::Int64=0)
    if ismissing(this.flux) error("No computed flux in this computation unit. To extract energy deposition, please use .run() method before.") end
    if type ∉ ["total","electrons","photons","positrons"] error("Unknown type of energy deposition.") end
    if type ∈ ["electrons","photons","positrons"] && type ∉ this.flux.get_particles() error("Energy deposition for the specified particle is not available.") end
    return energy_deposition(this.cross_sections,this.geometry,this.methods,this.sources,this.flux,type)
end

function get_charge_deposition(this::Computation_Unit,type::String,which_generations::Int64=0)
    if ismissing(this.flux) error("No computed flux in this computation unit. To extract charge deposition, please use .run() method before.") end
    if type ∉ ["total","electrons","photons","positrons"] error("Unknown type of charge deposition.") end
    if type ∈ ["electrons","photons","positrons"] && type ∉ this.flux.get_particles() error("Charge deposition for the specified particle is not available.") end
    return charge_deposition(this.cross_sections,this.geometry,this.methods,this.sources,this.flux,type)
end

function get_flux(this::Computation_Unit,type::String)
    if ismissing(this.flux) error("No computed flux in this computation unit. To extract flux, please use .run() method before.") end
    if type ∉ ["electrons","photons","positrons"] error("Unknown type of flux.") end
    if type ∈ ["electrons","photons","positrons"] && type ∉ this.flux.get_particles() error("Flux for the specified particle is not available.") end
    return flux(this.cross_sections,this.geometry,this.flux,type)
end

function get_voxels_position(this::Computation_Unit,axis::String)
    return this.geometry.get_voxels_position(axis)
end

function get_energies(this::Computation_Unit,particle::String)
    return this.cross_sections.get_energies(particle)
end