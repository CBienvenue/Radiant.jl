"""
    Computation_Unit

Structure used to consolidate the cross-sections, geometry, solvers and sources, execute transport calculations and extract its results.

# User-defined field(s)

- **Mandatory field(s)**
    - `cross_sections::Cross_Sections`: cross-section library.
    - `geometry::Geometry`: geometry.
    - `coupled_transport_solvers::Coupled_Transport_Solvers`: solvers.
    - `sources::Sources`: fixed sources.

- **Optional field(s) - with default values**
N/A

# System-defined field(s)
- `flux::Flux`: flux solution.

"""
mutable struct Computation_Unit

    # Variable(s)
    cross_sections            ::Union{Missing,Cross_Sections}
    geometry                  ::Union{Missing,Geometry}
    coupled_transport_solvers ::Union{Missing,Coupled_Transport_Solvers}
    sources                   ::Union{Missing,Fixed_Sources}
    flux                      ::Union{Missing,Flux}

    # Constructor(s)
    function Computation_Unit()

        this = new()

        this.cross_sections = missing
        this.geometry = missing
        this.coupled_transport_solvers = missing
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

"""
    set_cross_sections(this::Computation_Unit,cross_sections::Cross_Sections)

To set the cross-sections library for transport calculations.

# Input Argument(s)
- `this::Computation_Unit`: computation unit.
- `cross_sections::Cross_Sections`: cross-sections library.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cs = Cross_Sections()
julia> ... # Define cross-sections properties and generate multigroup cross-sections.
julia> cu = Computation_Unit()
julia> cu.set_cross_sections(cs)
```
"""
function set_cross_sections(this::Computation_Unit,cross_sections::Cross_Sections)
    this.cross_sections = cross_sections
end

"""
    set_geometry(this::Computation_Unit,geometry::Geometry)

To set the geometry for transport calculations.

# Input Argument(s)
- `this::Computation_Unit`: computation unit.
- `geometry::Geometry`: geometry.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> geo = Geometry()
julia> ... # Define geometry and its properties
julia> cu = Computation_Unit()
julia> cu.set_geometry(geo)
```
"""
function set_geometry(this::Computation_Unit,geometry::Geometry)
    this.geometry = geometry
end

"""
    set_methods(this::Computation_Unit,coupled_transport_solvers::Coupled_Transport_Solvers)

To set the discretization coupled_transport_solvers for transport calculations.

# Input Argument(s)
- `this::Computation_Unit`: computation unit.
- `coupled_transport_solvers::Coupled_Transport_Solvers`: collection of discretization method per particle.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ms = Coupled_Transport_Solvers()
julia> ... # Define all the discretization coupled_transport_solvers and their properties
julia> cu = Computation_Unit()
julia> cu.set_methods(ms)
```
"""
function set_methods(this::Computation_Unit,coupled_transport_solvers::Coupled_Transport_Solvers)
    this.coupled_transport_solvers = coupled_transport_solvers
end

"""
    set_sources(this::Computation_Unit,sources::Fixed_Sources)

To set the fixed sources for transport calculations.

# Input Argument(s)
- `this::Computation_Unit`: computation unit.
- `sources::Fixed_Sources`: collection of fixed sources.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> fs = Fixed_Sources()
julia> ... # Define all the fixed sources and their properties
julia> cu = Computation_Unit()
julia> cu.set_sources(fs)
```
"""
function set_sources(this::Computation_Unit,sources::Fixed_Sources)
    this.sources = sources
end

"""
    run(this::Computation_Unit)

To lauch the transport calculations to solve the transport equation and obtain the flux solution.

# Input Argument(s)
- `this::Computation_Unit`: computation unit.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cu = Computation_Unit()
julia> ... # Define the cross-sections, geometry, fixed sources and discretization coupled_transport_solvers
julia> cu.run()
```
"""
function run(this::Computation_Unit)
    is_CUDA = false
    this.flux = transport(this.cross_sections,this.geometry,this.coupled_transport_solvers,this.sources,is_CUDA)
end

function get_energy_deposition(this::Computation_Unit,type::String,which_generations::Int64=0)
    if ismissing(this.flux) error("No computed flux in this computation unit. To extract energy deposition, please use .run() method before.") end
    if type ∉ ["total","electrons","photons","positrons"] error("Unknown type of energy deposition.") end
    if type ∈ ["electrons","photons","positrons"] && type ∉ this.flux.get_particles() error("Energy deposition for the specified particle is not available.") end
    return energy_deposition(this.cross_sections,this.geometry,this.coupled_transport_solvers,this.sources,this.flux,type)
end

function get_charge_deposition(this::Computation_Unit,type::String,which_generations::Int64=0)
    if ismissing(this.flux) error("No computed flux in this computation unit. To extract charge deposition, please use .run() method before.") end
    if type ∉ ["total","electrons","photons","positrons"] error("Unknown type of charge deposition.") end
    if type ∈ ["electrons","photons","positrons"] && type ∉ this.flux.get_particles() error("Charge deposition for the specified particle is not available.") end
    return charge_deposition(this.cross_sections,this.geometry,this.coupled_transport_solvers,this.sources,this.flux,type)
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