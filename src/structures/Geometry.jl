"""
    Geometry

Structure used to define the geometry properties of the medium for transport calculations.

# Mandatory field(s)
- `name::String` : name (or identifier) of the Geometry structure.
- `dimension::Int64` : dimension of the geometry.
- `material_per_region::Array{Material}` : multidimensional array of the material per regions.
- `boundary_conditions::Dict{String,Int64}` : boundary conditions along each axis.
- `number_of_regions::Dict{String,Int64}` : number of regions along each axis.
- `voxels_per_region::Dict{String,Vector{Int64}}` : number of voxels inside each regions along each axis.
- `region_boundaries::Dict{String,Vector{Float64}}` : boundaries of each regions along each axis.

# Optional field(s) - with default values
- `type::String="cartesian"` : type of geometry.

"""
mutable struct Geometry

    # Variable(s)
    name                        ::Union{Missing,String}
    type                        ::Union{Missing,String}
    dimension                   ::Union{Missing,Int64}
    axis                        ::Union{Missing,Vector{String}}
    material_per_region         ::Union{Missing,Array{Material}}
    boundary_conditions         ::Dict{String,String}
    number_of_regions           ::Dict{String,Int64}
    voxels_per_region           ::Dict{String,Vector{Int64}}
    number_of_voxels            ::Dict{String,Int64}
    voxels_width                ::Dict{String,Vector{Float64}}
    voxels_position             ::Dict{String,Vector{Float64}}
    material_per_voxel          ::Union{Missing,Array{Int64}}
    volume_per_voxel            ::Union{Missing,Array{Float64}}
    region_boundaries           ::Dict{String,Vector{Float64}}
    is_build                    ::Bool

    # Constructor(s)
    function Geometry()

        this = new()

        this.name = missing
        this.type = "cartesian"
        this.dimension = missing
        this.axis = missing
        this.material_per_region = missing
        this.boundary_conditions = Dict{String,String}()
        this.number_of_regions = Dict{String,Int64}()
        this.voxels_per_region = Dict{String,Vector{Int64}}()
        this.region_boundaries = Dict{String,Vector{Float64}}()
        this.is_build = false
        this.number_of_voxels = Dict{String,Int64}()
        this.voxels_width = Dict{String,Vector{Float64}}()
        this.voxels_position = Dict{String,Vector{Float64}}()
        this.material_per_voxel = missing
        this.volume_per_voxel = missing

        return this
    end
end

# Method(s)
function is_ready_to_build(this::Geometry)
    if ismissing(this.type) error("The type of geometry has to be specified.") end
    if ismissing(this.dimension) error("The geometry dimension has to be specified.") end
    if this.type == "cartesian"
        if this.dimension ≥ 1
            if ~haskey(this.boundary_conditions,"X+") || ~haskey(this.boundary_conditions,"X-") error("Boundary conditions are not defined along the x-axis.") end
            if ~haskey(this.number_of_regions,"x") error("Number of regions are not defined along the x-axis.") end
            if ~haskey(this.voxels_per_region,"x") error("Number of voxels per region are not defined along the x-axis.") end
            if ~haskey(this.region_boundaries,"x") error("The region boundaries are not defined along the x-axis.") end
            if length(this.voxels_per_region["x"]) != this.number_of_regions["x"] error("The length of the voxel per region vector don`t fit the number of regions along the x-axis.") end
            if length(this.region_boundaries["x"]) != this.number_of_regions["x"] + 1 error("The length of the region boundaries vector don`t fit the number of regions along the x-axis.") end
        end
        if this.dimension ≥ 2
            if ~haskey(this.boundary_conditions,"Y+") || ~haskey(this.boundary_conditions,"Y-") error("Boundary conditions are not defined along the y-axis.") end
            if ~haskey(this.number_of_regions,"y") error("Number of regions are not defined along the y-axis.") end
            if ~haskey(this.voxels_per_region,"y") error("Number of voxels per region are not defined along the y-axis.") end
            if ~haskey(this.region_boundaries,"y") error("The region boundaries are not defined along the y-axis.") end
            if length(this.voxels_per_region["y"]) != this.number_of_regions["y"] error("The length of the voxel per region vector don`t fit the number of regions along the y-axis.") end
            if length(this.region_boundaries["y"]) != this.number_of_regions["y"] + 1 error("The length of the region boundaries vector don`t fit the number of regions along the y-axis.") end
        end
        if this.dimension ≥ 3
            if ~haskey(this.boundary_conditions,"Z+") || ~haskey(this.boundary_conditions,"Z-") error("Boundary conditions are not defined along the z-axis.") end
            if ~haskey(this.number_of_regions,"z") error("Number of regions are not defined along the z-axis.") end
            if ~haskey(this.voxels_per_region,"z") error("Number of voxels per region are not defined along the z-axis.") end
            if ~haskey(this.region_boundaries,"z") error("The region boundaries are not defined along the z-axis.") end
            if length(this.voxels_per_region["z"]) != this.number_of_regions["z"] error("The length of the voxel per region vector don`t fit the number of regions along the z-axis.") end
            if length(this.region_boundaries["z"]) != this.number_of_regions["z"] + 1 error("The length of the region boundaries vector don`t fit the number of regions along the z-axis.") end
        end
    end

end

"""
    build(this::Geometry,cs::Cross_Sections)

To build the geometry structure.

# Input Argument(s)
- `this::Geometry` : geometry.
- `cs::Cross_Sections` : cross-sections library.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> cs = Cross_Sections()
julia> ... # Defining the cross-sections library properties
julia> cs.build()
julia> geo = Geometry()
julia> ... # Defining the geometry properties
julia> geo.build(cs)
```
"""
function build(this::Geometry,cs::Cross_Sections)

    # Verification step
    is_ready_to_build(this)

    # Build geometry parameters
    geometry(this,cs)

end

"""
    set_type(this::Geometry,type::String)

To set the type of geometry.

# Input Argument(s)
- `this::Geometry` : geometry.
- `type::String` : type of geometry, which can takes the following value:
    -`type = "cartesian"` : Cartesian geometry.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> geo = Geometry()
julia> geo.set_type("cartesian")
```
"""
function set_type(this::Geometry,type::String)
    if lowercase(type) ∉ ["cartesian"] error("Unknown geometry type.") end
    this.type = lowercase(type)
end

"""
    set_dimension(this::Geometry,dimension::Int64)

To set the dimension of geometry.

# Input Argument(s)
- `this::Geometry` : geometry.
- `dimension::Int64` : dimension of the geometry, which can be either 1, 2 or 3.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> geo = Geometry()
julia> geo.set_dimension(2) # For 2D geometry
```
"""
function set_dimension(this::Geometry,dimension::Int64)
    if dimension == 1
        this.axis = ["x"]
    elseif dimension == 2
        this.axis = ["x","y"]
    elseif dimension == 3
        this.axis = ["x","y","z"]
    else
        error("Cartesian geometries are only available in 1D, 2D and 3D dimensions.")
    end 
    this.dimension = dimension
end

"""
    set_material_per_region(this::Geometry,material_per_region::Array{Material})

To set the material in each regions of the geometry.

# Input Argument(s)
- `this::Geometry` : geometry.
- `material_per_region::Array{Material}` : array containing the material for each regions. Its size should fit the number of regions per axis.

# Output Argument(s)
N/A

# Examples
```jldoctest
# Define material
julia> mat1 = Material(); mat2 = Material()
julia> ... # Define the material properties

# 1D geometry case
julia> geo1D = Geometry()
julia> geo1D.set_type("cartesian")
julia> geo1D.set_dimension(1)
julia> geo1D.set_number_of_regions("x",3)
julia> geo1D.set_material_per_region([mat1 mat2 mat1])

# 2D geometry case
julia> geo1D = Geometry()
julia> geo1D.set_type("cartesian")
julia> geo1D.set_dimension(2)
julia> geo1D.set_number_of_regions("x",2)
julia> geo1D.set_number_of_regions("y",3)
julia> geo1D.set_material_per_region([mat1 mat2 mat1 ; mat2 mat1 mat2])

```
"""
function set_material_per_region(this::Geometry,material_per_region::Array{Material})
    this.material_per_region = material_per_region
end

"""
    set_boundary_conditions(this::Geometry,boundary::String,boundary_condition::String)

To set the boundary conditions at the specified boundary.

# Input Argument(s)
- `this::Geometry` : geometry.
- `boundary::String` : boundary for which the boundary condition is applied, which can takes the following value:
    - `boundary = "x-"` : the lower bound along x-axis
    - `boundary = "x+"` : the upper bound along x-axis
    - `boundary = "y-"` : the lower bound along y-axis
    - `boundary = "y+"` : the upper bound along y-axis
    - `boundary = "z-"` : the lower bound along z-axis
    - `boundary = "z+"` : the upper bound along z-axis
- `boundary_condition::String` : boundary conditions, which can takes the following value:
    -`boundary = "void"` : void boundary conditions.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> geo = Geometry()
julia> geo.set_boundary_conditions("x-","void")
```
"""
function set_boundary_conditions(this::Geometry,boundary::String,boundary_condition::String)
    if uppercase(boundary) ∉ ["X-","X+","Y-","Y+","Z-","Z+"] error("Unknown boundary.") end
    if lowercase(boundary_condition) ∉ ["void"] error("Unkown boundary type.") end
    this.boundary_conditions[uppercase(boundary)] = boundary_condition
end

"""
    set_number_of_regions(this::Geometry,axis::String,number_of_regions::Int64)

To set the number of regions along a specified axis.

# Input Argument(s)
- `this::Geometry` : geometry.
- `axis::String` : axis along which the number of regions is specified, which can takes the following values:
    - `boundary = "x"` : along x-axis
    - `boundary = "y"` : along y-axis
    - `boundary = "z"` : along z-axis
- `number_of_regions::Int64` : number of regions.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> geo = Geometry()
julia> geo.set_number_of_regions("x",3)
```
"""
function set_number_of_regions(this::Geometry,axis::String,number_of_regions::Int64)
    if lowercase(axis) ∉ ["x","y","z"] error("Unknown axis.") end
    if number_of_regions ≤ 0 error("Number of regions should be at least one.") end
    this.number_of_regions[lowercase(axis)] = number_of_regions
end

"""
    set_voxels_per_region(this::Geometry,axis::String,voxels_per_region::Vector{Int64})

To set the number of voxels for each regions along a specified axis.

# Input Argument(s)
- `this::Geometry` : geometry.
- `axis::String` : axis along which the number of regions is specified, which can takes the following values:
    - `boundary = "x"` : along x-axis
    - `boundary = "y"` : along y-axis
    - `boundary = "z"` : along z-axis
- `voxels_per_region::Vector{Int64}` : vector with the number of voxels for each regions along the specified axis.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> geo = Geometry()
julia> geo.set_number_of_regions("x",3)
julia> geo.set_voxels_per_region("x",[10,5,2])
```
"""
function set_voxels_per_region(this::Geometry,axis::String,voxels_per_region::Vector{Int64})
    if lowercase(axis) ∉ ["x","y","z"] error("Unknown axis.") end
    if length(voxels_per_region) == 0 || any(x -> x ≤ 0, voxels_per_region) error("Number of voxels per region should be at least one.") end
    this.voxels_per_region[lowercase(axis)] = voxels_per_region
    this.number_of_voxels[lowercase(axis)] = sum(voxels_per_region)
end

"""
    set_region_boundaries(this::Geometry,axis::String,region_boundaries::Vector{Float64})

To set the boundaries of each regions along a specified axis.

# Input Argument(s)
- `this::Geometry` : geometry.
- `axis::String` : axis along which the number of regions is specified, which can takes the following values:
    - `boundary = "x"` : along x-axis
    - `boundary = "y"` : along y-axis
    - `boundary = "z"` : along z-axis
- `region_boundaries::Vector{Float64}` : vector with the regions boundaries along the specified axis in ascending order.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> geo = Geometry()
julia> geo.set_number_of_regions("x",3)
julia> geo.set_voxels_per_region("x",[10,5,2])
julia> geo.set_voxels_per_region("x",[0.0 0.3 0.5 1.0])
```
"""
function set_region_boundaries(this::Geometry,axis::String,region_boundaries::Vector{Float64})
    if lowercase(axis) ∉ ["x","y","z"] error("Unknown axis.") end
    if length(region_boundaries) ≤ 1 error("At least two boundaries are required.") end
    for i in range(1,length(region_boundaries)-1)
        if region_boundaries[i] ≥ region_boundaries[i+1] error("Boundaries should be provided in ascending order.") end
    end
    this.region_boundaries[lowercase(axis)] = region_boundaries
end

function get_type(this::Geometry)
    if ismissing(this.type) error("Unable to get geometry type. Missing data.") end
    return this.type
end

function get_dimension(this::Geometry)
    if ismissing(this.dimension) error("Unable to get geometry dimension. Missing data.") end
    return this.dimension
end

function get_axis(this::Geometry)
    if ismissing(this.axis) error("Unable to get geometry axis system. Missing data.") end
    return this.axis
end

function get_number_of_voxels(this::Geometry)
    N = this.get_dimension()
    axis = this.get_axis()
    Ns = ones(Int64,3)
    for n in range(1,N)
        if ~haskey(this.number_of_voxels,axis[n]) error("Unable to get the number of voxels along the ",axis[n],"-axis.") end
        Ns[n] = this.number_of_voxels[axis[n]]
    end
    return Ns
end

function get_voxels_width(this::Geometry)
    N = this.get_dimension()
    axis = this.get_axis()
    Δs = Vector{Vector{Float64}}(undef,3)
    for n in range(1,N)
        if ~haskey(this.voxels_width,axis[n]) error("Unable to get the number of voxels along the ",axis[n],"-axis.") end
        Δs[n] = this.voxels_width[axis[n]]
    end
    return Δs
end

function get_material_per_voxel(this::Geometry)
    if ismissing(this.material_per_voxel) error("Unable to get material per voxel array. Missing data.") end
    return this.material_per_voxel
end

function get_voxels_position(this::Geometry,axis::String)
    if ismissing(this.voxels_position) error("Unable to get voxel positions. Missing data.") end
    if axis ∉ this.get_axis() error("Unknown axis.") end
    return this.voxels_position[axis]
end