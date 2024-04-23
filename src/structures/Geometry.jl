
mutable struct Geometry

    # Variable(s)
    name                        ::Union{Missing,String}
    type                        ::Union{Missing,String}
    dimension                   ::Union{Missing,Int64}
    axis                        ::Union{Missing,Vector{String}}
    material_per_region         ::Union{Missing,Array{String}}
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
        this.type = missing
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
Base.propertynames(::Geometry) = 
(
    fieldnames(Geometry)...,
    :set_type,
    :set_dimension,
    :set_material_per_region,
    :set_boundary_conditions,
    :set_number_of_regions,
    :set_voxels_per_region,
    :set_region_boundaries,
    :get_type,
    :get_dimension,
    :get_axis,
    :get_number_of_voxels,
    :get_voxels_width,
    :get_material_per_voxel,
    :get_voxels_position,
    :build
)

function println(this::Geometry)
    entries = ["Name","Type","Dimension","Axis","Material per region","Boundary conditions","Number of regions","Voxels per region","Region boundaries [cm]"]
    values = [this.name,this.type,this.dimension,this.axis,this.material_per_region,this.boundary_conditions,this.number_of_regions,this.voxels_per_region,this.region_boundaries]
    N = length(entries); L = length.(entries); Lmax = maximum(L)
    println("Geometry")
    for n in range(1,N)
        println(string("   ",entries[n]," "^(Lmax-L[n])),"  :  ",values[n])
    end
end

function is_ready_to_build(this::Geometry)
    if ismissing(this.type) error("The type of geometry has to be specified.") end
    if ismissing(this.dimension) error("The geometry dimension has to be specified.") end
    if this.type == "Cartesian"
        if length(size(this.material_per_region)) != this.dimension error("Dimension of the material per region array don't fit the dimension.") end
        if this.dimension ≥ 1
            if ~haskey(this.boundary_conditions,"X+") || ~haskey(this.boundary_conditions,"X-") error("Boundary conditions are not defined along the x-axis.") end
            if ~haskey(this.number_of_regions,"X") error("Number of regions are not defined along the x-axis.") end
            if ~haskey(this.voxels_per_region,"X") error("Number of voxels per region are not defined along the x-axis.") end
            if ~haskey(this.region_boundaries,"X") error("The region boundaries are not defined along the x-axis.") end
            if this.number_of_regions["X"] != size(this.material_per_region,1) error("The size of the material per region array don't fit the number of regions along the x-axis.") end
            if length(this.voxels_per_region["X"]) != this.number_of_regions["X"] error("The length of the voxel per region vector don't fit the number of regions along the x-axis.") end
            if length(this.region_boundaries["X"]) != this.number_of_regions["X"] + 1 error("The length of the region boundaries vector don't fit the number of regions along the x-axis.") end
        end
        if this.dimension ≥ 2
            if ~haskey(this.boundary_conditions,"Y+") || ~haskey(this.boundary_conditions,"Y-") error("Boundary conditions are not defined along the y-axis.") end
            if ~haskey(this.number_of_regions,"Y") error("Number of regions are not defined along the y-axis.") end
            if ~haskey(this.voxels_per_region,"Y") error("Number of voxels per region are not defined along the y-axis.") end
            if ~haskey(this.region_boundaries,"Y") error("The region boundaries are not defined along the y-axis.") end
            if this.number_of_regions["Y"] != size(this.material_per_region,2) error("The size of the material per region array don't fit the number of regions along the y-axis.") end
            if length(this.voxels_per_region["Y"]) != this.number_of_regions["Y"] error("The length of the voxel per region vector don't fit the number of regions along the y-axis.") end
            if length(this.region_boundaries["Y"]) != this.number_of_regions["Y"] + 1 error("The length of the region boundaries vector don't fit the number of regions along the y-axis.") end
        end
        if this.dimension ≥ 3
            if ~haskey(this.boundary_conditions,"Z+") || ~haskey(this.boundary_conditions,"Z-") error("Boundary conditions are not defined along the z-axis.") end
            if ~haskey(this.number_of_regions,"Z") error("Number of regions are not defined along the z-axis.") end
            if ~haskey(this.voxels_per_region,"Z") error("Number of voxels per region are not defined along the z-axis.") end
            if ~haskey(this.region_boundaries,"Z") error("The region boundaries are not defined along the z-axis.") end
            if this.this.number_of_regions["Z"] != size(this.material_per_region,3) error("The size of the material per region array don't fit the number of regions along the z-axis.") end
            if length(this.voxels_per_region["Z"]) != this.number_of_regions["Z"] error("The length of the voxel per region vector don't fit the number of regions along the z-axis.") end
            if length(this.region_boundaries["Z"]) != this.number_of_regions["Z"] + 1 error("The length of the region boundaries vector don't fit the number of regions along the z-axis.") end
        end
    end

end

function build(this::Geometry,cs::Cross_Sections)

    # Verification step
    is_ready_to_build(this)

    # Build geometry parameters
    geometry(this,cs)

end

function set_type(this::Geometry,type::String)
    if lowercase(type) ∉ ["cartesian"] error("Unknown geometry type.") end
    this.type = lowercase(type)
end

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

function set_material_per_region(this::Geometry,material_per_region::Array{String})
    this.material_per_region = material_per_region
end

function set_boundary_conditions(this::Geometry,boundary::String,boundary_condition::String)
    if uppercase(boundary) ∉ ["X-","X+","Y-","Y+","Z-","Z+"] error("Unknown boundary.") end
    if lowercase(boundary_condition) ∉ ["void"] error("Unkown boundary type.") end
    this.boundary_conditions[uppercase(boundary)] = boundary_condition
end

function set_number_of_regions(this::Geometry,axis::String,number_of_regions::Int64)
    if lowercase(axis) ∉ ["x","y","z"] error("Unknown axis.") end
    if number_of_regions ≤ 0 error("Number of regions should be at least one.") end
    this.number_of_regions[lowercase(axis)] = number_of_regions
end

function set_voxels_per_region(this::Geometry,axis::String,voxels_per_region::Vector{Int64})
    if lowercase(axis) ∉ ["x","y","z"] error("Unknown axis.") end
    if length(voxels_per_region) == 0 || any(x -> x ≤ 0, voxels_per_region) error("Number of voxels per region should be at least one.") end
    this.voxels_per_region[lowercase(axis)] = voxels_per_region
    this.number_of_voxels[lowercase(axis)] = sum(voxels_per_region)
end

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