"""
    Discrete_Ordinates

Structure used to define the discretization method associated with the transport of a particle.

# Mandatory field(s)
- `name::String`: name (or identifier) of the Discrete_Ordinates structure.
- `particle::String`: particle for which the discretization methods is defined
- `solver_type::String`: type of solver for the transport calculations.
- `quadrature_type::String`: type of quadrature for the angular domain.
- `quadrature_order::Int64`: order of the quadrature for the angular domain.
- `legendre_order::Int64`: maximum order of the Legendre expansion for the differential cross-sections.
- `scheme_type::Dict{String,String}`: type of schemes for the spatial or energy discretization.
- `scheme_order::Dict{String,Int64}`: order of the expansion for the discretization schemes.

# Optional field(s) - with default values
- `angular_fokker_planck::String="finite-difference"`: type of discretization for the angular Fokker-Planck operation.
- `angular_boltzmann::String="galerkin-d"`: type of discretization for the Boltzmann operation.
- `convergence_criterion::Float64 = 1e-7`: convergence criterion of in-group iterations.
- `maximum_iteration::Int64 = 300`: maximum number of in-group iterations.
- `acceleration::Int64 = "none"`: acceleration method for the in-group iterations.

"""
mutable struct Discrete_Ordinates

    # Variable(s)
    name                       ::Union{Missing,String}
    particle                   ::Union{Missing,Particle}
    solver_type                ::Union{Missing,String}
    quadrature_type            ::Union{Missing,String}
    quadrature_order           ::Union{Missing,Int64}
    legendre_order             ::Union{Missing,Int64}
    angular_fokker_planck      ::Union{Missing,String}
    angular_boltzmann          ::Union{Missing,String}
    convergence_criterion      ::Float64
    maximum_iteration          ::Int64
    scheme_type                ::Dict{String,String}
    scheme_order               ::Dict{String,Int64}
    acceleration               ::String

    # Constructor(s)
    function Discrete_Ordinates()

        this = new()

        this.name = missing
        this.particle = missing
        this.solver_type = missing
        this.quadrature_type = missing
        this.quadrature_order = missing
        this.legendre_order = missing
        this.angular_fokker_planck = "finite-difference"
        this.angular_boltzmann = "galerkin-d"
        this.convergence_criterion = 1e-7 
        this.maximum_iteration = 300
        this.scheme_type = Dict{String,String}()
        this.scheme_order = Dict{String,Int64}()
        this.acceleration = "none"

        return this
    end
end

# Method(s)
"""
    set_particle(this::Discrete_Ordinates,particle::String)

To set the particle for which the transport discretization method is for. 

# Input Argument(s)
- `this::Discrete_Ordinates`: discretization method.
- `particle::String`: particle identifier, which is either:
    - `particle = "photons"`: photons.
    - `particle = "electrons"`: electrons.
    - `particle = "positrons"`: positrons.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = Discrete_Ordinates()
julia> m.set_particle("electrons") 
```
"""
function set_particle(this::Discrete_Ordinates,particle::Particle)
    this.particle = particle
end

"""
    set_solver_type(this::Discrete_Ordinates,solver_type::String)

To set the solver for the particle transport.

# Input Argument(s)
- `this::Discrete_Ordinates`: discretization method.
- `solver_type::String`: solver type, which can takes the following values:
    - `solver_type = "BTE"`: Boltzmann transport equation
    - `solver_type = "BFP"`: Boltzmann Fokker-Planck equation
    - `solver_type = "BCSD"`: Boltzmann-CSD equation
    - `solver_type = "FP"`: Fokker-Planck equation
    - `solver_type = "CSD"`: Continuous slowing-down only equation
    - `solver_type = "BFP-EF"`: Boltzmann Fokker-Planck without elastic

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = Discrete_Ordinates()
julia> m.set_particle("BFP")
```
"""
function set_solver_type(this::Discrete_Ordinates,solver_type::String)
    if uppercase(solver_type) ∉ ["BTE","BFP","BCSD","FP","CSD","BFP-EF"] error("Unkown type of solver.") end
    this.solver_type = uppercase(solver_type)
end

"""
    set_quadrature(this::Discrete_Ordinates,type::String,order::Int64)

To set the quadrature properties for the discretization of the angular domain.

# Input Argument(s)
- `this::Discrete_Ordinates`: discretization method.
- `type::String`: type of quadrature, which can takes the following values:
    - `type = "gauss-legendre"`: Gauss-Legendre quadrature (1D Cartesian geometry only)
    - `type = "gauss-lobatto"`: Gauss-Lobatto quadrature (1D Cartesian geometry only)
    - `type = "carlson"`: Carlson quadrature (2D or 3D Cartesian geometry only)
    - `type = "gauss-legendre-chebychev"`: product quadrature between Gauss-Legendre and Chebychev quadratures (2D or 3D Cartesian geometry only)
    - `type = "lebedev"`: Lebedev quadrature (2D or 3D Cartesian geometry only)
- `order::Int64`: order of the quadrature, which is any integer greater than 2.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = Discrete_Ordinates()
julia> m.set_quadrature("gauss-legendre",4)
```
"""
function set_quadrature(this::Discrete_Ordinates,type::String,order::Int64)
    if lowercase(type) ∉ ["gauss-legendre","gauss-lobatto","carlson","lebedev","gauss-legendre-chebychev"] error("Unknown quadrature type.") end
    if order ≤ 1 error("Quadrature order should be at least 2.") end
    this.quadrature_type = lowercase(type)
    this.quadrature_order = order
end

"""
    set_legendre_order(this::Discrete_Ordinates,legendre_order::Int64)

To set the maximum order of the Legendre expansion of the differential cross-sections.

# Input Argument(s)
- `this::Discrete_Ordinates`: discretization method.
- `legendre_order::Int64`: maximum order of the Legendre expansion of the differential cross-sections.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = Discrete_Ordinates()
julia> m.set_legendre_order(7)
```
"""
function set_legendre_order(this::Discrete_Ordinates,legendre_order::Int64)
    if legendre_order < 0 error("Legendre order should be at least 0.") end
    this.legendre_order = legendre_order
end

"""
    set_angular_fokker_planck(this::Discrete_Ordinates,angular_fokker_planck::String)

To set the discretization method for the angular Fokker-Planck term.

# Input Argument(s)
- `this::Discrete_Ordinates`: discretization method.
- `angular_fokker_planck::String`: discretization method for the angular Fokker-Planck term, which can takes the following values:
    - `angular_fokker_planck = "finite-difference"`: finite difference discretization.
    - `angular_fokker_planck = "galerkin"`: galerkin moment-based discretization.
    - `angular_fokker_planck = "differential-quadrature"`: finite difference discretization (1D Cartesian geometry only).
    
# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = Discrete_Ordinates()
julia> m.set_angular_fokker_planck("differential-quadrature")
```
"""
function set_angular_fokker_planck(this::Discrete_Ordinates,angular_fokker_planck::String)
    if angular_fokker_planck ∉ ["finite-difference","differential-quadrature","galerkin"] error("Unkown method to deal with the angular Fokker-Planck term.") end
    this.angular_fokker_planck = angular_fokker_planck
end

"""
    set_angular_boltzmann(this::Discrete_Ordinates,angular_boltzmann::String)

To set the angular discretization method for the Boltzmann operator.

# Input Argument(s)
- `this::Discrete_Ordinates`: discretization method.
- `angular_boltzmann::String`: angular discretization method for the Boltzmann operator, which can takes the following values:
    - `angular_boltzmann = "standard"`: standard discrete ordinates (SN) method.
    - `angular_boltzmann = "galerkin-m"`: Galerkin method by inversion of the discrete-to-moment M matrix.
    - `angular_boltzmann = "galerkin-d"`: Galerkin method by inversion of the moment-to-discrete D matrix.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = Discrete_Ordinates()
julia> m.set_angular_boltzmann("standard")
```
"""
function set_angular_boltzmann(this::Discrete_Ordinates,angular_boltzmann::String)
    if angular_boltzmann ∉ ["standard","galerkin-m","galerkin-d","galerkin"] error("Unkown method to deal with the Boltzmann kernel.") end
    if (angular_boltzmann == "galerkin") angular_boltzmann = "galerkin-d" end
    this.angular_boltzmann = angular_boltzmann
end

"""
    set_convergence_criterion(this::Discrete_Ordinates,convergence_criterion::Float64)

To set the convergence criterion for the in-group iterations.

# Input Argument(s)
- `this::Discrete_Ordinates`: discretization method.
- `convergence_criterion::Float64`: convergence criterion.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = Discrete_Ordinates()
julia> m.set_convergence_criterion(1e-5)
```
"""
function set_convergence_criterion(this::Discrete_Ordinates,convergence_criterion::Float64)
    if convergence_criterion ≤ 0 error("Convergence criterion has to be greater than 0.") end
    this.convergence_criterion = convergence_criterion
end

"""
    set_maximum_iteration(this::Discrete_Ordinates,maximum_iteration::Int64)

To set the maximum number of in-group iterations.

# Input Argument(s)
- `this::Discrete_Ordinates`: discretization method.
- `maximum_iteration::Int64`: maximum number of in-group iterations.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = Discrete_Ordinates()
julia> m.set_maximum_iteration(50)
```
"""
function set_maximum_iteration(this::Discrete_Ordinates,maximum_iteration::Int64)
    if maximum_iteration < 1 error("Maximum iteration has to be at least 1.") end
    this.maximum_iteration = maximum_iteration
end

"""
    set_scheme(this::Discrete_Ordinates,axis::String,scheme_type::String,scheme_order::Int64)

To set the type of discretization scheme for derivative along the specified spatial or energy axis.

# Input Argument(s)
- `this::Discrete_Ordinates`: discretization method.
- `axis::String`: variable of the derivative for which the scheme is applied, which can takes the following values:
    - `axis = "x"`: spatial x axis (discretization of the streaming term)
    - `axis = "y"`: spatial y axis (discretization of the streaming term)
    - `axis = "z"`: spatial z axis (discretization of the streaming term)
    - `axis = "E"`: spatial E axis (discretization of the continuous slowing-down term)
- `scheme_type::String`: type of scheme to be applied, which can takes the following values:
    - `scheme_type = "DD"`: diamond difference scheme (any order)
    - `scheme_type = "DG"`: discontinuous Galerkin scheme (any order)
    - `scheme_type = "AWD"`: adaptive weighted scheme (1st and 2nd order only)
- `scheme_order::Int64`: scheme order, which takes values greater than 1.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = Discrete_Ordinates()
julia> m.set_scheme("x","DD",1)
julia> m.set_scheme("E","DG",2)
```
"""
function set_scheme(this::Discrete_Ordinates,axis::String,scheme_type::String,scheme_order::Int64)
    if axis ∉ ["x","y","z","E"] error("Unknown axis.") end
    if uppercase(scheme_type) ∉ ["DD","DG","DG-","DG+","AWD"] error("Unknown type of scheme.") end
    if scheme_order ≤ 0 error("Scheme order should be at least of 1.") end
    this.scheme_type[axis] = scheme_type
    this.scheme_order[axis] = scheme_order
end

"""
    set_acceleration(this::Discrete_Ordinates,acceleration::String)

To set the acceleration method for the in-group iteration process.

# Input Argument(s)
- `this::Discrete_Ordinates`: discretization method.
- `acceleration::String`: acceleration method, which takes the following values
    - `acceleration = "none"`: none
    - `acceleration = "livolant"`: livolant acceleration method
# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = Discrete_Ordinates()
julia> m.set_acceleration("livolant")
```
"""
function set_acceleration(this::Discrete_Ordinates,acceleration::String)
    if lowercase(acceleration) ∉ ["none","livolant"] error("Unkown acceleration method.") end
    this.acceleration = acceleration
end

function get_legendre_order(this::Discrete_Ordinates)
    if ismissing(this.legendre_order) error("Unable to get Legendre order. Missing data.") end
    return this.legendre_order
end

function get_quadrature_order(this::Discrete_Ordinates)
    if ismissing(this.quadrature_order) error("Unable to get quadrature order. Missing data.") end
    return this.quadrature_order
end

function get_quadrature_type(this::Discrete_Ordinates)
    if ismissing(this.quadrature_type) error("Unable to get quadrature type. Missing data.") end
    return this.quadrature_type
end

function get_angular_boltzmann(this::Discrete_Ordinates)
    if ismissing(this.angular_boltzmann) error("Unable to get angular Boltzmann treatment type. Missing data.") end
    return this.angular_boltzmann
end

function get_angular_fokker_planck(this::Discrete_Ordinates)
    if ismissing(this.angular_fokker_planck) error("Unable to get angular Fokker-Planck treatment type. Missing data.") end
    return this.angular_fokker_planck
end

function get_particle(this::Discrete_Ordinates)
    if ismissing(this.particle) error("Unable to get particle. Missing data.") end
    return this.particle
end

function get_solver_type(this::Discrete_Ordinates)
    if ismissing(this.solver_type) error("Unable to get solver type. Missing data.") end
    if this.solver_type == "BTE"
        isCSD = false
        solver = 1
    elseif this.solver_type == "BFP"
        isCSD = true
        solver = 2
    elseif this.solver_type == "BCSD"
        isCSD = true
        solver = 3
    elseif this.solver_type == "FP"
        isCSD = true
        solver = 4
    elseif this.solver_type == "CSD"
        isCSD = true
        solver = 5
    elseif this.solver_type == "BFP-EF"
        isCSD = true
        solver = 6
    else
        error("Unknown type of solver.")
    end
    return solver, isCSD
end

function get_schemes(this::Discrete_Ordinates,geometry::Geometry,is_full_coupling::Bool)
    schemes = Vector{String}(undef,4)
    𝒪 = Vector{Int64}(undef,4)
    axis = geometry.get_axis()
    _, isCSD = this.get_solver_type()
    n = 1
    for x in ["x","y","z","E"]
        if x ∈ axis || (isCSD && x ∈ ["E"])
            if ~haskey(this.scheme_type,x) error("Scheme type is not defined along ",x,"-axis.") end
            if ~haskey(this.scheme_order,x) error("Scheme order is not defined along ",x,"-axis.") end
            schemes[n] = this.scheme_type[x]
            𝒪[n] = this.scheme_order[x]
        else
            schemes[n] = "DD"
            𝒪[n] = 1
        end
        n += 1
    end
    if is_full_coupling
        Nm = [𝒪[2]*𝒪[3]*𝒪[4],𝒪[1]*𝒪[3]*𝒪[4],𝒪[1]*𝒪[2]*𝒪[4],𝒪[1]*𝒪[2]*𝒪[3],prod(𝒪)]
    else
        Nm = [𝒪[2]+𝒪[3]+𝒪[4]+1,𝒪[1]+𝒪[3]+𝒪[4]+1,𝒪[1]+𝒪[2]+𝒪[4]+1,𝒪[1]+𝒪[2]+𝒪[3]+1,sum(𝒪)+1]
    end
    return schemes,𝒪,Nm
end

function get_convergence_criterion(this::Discrete_Ordinates)
    return this.convergence_criterion
end

function get_maximum_iteration(this::Discrete_Ordinates)
    return this.maximum_iteration
end

function get_acceleration(this::Discrete_Ordinates)
    return this.acceleration
end

function get_quadrature_dimension(this::Discrete_Ordinates,Ndims::Int64)
    if this.quadrature_type ∈ ["gauss-legendre","gauss-lobatto"]
        return 1
    elseif this.quadrature_type ∈ ["gauss-legendre-chebychev","lebedev","carlson"]
        if Ndims == 1
            return 3
        else
            return Ndims
        end
    else
        error("Unkown quadrature type.")
    end
end