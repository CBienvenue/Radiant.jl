"""
    SN

Structure used to define the discretization method associated with the transport of a particle.

# Mandatory field(s)
- `particle::Particle` : particle for which the discretization methods is defined
- `solver_type::String` : type of solver for the transport calculations.
- `quadrature_type::String` : type of quadrature for the angular domain.
- `quadrature_order::Int64` : order of the quadrature for the angular domain.
- `legendre_order::Int64` : maximum order of the Legendre expansion for the differential cross-sections.
- `scheme_type::Dict{String,String}` : type of schemes for the spatial or energy discretization.
- `scheme_order::Dict{String,Int64}` : order of the expansion for the discretization schemes.

# Optional field(s) - with default values
- `angular_fokker_planck::String = "finite-difference"` : type of discretization for the angular Fokker-Planck operation.
- `angular_boltzmann::String = "galerkin-d"` : type of discretization for the Boltzmann operation.
- `convergence_criterion::Float64 = 1e-7` : convergence criterion of in-group iterations.
- `maximum_iteration::Int64 = 300` : maximum number of in-group iterations.
- `acceleration::Int64 = "none"` : acceleration method for the in-group iterations.

"""
mutable struct SN

    # Variable(s)
    particle                   ::Union{Missing,Particle}
    solver_type                ::Union{Missing,String}
    quadrature_type            ::Union{Missing,String}
    quadrature_order           ::Union{Missing,Int64}
    quadrature_dimension       ::Int64
    legendre_order             ::Union{Missing,Int64}
    angular_fokker_planck      ::Union{Missing,String}
    angular_boltzmann          ::Union{Missing,String}
    convergence_criterion      ::Float64
    maximum_iteration          ::Int64
    scheme_type                ::Dict{String,String}
    scheme_order               ::Dict{String,Int64}
    acceleration               ::String
    isFC                       ::Bool

    # Constructor(s)
    function SN()
        this = new()
        this.particle = missing
        this.solver_type = missing
        this.quadrature_type = missing
        this.quadrature_order = missing
        this.quadrature_dimension = 0
        this.legendre_order = 64
        this.angular_fokker_planck = "finite-difference"
        this.angular_boltzmann = "galerkin-d"
        this.convergence_criterion = 1e-7 
        this.maximum_iteration = 300
        this.scheme_type = Dict{String,String}()
        this.scheme_order = Dict{String,Int64}()
        this.acceleration = "none"
        this.isFC = true
        return this
    end
end

# Backward-compatible alias: legacy user code using `Discrete_Ordinates()` keeps working.
const Discrete_Ordinates = SN

# Method(s)
"""
    set_particle(this::SN,particle::Particle)

To set the particle for which the transport discretization method is for. 

# Input Argument(s)
- `this::SN` : discretization method.
- `particle::Particle` : particle.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = SN()
julia> m.set_particle(electron) 
```
"""
function set_particle(this::SN,particle::Particle)
    this.particle = particle
end

"""
    set_solver_type(this::SN,solver_type::String)

To set the solver for the particle transport.

# Input Argument(s)
- `this::SN` : discretization method.
- `solver_type::String` : solver type, which can takes the following values:
    - `solver_type = "BTE"` : Boltzmann transport equation
    - `solver_type = "BFP"` : Boltzmann Fokker-Planck equation
    - `solver_type = "BCSD"` : Boltzmann-CSD equation
    - `solver_type = "FP"` : Fokker-Planck equation
    - `solver_type = "CSD"` : Continuous slowing-down only equation
    - `solver_type = "BFP-EF"` : Boltzmann Fokker-Planck without elastic

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = SN()
julia> m.set_particle("BFP")
```
"""
function set_solver_type(this::SN,solver_type::String)
    if uppercase(solver_type) ∉ ["BTE","BFP","BCSD","FP","CSD","BFP-EF"] error("Unkown type of solver.") end
    this.solver_type = uppercase(solver_type)
end

"""
    set_quadrature(this::SN,type::String,order::Int64,Qdims::Int64=0)

To set the quadrature properties for the discretization of the angular domain.

# Input Argument(s)
- `this::SN` : discretization method.
- `type::String` : type of quadrature, which can takes the following values:
    - `type = "gauss-legendre"` : Gauss-Legendre quadrature (1D Cartesian geometry only).
    - `type = "gauss-lobatto"` : Gauss-Lobatto quadrature (1D Cartesian geometry only).
    - `type = "carlson"` : Carlson quadrature (2D or 3D Cartesian geometry only).
    - `type = "gauss-legendre-chebychev"` : product quadrature between Gauss-Legendre and
      Chebychev quadratures (2D or 3D Cartesian geometry only).
    - `type = "lebedev"` : Lebedev quadrature (2D or 3D Cartesian geometry only).
- `order::Int64` : order of the quadrature, which is any integer greater than 2.
- `Qdims::Int64` : quadrature dimension, which can takes the following values:
    - `Qdims = 0` : default value, Qdims = 1 with 1D quadrature, Qdims = 3 with quadrature
      over the unit sphere in 1D or in 2D geometry or in 3D geometry.
    - `Qdims = 1` : 1D quadrature.
    - `Qdims = 2` : quadrature over half the unit-sphere.
    - `Qdims = 3` : quadrature over the unit sphere. 

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = SN()
julia> m.set_quadrature("gauss-legendre",4)
```
"""
function set_quadrature(this::SN,type::String,order::Int64,Qdims::Int64=0)
    if lowercase(type) ∉ ["gauss-legendre","gauss-lobatto","carlson","lebedev","gauss-legendre-chebychev"] error("Unknown quadrature type.") end
    if order ≤ 1 error("Quadrature order should be at least 2.") end
    if ~(0 ≤ Qdims ≤ 3) error("Quadrature dimension should be either 1, 2 or 3.") end
    this.quadrature_type = lowercase(type)
    this.quadrature_order = order
    this.quadrature_dimension = Qdims
end

"""
    set_legendre_order(this::SN,legendre_order::Int64)

To set the maximum order of the Legendre expansion of the differential cross-sections.

# Input Argument(s)
- `this::SN` : discretization method.
- `legendre_order::Int64` : maximum order of the Legendre expansion of the differential cross-sections.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = SN()
julia> m.set_legendre_order(7)
```
"""
function set_legendre_order(this::SN,legendre_order::Int64)
    if legendre_order < 0 error("Legendre order should be at least 0.") end
    this.legendre_order = legendre_order
end

"""
    set_angular_fokker_planck(this::SN,angular_fokker_planck::String)

To set the discretization method for the angular Fokker-Planck term.

# Input Argument(s)
- `this::SN` : discretization method.
- `angular_fokker_planck::String` : discretization method for the angular Fokker-Planck term, which can takes the following values:
    - `angular_fokker_planck = "finite-difference"` : finite difference discretization.
    - `angular_fokker_planck = "galerkin"` : galerkin moment-based discretization.
    - `angular_fokker_planck = "differential-quadrature"` : finite difference discretization (1D Cartesian geometry only).
    
# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = SN()
julia> m.set_angular_fokker_planck("differential-quadrature")
```
"""
function set_angular_fokker_planck(this::SN,angular_fokker_planck::String)
    if angular_fokker_planck ∉ ["finite-difference","differential-quadrature","galerkin"] error("Unkown method to deal with the angular Fokker-Planck term.") end
    this.angular_fokker_planck = angular_fokker_planck
end

"""
    set_angular_boltzmann(this::SN,angular_boltzmann::String)

To set the angular discretization method for the Boltzmann operator.

# Input Argument(s)
- `this::SN` : discretization method.
- `angular_boltzmann::String` : angular discretization method for the Boltzmann operator, which can takes the following values:
    - `angular_boltzmann = "standard"` : standard discrete ordinates (SN) method.
    - `angular_boltzmann = "galerkin-m"` : Galerkin method by inversion of the discrete-to-moment M matrix.
    - `angular_boltzmann = "galerkin-d"` : Galerkin method by inversion of the moment-to-discrete D matrix.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = SN()
julia> m.set_angular_boltzmann("standard")
```
"""
function set_angular_boltzmann(this::SN,angular_boltzmann::String)
    if angular_boltzmann ∉ ["standard","galerkin-m","galerkin-d","galerkin"] error("Unkown method to deal with the Boltzmann kernel.") end
    if (angular_boltzmann == "galerkin") angular_boltzmann = "galerkin-d" end
    this.angular_boltzmann = angular_boltzmann
end

"""
    set_convergence_criterion(this::SN,convergence_criterion::Float64)

To set the convergence criterion for the in-group iterations.

# Input Argument(s)
- `this::SN` : discretization method.
- `convergence_criterion::Float64` : convergence criterion.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = SN()
julia> m.set_convergence_criterion(1e-5)
```
"""
function set_convergence_criterion(this::SN,convergence_criterion::Float64)
    if convergence_criterion ≤ 0 error("Convergence criterion has to be greater than 0.") end
    this.convergence_criterion = convergence_criterion
end

"""
    set_maximum_iteration(this::SN,maximum_iteration::Int64)

To set the maximum number of in-group iterations.

# Input Argument(s)
- `this::SN` : discretization method.
- `maximum_iteration::Int64` : maximum number of in-group iterations.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = SN()
julia> m.set_maximum_iteration(50)
```
"""
function set_maximum_iteration(this::SN,maximum_iteration::Int64)
    if maximum_iteration < 1 error("Maximum iteration has to be at least 1.") end
    this.maximum_iteration = maximum_iteration
end

"""
    set_scheme(this::SN,axis::String,scheme_type::String,scheme_order::Int64)

To set the type of discretization scheme for derivative along the specified spatial or energy axis.

# Input Argument(s)
- `this::SN` : discretization method.
- `axis::String` : variable of the derivative for which the scheme is applied, which can takes the following values:
    - `axis = "x"` : spatial x axis (discretization of the streaming term)
    - `axis = "y"` : spatial y axis (discretization of the streaming term)
    - `axis = "z"` : spatial z axis (discretization of the streaming term)
    - `axis = "E"` : spatial E axis (discretization of the continuous slowing-down term)
- `scheme_type::String` : type of scheme to be applied, which can takes the following values:
    - `scheme_type = "DD"` : diamond difference scheme (any order)
    - `scheme_type = "DG"` : discontinuous Galerkin scheme (any order)
    - `scheme_type = "AWD"` : adaptive weighted scheme (1st and 2nd order only)
- `scheme_order::Int64` : scheme order, which takes values greater than 1.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = SN()
julia> m.set_scheme("x","DD",1)
julia> m.set_scheme("E","DG",2)
```
"""
function set_scheme(this::SN,axis::String,scheme_type::String,scheme_order::Int64)
    if axis ∉ ["x","y","z","E"] error("Unknown axis.") end
    if uppercase(scheme_type) ∉ ["DD","DG","DG-","DG+","AWD"] error("Unknown type of scheme.") end
    if scheme_order ≤ 0 error("Scheme order should be at least of 1.") end
    this.scheme_type[axis] = scheme_type
    this.scheme_order[axis] = scheme_order
end

"""
    set_acceleration(this::SN,acceleration::String)

To set the acceleration method for the in-group iteration process.

# Input Argument(s)
- `this::SN` : discretization method.
- `acceleration::String` : acceleration method, which takes the following values
    - `acceleration = "none"` : none
    - `acceleration = "livolant"` : livolant acceleration method
# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = SN()
julia> m.set_acceleration("livolant")
```
"""
function set_acceleration(this::SN,acceleration::String)
    if lowercase(acceleration) ∉ ["none","livolant"] error("Unkown acceleration method.") end
    this.acceleration = acceleration
end

"""
    get_legendre_order(this::SN)

Get the Legendre truncation order.

# Input Argument(s)
- `this::SN` : discretization method.

# Output Argument(s)
- `legendre_order::Int64` : Legendre truncation order.

"""
function get_legendre_order(this::SN)
    if ismissing(this.legendre_order) error("Unable to get Legendre order. Missing data.") end
    return this.legendre_order
end

"""
    get_quadrature_order(this::SN)

Get the quadrature order.

# Input Argument(s)
- `this::SN` : discretization method.

# Output Argument(s)
- `quadrature_order::Int64` : quadrature order.

"""
function get_quadrature_order(this::SN)
    if ismissing(this.quadrature_order) error("Unable to get quadrature order. Missing data.") end
    return this.quadrature_order
end

"""
    get_quadrature_type(this::SN)

Get the quadrature type.

# Input Argument(s)
- `this::SN` : discretization method.

# Output Argument(s)
- `quadrature_type::String` : quadrature type.

"""
function get_quadrature_type(this::SN)
    if ismissing(this.quadrature_type) error("Unable to get quadrature type. Missing data.") end
    return this.quadrature_type
end

"""
    get_angular_boltzmann(this::SN)

Get the type of angular discretization for the Boltzmann operator.

# Input Argument(s)
- `this::SN` : discretization method.

# Output Argument(s)
- `angular_boltzmann::String` : type of angular discretization for the Boltzmann operator.

"""
function get_angular_boltzmann(this::SN)
    if ismissing(this.angular_boltzmann) error("Unable to get angular Boltzmann treatment type. Missing data.") end
    return this.angular_boltzmann
end

"""
    get_angular_fokker_planck(this::SN)

Get the type of angular discretization for the Fokker-Planck operator.

# Input Argument(s)
- `this::SN` : discretization method.

# Output Argument(s)
- `angular_fokker_planck::String` : type of angular discretization for the Fokker-Planck
  operator.

"""
function get_angular_fokker_planck(this::SN)
    if ismissing(this.angular_fokker_planck) error("Unable to get angular Fokker-Planck treatment type. Missing data.") end
    return this.angular_fokker_planck
end

"""
    get_particle(this::SN)

Get the particle associated with the discretization methods.

# Input Argument(s)
- `this::SN` : discretization method.

# Output Argument(s)
- `particle::Particle` : particle.

"""
function get_particle(this::SN)
    if ismissing(this.particle) error("Unable to get particle. Missing data.") end
    return this.particle
end

"""
    get_solver_type(this::SN)

Get the type of solver for transport calculations.

# Input Argument(s)
- `this::SN` : discretization method.

# Output Argument(s)
- `solver::String` : type of solver for transport calculations.
- `isCSD::Bool` : indicate if continuous slowing-down term is used or not.

"""
function get_solver_type(this::SN)
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

"""
    get_schemes(this::SN,geometry::Geometry,isFC::Bool)

Get the space and/or energy schemes informations.

# Input Argument(s)
- `this::SN` : discretization method.
- `geometry::Geometry` : geometry.
- `isFC::Bool` : boolean indicating if the high-order moments are fully coupled.

# Output Argument(s)
- `schemes::Vector{String}` : scheme types.
- `𝒪::Vector{Int64}` : order of the schemes.
- `Nm::Vector{Int64}` : numbers of moments.

"""
function get_schemes(this::SN,geometry::Geometry,isFC::Bool)
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
    if isFC
        Nm = [𝒪[2]*𝒪[3]*𝒪[4],𝒪[1]*𝒪[3]*𝒪[4],𝒪[1]*𝒪[2]*𝒪[4],𝒪[1]*𝒪[2]*𝒪[3],prod(𝒪)]
    else
        Nm = [1+(𝒪[2]-1)+(𝒪[3]-1)+(𝒪[4]-1),1+(𝒪[1]-1)+(𝒪[3]-1)+(𝒪[4]-1),1+(𝒪[1]-1)+(𝒪[2]-1)+(𝒪[4]-1),1+(𝒪[1]-1)+(𝒪[2]-1)+(𝒪[3]-1),1+sum(𝒪.-1)]
    end
    return schemes,𝒪,Nm
end

"""
    get_convergence_criterion(this::SN)

Get the convergence criterion for in-group iteration convergence.

# Input Argument(s)
- `this::SN` : discretization method.

# Output Argument(s)
- `convergence_criterion::Float64` : convergence criterion.

"""
function get_convergence_criterion(this::SN)
    return this.convergence_criterion
end

"""
    get_maximum_iteration(this::SN)

Get the maximum number of iterations for in-group iteration convergence.

# Input Argument(s)
- `this::SN` : discretization method.

# Output Argument(s)
- `maximum_iteration::Int64` : maximum number of iterations.

"""
function get_maximum_iteration(this::SN)
    return this.maximum_iteration
end

"""
    get_acceleration(this::SN)

Get the acceleration method for in-group iteration convergence.

# Input Argument(s)
- `this::SN` : discretization method.

# Output Argument(s)
- `acceleration::String` : acceleration method.

"""
function get_acceleration(this::SN)
    return this.acceleration
end

"""
    get_quadrature_dimension(this::SN)

Get the quadrature dimension.

# Input Argument(s)
- `this::SN` : discretization method.

# Output Argument(s)
- `quadrature_dimension::Int64` : quadrature dimension.

"""
function get_quadrature_dimension(this::SN,Ndims::Int64)
    quadrature_type = this.get_quadrature_type()
    Qdims = this.get_quadrature_dimension()
    if quadrature_type ∈ ["gauss-legendre","gauss-lobatto"]
        if Qdims ∈ [2,3] error("Quadrature dimension should be 1 with quadrature over the direction cosine.") end
        return 1
    elseif quadrature_type ∈ ["gauss-legendre-chebychev","lebedev","carlson"]
        if Ndims ∈ [1,2]
            if Qdims == 1 error("Quadrature dimension should not be 1 with quadrature over the unit sphere.") end
            if Qdims == 2
                return 2
            else
                return 3
            end
        else
            if Qdims ∈ [1,2] error("Quadrature dimension should not be 1 or 2 in 3D geometry.") end
            return 3
        end
    else
        error("Unkown quadrature type.")
    end
end

"""
    set_is_full_coupling(this::SN,isFC::Bool)

Set, for multidimensional high-order schemes, if the high-order moments are fully coupled
or not. For example, with two linear schemes, the moments are either fully coupled
[00,10,01,11] or not [00,10,01].

# Input Argument(s)
- `this::SN` : discretization method.
- `isFC::Bool` : boolean indicating if the high-orde moments are fully coupled
  or not.

# Output Argument(s)
N/A

"""
function set_is_full_coupling(this::SN,isFC::Bool)
    this.isFC = isFC
end

"""
    get_is_full_coupling(this::SN)

Get, for multidimensional high-order schemes, if the high-order moments are fully coupled
or not. For example, with two linear schemes, the moments are either fully coupled
[00,10,01,11] or not [00,10,01].

# Input Argument(s)
- `this::SN` : discretization method.

# Output Argument(s)
- `isFC::Bool` : boolean indicating if the high-orde moments are fully coupled
  or not.

"""
function get_is_full_coupling(this::SN)
    return this.isFC
end

"""
    get_quadrature_dimension(this::SN)

Get the quadrature dimension.

# Input Argument(s)
- `this::SN` : discretization method.

# Output Argument(s)
- `quadrature_dimension::Int64` : quadrature dimension.

"""
function get_quadrature_dimension(this::SN)
    return this.quadrature_dimension
end