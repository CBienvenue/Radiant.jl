
mutable struct Method

    # Variable(s)
    name                       ::Union{Missing,String}
    particle                   ::Union{Missing,String}
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
    function Method()

        this = new()

        this.name = missing
        this.particle = missing
        this.solver_type = missing
        this.quadrature_type = missing
        this.quadrature_order = missing
        this.legendre_order = missing
        this.angular_fokker_planck = missing
        this.angular_boltzmann = missing
        this.convergence_criterion = 1e-7 
        this.maximum_iteration = 500
        this.scheme_type = Dict{String,String}()
        this.scheme_order = Dict{String,Int64}()
        this.acceleration = "none"

        return this
    end
end

# Method(s)
Base.propertynames(::Method) = 
(
    fieldnames(Method)...,
    :set_particle,
    :set_solver_type,
    :set_quadrature,
    :set_legendre_order,
    :set_angular_fokker_planck,
    :set_angular_boltzmann,
    :set_convergence_criterion,
    :set_maximum_iteration,
    :set_scheme,
    :set_acceleration,
    :get_legendre_order,
    :get_quadrature_order,
    :get_quadrature_type,
    :get_angular_boltzmann,
    :get_angular_fokker_planck,
    :get_particle,
    :get_solver_type,
    :get_schemes,
    :get_convergence_criterion,
    :get_maximum_iteration
    :get_acceleration
)

function println(this::Method)
    entries = ["Name","Particle","Solver type","Quadrature type","Quadrature order","Legendre order","Boltzmann angular","Fokker-Planck angular","Convergence criterion", "Maximum iterations","Schemes type","Scheme order"]
    values = [this.name,this.particle,this.solver_type,this.quadrature_type,this.quadrature_order,this.legendre_order,this.angular_boltzmann,this.angular_fokker_planck,this.convergence_criterion,this.maximum_iteration,this.scheme_type,this.scheme_order]
    N = length(entries); L = length.(entries); Lmax = maximum(L)
    println("Method")
    for n in range(1,N)
        println(string("   ",entries[n]," "^(Lmax-L[n])),"  :  ",values[n])
    end
end

function set_particle(this::Method,particle::String)
    if lowercase(particle) ∉ ["photons","electrons","positrons"] error("Unknown particle type") end
    this.particle = particle
end

function set_solver_type(this::Method,solver_type::String)
    if uppercase(solver_type) ∉ ["BTE","BFP","BCSD","FP","CSD","BFP-EF"] error("Unkown type of solver.") end
    this.solver_type = uppercase(solver_type)
end

function set_quadrature(this::Method,type::String,order::Int64)
    if lowercase(type) ∉ ["gauss-legendre","gauss-lobatto","carlson","lebedev","gauss-legendre-chebychev"] error("Unknown quadrature type.") end
    if order ≤ 1 error("Quadrature order should be at least 2.") end
    this.quadrature_type = lowercase(type)
    this.quadrature_order = order
end

function set_legendre_order(this::Method,legendre_order::Int64)
    if legendre_order < 0 error("Legendre order should be at least 0.") end
    this.legendre_order = legendre_order
end

function set_angular_fokker_planck(this::Method,angular_fokker_planck::String)
    if angular_fokker_planck ∉ ["finite-element","differential-quadrature","galerkin"] error("Unkown method to deal with the angular Fokker-Planck term.") end
    this.angular_fokker_planck = angular_fokker_planck
end

function set_angular_boltzmann(this::Method,angular_boltzmann::String)
    if angular_boltzmann ∉ ["standard","galerkin-m","galerkin-d"] error("Unkown method to deal with the Boltzmann kernel.") end
    this.angular_boltzmann = angular_boltzmann
end

function set_convergence_criterion(this::Method,convergence_criterion::Float64)
    if convergence_criterion ≤ 0 error("Convergence criterion has to be greater than 0.") end
    this.convergence_criterion = convergence_criterion
end

function set_maximum_iteration(this::Method,maximum_iteration::Int64)
    if maximum_iteration < 1 error("Maximum iteration has to be at least 1.") end
    this.maximum_iteration = maximum_iteration
end

function set_scheme(this::Method,axis::String,scheme_type::String,scheme_order::Int64)
    if axis ∉ ["x","y","z","E"] error("Unknown axis.") end
    if uppercase(scheme_type) ∉ ["DD","DG","DG-","DG+","AWD"] error("Unknown type of scheme.") end
    if scheme_order ≤ 0 error("Scheme order should be at least of 1.") end
    this.scheme_type[axis] = scheme_type
    this.scheme_order[axis] = scheme_order
end

function set_acceleration(this::Method,acceleration::String)
    if lowercase(acceleration) ∉ ["none","livolant","anderson"] error("Unkown acceleration method.") end
    this.acceleration = acceleration
end

function get_legendre_order(this::Method)
    if ismissing(this.legendre_order) error("Unable to get Legendre order. Missing data.") end
    return this.legendre_order
end

function get_quadrature_order(this::Method)
    if ismissing(this.quadrature_order) error("Unable to get quadrature order. Missing data.") end
    return this.quadrature_order
end

function get_quadrature_type(this::Method)
    if ismissing(this.quadrature_type) error("Unable to get quadrature type. Missing data.") end
    return this.quadrature_type
end

function get_angular_boltzmann(this::Method)
    if ismissing(this.angular_boltzmann) error("Unable to get angular Boltzmann treatment type. Missing data.") end
    return this.angular_boltzmann
end

function get_angular_fokker_planck(this::Method)
    if ismissing(this.angular_fokker_planck) error("Unable to get angular Fokker-Planck treatment type. Missing data.") end
    return this.angular_fokker_planck
end

function get_particle(this::Method)
    if ismissing(this.particle) error("Unable to get particle. Missing data.") end
    return this.particle
end

function get_solver_type(this::Method)
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

function get_schemes(this::Method,geometry::Geometry,is_full_coupling::Bool)
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

function get_convergence_criterion(this::Method)
    return this.convergence_criterion
end

function get_maximum_iteration(this::Method)
    return this.maximum_iteration
end

function get_acceleration(this::Method)
    return this.acceleration
end