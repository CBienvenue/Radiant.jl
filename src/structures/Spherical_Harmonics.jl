mutable struct Spherical_Harmonics

    # Variable(s)
    particle                   ::Union{Missing,Particle}
    solver_type                ::Union{Missing,String}
    legendre_order             ::Union{Missing,Int64}
    convergence_criterion      ::Float64
    maximum_iteration          ::Int64
    scheme_type                ::Dict{String,String}
    scheme_order               ::Dict{String,Int64}
    acceleration               ::String
    isFC                       ::Bool
    polynomial_basis           ::Union{Missing,String}
    angular_fokker_planck      ::Union{Missing,String}

    # Constructor(s)
    function Spherical_Harmonics()
        this = new()
        this.particle = missing
        this.solver_type = missing
        this.legendre_order = 64
        this.convergence_criterion = 1e-7 
        this.maximum_iteration = 300
        this.scheme_type = Dict{String,String}()
        this.scheme_order = Dict{String,Int64}()
        this.acceleration = "none"
        this.isFC = true
        this.polynomial_basis = missing
        this.angular_fokker_planck = "galerkin"
        return this
    end
end

function set_particle(this::Spherical_Harmonics,particle::Particle)
    this.particle = particle
end

function set_solver_type(this::Spherical_Harmonics,solver_type::String)
    if uppercase(solver_type) ∉ ["BTE","BFP","BCSD","FP","CSD","BFP-EF"] error("Unkown type of solver.") end
    this.solver_type = uppercase(solver_type)
end

function set_legendre_order(this::Spherical_Harmonics,legendre_order::Int64)
    if legendre_order < 0 error("Legendre order should be at least 0.") end
    this.legendre_order = legendre_order
end

function set_convergence_criterion(this::Spherical_Harmonics,convergence_criterion::Float64)
    if convergence_criterion ≤ 0 error("Convergence criterion has to be greater than 0.") end
    this.convergence_criterion = convergence_criterion
end

function set_maximum_iteration(this::Spherical_Harmonics,maximum_iteration::Int64)
    if maximum_iteration < 1 error("Maximum iteration has to be at least 1.") end
    this.maximum_iteration = maximum_iteration
end

function set_scheme(this::Spherical_Harmonics,axis::String,scheme_type::String,scheme_order::Int64)
    if axis ∉ ["x","y","z","E"] error("Unknown axis.") end
    if uppercase(scheme_type) ∉ ["DD","DG","DG-","DG+","AWD"] error("Unknown type of scheme.") end
    if scheme_order ≤ 0 error("Scheme order should be at least of 1.") end
    this.scheme_type[axis] = scheme_type
    this.scheme_order[axis] = scheme_order
end

function set_acceleration(this::Spherical_Harmonics,acceleration::String)
    if lowercase(acceleration) ∉ ["none","livolant"] error("Unkown acceleration method.") end
    this.acceleration = acceleration
end

function set_polynomial_basis(this::Spherical_Harmonics,basis::String)
    if lowercase(basis) ∉ ["legendre","spherical-harmonics","cartesian-harmonics"] error("Unknown polynomial basis.") end
    this.polynomial_basis = lowercase(basis)
end

function get_is_full_coupling(this::Spherical_Harmonics)
    return this.isFC
end

function get_schemes(this::Spherical_Harmonics,geometry::Geometry,isFC::Bool)
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

function get_solver_type(this::Spherical_Harmonics)
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

function get_legendre_order(this::Spherical_Harmonics)
    if ismissing(this.legendre_order) error("Unable to get Legendre order. Missing data.") end
    return this.legendre_order
end

function get_particle(this::Spherical_Harmonics)
    if ismissing(this.particle) error("Unable to get particle. Missing data.") end
    return this.particle
end

function get_acceleration(this::Spherical_Harmonics)
    return this.acceleration
end

function get_convergence_criterion(this::Spherical_Harmonics)
    return this.convergence_criterion
end

function get_maximum_iteration(this::Spherical_Harmonics)
    return this.maximum_iteration
end

function get_polynomial_basis(this::Spherical_Harmonics,Ndims::Int64)
    if ismissing(this.polynomial_basis)
        if Ndims == 1
            this.polynomial_basis = "legendre"
        else
            this.polynomial_basis = "spherical-harmonics"
        end
        return this.polynomial_basis
    else
        return this.polynomial_basis
    end
end
