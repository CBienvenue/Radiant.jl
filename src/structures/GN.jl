"""
    GN

Structure used to define the moment-based angular discretization method that subdivides
the angular domain and uses a polynomial expansion on each subdivision. It generalizes
the `DPN` solver.

# Mandatory field(s)
- `particle::Particle` : particle for which the discretization methods is defined.
- `solver_type::String` : type of solver for the transport calculations.
- `scheme_type::Dict{String,String}` : type of schemes for the spatial or energy discretization.
- `scheme_order::Dict{String,Int64}` : order of the expansion for the discretization schemes.

# Optional field(s) - with default values
- `legendre_order::Int64 = 16` : global Legendre order.
- `legendre_order_local::Int64 = 0` : per-subdivision local Legendre order.
- `subdivision::Int64 = 1` : number of angular subdivisions.
- `polynomial_basis::String` : polynomial basis (`"legendre"` in 1D,
  `"spherical-harmonics"` in 2D or 3D by default).
- `angular_fokker_planck::String = "galerkin"` : type of discretization for the angular
  Fokker-Planck operation.
- `convergence_criterion::Float64 = 1e-7` : convergence criterion of in-group iterations.
- `maximum_iteration::Int64 = 300` : maximum number of in-group iterations.
- `acceleration::String = "none"` : acceleration method for the in-group iterations.

"""
mutable struct GN

    # Variable(s)
    particle                   ::Union{Missing,Particle}
    solver_type                ::Union{Missing,String}
    legendre_order             ::Union{Missing,Int64}
    legendre_order_local ::Union{Missing,Int64}
    convergence_criterion      ::Float64
    maximum_iteration          ::Int64
    scheme_type                ::Dict{String,String}
    scheme_order               ::Dict{String,Int64}
    acceleration               ::String
    gmres_restart              ::Int64
    anderson_depth             ::Int64
    isFC                       ::Bool
    polynomial_basis           ::Union{Missing,String}
    angular_fokker_planck      ::Union{Missing,String}
    subdivision                ::Int64
    tiling                     ::String
    z_fold                     ::Bool

    # Constructor(s)
    function GN()
        this = new()
        this.particle = missing
        this.solver_type = missing
        this.legendre_order = 16
        this.legendre_order_local = 0
        this.convergence_criterion = 1e-7
        this.maximum_iteration = 300
        this.scheme_type = Dict{String,String}()
        this.scheme_order = Dict{String,Int64}()
        this.acceleration = "none"
        this.gmres_restart = 30
        this.anderson_depth = 3
        this.isFC = true
        this.polynomial_basis = missing
        this.angular_fokker_planck = "galerkin"
        this.subdivision = 1
        this.tiling = "symmetric"
        this.z_fold = true
        return this
    end
end

# Method(s)
"""
    set_particle(this::GN,particle::Particle)

To set the particle for which the transport discretization method is for.

# Input Argument(s)
- `this::GN` : discretization method.
- `particle::Particle` : particle.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = GN()
julia> m.set_particle(electron)
```
"""
function set_particle(this::GN,particle::Particle)
    this.particle = particle
end

"""
    set_solver_type(this::GN,solver_type::String)

To set the solver for the particle transport.

# Input Argument(s)
- `this::GN` : discretization method.
- `solver_type::String` : solver type, which can take the following values:
    - `solver_type = "BTE"` : Boltzmann transport equation.
    - `solver_type = "BFP"` : Boltzmann Fokker-Planck equation.
    - `solver_type = "BCSD"` : Boltzmann-CSD equation.
    - `solver_type = "FP"` : Fokker-Planck equation.
    - `solver_type = "CSD"` : Continuous slowing-down only equation.
    - `solver_type = "BFP-EF"` : Boltzmann Fokker-Planck without elastic.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = GN()
julia> m.set_solver_type("BFP")
```
"""
function set_solver_type(this::GN,solver_type::String)
    if uppercase(solver_type) ∉ ["BTE","BFP","BCSD","FP","CSD","BFP-EF"] error("Unkown type of solver.") end
    this.solver_type = uppercase(solver_type)
end

"""
    set_legendre_order(this::GN,legendre_order::Int64,legendre_order_local::Int64)

To set the global Legendre order and the per-subdivision local Legendre order.

# Input Argument(s)
- `this::GN` : discretization method.
- `legendre_order::Int64` : global Legendre order.
- `legendre_order_local::Int64` : per-subdivision local Legendre order.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = GN()
julia> m.set_legendre_order(16,0)
```
"""
function set_legendre_order(this::GN,legendre_order::Int64,legendre_order_local::Int64)
    if legendre_order < 0 || legendre_order_local < 0 error("Legendre order should be at least 0.") end
    this.legendre_order = legendre_order
    this.legendre_order_local = legendre_order_local
end

"""
    set_legendre_order(this::GN,legendre_order::Int64)

Set the global Legendre truncation order and couple the local patch order to it
(`legendre_order_local = legendre_order`), i.e. one fully-coupled patch per base domain.
This is the Double-PN convenience used by the [`DPN`](@ref) constructor.

# Input Argument(s)
- `this::GN` : discretization method.
- `legendre_order::Int64` : global (and local) Legendre truncation order.

# Output Argument(s)
N/A

"""
function set_legendre_order(this::GN,legendre_order::Int64)
    set_legendre_order(this,legendre_order,legendre_order)
end

"""
    set_convergence_criterion(this::GN,convergence_criterion::Float64)

To set the convergence criterion for the in-group iterations.

# Input Argument(s)
- `this::GN` : discretization method.
- `convergence_criterion::Float64` : convergence criterion.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = GN()
julia> m.set_convergence_criterion(1e-5)
```
"""
function set_convergence_criterion(this::GN,convergence_criterion::Float64)
    if convergence_criterion ≤ 0 error("Convergence criterion has to be greater than 0.") end
    this.convergence_criterion = convergence_criterion
end

"""
    set_maximum_iteration(this::GN,maximum_iteration::Int64)

To set the maximum number of in-group iterations.

# Input Argument(s)
- `this::GN` : discretization method.
- `maximum_iteration::Int64` : maximum number of in-group iterations.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = GN()
julia> m.set_maximum_iteration(50)
```
"""
function set_maximum_iteration(this::GN,maximum_iteration::Int64)
    if maximum_iteration < 1 error("Maximum iteration has to be at least 1.") end
    this.maximum_iteration = maximum_iteration
end

"""
    set_scheme(this::GN,axis::String,scheme_type::String,scheme_order::Int64)

To set the type of discretization scheme for derivative along the specified spatial or
energy axis.

# Input Argument(s)
- `this::GN` : discretization method.
- `axis::String` : variable of the derivative for which the scheme is applied, which can
  take the following values:
    - `axis = "x"` : spatial x axis (discretization of the streaming term).
    - `axis = "y"` : spatial y axis (discretization of the streaming term).
    - `axis = "z"` : spatial z axis (discretization of the streaming term).
    - `axis = "E"` : energy axis (discretization of the continuous slowing-down term).
- `scheme_type::String` : type of scheme to be applied, which can take the following values:
    - `scheme_type = "DD"` : diamond difference scheme (any order).
    - `scheme_type = "DG"` : discontinuous Galerkin scheme (any order).
    - `scheme_type = "AWD"` : adaptive weighted scheme (1st and 2nd order only).
- `scheme_order::Int64` : scheme order, which takes values greater than 1.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = GN()
julia> m.set_scheme("x","DD",1)
julia> m.set_scheme("E","DG",2)
```
"""
function set_scheme(this::GN,axis::String,scheme_type::String,scheme_order::Int64)
    if axis ∉ ["x","y","z","E"] error("Unknown axis.") end
    if uppercase(scheme_type) ∉ ["DD","DG","DG-","DG+","AWD"] error("Unknown type of scheme.") end
    if scheme_order ≤ 0 error("Scheme order should be at least of 1.") end
    this.scheme_type[axis] = scheme_type
    this.scheme_order[axis] = scheme_order
end

"""
    set_acceleration(this::GN,acceleration::String,parameter::Int64=0)

To set the acceleration method for the in-group iteration process.

# Input Argument(s)
- `this::GN` : discretization method.
- `acceleration::String` : acceleration method, which takes the following values:
    - `acceleration = "none"` : source iteration without acceleration.
    - `acceleration = "livolant"` : Livolant two-point extrapolation.
    - `acceleration = "anderson"` : depth-`m` Anderson acceleration.
    - `acceleration = "gmres"` : matrix-free restarted GMRES.
    - `acceleration = "bicgstab"` : matrix-free BiCGStab.
- `parameter::Int64` : optional tuning parameter for the chosen method. It is the Krylov
   subspace size before restart for `"gmres"` (default 30) and the memory depth for
   `"anderson"` (default 3); it is ignored by the other methods. A value of `0` keeps the
   current default.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = GN()
julia> m.set_acceleration("gmres")      # GMRES with the default restart of 30
julia> m.set_acceleration("gmres",50)   # GMRES restarted every 50 Krylov vectors
julia> m.set_acceleration("anderson",5) # Anderson acceleration with memory depth 5
```
"""
function set_acceleration(this::GN,acceleration::String,parameter::Int64=0)
    accel = lowercase(acceleration)
    if accel ∉ ["none","livolant","anderson","gmres","bicgstab"] error("Unkown acceleration method.") end
    this.acceleration = accel
    if parameter != 0
        if parameter < 1 error("Acceleration parameter has to be at least 1.") end
        if accel == "gmres"
            this.gmres_restart = parameter
        elseif accel == "anderson"
            this.anderson_depth = parameter
        end
    end
end

"""
    set_polynomial_basis(this::GN,basis::String)

To set the polynomial basis used for the angular discretization.

# Input Argument(s)
- `this::GN` : discretization method.
- `basis::String` : polynomial basis, which can take the following values:
    - `basis = "legendre"` : Legendre polynomials (1D only).
    - `basis = "spherical-harmonics"` : spherical harmonics (1D, 2D or 3D).

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = GN()
julia> m.set_polynomial_basis("legendre")
```
"""
function set_polynomial_basis(this::GN,basis::String)
    if lowercase(basis) ∉ ["legendre","spherical-harmonics"] error("Unknown polynomial basis.") end
    this.polynomial_basis = lowercase(basis)
end

"""
    set_is_full_coupling(this::GN,isFC::Bool)

Set, for multidimensional high-order schemes, if the high-order moments are fully coupled
or not. For example, with two linear schemes, the moments are either fully coupled
`[00,10,01,11]` or not `[00,10,01]`.

# Input Argument(s)
- `this::GN` : discretization method.
- `isFC::Bool` : boolean indicating if the high-order moments are fully coupled or not.

# Output Argument(s)
N/A

"""
function set_is_full_coupling(this::GN,isFC::Bool)
    this.isFC = isFC
end

"""
    set_angular_fokker_planck(this::GN,angular_fokker_planck::String)

To set the discretization method for the angular Fokker-Planck term.

# Input Argument(s)
- `this::GN` : discretization method.
- `angular_fokker_planck::String` : discretization method for the angular Fokker-Planck
  term, which can take the following values:
    - `angular_fokker_planck = "galerkin"` : Galerkin moment-based discretization
      (diagonal on full-range spherical harmonics, `ℳ[p,p] = -ℓ_p(ℓ_p+1)`).
    - `angular_fokker_planck = "finite-difference"` : finite-volume (TPFA)
      discretization of the Laplace-Beltrami operator on the GN patch graph,
      using the natural connectivity of the tiling (no Voronoi needed).

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = GN()
julia> m.set_angular_fokker_planck("finite-difference")
```
"""
function set_angular_fokker_planck(this::GN,angular_fokker_planck::String)
    if angular_fokker_planck ∉ ["galerkin","finite-difference"] error("Unkown method to deal with the angular Fokker-Planck term.") end
    this.angular_fokker_planck = angular_fokker_planck
end

"""
    set_subdivision(this::GN,subdivision::Int64)

To set the number of angular subdivisions used by the `GN` solver.

# Input Argument(s)
- `this::GN` : discretization method.
- `subdivision::Int64` : number of angular subdivisions (must be at least 1).

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = GN()
julia> m.set_subdivision(4)
```
"""
function set_subdivision(this::GN,subdivision::Int64)
    if subdivision < 1 error("Subdivision must be at least 1.") end
    this.subdivision = subdivision
end

"""
    set_tiling(this::GN,tiling::String)

To set the angular tiling strategy used by the `GN` solver with the spherical
harmonics basis. The tiling controls how each octant of the angular sphere is
partitioned into patches.

# Input Argument(s)
- `this::GN` : discretization method.
- `tiling::String` : tiling strategy, which can take the following values:
    - `tiling = "polar-anchored"` : default. Each octant is split into `Nv` polar
      bands in μ with `Nv+1-v` azimuthal patches per band (triangular tiling
      anchored at the polar axis, aligned with the spatial x-axis).
    - `tiling = "symmetric"` : each octant is barycentrically subdivided into
      `Nv²` spherical sub-triangles, treating the three axis-vertices (±x, ±y, ±z)
      equivalently. Restores symmetry under permutation of the spatial axes.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = GN()
julia> m.set_tiling("symmetric")
```
"""
function set_tiling(this::GN,tiling::String)
    if lowercase(tiling) ∉ ["polar-anchored","symmetric"] error("Unknown tiling strategy.") end
    this.tiling = lowercase(tiling)
end

"""
    get_is_full_coupling(this::GN)

Get, for multidimensional high-order schemes, if the high-order moments are fully coupled
or not.

# Input Argument(s)
- `this::GN` : discretization method.

# Output Argument(s)
- `isFC::Bool` : boolean indicating if the high-order moments are fully coupled or not.

"""
function get_is_full_coupling(this::GN)
    return this.isFC
end

"""
    get_schemes(this::GN,geometry::Geometry,isFC::Bool)

Get the space and/or energy schemes informations.

# Input Argument(s)
- `this::GN` : discretization method.
- `geometry::Geometry` : geometry.
- `isFC::Bool` : boolean indicating if the high-order moments are fully coupled.

# Output Argument(s)
- `schemes::Vector{String}` : scheme types.
- `𝒪::Vector{Int64}` : order of the schemes.
- `Nm::Vector{Int64}` : numbers of moments.

"""
function get_schemes(this::GN,geometry::Geometry,isFC::Bool)
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
    get_solver_type(this::GN)

Get the type of solver for transport calculations.

# Input Argument(s)
- `this::GN` : discretization method.

# Output Argument(s)
- `solver::Int64` : type of solver for transport calculations.
- `isCSD::Bool` : indicate if continuous slowing-down term is used or not.

"""
function get_solver_type(this::GN)
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
    get_legendre_order(this::GN)

Get the global Legendre truncation order.

# Input Argument(s)
- `this::GN` : discretization method.

# Output Argument(s)
- `legendre_order::Int64` : Legendre truncation order.

"""
function get_legendre_order(this::GN)
    if ismissing(this.legendre_order) error("Unable to get Legendre order. Missing data.") end
    return this.legendre_order
end

"""
    get_legendre_order_local(this::GN)

Get the per-subdivision local Legendre truncation order.

# Input Argument(s)
- `this::GN` : discretization method.

# Output Argument(s)
- `legendre_order_local::Int64` : per-subdivision local Legendre truncation order.

"""
function get_legendre_order_local(this::GN)
    if ismissing(this.legendre_order) || ismissing(this.legendre_order_local) error("Unable to get Legendre order. Missing data.") end
    return this.legendre_order_local
end

"""
    get_particle(this::GN)

Get the particle associated with the discretization methods.

# Input Argument(s)
- `this::GN` : discretization method.

# Output Argument(s)
- `particle::Particle` : particle.

"""
function get_particle(this::GN)
    if ismissing(this.particle) error("Unable to get particle. Missing data.") end
    return this.particle
end

"""
    get_acceleration(this::GN)

Get the acceleration method for in-group iteration convergence.

# Input Argument(s)
- `this::GN` : discretization method.

# Output Argument(s)
- `acceleration::String` : acceleration method.

"""
function get_acceleration(this::GN)
    return this.acceleration
end

"""
    get_gmres_restart(this::GN)

Get the GMRES restart parameter for in-group iteration convergence.

# Input Argument(s)
- `this::GN` : discretization method.

# Output Argument(s)
- `gmres_restart::Int64` : Krylov subspace size before restart.

"""
function get_gmres_restart(this::GN)
    return this.gmres_restart
end

"""
    get_anderson_depth(this::GN)

Get the Anderson acceleration memory depth for in-group iteration convergence.

# Input Argument(s)
- `this::GN` : discretization method.

# Output Argument(s)
- `anderson_depth::Int64` : Anderson memory depth.

"""
function get_anderson_depth(this::GN)
    return this.anderson_depth
end

"""
    get_convergence_criterion(this::GN)

Get the convergence criterion for in-group iteration convergence.

# Input Argument(s)
- `this::GN` : discretization method.

# Output Argument(s)
- `convergence_criterion::Float64` : convergence criterion.

"""
function get_convergence_criterion(this::GN)
    return this.convergence_criterion
end

"""
    get_maximum_iteration(this::GN)

Get the maximum number of iterations for in-group iteration convergence.

# Input Argument(s)
- `this::GN` : discretization method.

# Output Argument(s)
- `maximum_iteration::Int64` : maximum number of iterations.

"""
function get_maximum_iteration(this::GN)
    return this.maximum_iteration
end

"""
    get_polynomial_basis(this::GN,Ndims::Int64)

Get the polynomial basis used for the angular discretization. If the basis was not set
explicitly, the default basis is selected based on the geometry dimension (`"legendre"`
in 1D, `"spherical-harmonics"` otherwise).

# Input Argument(s)
- `this::GN` : discretization method.
- `Ndims::Int64` : geometry dimension.

# Output Argument(s)
- `polynomial_basis::String` : polynomial basis.

"""
function get_polynomial_basis(this::GN,Ndims::Int64)
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

"""
    get_angular_fokker_planck(this::GN)

Get the type of angular discretization for the Fokker-Planck operator.

# Input Argument(s)
- `this::GN` : discretization method.

# Output Argument(s)
- `angular_fokker_planck::String` : type of angular discretization for the Fokker-Planck
  operator.

"""
function get_angular_fokker_planck(this::GN)
    if ismissing(this.angular_fokker_planck) error("Unable to get angular Fokker-Planck treatment type. Missing data.") end
    return this.angular_fokker_planck
end

"""
    get_subdivision(this::GN)

Get the number of angular subdivisions.

# Input Argument(s)
- `this::GN` : discretization method.

# Output Argument(s)
- `subdivision::Int64` : number of angular subdivisions.

"""
function get_subdivision(this::GN)
    return this.subdivision
end

"""
    get_tiling(this::GN)

Get the angular tiling strategy.

# Input Argument(s)
- `this::GN` : discretization method.

# Output Argument(s)
- `tiling::String` : angular tiling strategy (`"polar-anchored"` or `"symmetric"`).

"""
function get_tiling(this::GN)
    return this.tiling
end

"""
    set_z_fold(this::GN,z_fold::Bool)

Enable or disable the angular-symmetry fold of the reduced-dimension angular domain
(`z_fold = true` by default, the DPN-equivalent treatment):

- In 1D with the spherical-harmonics basis, the azimuthally-symmetric flux is folded
  onto two half-spheres (`z_fold = true`) instead of the full octant tiling with
  azimuthal subdivision (`z_fold = false`).
- In 2D the `μ_z → -μ_z` symmetry folds the sphere onto four quadrants
  (`z_fold = true`, with subdivision `Nv = 1`) instead of eight octants
  (`z_fold = false`).
- In 3D it has no effect (eight octants either way).

With `z_fold = false` the full octant tiling is used in every dimension (the general
GN treatment). The 1D Legendre basis is always azimuthally collapsed and is unaffected.

# Input Argument(s)
- `this::GN` : discretization method.
- `z_fold::Bool` : whether to fold the angular domain by the problem's symmetry.

# Output Argument(s)
N/A

"""
function set_z_fold(this::GN,z_fold::Bool)
    this.z_fold = z_fold
end

"""
    get_z_fold(this::GN)

Get whether the 2D z-symmetry fold (four quadrants) is enabled.

# Input Argument(s)
- `this::GN` : discretization method.

# Output Argument(s)
- `z_fold::Bool` : whether the z-symmetry fold is enabled.

"""
function get_z_fold(this::GN)
    return this.z_fold
end
