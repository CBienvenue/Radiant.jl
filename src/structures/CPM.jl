"""
    CPM

Structure used to define the collision-probability method (CPM) for the transport of a
particle. The volumetric and surface fluxes are expanded in spherical harmonics of arbitrary
order; in one-dimensional Cartesian geometry the azimuthal symmetry reduces the expansion to
Legendre polynomials of degree up to `legendre_order` (see the TR-02 technical report). The
method assembles the collision, leakage and transmission probability matrices and solves the
resulting linear system directly, as an alternative to the `SN` and `GN` solvers.

# Mandatory field(s)
- `particle::Particle` : particle for which the discretization method is defined.
- `solver_type::String` : type of solver for the transport calculations (`"BTE"`).

# Optional field(s) - with default values
- `legendre_order::Int64 = 0` : Legendre order of the volumetric flux expansion.
- `surface_order::Int64 = 0` : order of the surface (half-range) flux expansion.
- `convergence_criterion::Float64 = 1e-7` : convergence criterion of the group iterations.
- `maximum_iteration::Int64 = 300` : maximum number of group iterations.

"""
mutable struct CPM

    # Variable(s)
    particle                   ::Union{Missing,Particle}
    solver_type                ::Union{Missing,String}
    mode                       ::String
    legendre_order             ::Int64
    surface_order              ::Int64
    convergence_criterion      ::Float64
    maximum_iteration          ::Int64
    acceleration               ::String
    gmres_restart              ::Int64
    anderson_depth             ::Int64
    isFC                       ::Bool

    # Constructor(s)
    function CPM()
        this = new()
        this.particle = missing
        this.solver_type = missing
        this.mode = "global"
        this.legendre_order = 0
        this.surface_order = 0
        this.convergence_criterion = 1e-7
        this.maximum_iteration = 300
        this.acceleration = "none"
        this.gmres_restart = 30
        this.anderson_depth = 3
        this.isFC = true
        return this
    end
end

# Method(s)
"""
    set_particle(this::CPM,particle::Particle)

To set the particle for which the transport discretization method is for.

# Input Argument(s)
- `this::CPM` : collision-probability method.
- `particle::Particle` : particle.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = CPM()
julia> m.set_particle(electron)
```
"""
function set_particle(this::CPM,particle::Particle)
    this.particle = particle
end

"""
    set_solver_type(this::CPM,solver_type::String)

To set the solver for the particle transport.

# Input Argument(s)
- `this::CPM` : collision-probability method.
- `solver_type::String` : solver type, which can take the following value:
    - `solver_type = "BTE"` : Boltzmann transport equation.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = CPM()
julia> m.set_solver_type("BTE")
```
"""
function set_solver_type(this::CPM,solver_type::String)
    if uppercase(solver_type) ∉ ["BTE"] error("The CPM solver only supports the BTE solver type.") end
    this.solver_type = uppercase(solver_type)
end

"""
    set_legendre_order(this::CPM,legendre_order::Int64)
    set_legendre_order(this::CPM,volume_order::Int64,surface_order::Int64)

To set the Legendre order of the volumetric and surface flux expansions. The single-argument
form sets both the volumetric and the surface order to the same value; the two-argument form
sets them independently.

# Input Argument(s)
- `this::CPM` : collision-probability method.
- `legendre_order::Int64` : Legendre order of both the volumetric and surface flux expansions.
- `volume_order::Int64` : Legendre order of the volumetric flux expansion.
- `surface_order::Int64` : Legendre order of the surface (half-range) flux expansion.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = CPM()
julia> m.set_legendre_order(2)      # volume = surface = 2
julia> m.set_legendre_order(2,4)    # volume = 2, surface = 4
```
"""
function set_legendre_order(this::CPM,volume_order::Int64,surface_order::Int64)
    if volume_order < 0 || surface_order < 0 error("Legendre order should be at least 0.") end
    this.legendre_order = volume_order
    this.surface_order = surface_order
end

function set_legendre_order(this::CPM,legendre_order::Int64)
    set_legendre_order(this,legendre_order,legendre_order)
end

"""
    set_surface_order(this::CPM,surface_order::Int64)

To set the order of the surface (half-range) flux expansion used by the boundary-coupling
probabilities.

# Input Argument(s)
- `this::CPM` : collision-probability method.
- `surface_order::Int64` : order of the surface flux expansion.

# Output Argument(s)
N/A

"""
function set_surface_order(this::CPM,surface_order::Int64)
    if surface_order < 0 error("Surface order should be at least 0.") end
    this.surface_order = surface_order
end

"""
    set_mode(this::CPM,mode::String)

To set the solution mode of the collision-probability method.

# Input Argument(s)
- `this::CPM` : collision-probability method.
- `mode::String` : solution mode, which can take the following values:
    - `mode = "global"` : finite-domain method. The collision-probability matrices couple all
      regions and surfaces directly and the resulting linear system is solved by a direct
      factorization. This is the default.
    - `mode = "sweeping"` : sweeping method. The unit cell is a single voxel and the surface
      (interface) fluxes are propagated cell by cell, in direct analogy with the discrete-
      ordinates sweep; the within-group scattering is resolved by source iteration.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> m = CPM()
julia> m.set_mode("sweeping")
```
"""
function set_mode(this::CPM,mode::String)
    if lowercase(mode) ∉ ["global","sweeping"] error("Unknown CPM mode (use \"global\" or \"sweeping\").") end
    this.mode = lowercase(mode)
end

"""
    get_mode(this::CPM)

Get the solution mode of the collision-probability method.

# Input Argument(s)
- `this::CPM` : collision-probability method.

# Output Argument(s)
- `mode::String` : solution mode (`"global"` or `"sweeping"`).

"""
function get_mode(this::CPM)
    return this.mode
end

"""
    set_convergence_criterion(this::CPM,convergence_criterion::Float64)

To set the convergence criterion for the group iterations.

# Input Argument(s)
- `this::CPM` : collision-probability method.
- `convergence_criterion::Float64` : convergence criterion.

# Output Argument(s)
N/A

"""
function set_convergence_criterion(this::CPM,convergence_criterion::Float64)
    if convergence_criterion ≤ 0 error("Convergence criterion has to be greater than 0.") end
    this.convergence_criterion = convergence_criterion
end

"""
    set_maximum_iteration(this::CPM,maximum_iteration::Int64)

To set the maximum number of group iterations.

# Input Argument(s)
- `this::CPM` : collision-probability method.
- `maximum_iteration::Int64` : maximum number of group iterations.

# Output Argument(s)
N/A

"""
function set_maximum_iteration(this::CPM,maximum_iteration::Int64)
    if maximum_iteration < 1 error("Maximum iteration has to be at least 1.") end
    this.maximum_iteration = maximum_iteration
end

"""
    get_particle(this::CPM)

Get the particle associated with the discretization method.

# Input Argument(s)
- `this::CPM` : collision-probability method.

# Output Argument(s)
- `particle::Particle` : particle.

"""
function get_particle(this::CPM)
    if ismissing(this.particle) error("Unable to get particle. Missing data.") end
    return this.particle
end

"""
    get_solver_type(this::CPM)

Get the type of solver for transport calculations.

# Input Argument(s)
- `this::CPM` : collision-probability method.

# Output Argument(s)
- `solver::Int64` : type of solver for transport calculations.
- `isCSD::Bool` : indicate if the continuous slowing-down term is used or not.

"""
function get_solver_type(this::CPM)
    if ismissing(this.solver_type) error("Unable to get solver type. Missing data.") end
    if this.solver_type == "BTE"
        return 1, false
    else
        error("The CPM solver only supports the BTE solver type.")
    end
end

"""
    get_legendre_order(this::CPM)

Get the Legendre order of the volumetric flux expansion.

# Input Argument(s)
- `this::CPM` : collision-probability method.

# Output Argument(s)
- `legendre_order::Int64` : Legendre order.

"""
function get_legendre_order(this::CPM)
    return this.legendre_order
end

"""
    get_surface_order(this::CPM)

Get the order of the surface (half-range) flux expansion.

# Input Argument(s)
- `this::CPM` : collision-probability method.

# Output Argument(s)
- `surface_order::Int64` : surface order.

"""
function get_surface_order(this::CPM)
    return this.surface_order
end

"""
    get_schemes(this::CPM,geometry::Geometry,isFC::Bool)

Get the space and/or energy schemes informations. The CPM assumes a flat (region-averaged)
flux, so a single spatial/energy moment is carried per region.

# Input Argument(s)
- `this::CPM` : collision-probability method.
- `geometry::Geometry` : geometry.
- `isFC::Bool` : boolean indicating if the high-order moments are fully coupled.

# Output Argument(s)
- `schemes::Vector{String}` : scheme types.
- `𝒪::Vector{Int64}` : order of the schemes.
- `Nm::Vector{Int64}` : numbers of moments.

"""
function get_schemes(this::CPM,geometry::Geometry,isFC::Bool)
    schemes = ["DD","DD","DD","DD"]
    𝒪 = [1,1,1,1]
    Nm = [1,1,1,1,1]
    return schemes,𝒪,Nm
end

"""
    get_is_full_coupling(this::CPM)

Get whether the high-order moments are fully coupled or not.

# Input Argument(s)
- `this::CPM` : collision-probability method.

# Output Argument(s)
- `isFC::Bool` : boolean indicating if the high-order moments are fully coupled or not.

"""
function get_is_full_coupling(this::CPM)
    return this.isFC
end

"""
    get_convergence_criterion(this::CPM)

Get the convergence criterion for the group iterations.

# Input Argument(s)
- `this::CPM` : collision-probability method.

# Output Argument(s)
- `convergence_criterion::Float64` : convergence criterion.

"""
function get_convergence_criterion(this::CPM)
    return this.convergence_criterion
end

"""
    get_maximum_iteration(this::CPM)

Get the maximum number of group iterations.

# Input Argument(s)
- `this::CPM` : collision-probability method.

# Output Argument(s)
- `maximum_iteration::Int64` : maximum number of group iterations.

"""
function get_maximum_iteration(this::CPM)
    return this.maximum_iteration
end

"""
    set_acceleration(this::CPM,acceleration::String,parameter::Int64=0)

To set the acceleration method for the in-group source iteration of the sweeping mode (unused by
the direct `"global"` solve).

# Input Argument(s)
- `this::CPM` : collision-probability method.
- `acceleration::String` : acceleration method, which takes the following values:
    - `acceleration = "none"` : plain source iteration.
    - `acceleration = "livolant"` : Livolant two-point extrapolation.
    - `acceleration = "anderson"` : depth-`m` Anderson acceleration.
    - `acceleration = "gmres"` : matrix-free restarted GMRES.
    - `acceleration = "bicgstab"` : matrix-free BiCGStab.
- `parameter::Int64` : optional tuning parameter (GMRES restart or Anderson depth; `0` keeps the
   current default).

# Output Argument(s)
N/A

"""
function set_acceleration(this::CPM,acceleration::String,parameter::Int64=0)
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
    get_acceleration(this::CPM)

Get the acceleration method (used only by the sweeping mode's in-group source iteration).

# Input Argument(s)
- `this::CPM` : collision-probability method.

# Output Argument(s)
- `acceleration::String` : acceleration method.

"""
function get_acceleration(this::CPM)
    return this.acceleration
end

"""
    get_gmres_restart(this::CPM)

Get the GMRES restart parameter (sweeping mode).
"""
function get_gmres_restart(this::CPM)
    return this.gmres_restart
end

"""
    get_anderson_depth(this::CPM)

Get the Anderson acceleration memory depth (sweeping mode).
"""
function get_anderson_depth(this::CPM)
    return this.anderson_depth
end
