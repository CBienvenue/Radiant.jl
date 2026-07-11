"""
    Surface_Source

Structure used to define a directionnal boundary source and its properties.

# Mandatory field(s)
- `name::String` : name (or identifier) of the Surface_Source structure.
- `particle::Particle` : type of particle emitted.
- `energy_group::Int64` : energy group index in which the particle are emitted.
- `direction::Vector{Float64}` : direction cosine.
- `location::String` : boundary at which the source is located.
- `boundaries::Vector{Float64}` : boundaries of the source along each axis [in cm].

# Optional field(s) - with default values
- `intensity::Float64=1.0` : intensity [# particles/cm⁽ᴺ⁻¹⁾, where N is the geometry dimension].
- `angular_moments::Vector{Float64}` : half-range moments of the incident angular flux
  (alternative to `direction` for a distributed angular source).
- `fcs::Bool=false` : whether the source is treated by the first-collision method
  (the uncollided flux is computed analytically and split from the transport solve);
  the default `false` keeps the ordinary boundary treatment.

"""
mutable struct Surface_Source

    # Variable(s)
    name                       ::Union{Missing,String}
    particle                   ::Union{Missing,Particle}
    intensity                  ::Union{Missing,Float64}
    energy_group               ::Union{Missing,Int64}
    direction                  ::Union{Missing,Vector{Float64}}
    angular_moments            ::Union{Missing,Vector{Float64}}
    location                   ::Union{Missing,String}
    boundaries                 ::Dict{String,Vector{Float64}}
    is_build                   ::Bool
    surface_sources            ::Union{Missing,Vector{Array{Float64}}}
    normalization_factor       ::Float64
    legendre_order             ::Int64
    fcs                        ::Bool
    uncollided_model           ::String
    ray_refinement             ::Int64
    uncollided_angular_order   ::Tuple{Int64,Int64}

    # Constructor(s)
    function Surface_Source()

        this = new()

        this.name = missing
        this.particle = missing
        this.intensity = 1.0
        this.energy_group = missing
        this.direction = missing
        this.angular_moments = missing
        this.location = missing
        this.boundaries = Dict{String,Vector{Float64}}()
        this.is_build = false
        this.normalization_factor = 0.0
        this.surface_sources = missing
        this.legendre_order = 64
        this.fcs = false
        this.uncollided_model = "goudsmit-saunderson"
        this.ray_refinement = 1
        this.uncollided_angular_order = (16,32)

        return this
    end
end

# Method(s)
"""
    set_particle(this::Surface_Source,particle::Particle)

To define the source particle.

# Input Argument(s)
- `this::Surface_Source` : surface source.
- `particle::String` : particle.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ss = Surface_Source()
julia> ss.set_particle(electron)
```
"""
function set_particle(this::Surface_Source,particle::Particle)
    this.particle = particle
end

"""
    set_intensity(this::Surface_Source,intensity::Real)

To define the intensity of the source.

# Input Argument(s)
- `this::Surface_Source` : surface source.
- `intensity::Float64` : intensity [# particles/cmᴺ, where N is the geometry dimension]

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ss = Surface_Source()
julia> ss.set_intensity(100)
```
"""
function set_intensity(this::Surface_Source,intensity::Real)
    if intensity < 0 error("The intensity should be greater than 0.") end
    this.intensity = intensity
end

"""
    set_energy_group(this::Surface_Source,energy_group::Int64)

To define the energy of the source by setting the energy group in which they are produced.

# Input Argument(s)
- `this::Surface_Source` : surface source.
- `energy_group::Int64` : energy group index.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ss = Surface_Source()
julia> ss.set_energy_group(1)
```
"""
function set_energy_group(this::Surface_Source,energy_group::Int64)
    if energy_group ≤ 0 error("The energy group number should be greater than 0.") end
    this.energy_group = energy_group
end

"""
    set_direction(this::Surface_Source,direction::Vector{Float64})

To set a direction of the source.

# Input Argument(s)
- `this::Surface_Source` : surface source.
- `direction::Vector{Float64}` : director cosines [μ,η,ξ]

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ss = Surface_Source()
julia> ss.set_direction([1.0,0.0,0.0])
```
"""
function set_direction(this::Surface_Source,direction::Vector{Float64})
    if length(direction) != 3 error("Three director cosines has to be provided.") end
    if abs(sum(direction.^2)-1) > 1e-3 error("The sum of the squared three director cosines is not equal to 1.") end
    if ~ismissing(this.angular_moments) error("A surface source is either monodirectional (set_direction) or distributed (set_angular_moments), not both.") end
    this.direction = direction
end

"""
    set_angular_moments(this::Surface_Source,angular_moments::Vector{Float64})

To set a distributed angular source, given by the half-range moments
`g[p+1] = ∫₀¹ R̄ₚ(μ̂) ψ_inc(μ̂) dμ̂` of the incident angular flux in the orthonormal
half-range Legendre basis `R̄ₚ(μ̂) = √(2p+1)·Pₚ(2μ̂-1)` (`μ̂ = |Ω⋅n̂|`, the basis of
`half_range_legendre_polynomials_up_to_L`). Alternative to `set_direction` and
mutually exclusive with it; available in 1D only.

# Input Argument(s)
- `this::Surface_Source` : surface source.
- `angular_moments::Vector{Float64}` : half-range moments (order 0 to L_src).

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ss = Surface_Source()
julia> ss.set_angular_moments([0.5,1/(2*sqrt(3))])  # ψ_inc(μ̂) = μ̂
```
"""
function set_angular_moments(this::Surface_Source,angular_moments::Vector{Float64})
    if length(angular_moments) == 0 error("At least the zeroth half-range moment has to be provided.") end
    if ~ismissing(this.direction) error("A surface source is either monodirectional (set_direction) or distributed (set_angular_moments), not both.") end
    this.angular_moments = angular_moments
end

"""
    set_fcs(this::Surface_Source,fcs::Bool)

To enable the first-collision treatment of the surface source. When `false`
(default), the source enters as an ordinary incoming boundary condition through
its truncated half-range moment expansion. When `true`, the uncollided flux is
computed analytically outside the solver, which then only transports the smooth
first-collision volume source (`first_collision_source!`). Available in 1D/2D/3D
Cartesian geometry with void boundaries, for the GN and SN solvers (BTE/BFP/BCSD).

# Input Argument(s)
- `this::Surface_Source` : surface source.
- `fcs::Bool` : whether to use the first-collision treatment.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ss = Surface_Source()
julia> ss.set_fcs(true)
```
"""
function set_fcs(this::Surface_Source,fcs::Bool)
    this.fcs = fcs
end

"""
    set_uncollided_model(this::Surface_Source,uncollided_model::String)

To set the model of the uncollided component in the first-collision treatment
(BFP solvers only; for BTE/BCSD both models coincide):
- `"goudsmit-saunderson"` (default) : the uncollided column carries the exact
  Fokker-Planck angular redistribution along its pathlength (order-ℓ moments
  damped as e^(-ℓ(ℓ+1)·∫T ds), mean-depth advance ⟨μ⟩ = μ₀·e^(-2∫T ds)); in 2D/3D
  each deposit is spread by a Gaussian kernel carrying the exact Lewis second
  moments (longitudinal straggling and transverse spread), while in 1D the
  column deposits at its mean depth.
- `"straight"` : the column travels straight with no Fokker-Planck broadening
  (δ-pure model), kept for comparison; it displaces the dose of FP-dominated
  problems toward depth.

See TR-09 for the derivation.

# Input Argument(s)
- `this::Surface_Source` : surface source.
- `uncollided_model::String` : model of the uncollided component.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ss = Surface_Source()
julia> ss.set_uncollided_model("straight")
```
"""
function set_uncollided_model(this::Surface_Source,uncollided_model::String)
    if lowercase(uncollided_model) ∉ ["goudsmit-saunderson","straight"] error("Unknown uncollided model.") end
    this.uncollided_model = lowercase(uncollided_model)
end

"""
    get_uncollided_model(this::Surface_Source)

Get the model of the uncollided component in the first-collision treatment.

# Input Argument(s)
- `this::Surface_Source` : surface source.

# Output Argument(s)
- `uncollided_model::String` : model (`"goudsmit-saunderson"` or `"straight"`).

"""
function get_uncollided_model(this::Surface_Source)
    return this.uncollided_model
end

"""
    set_ray_refinement(this::Surface_Source,ray_refinement::Int64)

To set the ray refinement of the multi-dimensional first-collision treatment: each
illuminated source face cell is subdivided into `ray_refinement^(Ndims-1)`
sub-elements (launch points). `1` (default) is exact for axis-aligned normal
beams; oblique beams and distributed sources converge as the refinement grows.

# Input Argument(s)
- `this::Surface_Source` : surface source.
- `ray_refinement::Int64` : number of sub-elements per face cell and tangent axis.

# Output Argument(s)
N/A

"""
function set_ray_refinement(this::Surface_Source,ray_refinement::Int64)
    if ray_refinement < 1 error("Ray refinement must be at least 1.") end
    this.ray_refinement = ray_refinement
end

"""
    get_ray_refinement(this::Surface_Source)

Get the ray refinement of the multi-dimensional first-collision treatment.

"""
function get_ray_refinement(this::Surface_Source)
    return this.ray_refinement
end

"""
    set_uncollided_angular_order(this::Surface_Source,Nμ::Int64,Nφ::Int64)

To set the incoming-hemisphere quadrature order used by the multi-dimensional
first-collision treatment of a distributed source: `Nμ` Gauss nodes in the polar
cosine and `Nφ` uniform azimuthal nodes about the face normal.

# Input Argument(s)
- `this::Surface_Source` : surface source.
- `Nμ::Int64`, `Nφ::Int64` : polar and azimuthal quadrature orders.

# Output Argument(s)
N/A

"""
function set_uncollided_angular_order(this::Surface_Source,Nμ::Int64,Nφ::Int64)
    if Nμ < 1 || Nφ < 1 error("Angular quadrature orders must be at least 1.") end
    this.uncollided_angular_order = (Nμ,Nφ)
end

"""
    get_uncollided_angular_order(this::Surface_Source)

Get the incoming-hemisphere quadrature order (`Nμ`,`Nφ`) of the multi-dimensional
first-collision treatment.

"""
function get_uncollided_angular_order(this::Surface_Source)
    return this.uncollided_angular_order
end

"""
    set_location(this::Surface_Source,location::String)

To set the location of the surface source.

# Input Argument(s)
- `this::Surface_Source` : surface source.
- `location::String` : boundary on which the surface source is, which can takes the following value:
    - `boundary = "x-"` : the lower bound along x-axis
    - `boundary = "x+"` : the upper bound along x-axis
    - `boundary = "y-"` : the lower bound along y-axis
    - `boundary = "y+"` : the upper bound along y-axis
    - `boundary = "z-"` : the lower bound along z-axis
    - `boundary = "z+"` : the upper bound along z-axis

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ss = Surface_Source()
julia> ss.set_location("x-")
```
"""
function set_location(this::Surface_Source,location::String)
    if uppercase(location) ∉ ["X-","X+","Y-","Y+","Z-","Z+"] error("Unknown location.") end
    this.location = location
end

"""
    set_boundaries(this::Surface_Source,axis::String,boundaries::Vector{Float64})

To define the boundaries of the source along the specified axis.

# Input Argument(s)
- `this::Surface_Source` : surface source.
- `axis::String` : axis along which the boundaries are defined.
- `boundaries::Vector{Float64}` : boundaries of the source in accending order [in cm]

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ss = Surface_Source()
julia> ss.set_boundaries("x",[1.0,3.0])
```
"""
function set_boundaries(this::Surface_Source,axis::String,boundaries::Vector{Float64})
    if lowercase(axis) ∉ ["x","y","z"] error("Unknown axis.") end
    if length(boundaries) != 2 error("The two boundary side should be provided.") end
    this.boundaries[axis] = boundaries
end

"""
    set_legendre_order(this::Surface_Source,legendre_order::Int64)

To define the Legendre order of the polynomial expansion of the boundary flux.

# Input Argument(s)
- `this::Surface_Source` : surface source.
- `legendre_order::Int64` : legendre order.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> ss = Surface_Source()
julia> ss.set_legendre_order(3)
```
"""
function set_legendre_order(this::Surface_Source,legendre_order::Int64)
    if legendre_order < 0 error("Legendre order should be 0 or greater.") end
    this.legendre_order = legendre_order
end

"""
    get_particle(this::Surface_Source)

Get the source particle.

# Input Argument(s)
- `this::Surface_Source` : surface source.

# Output Argument(s)
- `particle::Particle` : particle.

"""
function get_particle(this::Surface_Source)
    return this.particle
end

"""
    get_normalization_factor(this::Surface_Source)

Get the source normalization factor.

# Input Argument(s)
- `this::Surface_Source` : surface source.

# Output Argument(s)
- `normalization_factor::Float64` : source normalization factor.

"""
function get_normalization_factor(this::Surface_Source)
    return this.normalization_factor
end

"""
    get_legendre_order(this::Surface_Source)

To get the Legendre order of the polynomial expansion of the boundary flux.

# Input Argument(s)
- `this::Surface_Source` : surface source.

# Output Argument(s)
- `legendre_order::Int64` : legendre order.

"""
function get_legendre_order(this::Surface_Source)
    return this.legendre_order
end

"""
    get_angular_moments(this::Surface_Source)

Get the half-range moments of the incident angular flux (or `missing` if the
source is monodirectional).

# Input Argument(s)
- `this::Surface_Source` : surface source.

# Output Argument(s)
- `angular_moments::Union{Missing,Vector{Float64}}` : half-range moments.

"""
function get_angular_moments(this::Surface_Source)
    return this.angular_moments
end

"""
    get_fcs(this::Surface_Source)

Get whether the surface source is treated by the first-collision method.

# Input Argument(s)
- `this::Surface_Source` : surface source.

# Output Argument(s)
- `fcs::Bool` : whether the first-collision treatment is used.

"""
function get_fcs(this::Surface_Source)
    return this.fcs
end