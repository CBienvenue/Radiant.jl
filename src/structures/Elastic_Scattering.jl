"""
    IsotopeElasticDB

Stores ENDF elastic-scattering data for one target isotope.

# Field(s)
- `energies::Vector{Float64}` : incident energy grid [MeV].
- `mu_grid::Vector{Float64}` : cosine grid for tabulated angular distributions.
- `dcs::Vector{Vector{Float64}}` : tabulated differential cross-section values.
- `ltp::Int` : ENDF LAW=5 representation type.
- `identical::Bool` : whether projectile and target are identical particles.
- `Z1::Int` : projectile charge number.
- `Z2::Int` : target charge number.
- `spin_s::Float64` : spin factor used for identical-particle Coulomb terms.
- `m_mev::Float64` : projectile rest energy [MeV].
- `M_mev::Float64` : target rest energy [MeV].
- `mf6_nbt::Vector{Int}` : MF6 TAB2 interpolation breakpoints.
- `mf6_int::Vector{Int}` : MF6 TAB2 interpolation laws.
- `ba_tables_ltp1` : coefficient tables for LTP=1 data.
- `c_tables_ltp2` : coefficient tables for LTP=2 data.
- `mupni_tables` : tabulated cosine/probability data for tabulated LTP forms.
- `sigmaNI_E` : non-interference cross-section energy grid.
- `sigmaNI_S` : non-interference cross-section values.
- `sigmaNI_nbt` : non-interference interpolation breakpoints.
- `sigmaNI_int` : non-interference interpolation laws.
"""
struct IsotopeElasticDB
    energies::Vector{Float64}
    mu_grid::Vector{Float64}
    dcs::Vector{Vector{Float64}}
    ltp::Int
    identical::Bool
    Z1::Int
    Z2::Int
    spin_s::Float64
    m_mev::Float64
    M_mev::Float64
    mf6_nbt::Vector{Int}
    mf6_int::Vector{Int}
    ba_tables_ltp1::Union{Nothing,Dict{Float64,Tuple{Int,Vector{Float64},Vector{ComplexF64}}}}
    c_tables_ltp2::Union{Nothing,Dict{Float64,Tuple{Int,Vector{Float64}}}}
    mupni_tables::Union{Nothing,Dict{Float64,Vector{Float64}},Dict{Float64,Tuple{Vector{Float64},Vector{Float64}}}}
    sigmaNI_E::Union{Nothing,Vector{Float64}}
    sigmaNI_S::Union{Nothing,Vector{Float64}}
    sigmaNI_nbt::Union{Nothing,Vector{Int}}
    sigmaNI_int::Union{Nothing,Vector{Int}}
end

"""
    ElasticScatteringENDFDB

Stores the elastic-scattering ENDF database for one incoming particle type.

# Field(s)
- `dcs_source::String` : DCS representation used by cached evaluations.
- `coulomb_kinematics::String` : Coulomb kinematics convention.
- `screening::Bool` : whether Coulomb screening is enabled.
- `isotopes::Dict{Tuple{Int,Int},IsotopeElasticDB}` : isotope data indexed by `(Z, A)`.
"""
struct ElasticScatteringENDFDB
    dcs_source::String
    coulomb_kinematics::String
    screening::Bool
    isotopes::Dict{Tuple{Int,Int},IsotopeElasticDB}
end

"""
    ElasticDCSCache

Stores interpolated quantities used to evaluate elastic angular distributions at one
incident energy for one isotope.

# Field(s)
- `ltp::Int` : ENDF LAW=5 representation type.
- `E_in_mev::Float64` : incident kinetic energy [MeV].
- `k_coulomb::Float64` : Coulomb normalization factor.
- `eta_coulomb::Float64` : Coulomb Sommerfeld parameter.
- `delta_coulomb::Float64` : Coulomb phase-shift correction.
- `mu_grid::Vector{Float64}` : cosine grid for tabulated data.
- `sigma_grid::Vector{Float64}` : differential cross-section values on `mu_grid`.
- `pni_grid::Vector{Float64}` : probability values on `mu_grid`.
- `sigmaNI::Float64` : non-interference cross-section at this energy.
- `NL::Int` : highest Legendre order or number of tabulated cosines.
- `a::Vector{ComplexF64}` : complex coefficients for LTP=1 reconstruction.
- `b::Vector{Float64}` : real coefficients for LTP=1 reconstruction.
- `c::Vector{Float64}` : coefficients for LTP=2 reconstruction.
- `spin_s::Float64` : spin factor for identical-particle terms.
- `m_mev::Float64` : projectile rest energy [MeV].
- `M_mev::Float64` : target rest energy [MeV].
- `identical::Bool` : whether projectile and target are identical particles.
"""
struct ElasticDCSCache
    ltp::Int
    E_in_mev::Float64
    k_coulomb::Float64
    eta_coulomb::Float64
    delta_coulomb::Float64
    mu_grid::Vector{Float64}
    sigma_grid::Vector{Float64}
    pni_grid::Vector{Float64}
    sigmaNI::Float64
    NL::Int
    a::Vector{ComplexF64}
    b::Vector{Float64}
    c::Vector{Float64}
    spin_s::Float64
    m_mev::Float64
    M_mev::Float64
    identical::Bool
end

"""
    Elastic_Scattering

Structure used to define parameters for production of multigroup elastic-scattering
cross sections.

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{Type,Type},Vector{String}} = Dict((Proton,Proton) => ["S", "P"])` : interaction processes, where `"S"` is scattered-projectile production and `"P"` is recoil production.
- `scattering_model::String = "BFP"` : elastic-scattering transport model.
- `library::String = "TENDL2023"` : elastic ENDF library name.
- `data_path::String = "../../data"` : root directory for elastic ENDF data.
- `integrate_total_cs::String = "analytical"` : total-cross-section integration method.
- `integrate_fp::String = "tanh-sinh"` : Fokker-Planck integration method.
- `dcs_source::String = "coefficients"` : source used for DCS reconstruction.
- `coulomb_kinematics::String = "standard"` : Coulomb kinematics convention.
- `screening::Bool = true` : whether Coulomb screening is enabled.
"""
mutable struct Elastic_Scattering <: Interaction

    # Variable(s)
    name::String
    incoming_particle::Vector{Type}
    interaction_particles::Vector{Type}
    interaction_types::Dict{Tuple{Type,Type},Vector{String}}
    is_CSD::Bool
    is_AFP::Bool
    is_AFP_decomposition::Bool
    is_elastic::Bool
    is_subshells_dependant::Bool
    scattering_model::String
    library::String
    data_path::String
    endf_db::Dict{Type,ElasticScatteringENDFDB}
    integrate_total_cs::String
    integrate_fp::String
    dcs_source::String
    coulomb_kinematics::String
    screening::Bool
    dcs_cache::Dict{Tuple{Type,Int,Int,String,String,Bool},ElasticDCSCache}
    fp_tol::Float64
    fp_ts_h::Float64
    fp_ts_n::Int64
    fp_ts_max_iter::Int64
    fp_gk_split_points::Vector{Float64}

    # Constructor(s)
    function Elastic_Scattering()
        this = new()
        this.name = "Elastic_Scattering"
        this.interaction_types = Dict((Proton,Proton) => ["S","P"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.is_CSD = true
        this.is_AFP = true
        this.is_AFP_decomposition = false
        this.is_elastic = false
        this.is_subshells_dependant = false
        this.set_scattering_model("BFP")
        this.set_library("TENDL2023")
        this.dcs_cache = Dict{Tuple{Type,Int,Int,String,String,Bool},ElasticDCSCache}()
        this.set_data_path("../../data")
        this.set_endf_db(Dict{Type,ElasticScatteringENDFDB}())
        this.set_integrate_total_cs("analytical")
        this.set_integrate_fp("tanh-sinh")
        this.set_dcs_source("coefficients")
        this.set_coulomb_kinematics("standard")
        this.set_screening(true)
        this.set_fp_tol(1.0e-10)
        this.set_fp_tanh_sinh(0.03, 240, 12)
        this.set_fp_quadgk_split_points([0.99, 0.9999, 0.999999])
        return this
    end
end

const QUADGK_PKGID = Base.PkgId(Base.UUID("1fd47b50-473d-5c70-9696-f719f8f3bcdc"), "QuadGK")
const HAS_QUADGK = try
    Base.require(QUADGK_PKGID)
    true
catch
    false
end

"""
    require_quadgk()

Loads the QuadGK package object used by optional Gauss-Kronrod integrations.

# Output Argument(s)
- `quadgk` : loaded QuadGK module object.
"""
function require_quadgk()
    try
        return Base.require(QUADGK_PKGID)
    catch
        error("QuadGK is not available. Add it with: import Pkg; Pkg.add(\"QuadGK\")")
    end
end

const M_E_C2_MEV = 0.510999
const BARN_TO_CM2 = 1.0e-24

# Method(s)
"""
    initialize(this::Elastic_Scattering, particles, isotopes::Vector{Tuple{Int,Int}})

Initializes the elastic-scattering ENDF databases required by the incoming particles.
Already initialized particle databases are reused.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.
- `particles` : particle objects used in the simulation.
- `isotopes::Vector{Tuple{Int,Int}}` : target isotope pairs `(Z, A)` required by the materials.

# Output Argument(s)
N/A
"""
function initialize(this::Elastic_Scattering, particles, isotopes::Vector{Tuple{Int,Int}})
    endf_db = this.get_endf_db()
    for particle in particles
        ptype = get_type(particle)
        if ptype ∈ this.get_in_particles()
            if !haskey(endf_db, ptype)
                endf_db[ptype] = init_elastic_scattering_endf_db(
                    this.get_library(), isotopes, particle;
                    data_root=this.get_data_path(),
                    analytical=this.get_integrate_total_cs() == "analytical",
                    dcs_source=this.get_dcs_source(),
                    coulomb_kinematics=this.get_coulomb_kinematics(),
                    screening=this.get_screening())
            end
        end
    end
    this.set_endf_db(endf_db)
    return nothing
end

"""
    set_interaction_types(this::Elastic_Scattering,interaction_types::Dict{Tuple{DataType,DataType},Vector{String}})

To define the interaction types for elastic scattering processes.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}}` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `(Proton,Proton) => ["S"]` : scattering of incident protons
    - `(Proton,Proton) => ["P"]` : production of protons by incident protons on hydrogen (handled by target selection) 

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_scattering = Elastic_Scattering()
julia> elastic_scattering.set_interaction_types( Dict((Electron,Electron) => ["S"]) ) # No knock-on electrons.
```
"""
function set_interaction_types(this::Elastic_Scattering,interaction_types)
    this.interaction_types = interaction_types
end

"""
    scattering_model(this::Elastic_Scattering,scattering_model::String)

To define the solver for elastic scattering.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.
- `solver::String` : solver for elastic scattering, which can be:
    - `BFP` : Boltzmann Fokker-Planck solver.
    - `FP` : Fokker-Planck solver.
    - `BTE` : Boltzmann solver.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_scattering = Elastic_Scattering()
julia> elastic_scattering.set_scattering_model("FP")
```
"""
function set_scattering_model(this::Elastic_Scattering,scattering_model::String)
    if uppercase(scattering_model) ∉ ["BFP","FP","BTE"] error("Unknown scattering model (should be BFP, FP or BTE).") end
    this.scattering_model = uppercase(scattering_model)
end

"""
    set_library(this::Elastic_Scattering,library::String)

Define the elastic-scattering database library name (e.g., TENDL2023).

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.
- `library::String` : library name.

# Output Argument(s)
N/A
"""
function set_library(this::Elastic_Scattering,library::String)
    this.library = library
end

"""
    set_data_path(this::Elastic_Scattering,data_path::String)

Define the root directory containing the elastic-scattering ENDF data files.
Relative paths are resolved from the location of `Elastic_Scattering.jl`.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.
- `data_path::String` : data root path. The default `../../data` resolves to the
  package data directory.

# Output Argument(s)
N/A
"""
function set_data_path(this::Elastic_Scattering,data_path::String)
    path = isabspath(data_path) ? normpath(data_path) : normpath(joinpath(@__DIR__, data_path))
    this.data_path = path
    this.endf_db = Dict{Type,ElasticScatteringENDFDB}()
    empty!(this.dcs_cache)
end

"""
    set_integrate_total_cs(this::Elastic_Scattering,method::String)

Select the integration method for the total elastic cross section.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.
- `method::String` : "numerical-gl", "numerical-gk", "gauss-kronrod", "numerical-ts", or "analytical".

# Output Argument(s)
N/A
"""
function set_integrate_total_cs(this::Elastic_Scattering,method::String)
    method_lc = lowercase(method)
    if method_lc == "gauss-kronrod"
        method_lc = "numerical-gk"
    end
    if method_lc ∉ ["numerical-gl", "numerical-gk", "numerical-ts", "analytical"]
        error("Unknown integrate_total_cs method (should be numerical-gl, numerical-gk, numerical-ts, or analytical).")
    end
    if method_lc == "numerical-gk" && !HAS_QUADGK
        error("QuadGK is not available. Add it with: import Pkg; Pkg.add(\"QuadGK\")")
    end
    this.integrate_total_cs = method_lc
end

"""
    set_integrate_fp(this::Elastic_Scattering,method::String)

Select the integration method for Fokker-Planck terms (stopping power and momentum transfer).

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.
- `method::String` : "gauss-kronrod" or "tanh-sinh".

# Output Argument(s)
N/A
"""
function set_integrate_fp(this::Elastic_Scattering,method::String)
    method_lc = lowercase(method)
    if method_lc ∉ ["gauss-kronrod", "tanh-sinh"]
        error("Unknown integrate_fp method (should be gauss-kronrod or tanh-sinh).")
    end
    if method_lc == "gauss-kronrod" && !HAS_QUADGK
        error("QuadGK is not available. Add it with: import Pkg; Pkg.add(\"QuadGK\")")
    end
    this.integrate_fp = method_lc
end

"""
    get_integrate_fp(this::Elastic_Scattering)

Get the integration method for Fokker-Planck terms.
"""
function get_integrate_fp(this::Elastic_Scattering)
    return this.integrate_fp
end

"""
    set_dcs_source(this::Elastic_Scattering,source::String)

Set whether elastic DCS evaluations use the prebuilt database rows or the
coefficient-based reconstruction.
"""
function set_dcs_source(this::Elastic_Scattering,source::String)
    source_lc = lowercase(source)
    if source_lc ∉ ["database", "coefficients"]
        error("Unknown elastic DCS source (should be database or coefficients).")
    end
    this.dcs_source = source_lc
    empty!(this.dcs_cache)
end

"""
    get_dcs_source(this::Elastic_Scattering)

Gets the source used for elastic DCS evaluations.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.

# Output Argument(s)
- `source::String` : "database" or "coefficients".
"""
function get_dcs_source(this::Elastic_Scattering)
    return this.dcs_source
end

"""
    set_coulomb_kinematics(this::Elastic_Scattering,kinematics::String)

Set the Coulomb DCS kinematics convention ("standard" or "relativistic").
"""
function set_coulomb_kinematics(this::Elastic_Scattering,kinematics::String)
    kinematics_lc = lowercase(kinematics)
    if kinematics_lc ∉ ["standard", "relativistic"]
        error("Unknown Coulomb kinematics (should be standard or relativistic).")
    end
    this.coulomb_kinematics = kinematics_lc
    this.endf_db = Dict{Type,ElasticScatteringENDFDB}()
    empty!(this.dcs_cache)
end

"""
    get_coulomb_kinematics(this::Elastic_Scattering)

Gets the Coulomb DCS kinematics convention.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.

# Output Argument(s)
- `kinematics::String` : "standard" or "relativistic".
"""
function get_coulomb_kinematics(this::Elastic_Scattering)
    return this.coulomb_kinematics
end

"""
    set_screening(this::Elastic_Scattering,screening::Bool)

Set whether Coulomb DCS evaluations use Moliere screening.
"""
function set_screening(this::Elastic_Scattering,screening::Bool)
    this.screening = screening
    this.endf_db = Dict{Type,ElasticScatteringENDFDB}()
    empty!(this.dcs_cache)
end

"""
    get_screening(this::Elastic_Scattering)

Gets whether Coulomb screening is enabled.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.

# Output Argument(s)
- `screening::Bool` : screening flag.
"""
function get_screening(this::Elastic_Scattering)
    return this.screening
end

"""
    clear_dcs_cache(this::Elastic_Scattering)

Clears cached elastic DCS interpolation objects.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.

# Output Argument(s)
N/A
"""
function clear_dcs_cache(this::Elastic_Scattering)
    empty!(this.dcs_cache)
end

"""
    get_dcs_cache(this::Elastic_Scattering, ptype::Type, db::ElasticScatteringENDFDB,
    Z::Int, A::Int, E_in_mev::Float64; source::Union{Nothing,String}=nothing)

Gets or builds the cached DCS interpolation data for one particle, isotope, and incident
energy.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.
- `ptype::Type` : incoming particle type.
- `db::ElasticScatteringENDFDB` : ENDF database for the incoming particle type.
- `Z::Int` : target atomic number.
- `A::Int` : target isotope mass number.
- `E_in_mev::Float64` : incoming particle kinetic energy [MeV].
- `source::Union{Nothing,String}=nothing` : optional DCS source consistency check.

# Output Argument(s)
- `cache::ElasticDCSCache` : cached DCS interpolation data.
"""
function get_dcs_cache(this::Elastic_Scattering,ptype::Type,db::ElasticScatteringENDFDB,
                       Z::Int,A::Int,E_in_mev::Float64; source::Union{Nothing,String}=nothing)
    if !isnothing(source) && lowercase(source) != db.dcs_source
        error("Requested DCS source $(source), but the ENDF database was initialized with dcs_source=$(db.dcs_source).")
    end
    source_key = db.dcs_source
    coulomb_kinematics = db.coulomb_kinematics
    screening = db.screening
    key = (ptype, Z, A, source_key, coulomb_kinematics, screening)
    if !haskey(this.dcs_cache, key) || this.dcs_cache[key].E_in_mev != E_in_mev
        this.dcs_cache[key] = prepare_elastic_dsig_cache(db, Z, A, E_in_mev)
    end
    return this.dcs_cache[key]
end

"""
    set_fp_tol(this::Elastic_Scattering,tol::Float64)

Set the relative tolerance used for Fokker-Planck integrals.
"""
function set_fp_tol(this::Elastic_Scattering,tol::Float64)
    if tol <= 0.0
        error("Fokker-Planck integration tolerance must be positive.")
    end
    this.fp_tol = tol
end

"""
    get_fp_tol(this::Elastic_Scattering)

Gets the relative tolerance used for Fokker-Planck integrations.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.

# Output Argument(s)
- `tol::Float64` : relative integration tolerance.
"""
function get_fp_tol(this::Elastic_Scattering)
    return this.fp_tol
end

"""
    set_fp_tanh_sinh(this::Elastic_Scattering, h::Float64=0.03, n::Int=240, max_iter::Int=12)

Set tanh-sinh quadrature controls used for Fokker-Planck integrals.
"""
function set_fp_tanh_sinh(this::Elastic_Scattering, h::Float64=0.03, n::Int=240, max_iter::Int=12)
    if h <= 0.0
        error("Tanh-sinh step h must be positive.")
    end
    if n < 1
        error("Tanh-sinh node count n must be at least 1.")
    end
    if max_iter < 0
        error("Tanh-sinh max_iter must be nonnegative.")
    end
    this.fp_ts_h = h
    this.fp_ts_n = n
    this.fp_ts_max_iter = max_iter
end

"""
    get_fp_tanh_sinh(this::Elastic_Scattering)

Gets the tanh-sinh quadrature controls used for Fokker-Planck integrations.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.

# Output Argument(s)
- `controls::NamedTuple` : named tuple with `h`, `n`, and `max_iter`.
"""
function get_fp_tanh_sinh(this::Elastic_Scattering)
    return (h=this.fp_ts_h, n=this.fp_ts_n, max_iter=this.fp_ts_max_iter)
end

"""
    set_fp_quadgk_split_points(this::Elastic_Scattering,split_points::Vector{Float64})

Set the intermediate mu points used to split QuadGK Fokker-Planck integrals.
Only points strictly inside (0,1) are allowed.
"""
function set_fp_quadgk_split_points(this::Elastic_Scattering,split_points::Vector{Float64})
    if any(p -> !(0.0 < p < 1.0), split_points)
        error("QuadGK split points must lie strictly inside (0,1).")
    end
    this.fp_gk_split_points = sort(unique(copy(split_points)))
end

"""
    get_fp_quadgk_split_points(this::Elastic_Scattering)

Gets the intermediate split points used by QuadGK Fokker-Planck integrations.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.

# Output Argument(s)
- `split_points::Vector{Float64}` : sorted split points inside `(0, 1)`.
"""
function get_fp_quadgk_split_points(this::Elastic_Scattering)
    return copy(this.fp_gk_split_points)
end

"""
    integrate_fp_interval(this::Elastic_Scattering, integrand::Function,
    a::Float64, b::Float64)

Integrates a Fokker-Planck integrand over one cosine interval using the configured
quadrature method.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.
- `integrand::Function` : scalar integrand as a function of cosine.
- `a::Float64` : lower integration bound.
- `b::Float64` : upper integration bound.

# Output Argument(s)
- `val::Float64` : integral over `[a, b]`.
"""
function integrate_fp_interval(this::Elastic_Scattering,integrand::Function,a::Float64,b::Float64)
    if a >= b
        return 0.0
    end

    tol = get_fp_tol(this)
    method = this.get_integrate_fp()
    if method == "tanh-sinh"
        ts = get_fp_tanh_sinh(this)
        return tanh_sinh_integral(integrand, a, b; h=ts.h, n=ts.n, tol=tol, max_iter=ts.max_iter)
    elseif method == "gauss-kronrod"
        split_points = get_fp_quadgk_split_points(this)
        quadgk = require_quadgk()
        segments = Float64[a]
        for p in split_points
            if a < p < b
                push!(segments, p)
            end
        end
        push!(segments, b)

        val = 0.0
        for i in 1:(length(segments) - 1)
            left = segments[i]
            right = segments[i + 1]
            if left < right
                seg_val, _ = quadgk.quadgk(integrand, left, right; rtol=tol)
                val += seg_val
            end
        end
        return val
    else
        error("Unknown integrate_fp method (should be gauss-kronrod or tanh-sinh).")
    end
end

"""
    get_integrate_total_cs(this::Elastic_Scattering)

Get the integration method for the total elastic cross section.
"""
function get_integrate_total_cs(this::Elastic_Scattering)
    return this.integrate_total_cs
end

"""
    get_data_path(this::Elastic_Scattering)

Get the root directory used to load elastic-scattering ENDF data.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.

# Output Argument(s)
- `data_path::String` : resolved data root path.
"""
function get_data_path(this::Elastic_Scattering)
    return this.data_path
end

"""
    get_library(this::Elastic_Scattering)

Get the elastic-scattering database library name.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.

# Output Argument(s)
- `library::String` : library name.
"""
function get_library(this::Elastic_Scattering)
    return this.library
end

"""
    set_endf_db(this::Elastic_Scattering,endf_db::Dict{Type,ElasticScatteringENDFDB})

Set the elastic-scattering ENDF database cache by particle id.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.
- `endf_db::Dict{Type,ElasticScatteringENDFDB}` : ENDF DB cache.

# Output Argument(s)
N/A
"""
function set_endf_db(this::Elastic_Scattering,endf_db::Dict{Type,ElasticScatteringENDFDB})
    this.endf_db = endf_db
    empty!(this.dcs_cache)
end

"""
    get_endf_db(this::Elastic_Scattering)

Get the elastic-scattering ENDF database cache by particle id.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.

# Output Argument(s)
- `endf_db::Dict{Type,ElasticScatteringENDFDB}` : ENDF DB cache.
"""
function get_endf_db(this::Elastic_Scattering)
    return this.endf_db
end

"""
    in_distribution(this::Elastic_Scattering)

Describe the energy discretization method for the incoming particle in the elastic scattering
interaction.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function in_distribution(this::Elastic_Scattering)
    is_dirac = false
    N = 1
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    bounds(this::Elastic_Scattering,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String,
    Ec::Float64,particle::Particle,M_target::Float64)

Gives the integration energy bounds for the outgoing particle for elastic scattering
interaction. 

# Input Argument(s)
- `this::Annihilation` : annihilation structure.
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `Ei::Float64` : energy of the incoming particle.
- `type::String` : type of interaction.
- `Ec::Float64` : energy cutoff between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.
- `M_target::Float64` : target isotope mass [MeV/c^2].

# Output Argument(s)
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `isSkip::Bool` : define if the integration is skipped or not.

"""
function bounds(this::Elastic_Scattering,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String,Ec::Float64,particle::Particle,M_target::Float64)
    # Relativistic maximum energy loss in lab frame
    m_mev = get_mass(particle)
    M_mev = M_target
    E_in_mev = Ei * M_E_C2_MEV
    rel_tol = 2e-3
    mu_cm = abs(M_mev - m_mev) / m_mev < rel_tol ? 0.0 : -1.0
    E1 = E_in_mev + m_mev
    p1 = sqrt(max(E1^2 - m_mev^2, 0.0))
    s = m_mev^2 + M_mev^2 + 2.0 * M_mev * E1
    sqrt_s = sqrt(s)
    E1_star = (s + m_mev^2 - M_mev^2) / (2.0 * sqrt_s)
    p_star = sqrt(max(E1_star^2 - m_mev^2, 0.0))
    beta_cm = p1 / (E1 + M_mev)
    gamma_cm = 1.0 / sqrt(1.0 - beta_cm^2)
    E1_prime = gamma_cm * (E1_star + beta_cm * p_star * mu_cm)
    E_out_mev = max(E1_prime - m_mev, 0.0)
    Emax = max(E_in_mev - E_out_mev, 0.0) / M_E_C2_MEV
    # Scattered heavy particle
    if type == "S"
        Ef⁻ = min(Ef⁻,Ec)
        Ef⁺ = max(Ef⁺,0.0,Ei-Emax)
    # Knock-on nucleus
    elseif type == "P"
        Ef⁻ = min(Ef⁻,Ei-Emax)
        Ef⁺ = max(Ef⁺,Ei-Ec,0.0)
    else
        error("Unknown type.")
    end
    if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    return Ef⁻,Ef⁺,isSkip
end

"""
    dcs(this::Elastic_Scattering,Z::Int64,A::Int64,L::Int64,Ei::Float64,Ef⁻::Float64,Ef⁺::Float64,
    type::String,particle::Particle)

Gives the Legendre moments of the scattering cross-sections for elastic scattering
interaction. 

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.
- `Z::Int64` : atomic number.
- `A::Int64` : isotope mass number.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : incoming particle energy.
- `Ef⁻::Float64` : outgoing group upper bound.
- `Ef⁺::Float64` : outgoing group lower bound.
- `type::String` : interaction type ("S" or "P").
- `particle::Particle` : incoming particle.

# Output Argument(s)
- `σl::Vector{Float64}` : Legendre moments of the scattering cross-sections.

"""
function dcs(this::Elastic_Scattering,Z::Int64,A::Int64,L::Int64,Ei::Float64,Ef⁻::Float64,Ef⁺::Float64,type::String,particle::Particle)

    #----
    # Initialization
    #----
    σl = zeros(L+1)
    ptype = get_type(particle)
    endf_db = this.get_endf_db()
    db = endf_db[ptype]

    m_mev = get_mass(particle)
    M_mev = get_mass(Z, A)
    E_in_mev = Ei * M_E_C2_MEV
    cache = this.get_dcs_cache(ptype, db, Z, A, E_in_mev)

    mu_hi = mu_cm_from_Eout_relativistic(E_in_mev, Ef⁻ * M_E_C2_MEV, m_mev, M_mev)
    mu_lo = mu_cm_from_Eout_relativistic(E_in_mev, Ef⁺ * M_E_C2_MEV, m_mev, M_mev)

    Nq = max(16, L + 1)
    μ_nodes, w = quadrature(Nq, "gauss-legendre")
    Δ = 0.5 * (mu_hi - mu_lo)
    μ_mid = 0.5 * (mu_hi + mu_lo)

    for l in range(0, L)
        acc = 0.0
        for k in eachindex(μ_nodes)
            mu_cm = Δ * μ_nodes[k] + μ_mid
            mu_lab = mu_lab_from_mu_cm(mu_cm, E_in_mev, m_mev, M_mev)
            mu_cm_dsig = (type == "P" && is_proton(particle) && Z == 1 && A == 1) ? -mu_cm : mu_cm
            sigma = elastic_dsig(cache, mu_cm_dsig)
            acc += w[k] * sigma * legendre_polynomials(l, mu_lab)
        end
        σl[l+1] = 2.0 * pi * acc * Δ
    end
    return σl .* BARN_TO_CM2
end

"""
    tcs(this::Elastic_Scattering,Ei::Float64,Ec::Float64,particle::Particle,Z::Int64;
    A::Union{Nothing,Vector{Int64}}=nothing,
    atpercentA::Union{Nothing,Vector{Float64}}=nothing)

Gives the total cross-section for elastic scattering interaction. 

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure. 
- `Ei::Float64` : incoming particle energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.
- `Z::Int64` : atomic number.
- `A::Union{Nothing,Vector{Int64}}` : isotope mass numbers for this element.
- `atpercentA::Union{Nothing,Vector{Float64}}` : isotope atomic fractions.

# Output Argument(s)
- `σt::Float64` : total cross-section.

"""
function tcs(this::Elastic_Scattering,Ei::Float64,Ec::Float64,particle::Particle,Z::Int64;
    A::Union{Nothing,Vector{Int64}}=nothing,
    atpercentA::Union{Nothing,Vector{Float64}}=nothing)

    if isnothing(A) || isnothing(atpercentA)
        pairs = isotopic_composition(Z)
        A = [p[1] for p in pairs]
        atpercentA = [p[2] for p in pairs]
    end

    ptype = get_type(particle)
    endf_db = this.get_endf_db()
    db = endf_db[ptype]
    E_in_mev = Ei * M_E_C2_MEV
    Ec_mev = Ec * M_E_C2_MEV

    method = this.get_integrate_total_cs()

    local cache_source::Union{Nothing,String}
    local integrate_tcs
    if method == "numerical-gl"
        cache_source = nothing
        integrate_tcs = cache -> integrate_elastic_scattering_numerical_GL(cache, Ec_mev)
    elseif method == "numerical-gk"
        cache_source = nothing
        integrate_tcs = cache -> integrate_elastic_scattering_numerical_GK(cache, Ec_mev)
    elseif method == "numerical-ts"
        cache_source = nothing
        integrate_tcs = cache -> integrate_elastic_scattering_numerical_TS(cache, Ec_mev)
    elseif method == "analytical"
        cache_source = "coefficients"
        integrate_tcs = cache -> integrate_elastic_scattering_analytical(cache, Ec_mev)
    else
        error("Unknown integrate_total_cs method $(method).")
    end

    σt = 0.0
    for (Ai, atai) in zip(A, atpercentA)
        cache = this.get_dcs_cache(ptype, db, Z, Ai, E_in_mev; source=cache_source)
        σt += atai * integrate_tcs(cache)
    end

    return σt * BARN_TO_CM2
end


"""
    sp(this::Elastic_Scattering,Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,
    Ei::Float64,Ec::Float64,particle::Particle,
    A::Vector{Vector{Int64}},atpercentA::Vector{Vector{Float64}})

Gives the stopping power for elastic scattering interaction.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `atz::Vector{Float64}` : atomic fraction of the elements composing the material.
- `N_density::Float64` : number density of the material.
- `Ei::Float64` : incoming particle energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.
- `A::Vector{Vector{Int64}}` : isotope mass numbers per element.
- `atpercentA::Vector{Vector{Float64}}` : isotope atomic fractions per element.

# Output Argument(s)
- `S::Float64` : stopping power.

"""
function sp(this::Elastic_Scattering,Z::Vector{Int64},atz::Vector{Float64},N_density::Float64,Ei::Float64,Ec::Float64,particle::Particle,A::Vector{Vector{Int64}},atpercentA::Vector{Vector{Float64}})
    ptype = get_type(particle)
    endf_db = this.get_endf_db()
    db = endf_db[ptype]

    E_in_mev = Ei * M_E_C2_MEV
    Ec_mev = max(Ec, 0.0) * M_E_C2_MEV
    m_mev = get_mass(particle)

    S_mev_per_cm = 0.0
    for i in eachindex(Z)
        Zi = Z[i]
        S_elem = 0.0
        for (Ai, atai) in zip(A[i], atpercentA[i])
            iso = db.isotopes[(Zi, Ai)]
            M_mev = get_mass(Zi, Ai)
            cache = this.get_dcs_cache(ptype, db, Zi, Ai, E_in_mev)
            mu_b = mu_cm_from_Eout_relativistic(E_in_mev, Ec_mev, m_mev, M_mev)
            mu_lo = iso.identical ? 0.0 : -1.0
            mu_start = max(mu_b, mu_lo)
            mu_hi = 1.0
            if mu_start >= mu_hi
                continue
            end

            # For elastic kinematics, E_out(mu_cm) is affine in mu_cm, so the
            # energy-loss weight can be precomputed once per isotope.
            E1 = E_in_mev + m_mev
            p1 = sqrt(max(E1^2 - m_mev^2, 0.0))
            s = m_mev^2 + M_mev^2 + 2.0 * M_mev * E1
            sqrt_s = sqrt(s)
            E1_star = (s + m_mev^2 - M_mev^2) / (2.0 * sqrt_s)
            p_star = sqrt(max(E1_star^2 - m_mev^2, 0.0))
            beta_cm = p1 / (E1 + M_mev)
            gamma_cm = 1.0 / sqrt(1.0 - beta_cm^2)
            loss_const = E1 - gamma_cm * E1_star
            loss_slope = gamma_cm * beta_cm * p_star

            integrand = mu_cm -> begin
                sigma = elastic_dsig(cache, mu_cm)
                return sigma * max(loss_const - loss_slope * mu_cm, 0.0)
            end
            val = integrate_fp_interval(this, integrand, mu_start, mu_hi)
            S_elem += atai * 2.0 * pi * val * BARN_TO_CM2
        end

        S_mev_per_cm += atz[i] * N_density * S_elem 
    end

    # Keep internal stopping-power units consistent with other interactions (m_e c^2 / cm).
    return S_mev_per_cm  / M_E_C2_MEV
end


"""
    mt(this::Elastic_Scattering,Z::Int64,Ei::Float64,Ec::Float64,particle::Particle,
    A::Union{Nothing,Vector{Int64}}=nothing,
    atpercentA::Union{Nothing,Vector{Float64}}=nothing)

Gives the momentum transfer for elastic scattering interaction.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure. 

# Output Argument(s)
- `T::Float64` : momentum transfer.

"""
function mt(this::Elastic_Scattering,Z::Int64,Ei::Float64,Ec::Float64,particle::Particle,A::Vector{Int64},atpercentA::Vector{Float64})

    ptype = get_type(particle)
    endf_db = this.get_endf_db()
    db = endf_db[ptype]

    E_in_mev = Ei * M_E_C2_MEV
    E_out_mev = max(Ec, 0.0) * M_E_C2_MEV
    m_mev = get_mass(particle)

    T = 0.0
    for (Ai, atai) in zip(A, atpercentA)
        iso = db.isotopes[(Z, Ai)]
        M_mev = get_mass(Z, Ai)
        cache = this.get_dcs_cache(ptype, db, Z, Ai, E_in_mev)
        mu_b = mu_cm_from_Eout_relativistic(E_in_mev, E_out_mev, m_mev, M_mev)
        mu_lo = iso.identical ? 0.0 : -1.0
        mu_start = max(mu_b, mu_lo)
        mu_hi = 1.0
        if mu_start >= mu_hi
            continue
        end

        integrand = mu_cm -> begin
            mu_lab = mu_lab_from_mu_cm(mu_cm, E_in_mev, m_mev, M_mev)
            sigma = elastic_dsig(cache, mu_cm)
            return sigma * (1.0 - mu_lab)
        end
        T += atai * integrate_fp_interval(this, integrand, mu_start, mu_hi)
    end

    return 2.0 * pi * T * BARN_TO_CM2
end

