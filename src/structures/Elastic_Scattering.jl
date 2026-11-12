"""
    IsotopeElasticDB

Stores ENDF elastic-scattering data for one target isotope.

# Field(s)
- `energies::Vector{Float64}` : incident energies of the ENDF tables [MeV].
- `mu_grid::Vector{Float64}` : cosine grid common to the LTP>2 tables (empty otherwise).
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
- `sigmaNI_E` : energy grid of the nuclear-plus-interference cross-section σ_NI (MF=3).
- `sigmaNI_S` : nuclear-plus-interference cross-section σ_NI [b].
- `sigmaNI_nbt` : σ_NI interpolation breakpoints.
- `sigmaNI_int` : σ_NI interpolation laws.
"""
struct IsotopeElasticDB
    energies::Vector{Float64}
    mu_grid::Vector{Float64}
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
- `coulomb_kinematics::String` : Coulomb kinematics convention.
- `screening::Bool` : whether Coulomb screening is enabled.
- `isotopes::Dict{Tuple{Int,Int},IsotopeElasticDB}` : isotope data indexed by `(Z, A)`.
"""
struct ElasticScatteringENDFDB
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
- `k_coulomb::Float64` : center-of-mass wave number [b^-1/2].
- `eta_coulomb::Float64` : Coulomb Sommerfeld parameter.
- `delta_coulomb::Float64` : Molière screening shift δ in 1 - μ + δ.
- `mu_grid::Vector{Float64}` : cosine grid of a tabulated cross-section.
- `sigma_grid::Vector{Float64}` : differential cross-section on `mu_grid` [b/sr]
  (Coulomb + σ_NI·P_NI for LTP>2), interpolated linearly by `elastic_dsig_raw`.
- `NL::Int` : highest Legendre order of the LTP=1/2 coefficients.
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
- `coulomb_kinematics::String = "standard"` : Coulomb kinematics convention.
- `screening::Bool = true` : whether Coulomb screening is enabled.
- `coulomb_only::Bool = false` : whether only the Coulomb (Rutherford) cross section is used, ignoring the ENDF library.
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
    coulomb_kinematics::String
    screening::Bool
    coulomb_only::Bool
    dcs_cache::Dict{Tuple{Type,Int,Int,String,Bool},ElasticDCSCache}

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
        this.dcs_cache = Dict{Tuple{Type,Int,Int,String,Bool},ElasticDCSCache}()
        this.set_scattering_model("BFP")
        this.set_library("TENDL2023")
        this.set_data_path("../../data")
        this.set_endf_db(Dict{Type,ElasticScatteringENDFDB}())
        this.set_coulomb_kinematics("standard")
        this.set_screening(true)
        this.set_coulomb_only(false)
        return this
    end
end

const QUADGK_PKGID = Base.PkgId(Base.UUID("1fd47b50-473d-5c70-9696-f719f8f3bcdc"), "QuadGK")

"""
    require_quadgk()

Loads the QuadGK package object used by `integrate_elastic_scattering_numerical_GK`,
kept for debugging.

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
            if haskey(endf_db, ptype)
                continue
            elseif this.get_coulomb_only()
                endf_db[ptype] = init_elastic_scattering_coulomb_db(isotopes, particle;
                    coulomb_kinematics=this.get_coulomb_kinematics(),
                    screening=this.get_screening())
            else
                endf_db[ptype] = init_elastic_scattering_endf_db(
                    this.get_library(), isotopes, particle;
                    data_root=this.get_data_path(),
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
julia> elastic_scattering.set_interaction_types( Dict((Proton,Proton) => ["S"]) ) # No recoil protons.
```
"""
function set_interaction_types(this::Elastic_Scattering,interaction_types)
    this.interaction_types = interaction_types
end

"""
    set_scattering_model(this::Elastic_Scattering,scattering_model::String)

To define the solver for elastic scattering.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.
- `scattering_model::String` : solver for elastic scattering, which can be:
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
    this.endf_db = Dict{Type,ElasticScatteringENDFDB}()
    empty!(this.dcs_cache)
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
    set_coulomb_only(this::Elastic_Scattering,coulomb_only::Bool)

Set whether elastic scattering uses only the Coulomb (Rutherford) cross section. When
enabled, the ENDF library is not read and the nuclear and Coulomb-nuclear interference
terms are dropped from every calculation (differential, total, stopping power and
momentum transfer). The Coulomb kinematics and screening settings still apply.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.
- `coulomb_only::Bool` : Coulomb-only flag (default `false`).

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_scattering = Elastic_Scattering()
julia> elastic_scattering.set_coulomb_only(true)
```
"""
function set_coulomb_only(this::Elastic_Scattering,coulomb_only::Bool)
    this.coulomb_only = coulomb_only
    this.endf_db = Dict{Type,ElasticScatteringENDFDB}()
    empty!(this.dcs_cache)
end

"""
    get_coulomb_only(this::Elastic_Scattering)

Gets whether elastic scattering uses only the Coulomb (Rutherford) cross section.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.

# Output Argument(s)
- `coulomb_only::Bool` : Coulomb-only flag.
"""
function get_coulomb_only(this::Elastic_Scattering)
    return this.coulomb_only
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
    Z::Int, A::Int, E_in_mev::Float64)

Gets or builds the cached DCS interpolation data for one particle, isotope, and incident
energy.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.
- `ptype::Type` : incoming particle type.
- `db::ElasticScatteringENDFDB` : ENDF database for the incoming particle type.
- `Z::Int` : target atomic number.
- `A::Int` : target isotope mass number.
- `E_in_mev::Float64` : incoming particle kinetic energy [MeV].

# Output Argument(s)
- `cache::ElasticDCSCache` : cached DCS interpolation data.
"""
function get_dcs_cache(this::Elastic_Scattering,ptype::Type,db::ElasticScatteringENDFDB,
                       Z::Int,A::Int,E_in_mev::Float64)
    key = (ptype, Z, A, db.coulomb_kinematics, db.screening)
    if !haskey(this.dcs_cache, key) || this.dcs_cache[key].E_in_mev != E_in_mev
        this.dcs_cache[key] = prepare_elastic_dsig_cache(db, Z, A, E_in_mev)
    end
    return this.dcs_cache[key]
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
    dcs_all_groups!(this::Elastic_Scattering,𝓕::Array{Float64,2},𝓕ₑ::Vector{Float64},
    Z::Int64,A::Int64,L::Int64,Ei::Float64,Eout::Vector{Float64},Ec::Float64,
    particle::Particle,scale::Float64;recoil::Bool=false)

Adds to the feed functions the Legendre moments and the energy moment of every outgoing
group at once, for the scattered projectile, or with `recoil` for the recoiling target of
identical particles (the proton of p + H-1). The catastrophic range, from backscattering
to the cutoff `Ec`, is walked once by `elastic_group_moments!`, cut at the group
boundaries and wherever the cross-section is not smooth; only its positive part is kept,
as in `elastic_dsig`. The moments of order zero of the scattered projectile over all
groups, plus the part below the last group, add up to `tcs` to round-off. Nothing is
added below the last group: that part is left to the absorption cross-section.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.
- `𝓕::Array{Float64,2}` : feed function to add to [cm²], size (Ng+1, L+1).
- `𝓕ₑ::Vector{Float64}` : energy-weighted feed function to add to [cm² × mₑc²].
- `Z::Int64` : atomic number.
- `A::Int64` : isotope mass number.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : incoming particle energy [mₑc²].
- `Eout::Vector{Float64}` : outgoing group boundaries [mₑc²], decreasing.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions [mₑc²].
- `particle::Particle` : incoming particle.
- `scale::Float64` : factor applied to every contribution (atomic fraction of the
  element times that of the isotope).
- `recoil::Bool` : boolean to give the recoiling target instead of the projectile.

# Output Argument(s)
N/A

"""
function dcs_all_groups!(this::Elastic_Scattering,𝓕::Array{Float64,2},𝓕ₑ::Vector{Float64},Z::Int64,A::Int64,L::Int64,Ei::Float64,Eout::Vector{Float64},Ec::Float64,particle::Particle,scale::Float64;recoil::Bool=false)
    ptype = get_type(particle)
    db = this.get_endf_db()[ptype]
    E_in_mev = Ei * M_E_C2_MEV
    cache = this.get_dcs_cache(ptype, db, Z, A, E_in_mev)
    recoil && !cache.identical && error("Recoil feed functions are only given for identical particles.")

    # Same range as the total cross-section
    mu_lo = cache.identical ? 0.0 : -1.0
    mu_c = mu_cm_from_Eout_relativistic(E_in_mev, max(Ec, 0.0) * M_E_C2_MEV, cache.m_mev, cache.M_mev)
    mu_c <= mu_lo && return nothing

    Ng = size(𝓕, 1) - 1
    σl = zeros(Ng, L+1)
    σe = zeros(Ng)
    elastic_group_moments!(σl, σe, cache, Eout[1:Ng+1] .* M_E_C2_MEV, mu_lo, mu_c, L; recoil=recoil)
    𝓕[1:Ng,:] .+= σl .* (BARN_TO_CM2 * scale)
    𝓕ₑ[1:Ng] .+= σe .* (BARN_TO_CM2 * scale / M_E_C2_MEV)
    return nothing
end

"""
    tcs(this::Elastic_Scattering,Ei::Float64,Ec::Float64,particle::Particle,Z::Int64;
    A::Union{Nothing,Vector{Int64}}=nothing,
    atpercentA::Union{Nothing,Vector{Float64}}=nothing,
    Eout::Vector{Float64})

Gives the part of the total cross-section for elastic scattering interaction that is not
in the feed functions: the catastrophic collisions whose scattered projectile ends below
the last outgoing group (`integrate_elastic_scattering_below`). The multigroup
calculation adds to it the moments of order zero of the feed functions, summed over the
groups, to obtain the total cross-section; the absorption cross-section is then exactly
this part.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.
- `Ei::Float64` : incoming particle energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.
- `Z::Int64` : atomic number.
- `A::Union{Nothing,Vector{Int64}}` : isotope mass numbers for this element.
- `atpercentA::Union{Nothing,Vector{Float64}}` : isotope atomic fractions.
- `Eout::Vector{Float64}` : outgoing group boundaries, decreasing.

# Output Argument(s)
- `σt::Float64` : total cross-section below the last outgoing group.

"""
function tcs(this::Elastic_Scattering,Ei::Float64,Ec::Float64,particle::Particle,Z::Int64;
    A::Union{Nothing,Vector{Int64}}=nothing,
    atpercentA::Union{Nothing,Vector{Float64}}=nothing,
    Eout::Vector{Float64})

    if isnothing(A) || isnothing(atpercentA)
        pairs = isotopic_composition(Z)
        A = [p[1] for p in pairs]
        atpercentA = [p[2] for p in pairs]
    end

    ptype = get_type(particle)
    db = this.get_endf_db()[ptype]
    E_in_mev = Ei * M_E_C2_MEV
    Ec_mev = Ec * M_E_C2_MEV
    E_last_mev = Eout[end] * M_E_C2_MEV

    σt = 0.0
    for (Ai, atai) in zip(A, atpercentA)
        cache = this.get_dcs_cache(ptype, db, Z, Ai, E_in_mev)
        σt += atai * integrate_elastic_scattering_below(cache, Ec_mev, E_last_mev)
    end
    return σt * BARN_TO_CM2
end


"""
    sp(this::Elastic_Scattering,Z::Vector{Int64},atz::Vector{Float64},N_density::Float64,
    Ei::Float64,Ec::Float64,particle::Particle,
    A::Vector{Vector{Int64}},atpercentA::Vector{Vector{Float64}})

Gives the stopping power of the soft elastic collisions, energy loss below the cutoff
(center-of-mass cosines from μ_c to 1), 2π ∫ dσ/dΩ (E_in - E_out) dμ. The energy lost is
e1 (1 - μ_cm) exactly, and the integral is `elastic_soft_integrals`.

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
    db = this.get_endf_db()[ptype]
    E_in_mev = Ei * M_E_C2_MEV
    Ec_mev = max(Ec, 0.0) * M_E_C2_MEV

    S_mev_per_cm = 0.0
    for i in eachindex(Z)
        S_elem = 0.0
        for (Ai, atai) in zip(A[i], atpercentA[i])
            cache = this.get_dcs_cache(ptype, db, Z[i], Ai, E_in_mev)
            mu_c = max(mu_cm_from_Eout_relativistic(E_in_mev, Ec_mev, cache.m_mev, cache.M_mev),
                       cache.identical ? 0.0 : -1.0)
            # The energy lost is e1 (1 - μ_cm) exactly
            _, e1 = scattered_energy_affine(cache)
            I_cm, _ = elastic_soft_integrals(cache, mu_c)
            S_elem += atai * 2.0 * pi * e1 * I_cm * BARN_TO_CM2
        end
        S_mev_per_cm += atz[i] * N_density * S_elem
    end

    # Keep internal stopping-power units consistent with other interactions (m_e c^2 / cm).
    return S_mev_per_cm / M_E_C2_MEV
end


"""
    mt(this::Elastic_Scattering,Z::Int64,Ei::Float64,Ec::Float64,particle::Particle,
    A::Vector{Int64},atpercentA::Vector{Float64})

Gives the momentum transfer of the soft elastic collisions (center-of-mass cosines from
μ_c to 1), π ∫ dσ/dΩ (1 - μ_lab) dμ, half the transport cross-section, from
`elastic_soft_integrals`.

# Input Argument(s)
- `this::Elastic_Scattering` : elastic scattering structure.
- `Z::Int64` : atomic number.
- `Ei::Float64` : incoming particle energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.
- `A::Vector{Int64}` : isotope mass numbers for this element.
- `atpercentA::Vector{Float64}` : isotope atomic fractions.

# Output Argument(s)
- `T::Float64` : momentum transfer.

"""
function mt(this::Elastic_Scattering,Z::Int64,Ei::Float64,Ec::Float64,particle::Particle,A::Vector{Int64},atpercentA::Vector{Float64})
    ptype = get_type(particle)
    db = this.get_endf_db()[ptype]
    E_in_mev = Ei * M_E_C2_MEV
    E_out_mev = max(Ec, 0.0) * M_E_C2_MEV

    T = 0.0
    for (Ai, atai) in zip(A, atpercentA)
        cache = this.get_dcs_cache(ptype, db, Z, Ai, E_in_mev)
        mu_c = max(mu_cm_from_Eout_relativistic(E_in_mev, E_out_mev, cache.m_mev, cache.M_mev),
                   cache.identical ? 0.0 : -1.0)
        _, I_lab = elastic_soft_integrals(cache, mu_c)
        T += atai * I_lab
    end

    # Fokker-Planck momentum transfer is half the transport cross section:
    # T = 1/2 ∫ (1-μ) σ dΩ = π ∫ (1-μ) σ dμ.
    return pi * T * BARN_TO_CM2
end

