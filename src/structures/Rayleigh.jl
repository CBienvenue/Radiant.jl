"""
    Rayleigh

Structure used to define parameters for production of multigroup Rayleigh cross-sections.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{String,String},Vector{String}} = Dict(("photons","photons") => ["S"])`: Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `("photons","photons") => ["S"]`: elastic scattering of photons

"""
mutable struct Rayleigh <: Interaction

    # Variable(s)
    name::String
    incoming_particle::Vector{String}
    interaction_particles::Vector{String}
    interaction_types::Dict{Tuple{String,String},Vector{String}}
    is_CSD::Bool
    is_AFP::Bool
    is_elastic::Bool
    is_preload_data::Bool
    is_subshells_dependant::Bool
    is_ETC::Bool

    rayleigh_cross_sections::Function

    # Constructor(s)
    function Rayleigh()
        this = new()
        this.name = "rayleigh"
        this.interaction_types = Dict(("photons","photons") => ["S"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.is_CSD = false
        this.is_AFP = false
        this.is_elastic = true
        this.is_preload_data = true
        this.is_subshells_dependant = false
        this.is_ETC = true
        return this
    end

end

# Method(s)
function in_distribution(this::Rayleigh)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

function out_distribution(this::Rayleigh)
    is_dirac = true
    N = 1
    quadrature = "dirac"
    return is_dirac, N, quadrature
end

function bounds(this::Rayleigh,Ef⁻::Float64,Ef⁺::Float64,gi::Int64,gf::Int64)
    if (gf != gi) isSkip = true else isSkip = false end
    return Ef⁻,Ef⁺,isSkip
end

function dcs(this::Rayleigh,L::Int64,Ei::Float64,Z::Int64,iz::Int64)

    # Initialization
    rₑ = 2.81794092E-13 # (in cm)
    σℓ = zeros(L+1)

    # Angular quadrature
    μ,w = quadrature(32,"gauss-legendre")
    for n in range(1,32)
        σs = this.rayleigh_cross_sections(iz,Z,Ei,μ[n])
        Pℓμ = legendre_polynomials(L,μ[n])
        for ℓ in range(0,L) σℓ[ℓ+1] += w[n] * σs * Pℓμ[ℓ+1] end
    end
    return σℓ
end

function tcs(this::Rayleigh,Ei::Float64,Z::Int64,iz::Int64)
    σt = dcs(this,0,Ei,Z,iz)[1]
    return σt
end

function scattering_centers(this::Rayleigh,Z::Union{Int64,Vector{Int64}},ρ::Float64)
    return nuclei_density.(Z,ρ)
end

function preload_data(this::Rayleigh,Z::Vector{Int64})
    this.preload_rayleigh_cross_sections(Z)
end

function preload_rayleigh_cross_sections(this::Rayleigh,Z::Vector{Int64})

    path = joinpath(find_package_root(), "data", "rayleigh_factors_JENDL5.jld2")
    data = load(path)
    Nz = length(Z)
    x = Vector{Vector{Float64}}(undef,Nz)
    F = Vector{Vector{Float64}}(undef,Nz)
    E_real = Vector{Vector{Float64}}(undef,Nz)
    f_real = Vector{Vector{Float64}}(undef,Nz)
    E_imag = Vector{Vector{Float64}}(undef,Nz)
    f_imag = Vector{Vector{Float64}}(undef,Nz)
    for iz in range(1,Nz)
        x[iz] = data["x"][Z[iz]]
        F[iz] = data["F"][Z[iz]]
        E_real[iz] = data["E_real"][Z[iz]]
        f_real[iz] = data["f_real"][Z[iz]]
        E_imag[iz] = data["E_imag"][Z[iz]]
        f_imag[iz] = data["f_imag"][Z[iz]]
    end

    # Return the interpolation function
    this.rayleigh_cross_sections = function rayleigh_cross_sections(iz::Int64,Z::Int64,Ei::Float64,μ::Float64)
        rₑ = 2.81794092E-13 # (in cm)
        hc = 1/20.60744 # (hc in mₑc² × Å)
        xi = 2*Ei/(hc)*sqrt((1-μ)/2)
        Fi = linear_interpolation(xi,x[iz],F[iz])
        if Ei < E_real[iz][end] fi_real = linear_interpolation(Ei,E_real[iz],f_real[iz]) else fi_real = 0 end
        if Ei < E_imag[iz][end] fi_imag = linear_interpolation(Ei,E_imag[iz],f_imag[iz]) else fi_imag = 0 end
        return π*rₑ^2 * (1+μ^2) * ((Fi + fi_real)^2 + fi_imag^2)
    end
end