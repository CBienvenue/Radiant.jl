"""
    Pair_Production

Structure used to define parameters for production of multigroup pair production cross-sections.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}} = Dict((Photon,Photon) => ["A"],(Photon,Electron) => ["P"],(Photon,Positron) => ["P"])` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `(Photon,Photon) => ["A"]` : absorption of incoming photon.
    - `(Photon,Electron) => ["P"]` : produced electron.
    - `(Photon,Positron) => ["P"]` : produced positron.
- `angular_scattering_type::String=modified_dipole` : type of angular scattering, which can takes the following values:
    - `angular_scattering_type = modified_dipole` : modified dipôle distribution, based on Poskus (2019) shape functions.
    - `angular_scattering_type = sommerfield` : Sommerfield distribution.

"""
mutable struct Pair_Production <: Interaction

    # Variable(s)
    name::String
    incoming_particle::Vector{Type}
    interaction_particles::Vector{Type}
    interaction_types::Dict{Tuple{Type,Type},Vector{String}}
    is_CSD::Bool
    is_AFP::Bool
    is_AFP_decomposition::Bool
    is_elastic::Bool
    is_preload_data::Bool
    is_subshells_dependant::Bool
    angular_scattering_type::String
    Cℓk::Array{Float64}
    normalization_factor::Function
    angular_distribution::Function
    is_triplet_contribution::Bool
    scattering_model::String

    # Constructor(s)
    function Pair_Production()
        this = new()
        this.name = "pair_production"
        this.interaction_types = Dict((Photon,Photon) => ["A"],(Photon,Electron) => ["P"],(Photon,Positron) => ["P"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.is_CSD = false
        this.is_AFP = false
        this.is_AFP_decomposition = false
        this.is_elastic = false
        this.is_preload_data = true
        this.is_subshells_dependant = false
        this.set_angular_scattering_type("modified_dipole")
        this.scattering_model = "BTE"
        return this
    end
end

# Method(s)
"""
    set_interaction_types(this::Pair_Production,interaction_types::Dict{Tuple{DataType,DataType},Vector{String}})

To define the interaction types for pair production processes.

# Input Argument(s)
- `this::Pair_Production` : pair production structure.
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}}` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `(Photon,Photon) => ["A"]` : absorption of incoming photon.
    - `(Photon,Electron) => ["P"]` : produced electron.
    - `(Photon,Positron) => ["P"]` : produced positron.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> pair_production = Pair_Production()
julia> pair_production.set_interaction_types( Dict((Electron,Electron) => ["S"]) ) # Only electron scattering, with photon absorption.
```
"""
function set_interaction_types(this::Pair_Production,interaction_types)
    this.interaction_types = interaction_types
end

"""
    set_angular_scattering_type(this::Pair_Production,angular_scattering_type::String)

To define the pair_production photons angular distribution.

# Input Argument(s)
- `this::Pair_Production` : pair_production structure.
- `angular_scattering_type::String` : angular scattering type.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> pair_production = Pair_Production()
julia> pair_production.set_angular_scattering_type("sommerfield")
```
"""
function set_angular_scattering_type(this::Pair_Production,angular_scattering_type::String)
    if lowercase(angular_scattering_type) ∉ ["modified_dipole","sommerfield"] error("Undefined angular distribution.") end
    this.angular_scattering_type = lowercase(angular_scattering_type)
end

"""
    in_distribution(this::Pair_Production)

Describe the energy discretization method for the incoming particle in the pair production
interaction.

# Input Argument(s)
- `this::Pair_Production` : pair production structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function in_distribution(this::Pair_Production)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    out_distribution(this::Pair_Production)

Describe the energy discretization method for the outgoing particle in the pair production
interaction.

# Input Argument(s)
- `this::Pair_Production` : pair production structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function out_distribution(this::Pair_Production)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    bounds(this::Pair_Production,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String)

Gives the integration energy bounds for the outgoing particle for pair production. 

# Input Argument(s)
- `this::Pair_Production` : pair production structure.
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `Ei::Float64` : energy of the incoming particle.
- `type::String` : type of interaction.

# Output Argument(s)
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `isSkip::Bool` : define if the integration is skipped or not.

"""
function bounds(this::Pair_Production,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String)
    # Electron/positron production
    if type == "P" || type == "A"
        Ef⁻ = min(Ef⁻,Ei-2)
        if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    else
        error("Unknown type of method for pair production.")
    end
    return Ef⁻,Ef⁺,isSkip
end

"""
    dcs(this::Pair_Production,L::Int64,Ei::Float64,Ef::Float64,Z::Int64,particle::Particle,
    type::String,iz::Int64)

Gives the Legendre moments of the scattering cross-sections for pair production. 

# Input Argument(s)
- `this::Pair_Production` : pair production structure.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : incoming particle energy.
- `Ef::Float64` : outgoing particle energy.
- `Z::Int64` : atomic number.
- `type::String` : type of interaction.
- `iz::Int64` : index of the element in the material.
- `particles::Vector{Particle}` : list of particles.

# Output Argument(s)
- `σℓ::Vector{Float64}` : Legendre moments of the scattering cross-sections.

"""
function dcs(this::Pair_Production,L::Int64,Ei::Float64,Ef::Float64,Z::Int64,type::String,iz::Int64,particles::Vector{Particle})

    # Initialization
    mₑc² = 0.510999
    rₑ = 2.81794092e-13 # (in cm)
    β² = Ei*(Ei+2)/(Ei+1)^2
    α = 1/137
    a = α*Z
    σs = 0.0
    σℓ = zeros(L+1)
    rs, n∞ = baro_coefficient(Z)

    # Electron/Positron production
    if type == "P"
        if 0 < Ef < Ei-2

            # High-energy Coulomb correction
            if Ei ≥ 50/mₑc²
                fc = a^2*(1/(1+a^2) + 0.202059 - 0.03693*a^2 + 0.00835*a^4 - 0.00201*a^6 + 0.00049*a^8 - 0.00012*a^10 + 0.00003*a^12)
            else
                fc = 0
            end

            # Normalization factor
            A = this.normalization_factor(iz,Ei)

            ϵ = (Ef+1)/Ei
            ϵ₀ = 1/Ei
            b = rs/2 * ϵ₀/(ϵ*(1-ϵ))
            g1 = 7/3 - 2*log(1+b^2) - 6*b*atan(1/b) - b^2*(4 - 4*b*atan(1/b) - 3*log(1+1/b^2))
            g2 = 11/6 - 2*log(1+b^2) - 3*b*atan(1/b) + b^2/2*(4 - 4*b*atan(1/b) - 3*log(1+1/b^2))
            g0 = 4*log(rs) - 4*fc
            ϕ₁ = g1 + g0
            ϕ₂ = g2 + g0
            σs = A * α * rₑ^2 * Z*(Z+0) * 1/Ei * (2*(1/2-ϵ)^2*ϕ₁+ϕ₂)

            # If no explicit positron transport, generate two electrons
            if ~any(is_positron.(particles)) σs *= 2 end
 
            # Sommerfield angular distribution
            β = sqrt(β²)
            if this.angular_scattering_type == "sommerfield"
                Wℓ = zeros(L+1)
                @inbounds for ℓ in range(0,L)
                    if ℓ == 0
                        Wℓ[ℓ+1] = 1
                    elseif ℓ == 1
                        Wℓ[ℓ+1] = (2*β + (1-β^2)*log((1-β)/(1+β)))/(2*β^2)
                    else
                        Wℓ[ℓ+1] = ((2*ℓ-1)*Wℓ[ℓ] - ℓ*β*Wℓ[ℓ-1])/((ℓ-1)*β)
                    end
                    σℓ[ℓ+1] += σs * Wℓ[ℓ+1]
                end

            # Dipôle angular distribution
            elseif this.angular_scattering_type == "modified_dipole"
                A,B,C = this.angular_distribution(iz,Z,Ei,Ef/Ei)
                αi = [C^2,-2*C,1] .* (A-B)
                𝒢a = zeros(L+3)
                𝒢b = zeros(L+3)
                @inbounds for i in range(0,L+2)
                    𝒢a[i+1] = 𝒢₃(i,-2,1,-C,0,1,1)-𝒢₃(i,-2,1,-C,0,1,-1)
                    𝒢b[i+1] = 𝒢₃(i,-4,1,-C,0,1,1)-𝒢₃(i,-4,1,-C,0,1,-1)
                end
                @inbounds for ℓ in range(0,L)
                    for k in range(0,div(ℓ,2))
                        σℓk = 0.0
                        σℓk += (A+B)*𝒢a[ℓ-2*k+1]
                        for i in range(0,2)
                            σℓk += αi[i+1] * 𝒢b[ℓ-2*k+i+1]
                        end
                        σℓ[ℓ+1] += this.Cℓk[ℓ+1,k+1] * σℓk
                    end
                    σℓ[ℓ+1] *= 3/(4*(2*A+B)) * (1-C^2)/(2^ℓ) * σs
                end
                # Correction to deal with high-order Legendre moments
                for ℓ in range(1,L)
                    if abs(σℓ[1]) < abs(σℓ[ℓ+1])
                        σℓ[ℓ+1:end] .= 0.0
                        break
                    end
                end
            else
                error("Unkown angular distribution.")
            end
        end
    else
        error("Unknown interaction.")
    end
    return σℓ
end

"""
    tcs(this::Pair_Production,Ei::Float64,Z::Int64,iz::Int64,Eout::Vector{Float64},
    type::String)

Gives the total cross-section for pair production. 

# Input Argument(s)
- `this::Pair_Production` : pair production structure. 
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.
- `iz::Int64` : index of the element in the material.
- `Eout::Vector{Float64}` : outgoing energy boundaries.
- `type::String` : type of interaction.

# Output Argument(s)
- `σt::Float64` : total cross-section.

"""
function tcs(this::Pair_Production,Ei::Float64,Z::Int64,iz::Int64,Eout::Vector{Float64})

    # Initialization
    mₑc² = 0.510999
    rₑ = 2.81794092e-13 # (in cm)
    α = 1/137
    a = α*Z
    σt = 0.0
    rs, n∞ = baro_coefficient(Z)

    if Ei-2 > 0

        # High-energy Coulomb correction
        if Ei ≥ 50/mₑc²
            fc = a^2*(1/(1+a^2) + 0.202059 - 0.03693*a^2 + 0.00835*a^4 - 0.00201*a^6 + 0.00049*a^8 - 0.00012*a^10 + 0.00003*a^12)
        else
            fc = 0
        end

        # Normalization factor
        A = this.normalization_factor(iz,Ei)

        # Compute total cross-sections
        Ngf = length(Eout)-1
        is_dirac, Np, q_type = out_distribution(this)
        if is_dirac Np = 1; u = [0]; w = [2] else u,w = quadrature(Np,q_type) end
        for gf in range(1,Ngf+1)
            Ef⁻ = Eout[gf]
            if (gf != Ngf+1) Ef⁺ = Eout[gf+1] else Ef⁺ = 0.0 end
            Ef⁻,Ef⁺,isSkip = bounds(this,Ef⁻,Ef⁺,Ei,"A")
            if isSkip continue end
            ΔEf = Ef⁻ - Ef⁺
            for n in range(1,Np)
                Ef = (u[n]*ΔEf + (Ef⁻+Ef⁺))/2

                ϵ = (Ef+1)/Ei
                ϵ₀ = 1/Ei
                b = rs/2 * ϵ₀/(ϵ*(1-ϵ))
                g1 = 7/3 - 2*log(1+b^2) - 6*b*atan(1/b) - b^2*(4 - 4*b*atan(1/b) - 3*log(1+1/b^2))
                g2 = 11/6 - 2*log(1+b^2) - 3*b*atan(1/b) + b^2/2*(4 - 4*b*atan(1/b) - 3*log(1+1/b^2))
                g0 = 4*log(rs) - 4*fc
                ϕ₁ = g1 + g0
                ϕ₂ = g2 + g0
                σt += ΔEf/2 * w[n] * A * α * rₑ^2 * Z*(Z+0) * 1/Ei * (2*(1/2-ϵ)^2*ϕ₁+ϕ₂)
            end
        end

    end
    return σt
end

"""
    preload_data(this::Pair_Production,Z::Vector{Int64},Emax::Float64,Emin::Float64,E_out::Vector{Float64},L::Int64)

Preload data for multigroup pair production calculations.

# Input Argument(s)
- `this::Pair_Production` : pair production structure. 
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `Emax::Float64` : maximum energy of the incoming particle.
- `Emin::Float64` : minimum energy of the incoming particle.
- `Eout::Vector{Float64}` : outgoing energy boundaries.
- `L::Int64` : Legendre truncation order.

# Output Argument(s)
N/A

"""
function preload_data(this::Pair_Production,Z::Vector{Int64},Emax::Float64,Emin::Float64,E_out::Vector{Float64},L::Int64)
    this.preload_normalization_factor(Z,Emax,Emin,E_out)
    if this.angular_scattering_type == "modified_dipole"

        # Precompute angular integration factors
        this.Cℓk = zeros(L+1,div(L,2)+1)
        for ℓ in range(0,L), k in range(0,div(L,2))
            this.Cℓk[ℓ+1,k+1] = (-1)^k * exp( sum(log.(1:2*ℓ-2*k)) - sum(log.(1:k)) - sum(log.(1:ℓ-k)) - sum(log.(1:ℓ-2*k)) )
        end

        # Preload angular distribution from Poskus (2019)
        this.preload_angular_distribution(Z)

    end
end

"""
    preload_data(this::Pair_Production,Z::Vector{Int64},Emax::Float64,Emin::Float64,E_out::Vector{Float64})

Preload renormalization factor for pair production cross-sections

# Input Argument(s)
- `this::Pair_Production` : pair production structure. 
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `Emax::Float64` : maximum energy of the incoming particle.
- `Emin::Float64` : minimum energy of the incoming particle.
- `Eout::Vector{Float64}` : outgoing energy boundaries.
- `L::Int64` : Legendre truncation order.

# Output Argument(s)
N/A

"""
function preload_normalization_factor(this::Pair_Production,Z::Vector{Int64},Emax::Float64,Emin::Float64,E_out::Vector{Float64})
    path = joinpath(find_package_root(), "data", "pair_production_JENDL5.jld2")
    data = load(path)
    Nz = length(Z)
    E = Vector{Vector{Float64}}(undef,Nz)
    A = Vector{Vector{Float64}}(undef,Nz)
    spline_A = Vector{Function}(undef,Nz)
    this.normalization_factor = function normalization_factor_equal_to_one(iz::Int64,Ei::Float64) return 1 end
    for iz in range(1,Nz)
        E[iz] = data["E"][Z[iz]]
        σt_exp = data["σ"][Z[iz]]

        if ~all(E[iz][i] < E[iz][i+1] for i in 1:length(E[iz])-1) error("Incorrect input data.") end

        # Define the interval corresponding to the required interpolation data for calculations
        i⁻ = max(searchsortedfirst(E[iz],Emin) - 1,1)
        i⁺ = searchsortedfirst(E[iz],Emax)
        E[iz] = E[iz][i⁻:i⁺]
        σt_exp = σt_exp[i⁻:i⁺]

        Ng = length(E[iz])
        σt = zeros(Ng)
        for i in range(1,Ng)
            σt[i] = tcs(this,E[iz][i],Z[iz],iz,E_out)
        end
        A[iz] = zeros(Ng)
        for i in range(1,Ng)
            if (σt[i] != 0) A[iz][i] = σt_exp[i]/σt[i] else A[iz][i] = 0.0 end
        end
        if length(E[iz]) != 1
            spline_A[iz] = cubic_hermite_spline(E[iz],A[iz])
        end
    end
    this.normalization_factor = function normalization_factor(iz::Int64,Ei::Float64)
        return spline_A[iz](Ei)
    end
end

"""
    preload_angular_distribution(this::Pair_Production,Z::Vector{Int64})

Preload angular distribution of produced leptons.

# Input Argument(s)
- `this::Pair_Production` : pair production structure. 
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.

# Output Argument(s)
N/A

"""
function preload_angular_distribution(this::Pair_Production,Z::Vector{Int64})

    # Extract vectors
    path = joinpath(find_package_root(), "data", "bremsstrahlung_photons_distribution_poskus_2019.jld2")
    data = load(path)

    # Extract scaled cross-sections
    Nz = length(Z)
    A = Vector{Array{Float64}}(undef,Nz)
    B = Vector{Array{Float64}}(undef,Nz)
    C = Vector{Array{Float64}}(undef,Nz)
    E = Vector{Vector{Float64}}(undef,Nz)
    r = Vector{Array{Float64}}(undef,Nz)
    for iz in range(1,Nz)
        A[iz] = data["A"][Z[iz]]
        B[iz] = data["B"][Z[iz]]
        C[iz] = data["C"][Z[iz]]
        E[iz] = data["E"][Z[iz]]
        r[iz] = data["r"][Z[iz]] 
    end

    # Return the interpolation function
    this.angular_distribution = function angular_distribution(iz::Int64,Z::Int64,Ei::Float64,ri::Float64)
        mₑc² = 0.510999
        β = sqrt(Ei*(Ei+2))/(Ei+1)
        if Ei ≥ 3/mₑc²
            Ai = 1.0
            Bi = 0.0
            Ci = β
        else
            # Find index vector E
            i = searchsortedfirst(E[iz],Ei)

            # Find index vector r
            i⁻ = searchsortedfirst(r[iz][i-1,:],ri)
            i⁺ = searchsortedfirst(r[iz][i,:],ri)

            # Interpolation of parameter
            if i⁻ < 13
                A⁻ = cubic_hermite_spline(r[iz][i-1,i⁻-1:i⁻],A[iz][i-1,i⁻-1:i⁻])(ri)
                B⁻ = cubic_hermite_spline(r[iz][i-1,i⁻-1:i⁻],B[iz][i-1,i⁻-1:i⁻])(ri)
                C⁻ = cubic_hermite_spline(r[iz][i-1,i⁻-1:i⁻],C[iz][i-1,i⁻-1:i⁻])(ri)
            else
                A⁻ = A[iz][i-1,end]
                B⁻ = B[iz][i-1,end]
                C⁻ = C[iz][i-1,end]
            end
            if i⁺ < 13
                A⁺ = cubic_hermite_spline(r[iz][i,i⁺-1:i⁺],A[iz][i,i⁺-1:i⁺])(ri)
                B⁺ = cubic_hermite_spline(r[iz][i,i⁺-1:i⁺],B[iz][i,i⁺-1:i⁺])(ri)
                C⁺ = cubic_hermite_spline(r[iz][i,i⁺-1:i⁺],C[iz][i,i⁺-1:i⁺])(ri)
            else
                A⁺ = A[iz][i,end]
                B⁺ = B[iz][i,end]
                C⁺ = C[iz][i,end]
            end
            Ai = cubic_hermite_spline(E[iz][i-1:i],[A⁻,A⁺])(Ei)
            Bi = cubic_hermite_spline(E[iz][i-1:i],[B⁻,B⁺])(Ei)
            Ci = cubic_hermite_spline(E[iz][i-1:i],[C⁻,C⁺])(Ei)
        end
        return Ai,Bi,Ci
    end
end