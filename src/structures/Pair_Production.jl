"""
    Pair_Production

Structure used to define parameters for production of multigroup pair production cross-sections.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{String,String},Vector{String}} = Dict(("photons","photons") => ["A"],("photons","electrons") => ["P"],("photons","positrons") => ["P"])`: Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `("photons","photons") => ["A"]`: absorption of incoming photon.
    - `("photons","electrons") => ["P"]`: produced electron.
    - `("photons","positrons") => ["P"]`: produced positron.
- `angular_scattering_type::String=modified_dipole`: type of angular scattering, which can takes the following values:
    - `angular_scattering_type = modified_dipole`: modified dipôle distribution, based on Poskus (2019) shape functions.
    - `angular_scattering_type = sommerfield`: Sommerfield distribution.

"""
mutable struct Pair_Production <: Interaction

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
    angular_scattering_type::String
    Cℓk::Array{Float64}
    normalization_factor::Function
    angular_distribution::Function
    is_triplet_contribution::Bool

    # Constructor(s)
    function Pair_Production()
        this = new()
        this.name = "pair_production"
        this.interaction_types = Dict(("photons","photons") => ["A"],("photons","electrons") => ["P"],("photons","positrons") => ["P"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.is_CSD = false
        this.is_AFP = false
        this.is_elastic = false
        this.is_preload_data = true
        this.is_subshells_dependant = false
        this.angular_scattering_type = "modified_dipole"
        return this
    end
end

# Method(s)
"""
    set_interaction_types(this::Pair_Production,interaction_types::Dict{Tuple{String,String},Vector{String}})

To define the interaction types for pair production processes.

# Input Argument(s)
- `this::Pair_Production`: pair production structure.
- `interaction_types::Dict{Tuple{String,String},Vector{String}}`: Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `("photons","photons") => ["A"]`: absorption of incoming photon.
    - `("photons","electrons") => ["P"]`: produced electron.
    - `("photons","positrons") => ["P"]`: produced positron.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> pair_production = Pair_Production()
julia> pair_production.set_interaction_types( Dict(("electrons","electrons") => ["S"]) ) # Only electron scattering, with photon absorption.
```
"""
function set_interaction_types(this::Pair_Production,interaction_types::Dict{Tuple{String,String},Vector{String}})
    this.interaction_types = interaction_types
end

"""
    set_angular_scattering_type(this::Pair_Production,angular_scattering_type::String)

To define the pair_production photons angular distribution.

# Input Argument(s)
- `this::Pair_Production`: pair_production structure.
- `angular_scattering_type::String`: angular scattering type.

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

function in_distribution(this::Pair_Production)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

function out_distribution(this::Pair_Production)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

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

function dcs(this::Pair_Production,L::Int64,Ei::Float64,Ef::Float64,Z::Int64,particle::String,type::String,iz::Int64,particles::Vector{String},Ecutoff::Float64)

    # Initialization
    β² = Ei*(Ei+2)/(Ei+1)^2
    α = 1/137
    a = α*Z
    σs = 0.0
    σℓ = zeros(L+1)
    rs, n∞ = baro_coefficient(Z)
    η = 0

    #=
    if Ei ≥ 4 && this.is_triplet_contribution
        a = α*Z
        ν = (0.2840-0.1909*a)*log(4/Ei) + (0.1095+0.2206*a)*(log(4/Ei))^2 + (0.02888-0.04269*a)*(log(4/Ei))^3 + (0.002527+0.002623*a)*(log(4/Ei))^4
        η = (1-exp(-ν))*n∞
    else
        η = 0
    end
    =#

    # Electron/Positron production
    if type == "P"
        if 0 < Ef < Ei-2

            # High-energy Coulomb correction
            fc = a^2*(1/(1+a^2) + 0.202059 - 0.03693*a^2 + 0.00835*a^4 - 0.00201*a^6 + 0.00049*a^8 - 0.00012*a^10 + 0.00003*a^12)
            #F₀ = (-0.1774-12.10*a+11.18*a^2)*sqrt(2/Ei) + (8.523+73.26*a-44.41*a^2)*(2/Ei) - (13.52+121.1*a-96.41*a^2)*(2/Ei)^(3/2) + (8.946+62.05*a-63.41*a^2)*(2/Ei)^2

            # Normalization factor
            A = this.normalization_factor(iz,Ei)

            # Screening functions
            ϵ = (Ef+1)/Ei
            ϵ₀ = 1/Ei
            b = rs/2 * ϵ₀/(ϵ*(1-ϵ))
            g1 = 7/3 - 2*log(1+b^2) - 6*b*atan(1/b) - b^2*(4 - 4*b*atan(1/b) - 3*log(1+1/b^2))
            g2 = 11/6 - 2*log(1+b^2) - 3*b*atan(1/b) + b^2/2*(4 - 4*b*atan(1/b) - 3*log(1+1/b^2))
            g0 = 4*log(rs) - 4*fc
            ϕ₁ = max(g1 + g0,0)
            ϕ₂ = max(g2 + g0,0)

            # Cross-sections
            σs = A * Z*(Z+η) * (2*(1/2-ϵ)^2*ϕ₁+ϕ₂)

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

function tcs(this::Pair_Production,Ei::Float64,Z::Int64,iz::Int64,Eout::Vector{Float64},type::String)

    # Initialization
    β² = Ei*(Ei+2)/(Ei+1)^2
    α = 1/137
    a = α*Z
    σt = 0.0
    rs, n∞ = baro_coefficient(Z)
    η = 0

    #=
    if Ei ≥ 4 && this.is_triplet_contribution
        a = α*Z
        ν = (0.2840-0.1909*a)*log(4/Ei) + (0.1095+0.2206*a)*(log(4/Ei))^2 + (0.02888-0.04269*a)*(log(4/Ei))^3 + (0.002527+0.002623*a)*(log(4/Ei))^4
        η = (1-exp(-ν))*n∞
    else
        η = 0
    end
    =#

    if Ei > 2

        # High-energy Coulomb correction
        fc = a^2*(1/(1+a^2) + 0.202059 - 0.03693*a^2 + 0.00835*a^4 - 0.00201*a^6 + 0.00049*a^8 - 0.00012*a^10 + 0.00003*a^12)
        #F₀ = (-0.1774-12.10*a+11.18*a^2)*sqrt(2/Ei) + (8.523+73.26*a-44.41*a^2)*(2/Ei) - (13.52+121.1*a-96.41*a^2)*(2/Ei)^(3/2) + (8.946+62.05*a-63.41*a^2)*(2/Ei)^2

        # Normalization factor
        A = this.normalization_factor(iz,Ei)

        # Compute total cross-sections
        Ngf = length(Eout)-1
        is_dirac, Np, q_type = out_distribution(this)
        if is_dirac Np = 1; u = [0]; w = [2] else u,w = quadrature(Np,q_type) end
        for gf in range(1,Ngf+1)
            Ef⁻ = Eout[gf]
            if (gf != Ngf+1) Ef⁺ = Eout[gf+1] else Ef⁺ = 0.0 end
            Ef⁻,Ef⁺,isSkip = bounds(this,Ef⁻,Ef⁺,Ei,type)
            if isSkip continue end
            ΔEf = Ef⁻ - Ef⁺
            for n in range(1,Np)
                Ef = (u[n]*ΔEf + (Ef⁻+Ef⁺))/2

                # Screening functions
                ϵ = (Ef+1)/Ei
                ϵ₀ = 1/Ei
                b = rs/2 * ϵ₀/(ϵ*(1-ϵ))
                g1 = 7/3 - 2*log(1+b^2) - 6*b*atan(1/b) - b^2*(4 - 4*b*atan(1/b) - 3*log(1+1/b^2))
                g2 = 11/6 - 2*log(1+b^2) - 3*b*atan(1/b) + b^2/2*(4 - 4*b*atan(1/b) - 3*log(1+1/b^2))
                g0 = 4*log(rs) - 4*fc
                ϕ₁ = max(g1 + g0,0)
                ϕ₂ = max(g2 + g0,0)
                
                # Cross-sections
                σt += ΔEf/2 * w[n] * A * Z*(Z+η) * (2*(1/2-ϵ)^2*ϕ₁+ϕ₂)
            end
        end

    end
    return σt
end

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

function preload_normalization_factor(this::Pair_Production,Z::Vector{Int64},Emax::Float64,Emin::Float64,E_out::Vector{Float64})
    data = load("./data/pair_production_JENDL5.jld2")
    Nz = length(Z)
    E = Vector{Vector{Float64}}(undef,Nz)
    A = Vector{Vector{Float64}}(undef,Nz)
    spline_A = Vector{Function}(undef,Nz)
    this.normalization_factor = function normalization_factor_equal_to_one(iz::Int64,Ei::Float64) return 1 end
    this.is_triplet_contribution = false
    for iz in range(1,Nz)
        E[iz] = data["E"][Z[iz]]
        σt_exp = data["σ"][Z[iz]]

        # Define the interval corresponding to the required interpolation data for calculations
        i⁻ = max(searchsortedfirst(E[iz],Emin) - 1,1)
        i⁺ = searchsortedfirst(E[iz],Emax)
        E[iz] = E[iz][i⁻:i⁺]
        σt_exp = σt_exp[i⁻:i⁺]

        Ng = length(E[iz])
        σt = zeros(Ng)
        for i in range(1,Ng)
            σt[i] = tcs(this,E[iz][i],Z[iz],iz,E_out,"A")
        end
        A[iz] = zeros(Ng)
        for i in range(1,Ng)
            if (σt[i] != 0) A[iz][i] = σt_exp[i]/σt[i] else A[iz][i] = 0.0 end
        end
        spline_A[iz] = cubic_hermite_spline(E[iz],A[iz])
    end
    this.normalization_factor = function normalization_factor(iz::Int64,Ei::Float64)
        return spline_A[iz](Ei)
    end
    this.is_triplet_contribution = true
end

function preload_angular_distribution(this::Pair_Production,Z::Vector{Int64})

    # Extract vectors
    file = "./data/bremsstrahlung_photons_distribution_poskus_2019.jld2"
    data = load(file)

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
            Ai = cubic_hermite_spline(E[iz][i-1:i],[A⁺,A⁻])(Ei)
            Bi = cubic_hermite_spline(E[iz][i-1:i],[B⁺,B⁻])(Ei)
            Ci = cubic_hermite_spline(E[iz][i-1:i],[C⁺,C⁻])(Ei)
        end
        return Ai,Bi,Ci
    end
end