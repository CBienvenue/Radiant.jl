"""
    Bremsstrahlung

Structure used to define parameters for production of multigroup bremsstrahlung cross-sections.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}} = Dict((Electron,Electron) => ["S"],(Electron,Photon) => ["P"],(Positron,Positron) => ["S"],(Positron,Photon) => ["P"])` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `(Electron,Electron) => ["S"]` : scattering of incident electron following Bremsstrahlung interaction.
    - `(Electron,Photon) => ["P"]` : produced photon following Bremsstrahlung interaction by incident electron.
    - `(Positron,Positron) => ["S"]` : scattering of incident positron following Bremsstrahlung interaction.
    - `(Positron,Photon) => ["P"]` : produced photon following Bremsstrahlung interaction by incident positron.
- `angular_scattering_type::String=modified_dipole` : type of angular scattering, which can takes the following values:
    - `angular_scattering_type = modified_dipole` : modified dipôle distribution, based on Poskus (2019) shape functions.
    - `angular_scattering_type = sommerfield` : Sommerfield distribution.

"""
mutable struct Bremsstrahlung <: Interaction

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
    Cℓk::Array{Float64}
    angular_scattering_type::String
    bremsstrahlung_cross_sections::Function
    bremsstrahlung_stopping_powers::Function
    bremsstrahlung_photons_distribution::Function
    scattering_model::String

    # Constructor(s)
    function Bremsstrahlung()
        this = new()
        this.name = "bremsstrahlung"
        this.interaction_types = Dict((Electron,Electron) => ["S"],(Electron,Photon) => ["P"],(Positron,Positron) => ["S"],(Positron,Photon) => ["P"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.is_CSD = true
        this.is_AFP = true
        this.is_AFP_decomposition = false
        this.is_elastic = false
        this.is_preload_data = true
        this.is_subshells_dependant = false
        this.set_angular_scattering_type("modified_dipole")
        this.set_scattering_model("BFP")
        return this
    end

end

# Method(s)
"""
    set_interaction_types(this::Bremsstrahlung,interaction_types::Dict{Tuple{DataType,DataType},Vector{String}})

To define the interaction types for bremsstrahlung processes.

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure.
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}}` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `(Electron,Electron) => ["S"]` : scattering of incident electron following Bremsstrahlung interaction.
    - `(Electron,Photon) => ["P"]` : produced photon following Bremsstrahlung interaction by incident electron.
    - `(Positron,Positron) => ["S"]` : scattering of incident positron following Bremsstrahlung interaction.
    - `(Positron,Photon) => ["P"]` : produced photon following Bremsstrahlung interaction by incident positron.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> bremsstrahlung = Bremsstrahlung()
julia> bremsstrahlung.set_interaction_types( Dict((Electron,Electron) => ["S"]) ) # Only electron scattering, with photon absorption.
```
"""
function set_interaction_types(this::Bremsstrahlung,interaction_types)
    this.interaction_types = interaction_types
end

"""
    set_angular_scattering_type(this::Bremsstrahlung,angular_scattering_type::String)

To define the bremsstrahlung photons angular distribution.

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure.
- `angular_scattering_type::String` : angular scattering type.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> bremsstrahlung = Bremsstrahlung()
julia> bremsstrahlung.set_angular_scattering_type("sommerfield")
```
"""
function set_angular_scattering_type(this::Bremsstrahlung,angular_scattering_type::String)
    if lowercase(angular_scattering_type) ∉ ["modified_dipole","sommerfield"] error("Undefined angular distribution.") end
    this.angular_scattering_type = lowercase(angular_scattering_type)
end

"""
    scattering_model(this::Bremsstrahlung,scattering_model::String)

To define the solver for bremsstrahlung scattering.

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure.
- `solver::String` : solver for bremsstrahlung scattering, which can be:
    - `BFP` : Boltzmann Fokker-Planck solver.
    - `FP` : Fokker-Planck solver.
    - `BTE` : Boltzmann solver.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> bremsstrahlung = Bremsstrahlung()
julia> bremsstrahlung.set_scattering_model("FP")
```
"""
function set_scattering_model(this::Bremsstrahlung,scattering_model::String)
    if uppercase(scattering_model) ∉ ["BFP","FP","BTE"] error("Unknown scattering model (should be BFP, FP or BTE).") end
    this.scattering_model = uppercase(scattering_model)
end

"""
    in_distribution(this::Bremsstrahlung)

Describe the energy discretization method for the incoming particle in the bremsstrahlung
interaction.

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function in_distribution(this::Bremsstrahlung)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    out_distribution(this::Bremsstrahlung)

Describe the energy discretization method for the outgoing particle in the bremsstrahlung
interaction.

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function out_distribution(this::Bremsstrahlung)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    bounds(this::Bremsstrahlung,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String,
    Ec::Float64)

Gives the integration energy bounds for the outgoing particle for bremsstrahlung. 

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure.
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `Ei::Float64` : energy of the incoming particle.
- `type::String` : type of interaction.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.

# Output Argument(s)
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `isSkip::Bool` : define if the integration is skipped or not.

"""
function bounds(this::Bremsstrahlung,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String,Ec::Float64)
    # Scattered electron or positron
    if type == "S" 
        Ef⁻ = min(Ef⁻,Ec)
        if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    # Produced photon
    elseif type == "P" 
        Ef⁻ = min(Ef⁻,Ei)
        if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    elseif type == "Pₐ"
        Ef⁻ = min(Ef⁻,Ei-Ec)
        if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    else
        error("Unknown type of method for Bremsstrahlung scattering.")
    end
    return Ef⁻,Ef⁺,isSkip
end

"""
    dcs(this::Bremsstrahlung,L::Int64,Ei::Float64,Ef::Float64,Z::Int64,particle::Particle,
    type::String,iz::Int64)

Gives the Legendre moments of the scattering cross-sections for bremsstrahlung. 

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : incoming particle energy.
- `Ef::Float64` : outgoing particle energy.
- `Z::Int64` : atomic number.
- `particle::Particle` : incoming particle.
- `type::String` : type of interaction.
- `iz::Int64` : index of the element in the material.

# Output Argument(s)
- `σℓ::Vector{Float64}` : Legendre moments of the scattering cross-sections.

"""
function dcs(this::Bremsstrahlung,L::Int64,Ei::Float64,Ef::Float64,Z::Int64,particle::Particle,type::String,iz::Int64)

    # Inititalisation
    β² = Ei*(Ei+2)/(Ei+1)^2
    β = sqrt(β²)
    σs = 0.0
    σℓ = zeros(L+1)

    # Correction for positrons
    Fp = 1
    if is_positron(particle)
        t = log(1+1e6/Z^2*Ei)
        Fp = 1 - exp(-1.2359e-1*t + 6.1274e-2*t^2-3.1516e-2*t^3+7.7446e-3*t^4-1.0595e-3*t^5+7.0568e-5*t^6-1.8080e-6*t^7)
    end

    # Scattered electron
    if type == "S" 
        Eγ = Ei-Ef
        if Ei ≥ Ef
            σs = Fp * this.bremsstrahlung_cross_sections(iz,Z,Ei,Eγ)
        end
        
        # Compute the Legendre moments of the flux, μ = 1.0 (Forward peaked)
        for ℓ in range(0,L) σℓ[ℓ+1] += σs end

    # Produced photons
    elseif type == "P"

        # Compute the differential scattering cross section
        Eγ = Ef
        if Ei ≥ Ef
            σs = Fp * this.bremsstrahlung_cross_sections(iz,Z,Ei,Eγ)
        end

        # Sommerfield angular distribution
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
            A,B,C = this.bremsstrahlung_photons_distribution(iz,Z,Ei,Eγ/Ei)
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
    else
        error("Unknown interaction.")
    end
    return σℓ
end

"""
    tcs(this::Bremsstrahlung,Ei::Float64,Z::Int64,Ec::Float64,iz::Int64,particle::Particle,
    type::String,Eout::Vector{Float64})

Gives the total cross-section for bremsstrahlung. 

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure. 
- `Ei::Float64` : incoming particle energy.
- `Z::Int64` : atomic number.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `iz::Int64` : index of the element in the material.
- `particle::Particle` : incoming particle.
- `type::String` : type of interaction.
- `Eout::Vector{Float64}` : outgoing energy boundaries.

# Output Argument(s)
- `σt::Float64` : total cross-section.

"""
function tcs(this::Bremsstrahlung,Ei::Float64,Z::Int64,Ec::Float64,iz::Int64,particle::Particle,Eout::Vector{Float64})

    # Inititalization
    β² = 1-1/(Ei+1)^2

    # Correction for positrons
    Fp = 1
    if is_positron(particle)
        t = log(1+1e6/Z^2*Ei)
        Fp = 1 - exp(-1.2359e-1*t + 6.1274e-2*t^2-3.1516e-2*t^3+7.7446e-3*t^4-1.0595e-3*t^5+7.0568e-5*t^6-1.8080e-6*t^7)
    end

    # Compute total cross section
    σt = 0.0
    Ngf = length(Eout)-1
    is_dirac, Np, q_type = out_distribution(this)
    if is_dirac Np = 1; u = [0]; w = [2] else u,w = quadrature(Np,q_type) end
    for gf in range(1,Ngf+1)
        Ef⁻ = Eout[gf]
        if (gf != Ngf+1) Ef⁺ = Eout[gf+1] else Ef⁺ = 0.0 end
        Ef⁻,Ef⁺,isSkip = bounds(this,Ef⁻,Ef⁺,Ei,"S",Ec)
        if isSkip continue end
        ΔEf = Ef⁻ - Ef⁺
        for n in range(1,Np)
            Ef = (u[n]*ΔEf + (Ef⁻+Ef⁺))/2
            Eγ = Ei-Ef
            if Ei ≥ Ef && ΔEf ≥ 0
                σt += ΔEf/2 * w[n] * Fp * this.bremsstrahlung_cross_sections(iz,Z,Ei,Eγ)
            end
        end
    end
    return σt
end

"""
    sp(this::Bremsstrahlung,Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,Ei::Float64,
    Ec::Float64,Eout::Vector{Float64},particle::Particle)

Gives the stopping power for bremsstrahlung.

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure. 
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `ωz::Vector{Float64}` : weight fraction of the elements composing the material.
- `ρ::Float64` : material density.
- `Ei::Float64` : incoming particle energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `Eout::Vector{Float64}` : energy boundaries associated with outgoing particles.
- `particle::Particle` : incoming particle.

# Output Argument(s)
- `S::Float64` : stopping power.

"""
function sp(this::Bremsstrahlung,Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,Ei::Float64,Ec::Float64,Eout::Vector{Float64},particle::Particle)

    # Initialization
    rₑ = 2.81794092e-13 # (in cm)
    β² = Ei*(Ei+2)/(Ei+1)^2
    α = 1/137
    𝒩ₙ = nuclei_density.(Z,ρ)
    Nz = length(Z)

    # Compute the total stopping power 
    St = 0.0
    for iz in range(1,Nz)
        Fp = 1
        if is_positron(particle)
            t = log(1+1e6/Z[iz]^2*Ei)
            Fp = 1 - exp(-1.2359e-1*t+6.1274e-2*t^2-3.1516e-2*t^3+7.7446e-3*t^4-1.0595e-3*t^5+7.0568e-5*t^6-1.8080e-6*t^7)
        end
        St += ωz[iz] * 𝒩ₙ[iz] * Fp * this.bremsstrahlung_stopping_powers(iz,Z[iz],Ei)
    end

    # Compute the catastrophic stopping power
    Sc = 0.0
    Ngf = length(Eout)-1
    is_dirac, Np, q_type = out_distribution(this)
    if is_dirac Np = 1; u = [0]; w = [2] else u,w = quadrature(Np,q_type) end
    for iz in range(1,Nz)
        Fp = 1
        if is_positron(particle)
            t = log(1+1e6/Z[iz]^2*Ei)
            Fp = 1 - exp(-1.2359e-1*t+6.1274e-2*t^2-3.1516e-2*t^3+7.7446e-3*t^4-1.0595e-3*t^5+7.0568e-5*t^6-1.8080e-6*t^7)
        end
        for gf in range(1,Ngf+1)
            Ef⁻ = Eout[gf]
            if (gf != Ngf+1) Ef⁺ = Eout[gf+1] else Ef⁺ = 0.0 end
            Ef⁻,Ef⁺,isSkip = bounds(this,Ef⁻,Ef⁺,Ei,"S",Ec)
            if isSkip continue end
            ΔEf = Ef⁻ - Ef⁺
            for n in range(1,Np)
                Ef = (u[n]*ΔEf + (Ef⁻+Ef⁺))/2
                Eγ = Ei-Ef
                if Ei ≥ Ef && ΔEf ≥ 0
                    Sc += ωz[iz] * ΔEf/2 * w[n] * 𝒩ₙ[iz] * Eγ * Fp * this.bremsstrahlung_cross_sections(iz,Z[iz],Ei,Eγ)
                end
            end
        end
    end

    # Compute the soft stopping powers
    S = St-Sc
    return S
end

"""
    mt(this::Bremsstrahlung)

Gives the momentum transfer for bremsstrahlung.

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure. 

# Output Argument(s)
- `T::Float64` : momentum transfer.

"""
function mt(this::Bremsstrahlung)
    T = 0.0 # Because μ = 1.0
    return T
end

"""
    preload_data(this::Bremsstrahlung,Z::Vector{Int64},Emax::Float64,Emin::Float64,L::Int64)

Preload data for multigroup bremsstrahlung calculations.

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure. 
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `Emax::Float64` : maximum energy of the incoming particle.
- `Emin::Float64` : minimum energy of the incoming particle.
- `L::Int64` : Legendre truncation order.

# Output Argument(s)
N/A

"""
function preload_data(this::Bremsstrahlung,Z::Vector{Int64},Emax::Float64,Emin::Float64,L::Int64)

    # Preload cross-sections and stopping powers functions
    this.preload_bremsstrahlung_cross_sections(Z,Emax,Emin)
    this.preload_bremsstrahlung_stopping_powers(Z,Emax,Emin)

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
    preload_bremsstrahlung_cross_sections(this::Bremsstrahlung,Z::Vector{Int64},
    Emax::Float64,Emin::Float64)

Preload data for multigroup bremsstrahlung cross-sections.

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure. 
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `Emax::Float64` : maximum energy of the incoming particle.
- `Emin::Float64` : minimum energy of the incoming particle.

# Output Argument(s)
N/A

"""
function preload_bremsstrahlung_cross_sections(this::Bremsstrahlung,Z::Vector{Int64},Emax::Float64,Emin::Float64)

    # Initialize
    χ = Vector{Array{Float64}}()

    # Extract vectors
    path = joinpath(find_package_root(), "data", "bremsstrahlung_data_seltzer_berger_1986.jld2")
    data = load(path)
    E = data["incident_electron_energy"] # (in mₑc²)
    r = data["radiative_energy_fraction"]
    if ( Emin < E[1] || Emax > E[end]) error("Energy less than 1 keV or more than 10 GeV not tabulated.") end

    # Define the interval corresponding to the required interpolation data for calculations
    i⁻ = max(searchsortedfirst(E,Emin) - 1,1)
    i⁺ = searchsortedfirst(E,Emax)
    E = E[i⁻:i⁺]

    # Extract scaled cross-sections
    Nz = length(Z)
    for iz in range(1,Nz)
        if ~haskey(data["scaled_cross_sections"],Z[iz]) error(string("Undefined scaled Bremsstrahlung cross-sections for Z = ",Z[iz])) end
        χi = data["scaled_cross_sections"][Z[iz]]
        push!(χ,χi[i⁻:i⁺,:])
    end

    # Return the interpolation function
    this.bremsstrahlung_cross_sections = function bremsstrahlung_cross_sections(iz::Int64,Z::Int64,Ei::Float64,Eγ::Float64)
        ri = Eγ/Ei
        β² = Ei*(Ei+2)/(Ei+1)^2
        if E[1] ≤ Ei ≤ E[end] && r[1] ≤ ri ≤ r[end]

            # Find index in E and r
            i = searchsortedfirst(E,Ei)
            j = searchsortedfirst(r,ri)
            if (i == 1 && E[1] ≈ Ei) i=2; Ei=E[1] end
            if (j == 1 && r[1] ≈ ri) j=2; ri=r[1] end
            if (i == length(E)+1 && E[end] ≈ Ei) i=length(E); Ei=E[end] end
            if (j == length(r)+1 && r[end] ≈ ri) j=length(r); ri=r[end] end
            if (i == 1 || i > length(E)) error("Interpolation value is outside the interpolation vector.") end
            if (j == 1 || j > length(r)) error("Interpolation value is outside the interpolation vector.") end
            χij⁻ = cubic_hermite_spline(E[i-1:i],χ[iz][i-1:i,j-1])(Ei)
            χij⁺ = cubic_hermite_spline(E[i-1:i],χ[iz][i-1:i,j])(Ei)
            χij = cubic_hermite_spline(r[j-1:j],[χij⁻,χij⁺])(ri)

            σs = χij/1e27/Eγ * Z^2/β²
        else
            σs = 0
        end
        return σs
    end
end

"""
    preload_bremsstrahlung_stopping_powers(this::Bremsstrahlung,Z::Vector{Int64},
    Emax::Float64,Emin::Float64)

Preload data for multigroup bremsstrahlung stopping powers.

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure. 
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `Emax::Float64` : maximum energy of the incoming particle.
- `Emin::Float64` : minimum energy of the incoming particle.

# Output Argument(s)
N/A

"""
function preload_bremsstrahlung_stopping_powers(this::Bremsstrahlung,Z::Vector{Int64},Emax::Float64,Emin::Float64)

    # Initialize
    ϕ = Vector{Vector{Float64}}()

    # Extract vectors
    path = joinpath(find_package_root(), "data", "bremsstrahlung_data_seltzer_berger_1986.jld2")
    data = load(path)
    E = data["incident_electron_energy"] # (in mₑc²)
    if ( Emin < E[1] || Emax > E[end]) error("Energy less than 1 keV or more than 10 GeV not tabulated.") end

    # Define the interval corresponding to the required interpolation data for calculations
    i⁻ = max(searchsortedfirst(E,Emin) - 1,1)
    i⁺ = searchsortedfirst(E,Emax)
    E = E[i⁻:i⁺]

    # Extract scaled cross-sections
    Nz = length(Z)
    for iz in range(1,Nz)
        if ~haskey(data["scaled_stopping_powers"],Z[iz]) error(string("Undefined radiative stopping power for Z = ",Z[iz])) end
        ϕi = data["scaled_stopping_powers"][Z[iz]]
        push!(ϕ,ϕi[i⁻:i⁺])
    end

    # Cubic hermite spline
    ϕ_spline = Vector{Function}(undef,Nz)
    for iz in range(1,Nz)
        ϕ_spline[iz] = cubic_hermite_spline(E,ϕ[iz])
    end

    # Return the interpolation function
    this.bremsstrahlung_stopping_powers = function bremsstrahlung_stopping_powers(iz::Int64,Z::Int64,Ei::Float64)
        α = 1/137
        rₑ = 2.81794092e-13 # (in cm)
        if E[1] ≤ Ei ≤ E[end]
            s = ϕ_spline[iz](Ei) * α * rₑ^2 * Z^2 * (Ei+1)
        else
            s = 0
        end
        return s
    end

end

"""
    preload_angular_distribution(this::Bremsstrahlung,Z::Vector{Int64})

Preload data for bremsstrahlung photon angular distribution.

# Input Argument(s)
- `this::Bremsstrahlung` : bremsstrahlung structure. 
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.

# Output Argument(s)
N/A

"""
function preload_angular_distribution(this::Bremsstrahlung,Z::Vector{Int64})

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
    this.bremsstrahlung_photons_distribution = function bremsstrahlung_photons_distribution(iz::Int64,Z::Int64,Ei::Float64,ri::Float64)
        mₑc² = 0.510999
        if Ei ≥ 3/mₑc²
            β = sqrt(Ei*(Ei+2))/(Ei+1)
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