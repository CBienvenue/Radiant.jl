"""
    ECUT_DEFAULTS

Energy above which the corrected Bethe formula is used, per particle [in MeV]. `bethe`
reads this unless the caller overrides it, and `Inelastic_Collision` seeds `set_ecut` from
it, so the defaults are written down once.

"""
const ECUT_DEFAULTS = Dict{Type,Float64}(
    Electron => 0.001, Positron => 0.001,   # 1 keV
    Muon     => 0.15,  Antimuon => 0.15,    # 150 keV
    Proton   => 0.75,  Antiproton => 0.75,  # 750 keV
    Alpha    => 5.0                         # 5 MeV
)

"""
    EMAX_DEFAULTS

Energy of the stopping-power maximum used by the analytic low-energy extension, per
particle [in MeV]. Read only by that extension, and seeded into `set_emax` the same way.

"""
const EMAX_DEFAULTS = Dict{Type,Float64}(
    Electron => 0.0001, Positron => 0.0001,  # 100 eV
    Muon     => 0.02,   Antimuon => 0.02,    # 20 keV
    Proton   => 0.05,   Antiproton => 0.05,  # 50 keV
    Alpha    => 0.8                          # 800 keV
)

"""
    bethe(Z::Vector{Int64},ωz::Vector{Float64},atz::Vector{Float64},ρ::Float64,
    N_density::Float64,Ei::Float64,particle::Particle,density_correction::String,
    state_of_matter::String,I_eff::Float64;
    atomic_weights::Union{Nothing,Vector{Float64}}=nothing,user_ecut::Float64=NaN,
    user_emax::Float64=NaN,is_emax_fit::Bool=false,low_energy_model::String="analytic",
    low_energy_source::String="srim",is_bragg_additivity::Bool=false,
    is_bethe_correction::Bool=true,is_density_correction::Bool=true,
    is_shell_correction::Bool=true,is_lindhard_sorensen::Bool=true,is_barkas::Bool=true)

Gives the collisional stopping power using Bethe corrected formula for electron, positron
and heavy ions (proton, muon, alpha particle).

# Input Argument(s)
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `ωz::Vector{Float64}` : weight fraction of the elements composing the material.
- `atz::Vector{Float64}` : atomic fraction of the elements composing the material.
- `ρ::Float64` : material density.
- `N_density::Float64` : atom density of the material [in cm⁻³].
- `Ei::Float64` : incoming particle energy.
- `particle::Particle` : incoming particle.
- `density_correction::String` : density-effect correction type.
- `state_of_matter::String` : material state.
- `I_eff::Float64` : user-defined mean excitation energy override [in mₑc²]; NaN ⟹ table/additivity.
- `atomic_weights::Union{Nothing,Vector{Float64}}` : element atomic weights [u].
- `is_bethe_correction::Bool` : boolean to enable or not the Bethe correction.
- `is_density_correction::Bool` : boolean to enable or not the density correction.
- `is_shell_correction::Bool` : boolean to enable or not the shell correction.
- `is_lindhard_sorensen::Bool` : boolean to enable or not the Lindhard-Sorensen correction.
- `is_barkas::Bool` : boolean to enable or not the Barkas correction.
- `user_ecut::Float64` : energy above which the corrected Bethe formula is used [in MeV];
  NaN ⟹ the per-particle default.
- `user_emax::Float64` : energy of the stopping-power maximum of the analytic low-energy
  extension [in MeV]; NaN ⟹ the per-particle default.
- `is_emax_fit::Bool` : take that energy from the Z-dependent fit
  `emax_stopping_maximum` instead, which overrides `user_emax`.
- `low_energy_model::String` : model used below `Ecut`, `"analytic"` or `"icru49"`.
- `low_energy_source::String` : coefficient source of the fitted low-energy model.
- `is_bragg_additivity::Bool` : boolean to build a compound from its elements rather than
  from its own composition.

# Output Argument(s)
- `S::Float64` : stopping power.

# Reference(s)
- Salvat (2022), Bethe stopping-power formula and its corrections.
- Salvat and Andreo (2023), SBETHE: Stopping powers of materials for swift charged
  particles from the corrected Bethe formula.

"""
function bethe(Z::Vector{Int64},ωz::Vector{Float64},atz::Vector{Float64},ρ::Float64,N_density::Float64,Ei::Float64,particle::Particle,density_correction::String="fano",state_of_matter::String="solid",I_eff::Float64=NaN; atomic_weights::Union{Nothing,Vector{Float64}}=nothing,is_bethe_correction::Bool=true,is_density_correction::Bool=true,is_shell_correction::Bool=true,is_lindhard_sorensen::Bool=true,is_barkas::Bool=true,user_ecut::Float64=NaN,user_emax::Float64=NaN,is_emax_fit::Bool=false,low_energy_model::String="analytic",low_energy_source::String="srim",is_bragg_additivity::Bool=false)

    #----
    # Set the minimum energy cutoff over which the corrected Bethe formula is valid.
    #----
    mₑ = 0.51099895069
    particle_type = get_type(particle)
    if !haskey(ECUT_DEFAULTS,particle_type) error("Unknown particle.") end
    Ecut = ECUT_DEFAULTS[particle_type]/mₑ
    Emax = EMAX_DEFAULTS[particle_type]/mₑ

    # A user-supplied cutoff overrides the default. Emax is read only by the analytic
    # low-energy extension and means nothing at or above Ecut, so nothing there changes.
    if !isnan(user_ecut) Ecut = user_ecut/mₑ end
    if is_emax_fit
        Emax = emax_stopping_maximum(Z,atz,particle)/mₑ
    elseif !isnan(user_emax)
        Emax = user_emax/mₑ
    end

    #----
    # Bragg additivity
    #----
    # Build the compound by evaluating each element on its own and combining with mass
    # fractions. Mass stopping power is additive in mass fractions; this is a different
    # question from the Eq. (53) electron weighting inside the bracket.
    #
    # ρ is the compound's density, not the element's: the density effect in a compound is
    # set by the compound's electron density, not by what graphite would do on its own.
    # That also means no elemental density table is needed, and since every element is
    # evaluated at the same density the mass-stopping-power conversion cancels and the
    # combination is a plain sum over the linear stopping powers bethe returns.
    #
    # N_density is per element at the compound's mass density. I_eff is NaN so each
    # element carries its own tabulated mean excitation energy, which is the whole point.
    # The recursion passes is_bragg_additivity = false, or it would not terminate.
    if is_bragg_additivity && length(Z) > 1
        element_atomic_weights = isnothing(atomic_weights) ? standard_atomic_weight.(Z) : atomic_weights
        S = 0.0
        for j in eachindex(Z)
            S += ωz[j] * bethe([Z[j]],[1.0],[1.0],ρ,nuclei_density(element_atomic_weights[j],ρ),Ei,particle,density_correction,state_of_matter,NaN;
                atomic_weights=[element_atomic_weights[j]],
                is_bethe_correction=is_bethe_correction,
                is_density_correction=is_density_correction,
                is_shell_correction=is_shell_correction,
                is_lindhard_sorensen=is_lindhard_sorensen,
                is_barkas=is_barkas,
                user_ecut=user_ecut,
                user_emax=user_emax,
                is_emax_fit=is_emax_fit,
                low_energy_model=low_energy_model,
                low_energy_source=low_energy_source,
                is_bragg_additivity=false)
        end
        return S
    end

    #----
    # Corrected Bethe formula
    #----
    if Ei ≥ Ecut

        #----
        # Initialization
        #----
        α = 1/137.035999177
        rₑ = 2.81794092e-13 # (in cm)
        ratio_mass = get_mass(particle)/mₑ
        γ = (Ei + ratio_mass)/ratio_mass
        β² = (γ^2-1)/γ^2
        charge = get_charge(particle)
        is_heavy = get_mass(particle) > 10*mₑ # (Is particle much heavier than the electron?)

        #----
        # Effective electron density and ionization energy
        #----
        𝒩ₑ_eff = N_density * sum(atz .* Z) # (in cm⁻³)
        I = effective_mean_excitation_energy(Z,ωz,I_eff) # (in mₑc²)

        #----
        # Bethe correction factor
        #----
        f = 0
        if is_bethe_correction
            if is_electron(particle)
                f = (2*γ^2-1)/γ^2 + ((γ-1)/γ)^2/8 - (4-((γ-1)/γ)^2)*log(2) - log(γ+1)
            elseif is_positron(particle)
                f = (γ^2-1)/(12*γ^2) * (1 - 14/(γ+1)-10/(γ+1)^2-4/(γ+1)^3) - log(2) - log(γ+1)
            elseif is_heavy
                R = 1/(1+(1/ratio_mass)^2+2*γ/ratio_mass)
                f = log(R) + (1/ratio_mass * (γ^2-1)/γ * R)^2
            else
                error("Bethe correction factor is not defined for the given particle.")
            end
        end

        #----
        # Density effect
        #----
        δF = 0
        if is_density_correction
            δF = fermi_density_effect(Z,atz,𝒩ₑ_eff,Ei,state_of_matter,density_correction,ratio_mass,I)
        end

        #----
        # Shell correction
        #----
        Cz = 0
        if is_shell_correction
            Cz = shell_correction(Z,atz,γ,particle)
        end

        #----
        # Lindhard-Sorensen correction
        #----
        ΔL_LS = 0
        if is_lindhard_sorensen && is_heavy
            if charge == 2
                A = 90.59
            elseif charge == 1
                A = 180.20
            elseif charge == -1
                A = -178.34
            elseif charge == -2
                A = -88.73
            else
                error("Lindhard-Sorensen correction not valid for particle with charge |Z| > 2.")
            end
            η = charge*α/sqrt(β²)
            L_Bloch = 0; L_Bloch⁻ = 0
            ϵ = Inf
            N = 1
            while ϵ > 1e-5
                L_Bloch += 1/(N*(N^2+η^2))
                ϵ = abs((L_Bloch - L_Bloch⁻)/max(L_Bloch⁻,1e-5))
                L_Bloch⁻ = copy(L_Bloch)
                N += 1
            end
            L_Bloch = - η^2/charge^2 * L_Bloch
            ΔL_LS = ((1+A)/(1+1.92*(γ-1)^1.41)-A) * charge^2 * L_Bloch
        end

        #----
        # Barkas correction
        #----
        ΔL_B = 0
        if is_barkas && is_heavy
            ΔL_B = barkas_correction(Z,atz,𝒩ₑ_eff,γ,particle,I)
        end

        #----
        # To take into account electron capture by positively charged ions
        #----
        Zc = charge
        if is_proton(particle) || is_alpha(particle)
            Zc *= (1-exp(-sqrt(β²)/(α*charge)))
        end

        #----
        # Total stopping power
        #----
        S = 4*π*rₑ^2/β² * Zc^2 * 𝒩ₑ_eff * ( log(2*β²/I) + log(γ^2) - β² + f/2 - Cz - δF/2 + ΔL_LS + ΔL_B )

        #----
        # Handover window, from Ecut up to 2*Ecut
        #----
        # The fitted form is not tied to the Bethe branch at Ecut, so without blending the
        # handover leaves a discontinuity of up to 13.6%, which multigroup.jl would carry
        # straight into the CSD operator. The window sits above Ecut, where S is already in
        # hand: the Bethe branch is never evaluated below its own cutoff. E_top is a
        # multiple of Ecut, not a constant, so set_ecut moves both ends together and the
        # window can never invert. The guard is the same predicate the branch below Ecut
        # uses, or the two would disagree about a material.
        E_MeV = Ei*mₑ
        Ecut_MeV = Ecut*mₑ
        E_top = 2*Ecut_MeV
        if low_energy_model != "analytic" && E_MeV < E_top && has_low_energy_data(Z,particle)
            w = blend_weight(E_MeV,Ecut_MeV,E_top)
            S_fit = icru49_stopping(Z,ωz,E_MeV,particle,state_of_matter,low_energy_source) * ρ/mₑ
            S = (1-w)*S_fit + w*S
        end

    #----
    # Low-energy branch
    #----
    else
        S = low_energy_stopping(Z,ωz,atz,ρ,N_density,Ei,particle,density_correction,state_of_matter,I_eff,Ecut,Emax;
            atomic_weights=atomic_weights,
            is_bethe_correction=is_bethe_correction,
            is_density_correction=is_density_correction,
            is_shell_correction=is_shell_correction,
            is_lindhard_sorensen=is_lindhard_sorensen,
            is_barkas=is_barkas,
            user_ecut=user_ecut,
            low_energy_model=low_energy_model,
            low_energy_source=low_energy_source)
    end
    return S
end

"""
    shell_correction(Z::Vector{Int64},atz::Vector{Float64},γi::Float64,particle::Particle)

Gives the shell correction to the Bethe formula.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `atz::Vector{Float64}` : atomic fraction of the elements composing the material.
- `γi::Float64` : total energy of the incoming particle in unit of its rest energy.
- `particle::Particle` : incoming particle.

# Output Argument(s)
- `Cz::Float64` : shell correction.

# Reference(s)
- Salvat (2022), Bethe stopping-power formula and its corrections.
- Salvat and Andreo (2023), SBETHE: Stopping powers of materials for swift charged
  particles from the corrected Bethe formula.

"""
function shell_correction(Z::Vector{Int64},atz::Vector{Float64},γi::Float64,particle::Particle)

    #----
    # Initialization
    #----
    mₑ = 0.51099895069
    is_heavy = get_mass(particle) > 10*mₑ # (Is particle much heavier than the electron?)
    Nz = length(Z)
    S₀ = zeros(Nz)
    Ec = zeros(Nz)
    pn = zeros(Nz,6)
    γ = zeros(Nz,153)
    Cz_prime = zeros(Nz,153)
    spline_Cz = Vector{Function}(undef,Nz)
    if is_electron(particle)
        particle_name = "electron"
    elseif is_positron(particle)
        particle_name = "positron"
    elseif is_heavy
        particle_name = "heavy_particle"
    else
        error("Undefined particle.")
    end

    #----
    # Read shell correction data
    #----

    # Extract data
    data = fast_load("shell_corrections_salvat_2023.jld2")

    # Extract parameters
    for iz in range(1,Nz)
        datai = data[particle_name][Z[iz]]
        for n in range(1,153)
            Cz_prime[iz,n] = datai["modified_shell_corrections"][n]
            γ[iz,n] = datai["gamma"][n]
        end
        if γi < minimum(γ[iz,:]) || γi > maximum(γ[iz,:])
            error("Shell-correction table for Z = $(Z[iz]) does not cover γ = $γi " *
                  "(tabulated range $(minimum(γ[iz,:])) to $(maximum(γ[iz,:]))). " *
                  "If Ecut was lowered with set_ecut, raise it.")
        end
        Ec[iz] = datai["cutoff_energy"]
        for n in range(1,6)
            pn[iz,n] = datai["fit_parameters"][n]
        end
        S₀[iz] = datai["integrated_gos"]
        spline_Cz[iz] = cubic_hermite_spline(γ[iz,:],Cz_prime[iz,:])
    end

    #----
    # Compute the shell correction
    #----
    Cz = 0.0
    β² = (γi^2-1)/γi^2
    Nz = length(Z)

    # Stoichiometric electron fraction of each element, SBETHE Eq. (53):
    # C(γ)/Z = (1/Z) Σⱼ nⱼZⱼ Cⱼ(γ)/Zⱼ with Z = Σⱼ nⱼZⱼ. The mixing rule takes atom
    # fractions, not mass fractions.
    g = atz .* Z ./ sum(atz .* Z)
    for iz in range(1,Nz)
        if γi < maximum(γ[iz,:])
            Czprime = spline_Cz[iz](γi)
        elseif γi < 2
            Czprime = sum(pn[iz,:] .* (γi-1).^(1:6))
        else
            Czprime = sum(pn[iz,:] .* 2 .^(1:6))
        end
        Cz += g[iz] * (Czprime + (Z[iz]-S₀[iz])/(2*Z[iz])*(log(β²/(1-β²))-β²))
    end
    return Cz
end

"""
    barkas_correction(Z::Vector{Int64},atz::Vector{Float64},𝒩ₑ_eff::Float64,
    γ::Float64,particle::Particle,I::Float64)

Gives the Barkas correction to the Bethe formula.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `atz::Vector{Float64}` : atomic fraction of the elements composing the material.
- `𝒩ₑ_eff::Float64` : effective electron density of the material [in cm⁻³], already
  built by the caller, so that the Bethe bracket and the plasma energy cannot drift
  apart.
- `γ::Float64` : total energy of the incoming particle in unit of its rest energy.
- `particle::Particle` : incoming particle.
- `I::Float64` : mean excitation energy of the material [in mₑc²], already resolved by the
  caller. Passing it in is what makes the Barkas term honour a user-set value; resolving it
  here would silently ignore `set_mean_excitation_energy`.

# Output Argument(s)
- `ΔLB::Float64` : Barkas correction.

# Reference(s)
- Salvat (2022), Bethe stopping-power formula and its corrections.
- Salvat and Andreo (2023), SBETHE: Stopping powers of materials for swift charged
  particles from the corrected Bethe formula.

"""
function barkas_correction(Z::Vector{Int64},atz::Vector{Float64},𝒩ₑ_eff::Float64,γ::Float64,particle::Particle,I::Float64)

    #----
    # Initialization
    #----
    α = 1/137.035999177
    β² = (γ^2-1)/γ^2
    Nz = length(Z)

    # SBETHE Eq. (54): CB = max{1, Z̄/10} with Z̄ = Z(Σⱼnⱼ)⁻¹, the mean atomic number
    # per atom. Atom fractions sum to unity, so Z̄ is their Z-weighted sum, which is also
    # the number of electrons on an average atom of the material.
    Zbar = sum(atz .* Z)
    CB = max(1,Zbar/10)
    charge = get_charge(particle)
    Ωp = plasma_energy(𝒩ₑ_eff)

    #----
    # Extract the subshell of all atoms of the compound
    #----
    Zi = zeros(0); Ui = zeros(0); Wi = zeros(0)
    for iz in range(1,Nz)
        _,Ze,Ue,_,_,_ = electron_subshells(Z[iz])
        We = resonance_energy(Z[iz],Ωp,I)
        append!(Zi,(atz[iz]/Zbar).*Ze)
        append!(Ui,Ue)
        append!(Wi,We)
    end

    #----
    # Compute the Barkas correction
    #----
    ΔLB = 0
    for δi in range(1,length(Zi))
        ξ = 0.5616 * CB * Wi[δi]/(γ*β²)
        ΔLB += charge*α/(γ^2*β²^(3/2)) * Zi[δi] * Wi[δi] * (I₁_barkas(ξ) + I₂_barkas(ξ)/γ^2)
    end
    return ΔLB
end

"""
    I₁_barkas(x::Float64)

Approximation of integral I₁ (Eq. 96a) from Salvat (2022).

# Input Argument(s)
- `x::Float64` : evaluation point.

# Output Argument(s)
- `I₁_barkas::Float64` : evaluated integral.

# Reference(s)
- Salvat (2022), Bethe stopping-power formula and its corrections.
- Salvat and Andreo (2023), SBETHE: Stopping powers of materials for swift charged
  particles from the corrected Bethe formula.

"""
function I₁_barkas(x::Float64)
    if x < 1.0e-10
        return 1.5 * π * log(0.375 / x)
    elseif x < 0.25
        p = [17.96676, 48.25194, 46.61380, 22.69807, 5.40918, 0.5116173]
        xl = log(x)
        return 1.5 * π * log(0.375 / x) - 2.0 * π * x^2 * (p[1] + xl * (p[2] + xl * (p[3] + xl * (p[4] + xl * (p[5] + xl * p[6])))))
    elseif x < 2.0
        p = [-3.062544, 17.87672, -41.00572, 46.52797, -27.28692, 8.070905, -0.9587701]
        xr = sqrt(1.0 / x)
        return (p[1] + xr * (p[2] + xr * (p[3] + xr * (p[4] + xr * (p[5] + xr * (p[6] + xr * p[7])))))) * xr^1.25
    elseif x < 15.01
        p = [-0.001018033, 0.1080448, 1.568923, -2.77911, 8.07993, -7.244712]
        xr = 1.0 / x
        return (p[1] + xr * (p[2] + xr * (p[3] + xr * (p[4] + xr * (p[5] + xr * p[6]))))) * exp(-2.0 * x)
    elseif x < 50.01
        p = [-27.98361, 3.964525, -0.5706711, 0.02652851, -0.0006442942, 8.003948e-6, -4.016094e-8]
        return exp(p[1] + x * (p[2] + x * (p[3] + x * (p[4] + x * (p[5] + x * (p[6] + x * p[7]))))))
    else
        return 0.0
    end
end

"""
    I₂_barkas(x::Float64)

Approximation of integral I₂ (Eq. 96b) from Salvat (2022).

# Input Argument(s)
- `x::Float64` : evaluation point.

# Output Argument(s)
- `I₂_barkas::Float64` : evaluated integral.

# Reference(s)
- Salvat (2022), Bethe stopping-power formula and its corrections.
- Salvat and Andreo (2023), SBETHE: Stopping powers of materials for swift charged
  particles from the corrected Bethe formula.

"""
function I₂_barkas(x::Float64)
    if x < 1.0e-9
        return 2.17759
    elseif x < 0.25
        p = [2.17759, 0.5689823, 1.038828, 0.0002877808]
        xl = log(x)^2
        return p[1] - 2.0 * π * x^2 * (p[2] + xl * (p[3] + xl * p[4]))
    elseif x < 2.0
        p = [3.768512, -21.35702, 46.01060, -47.14780, 25.51946, -7.102352, 0.804296]
        xr = sqrt(1.0 / x)
        return (p[1] + xr * (p[2] + xr * (p[3] + xr * (p[4] + xr * (p[5] + xr * (p[6] + xr * p[7])))))) * xr^2 * exp(-1.5 * x)
    elseif x < 15.01
        p = [7.431717e-5, 0.06662051, 2.14271, -7.407167, 23.27532, -34.97742, 19.50319]
        xr = 1.0 / x
        return (p[1] + xr * (p[2] + xr * (p[3] + xr * (p[4] + xr * (p[5] + xr * (p[6] + xr * p[7])))))) * exp(-2.0 * x)
    elseif x < 50.01
        p = [-28.08307, 3.977364, -0.5714987, 0.02655846, -0.0006449171, 8.010890e-6, -4.019306e-8]
        return exp(p[1] + x * (p[2] + x * (p[3] + x * (p[4] + x * (p[5] + x * (p[6] + x * p[7]))))))
    else
        return 0.0
    end
end