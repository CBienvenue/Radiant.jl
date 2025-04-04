"""
    bethe(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,Ei::Float64,particle::Particle,
    is_bethe_correction::Bool=true,is_density_correction::Bool=true,
    is_shell_correction::Bool=true,is_lindhard_sorensen::Bool=true,is_barkas::Bool=true)

Gives the collisional stopping power using Bethe corrected formula for electron, positron
and heavy ions (proton, muon, alpha particle).

# Input Argument(s)
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `ωz::Vector{Float64}` : weight fraction of the elements composing the material.
- `ρ::Float64` : material density.
- `Ei::Float64` : incoming particle energy.
- `particle::Particle` : incoming particle.
- `is_bethe_correction::Bool` : boolean to enable or not the Bethe correction.
- `is_density_correction::Bool` : boolean to enable or not the density correction.
- `is_shell_correction::Bool` : boolean to enable or not the shell correction.
- `is_lindhard_sorensen::Bool` : boolean to enable or not the Lindhard-Sorensen correction.
- `is_barkas::Bool` : boolean to enable or not the Barkas correction.

# Output Argument(s)
- `S::Float64` : stopping power.

# Reference(s)
- Salvat (2022), Bethe stopping-power formula and its corrections.
- Salvat and Andreo (2023), SBETHE: Stopping powers of materials for swift charged
  particles from the corrected Bethe formula.

"""
function bethe(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,Ei::Float64,particle::Particle,is_bethe_correction::Bool=true,is_density_correction::Bool=true,is_shell_correction::Bool=true,is_lindhard_sorensen::Bool=true,is_barkas::Bool=true)

    #----
    # Set the minimum energy cutoff over which the corrected Bethe formula is valid. 
    #----
    mₑ = 0.51099895069
    if is_electron(particle) || is_positron(particle)
        Ecut = 0.001/mₑ  # 1 keV
        Emax = 0.0001/mₑ # 100 eV
    elseif is_muon(particle) || is_antimuon(particle)
        Ecut = 0.15/mₑ # 150 keV
        Emax = 0.02/mₑ # 20 keV
    elseif is_proton(particle) || is_antiproton(particle)
        Ecut = 0.75/mₑ # 750 keV
        Emax = 0.05/mₑ # 50 keV
    elseif is_alpha(particle)
        Ecut = 5/mₑ   # 5 MeV
        Emax = 0.8/mₑ # 800 keV
    else
        error("Unknown particle.")
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
        # Effective density and ionization energy
        #----
        𝒩ₑ_eff = sum(ωz.*Z.*nuclei_density.(Z,ρ)) # (in cm⁻³)
        I = effective_mean_excitation_energy(Z,ωz) # (in mₑc²)

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
            δF = fermi_density_effect(Z,ωz,ρ,Ei,"solid","fano",ratio_mass)
        end

        #----
        # Shell correction
        #----
        Cz = 0
        if is_shell_correction 
            Cz = shell_correction(Z,ωz,γ,particle)
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
            ΔL_B = barkas_correction(Z,ωz,ρ,γ,particle)
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
    # Analytical extension to low-energies
    #----
    else

        #----
        # Estimate parameters A and B
        #----
        h = 0.00005
        Sc = bethe(Z,ωz,ρ,Ecut,particle)
        Sc_plus_h = bethe(Z,ωz,ρ,Ecut+h,particle)
        dSdE = (Sc_plus_h - Sc)/h
        A = (3*log(Sc)*Sc - 2*dSdE*Ecut*log(Ecut/Emax))/(3*Sc)
        B = -2*Ecut*dSdE/(3*Sc*sqrt(log(Ecut/Emax)))

        #----
        # Compute total cross-sections using analytical formula
        #----
        t = Ei/Emax
        if t ≥ 1
            S = exp(A-B*(log(t))^1.5)
        else
            S = 1.5*exp(A)*sqrt(t)*(1-t/3)
        end
    end
    return S
end

"""
    shell_correction(Z::Vector{Int64},ωz::Vector{Float64},γi::Float64,particle::Particle)

Gives the shell correction to the Bethe formula.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `ωz::Vector{Float64}` : weight fraction of the elements composing the material.
- `ρ::Float64` : material density.
- `γi::Float64` : total energy of the incoming particle in unit of its rest energy.
- `particle::Particle` : incoming particle.

# Output Argument(s)
- `Cz::Float64` : shell correction.

# Reference(s)
- Salvat (2022), Bethe stopping-power formula and its corrections.
- Salvat and Andreo (2023), SBETHE: Stopping powers of materials for swift charged
  particles from the corrected Bethe formula.

"""
function shell_correction(Z::Vector{Int64},ωz::Vector{Float64},γi::Float64,particle::Particle)

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
    if ~haskey(cache_files[],"shell_corrections_salvat_2023")
        cache_files[]["shell_corrections_salvat_2023"] = load(joinpath(find_package_root(), "data", "shell_corrections_salvat_2023.jld2"))
    end
    data = cache_files[]["shell_corrections_salvat_2023"]
    for iz in range(1,Nz)
        datai = data[particle_name][Z[iz]]
        for n in range(1,153)
            Cz_prime[iz,n] = datai["modified_shell_corrections"][n]
            γ[iz,n] = datai["gamma"][n]
        end
        if γi < minimum(γ[iz,:]) || γi > maximum(γ[iz,:]) error("Energy is out of range.") end
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
    Zeff = sum(ωz.*Z) 
    for iz in range(1,Nz)
        if γi < maximum(γ[iz,:])
            Czprime = spline_Cz[iz](γi)
        elseif γi < 2
            Czprime = sum(pn[iz,:] .* (γi-1).^(1:6))
        else
            Czprime = sum(pn[iz,:] .* 2 .^(1:6))
        end
        Cz += Z[iz]/Zeff * (Czprime + (Z[iz]-S₀[iz])/(2*Z[iz])*(log(β²/(1-β²))-β²))
    end 
    return Cz
end

"""
    barkas_correction(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,γ::Float64,
    particle::Particle)

Gives the Barkas correction to the Bethe formula.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `ωz::Vector{Float64}` : weight fraction of the elements composing the material.
- `ρ::Float64` : material density.
- `γ::Float64` : total energy of the incoming particle in unit of its rest energy.
- `particle::Particle` : incoming particle.

# Output Argument(s)
- `ΔLB::Float64` : Barkas correction.

# Reference(s)
- Salvat (2022), Bethe stopping-power formula and its corrections.
- Salvat and Andreo (2023), SBETHE: Stopping powers of materials for swift charged
  particles from the corrected Bethe formula.

"""
function barkas_correction(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,γ::Float64,particle::Particle)

    #----
    # Initialization
    #----
    α = 1/137.035999177
    β² = (γ^2-1)/γ^2
    Nz = length(Z)
    Zeff = sum(ωz.*Z)
    CB = max(1,Zeff/10)
    charge = get_charge(particle)
    Ωp = plasma_energy(Z,ωz,ρ)
    I = effective_mean_excitation_energy(Z,ωz)

    #----
    # Extract the subshell of all atoms of the compound
    #----
    Zi = zeros(0); Ui = zeros(0); Wi = zeros(0)
    for iz in range(1,Nz)
        _,Ze,Ue,_,_,_ = electron_subshells(Z[iz])
        We = resonance_energy(Z[iz],Ωp,I)
        append!(Zi,ωz[iz].*Ze)
        append!(Ui,Ue)
        append!(Wi,We)
    end

    #----
    # Compute the Barkas correction
    #----
    ΔLB = 0
    for iz in range(1,Nz)
        for δi in range(1,length(Zi))
            ξ = 0.5616 * CB * Wi[δi]/(γ*β²)
            ΔLB += ωz[iz] * charge*α/(γ^2*β²^(3/2))/Z[iz] * Zi[δi] * Wi[δi] * (I₁(ξ) + I₂(ξ)/γ^2)
        end
    end
    return ΔLB
end

"""
    I₁(x::Float64)

Approximation of integral I₁ (Eq. 96a) from Salvat (2022).

# Input Argument(s)
- `x::Float64` : evaluation point.

# Output Argument(s)
- `I₁::Float64` : evaluated integral.

# Reference(s)
- Salvat (2022), Bethe stopping-power formula and its corrections.
- Salvat and Andreo (2023), SBETHE: Stopping powers of materials for swift charged
  particles from the corrected Bethe formula.

"""
function I₁(x::Float64)
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
    I₂(x::Float64)

Approximation of integral I₂ (Eq. 96b) from Salvat (2022).

# Input Argument(s)
- `x::Float64` : evaluation point.

# Output Argument(s)
- `I₂::Float64` : evaluated integral.

# Reference(s)
- Salvat (2022), Bethe stopping-power formula and its corrections.
- Salvat and Andreo (2023), SBETHE: Stopping powers of materials for swift charged
  particles from the corrected Bethe formula.

"""
function I₂(x::Float64)
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