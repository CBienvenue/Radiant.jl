"""
    fermi_density_effect(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,Ei::Float64,
    type::String)

Compute the Fermi density effect correction.

# Input Argument(s)
- 'Z::Vector{Int64}': atomic number of the element(s) composing the material.
- 'ωz::Vector{Float64}': weight fraction of the element(s) composing the material.
- 'ρ::Float64': density of the material [in g/cm³].
- 'Ei::Float64': energy of the incoming particle [in mₑc²].
- 'type::String': type of calculation.

# Output Argument(s)
- 'δF::Float64': Fermi density effect correction.

# Author(s)
Charles Bienvenue

# Reference(s)
- Salvat (2019), PENELOPE-2018: A Code System for Monte Carlo Simulation of Electron and
  Photon Transport.
- Lorence (1989) : Physics Guide to CEPXS: A Multigroup Coupled Electron-Photon
  Cross-Section Generating Code.

"""
function fermi_density_effect(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,Ei::Float64,state_of_matter::String,type::String)

if type == "fano"

    # Initialization
    β² = Ei*(Ei+2)/(Ei+1)^2
    Ωp = plasma_energy(Z,ωz,ρ)
    Zeff = sum(ωz.*Z)
    I = effective_mean_excitation_energy(Z,ωz)

    # Extract the subshell of all atoms of the compound
    Zi = zeros(0); Ui = zeros(0); Wi = zeros(0); Nz = length(Z)
    for i in range(1,Nz)
        _,Ze,Ue,_,_,_ = electron_subshells(Z[i])
        We = resonance_energy(Z[i],Ωp,I)
        append!(Zi,ωz[i].*Ze)
        append!(Ui,Ue)
        append!(Wi,We)
    end

    if (1-β²) < 1/Zeff * Ωp^2 * sum(Zi./Wi.^2)

        # Find the root L
        f(L) = Ωp^2/Zeff * sum(Zi./(Wi.^2 .+ L^2)) - (1-β²)
        dfdL(L) = -2*L*Ωp^2/Zeff * sum(Zi./((Wi.^2 .+ L.^2).^2))
        L = newton_bissection(f,dfdL,10*Ωp/sqrt(1-β²),0)

        # Compute the Fermi density effect
        δF = 1/Zeff * sum(Zi.* log.(1 .+ L^2 ./Wi.^2)) - L^2/Ωp^2 * (1-β²)

    else
        # No root for L
        δF = 0
    end

elseif type == "sternheimer"

    mₑc² = 0.510999
    I = effective_mean_excitation_energy(Z,ωz)
    Ωp = plasma_energy(Z,ωz,ρ)

    C = -2*log(I/Ωp) - 1
    X = log10(sqrt(2*Ei+Ei^2))

    if state_of_matter ∈ ["solid","liquid"]
        if I ≥ 0.0001/mₑc² && -C ≥ 5.215
            X0 = -0.326*C-1.5
        elseif I < 0.0001/mₑc² && -C ≥ 3.681
            X0 = -0.326*C-1
        else
            X0 = 0.2
        end
        if I ≥ 0.0001/mₑc²
            X1 = 3
        else
            X1 = 2
        end
    elseif state_of_matter ∈ ["gaz"]
        if 13.804 ≤ -C
            X0 = -0.326*C-2.5
        elseif 11.5 ≤ -C < 13.804
            X0 = 2
        elseif 11 ≤ -C < 11.5
            X0 = 1.9
        elseif 10.5 ≤ -C < 11
            X0 = 1.8
        elseif 10 ≤ -C < 10.5
            X0 = 1.7
        else
            X0 = 1.6
        end
        if 12.25 ≤ -C
            X1 = 5
        else
            X1 = 4
        end
    else
        error("The state of matter is either solid, liquid or gaz.")
    end
    B = (-C - 4.606*X0)/(X1-X0)^3
    if X ≤ X0
        D = 0.0
    elseif X0 < X ≤ X1 
        D = 4.606*X + C + B*(X1-X)^3
    else
        D = 4.606*X + C
    end
    δF = D

else
    error("Unknown type of Fermi density effect calculations")
end

return δF
end