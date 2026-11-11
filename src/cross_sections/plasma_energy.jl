"""
    plasma_energy(𝒩ₑ_eff::Float64)

Calculate the plasma energy of a free-electron gas with the electron density of the medium.

# Input Argument(s)
- `𝒩ₑ_eff::Float64`: effective electron density in the medium [in cm⁻³].

# Output Argument(s)
- `Ωp::Float64`: plasma energy [in mₑc²].

# Reference(s)
- Salvat (2019), PENELOPE-2018: A Code System for Monte Carlo Simulation of Electron and
  Photon Transport.

"""
function plasma_energy(𝒩ₑ_eff::Float64)
    ħc = 3.86189603E-11                          # (mₑc²⋅cm)
    rₑ = 2.81794092E-13                          # (cm)
    Ωp = sqrt(4*π*𝒩ₑ_eff*ħc^2*rₑ)
    return Ωp
end
