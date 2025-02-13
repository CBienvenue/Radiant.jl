"""
    plasma_energy(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64)

Calculate the plasma energy of a free-electron gas with the electron density of the medium.

# Input Argument(s)
- `Z::Vector{Int64}`: atomic number of the element(s) composing the material.
- `ωz::Vector{Float64}`: weight fraction of the element(s) composing the material.
- `ρ::Float64`: density of the material [in g/cm³].

# Output Argument(s)
- `Ωp::Float64`: plasma energy [in mₑc²].

# Reference(s)
- Salvat (2019), PENELOPE-2018: A Code System for Monte Carlo Simulation of Electron and
  Photon Transport.

"""
function plasma_energy(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64)

Nₐ = 6.022094e23
𝒩ₑ_eff = sum(ωz.*Z./atomic_weight.(Z))*Nₐ*ρ # (cm⁻³)
ħc = 3.86189603E-11                          # (mₑc²⋅cm)
rₑ = 2.81794092E-13                          # (cm)
Ωp = sqrt(4*π*𝒩ₑ_eff*ħc^2*rₑ)               # (mₑc²)

return Ωp
end