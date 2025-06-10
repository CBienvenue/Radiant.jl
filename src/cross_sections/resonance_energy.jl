"""
    resonance_energy(Z::Int64)

Compute the resonance energy of a set of subshells.

# Input Argument(s)
- `Z::Int64`: atomic number of the element.
- `Ωp::Float64`: plasma energy of the medium [in mₑc²].

# Output Argument(s)
- `Wi::Vector{Float64}`: resonance energy per subshells [in mₑc²].

# Reference(s)
- Salvat (2019), PENELOPE-2018: A Code System for Monte Carlo Simulation of Electron and
  Photon Transport.

"""
function resonance_energy(Z::Int64,Ωp::Float64,I::Float64)

Nshells,Zi,Ui,_,_,_ = electron_subshells(Z)

# Find parameter a
f(a) = sum(Zi.*log.(sqrt.((a.*Ui).^2 + 2/3*Zi./Z.*Ωp^2))) - Z*log(I)
dfda(a) = sum(Zi.*Ui.^2 .*a./(Ui.^2 .*a.^2 .+ 2/3 .* Zi./Z.*Ωp^2))
a = newton_bisection(f,dfda,0,100,1e-4,2000)

# Compute the resonance energies
Wi = zeros(Nshells)
for i in range(1,Nshells)
    Wi[i] = sqrt((a*Ui[i])^2 + 2/3*Zi[i]/Z*Ωp^2)
end

return Wi
end