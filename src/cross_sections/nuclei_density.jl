"""
    nuclei_density(Z::Int64,ρ::Float64)

Compute the nuclei density for atomic number Z ∈ {1,100} in a material of density ρ.

# Input Argument(s)
- 'Z::Int64': atomic number of the element.
- 'ρ::Float64': density of the material [in g/cm³].

# Output Argument(s)
- '𝒩::Float64': nuclei density [in cm⁻³].

# Reference(s)
- Hébert (2016) : Applied Reactor Physics.

"""
function nuclei_density(Z::Int64,ρ::Float64)

Nₐ = 6.022094e23 # (Avogadro number in u/g)
𝒩 = ρ*Nₐ/atomic_weight(Z)

return 𝒩
end