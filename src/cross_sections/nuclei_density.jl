Nₐ = 6.022094e23 # (Avogadro number in u/g)

"""
    nuclei_density(M::Float64,ρ::Float64)

Compute the nuclei density for a material with molar mass M and density ρ.

# Input Argument(s)
- `M::Float64`: molar mass of the material [in g/mol].
- `ρ::Float64`: density of the material [in g/cm³].

# Output Argument(s)
- `𝒩::Float64`: nuclei density [in cm⁻³].
"""
function nuclei_density(M::Float64,ρ::Float64)

𝒩 = ρ*Nₐ/M

return 𝒩
end

"""
    nuclei_density(Z::Int64,ρ::Float64)

Compute the nuclei density for atomic number Z ∈ {1,100} in a material of density ρ.

# Input Argument(s)
- `Z::Int64`: atomic number of the element.
- `ρ::Float64`: density of the material [in g/cm³].

# Output Argument(s)
- `𝒩::Float64`: nuclei density [in cm⁻³].

# Reference(s)
- Hébert (2016), Applied Reactor Physics.

"""
function nuclei_density(Z::Int64,ρ::Float64)

𝒩 = nuclei_density(standard_atomic_weight(Z),ρ)

return 𝒩
end