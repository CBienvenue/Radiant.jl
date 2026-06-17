const A0_M = 5.29177210544e-11

"""
    moliere_screening(Z::Int64,Ei::Float64,is_seltzer_correction::Bool=true)

Gives the Molière screening factor.

# Input Argument(s)
- `Z::Int64` : atomic number.
- `Ei::Float64` : kinetic energy of the incoming particle.
- `is_seltzer_correction::Bool=true` : boolean to enable or not the Seltzer correction.

# Output Argument(s)
- `η::Float64` : Molière screening factor.

# Reference(s)
- Lorence et al. (1989), Physics guide to CEPXS: a multigroup coupled electron-photon
  cross-section generating code.
- Seltzer (1988), An Overview of ETRAN Monte Carlo Methods.

"""
function moliere_screening(Z::Int64,Ei::Float64,is_seltzer_correction::Bool=true)
    α = 1/137.035999177 
    β² = Ei*(Ei+2)/(Ei+1)^2
    if is_seltzer_correction
        g = sqrt(Ei/(Ei+1))
    else
        g = 1
    end
    η = Z^(2/3) * α^2 * (1.13 + 3.76 * (Z*α)^2/β² * g ) / (4 * (9*π^2/128)^(2/3) * Ei * (Ei+2) )
    return η
end

"""
    moliere_chi_a2(E_in_eV::Float64, m_eV::Float64, z::Int, Z::Int)

Gives the squared Moliere screening angle for Coulomb elastic scattering.

# Input Argument(s)
- `E_in_eV::Float64` : projectile kinetic energy in eV.
- `m_eV::Float64` : projectile rest mass energy in eV.
- `z::Int` : projectile charge number.
- `Z::Int` : target atomic number.

# Output Argument(s)
- `χa²::Float64` : squared Moliere screening angle in radian squared.

# Reference(s)
- Bethe (1953), Moliere theory of multiple scattering.
"""
function moliere_chi_a2(E_in_eV::Float64, m_eV::Float64, z::Int, Z::Int)::Float64
    gamma = 1.0 + E_in_eV / m_eV
    beta = sqrt(1.0 - 1.0 / (gamma^2))
    alpha = (z * Z * ALPHA_FS) / beta
    a = (9.0 * π^2 / (128.0 * Z))^(1.0 / 3.0) * A0_M
    pc_eV = sqrt((E_in_eV + m_eV)^2 - (m_eV^2))
    chi_0 = (HBAR_eVs * C_m_per_s) / (pc_eV * a)
    chi_a2 = chi_0^2 * (1.13 + 3.76 * alpha^2)
    return chi_a2
end