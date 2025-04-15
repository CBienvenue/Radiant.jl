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