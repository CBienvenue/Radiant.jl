"""
    ratio_positron_electron_bremsstrahlung(Z::Int64,Ei::Float64)

Gives the ratio of the positron and electron Bremsstrahlung cross-sections.

# Input Argument(s)
- `Z::Int64` : atomic number.
- `Ei::Float64` : incoming particle energy.

# Output Argument(s)
- `r::Float64` : ratio of the positron and electron Bremsstrahlung cross-sections.

# Reference(s)
- Kim et al. (1986), Ratio of positron to electron bremsstrahlung energy loss: An
  approximate scaling law.
- Salvat and Andreo (2023), SBETHE: Stopping powers of materials for swift charged
  particles from the corrected Bethe formula

"""
function ratio_positron_electron_bremsstrahlung(Z::Int64,Ei::Float64)
    t = log(1+1e6/Z^2*Ei)
    r = 1 - exp(-1.2359e-1*t + 6.1274e-2*t^2-3.1516e-2*t^3+7.7446e-3*t^4-1.0595e-3*t^5+7.0568e-5*t^6-1.8080e-6*t^7)
    return r
end