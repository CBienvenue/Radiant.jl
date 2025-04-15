"""
    high_energy_Coulomb_correction(Z::Int64,is_approx::Bool=true)

Gives the high-energy Coulomb correction of Davies et al.

# Input Argument(s)
- `Z::Int64` : atomic numbers.
- `is_approx` : boolean to enable or not the approximation of this correction.

# Output Argument(s)
- `S::Float64` : stopping power.

# Reference(s)
- Davies et al. (1954), Theory of Bremsstrahlung and Pair Production. II. Integral Cross
  Section for Pair Production.
- Baró et al. (1994), Analytical cross-sections for Monte-Carlo simulation of photon
  transport.

"""
function high_energy_Coulomb_correction(Z::Int64,is_approx::Bool=true)

    # Initialization
    α = 1/137.035999177
    a = α*Z

    # Approximate form
    if is_approx
        fc = a^2*(1/(1+a^2) + 0.202059 - 0.03693*a^2 + 0.00835*a^4 - 0.00201*a^6 + 0.00049*a^8 - 0.00012*a^10 + 0.00003*a^12)

    # Analytical form (as an infinite summation)
    else
        fc = 0; fc⁻ = 0
        ϵ = Inf
        n = 1
        while ϵ > 1e-5
            fc += 1/(n*(n^2+a^2))
            ϵ = abs((fc - fc⁻)/max(fc⁻,1e-5))
            fc⁻ = copy(fc)
            n += 1
        end
        fc *= a^2
    end
    return fc
end