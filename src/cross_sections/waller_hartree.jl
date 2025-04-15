"""
    waller_hartree(Z::Int64,Ei::Float64,Ef::Float64)

Gives the Waller-Hartree differential cross-section.

# Input Argument(s)
- `Z::Int64` : number of electron(s).
- `Ei::Float64` : incoming photon energy.
- `Ef::Float64` : outgoing photon energy.

# Output Argument(s)
- `σs::Float64` : Waller-Hartree differential cross-section.

# Reference(s)
- Brusa et al. (1996), Fast sampling algorithm for the simulation of photon Compton
  scattering.
- MacFarlane et al. (2019), The NJOY Nuclear Data Processing System, Version 2016.

"""
function waller_hartree(Z::Int64,Ei::Float64,Ef::Float64)

    # Extract data from cache
    data = fast_load("compton_factors_EPDL97.jld2")
    x = data["x"][Z]
    S = data["F"][Z]

    # Compute spline to interpolate Waller-Hartree factor
    S_spline = fast_spline(x,S,Z,"waller_hartree_factor")

    # Compute Waller-Hartree incoherent scattering factor
    hc = 1/20.60744 # (hc in mₑc² × Å)
    μ = max(min(1 + 1/Ei - 1/Ef,1),-1)
    xi = 2*Ei/(hc)*sqrt((1-μ)/2) * sqrt((1+(Ei^2+2*Ei)*(1-μ)/2))/(1+Ei*(1-μ))
    S = S_spline(xi)

    # Compute the scattering cross-section
    σs = S * klein_nishina(Ei,Ef)
    return σs
end