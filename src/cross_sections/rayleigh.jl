"""
    rayleigh(Z::Int64,Ei::Float64,μ::Float64)

Gives the Rayleigh differential cross-section.

# Input Argument(s)
- `Z::Int64` : atomic number.
- `Ei::Float64` : incoming particle energy.
- `μ::Float64` : direction cosine (μ = cos(θ)).

# Output Argument(s)
- `σs::Float64` : Rayleigh differential cross-section.

# Reference(s)
- Cullen et al. (1997), EPDL97, the Evaluated Photon Data Library, '97 Version. 

"""
function rayleigh(Z::Int64,Ei::Float64,μ::Float64)

    #----
    # Extract data from cache
    #----
    data = fast_load("rayleigh_factors_EPDL97.jld2")
    x= data["x"][Z]
    F= data["F"][Z]
    E_real= data["E_real"][Z]
    f_real= data["f_real"][Z]
    E_imag= data["E_imag"][Z]
    f_imag= data["f_imag"][Z]
    rayleigh_scattering_factor = fast_spline(x,F,Z,"rayleigh_scattering_factor")
    rayleigh_anomalous_factor_real = fast_spline(E_real,f_real,Z,"rayleigh_anomalous_factor_real")
    rayleigh_anomalous_factor_imag = fast_spline(E_imag,f_imag,Z,"rayleigh_anomalous_factor_imag")

    #----
    # Compute the Rayleigh cross-section
    #----
    rₑ = 2.81794092E-13 # (in cm)
    hc = 1/20.60744 # (hc in mₑc² × Å)
    xi = 2*Ei/(hc)*sqrt((1-μ)/2)
    Fi = rayleigh_scattering_factor(xi)
    if Ei < E_real[end] fi_real = rayleigh_anomalous_factor_real(Ei) else fi_real = 0 end
    if Ei < E_imag[end] fi_imag = rayleigh_anomalous_factor_imag(Ei) else fi_imag = 0 end
    σs =  π*rₑ^2 * (1+μ^2) * ((Fi + fi_real)^2 + fi_imag^2)

    return σs
end