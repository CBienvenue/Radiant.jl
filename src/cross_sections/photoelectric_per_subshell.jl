"""
    photoelectric_per_subshell(Z::Int64,Ei::Float64,subshell::String)

Gives the photoelectric absorption cross-section based on the Biggs and Lighthill
interpolation.

# Input Argument(s)
- `Zi::Real` : number of electron(s).
- `Ei::Float64` : incoming particle energy.
- `δi::Int64` : subshell index.

# Output Argument(s)
- `σa::Float64` : photoelectric differential cross-sections.

# Reference(s)
- Cullen et al. (1997), EPDL97, the Evaluated Photon Data Library, '97 Version. 

"""
function photoelectric_per_subshell(Z::Int64,Ei::Float64,δi::Int64)

    #----
    # Extract data
    #----
    data = fast_load("photoelectric_EPDL97.jld2")
    E = data["E"][Z][δi]
    σ = data["σ"][Z][δi]
    photoelectric_spline = fast_spline(E,σ,[Z,δi],"photoelectric_per_subshell")

    #----
    # Compute absorption cross-section
    #----
    σa = 0.0
    if (Ei > E[1]) σa = photoelectric_spline(Ei) end
    return σa
end