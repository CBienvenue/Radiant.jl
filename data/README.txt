DATA STRUCTURES
    Based on JLD2.jl package, to save and load data structure in format HDF5.

-------------------------------------------------------------------------------------------
    FILE: shell_corrections_salvat_2023.jld2
-------------------------------------------------------------------------------------------
REFERENCE: 
    - Salvat, F., & Andreo, P. (2023). SBETHE: Stopping powers of materials for swift 
      charged particles from the corrected Bethe formula. Computer Physics Communications,
      287, 108697.
-------------------------------------------------------------------------------------------
DATA STRUCTURE:
├── [particle]
│   ├── [Z]
│   │   ├── ["gamma"] -> γ
│   │   ├── ["modified_shell_corrections"] -> C'/Z
│   │   ├── ["cutoff_energy"] -> Ec
│   │   ├── ["fit_parameters"] -> p_n
│   │   ├── ["integrated_gos"] -> S_0
where:
- particle::String = "electron" or "positron" : type of particle. 
- Z::Int64 = 1 or 2 or ... or 99 : atomic number.
- γ::Vector{Float64} (size 153) : total energy of the projectile in units of its rest
  energy.
- C'/Z::Vector{Float64} (size 153) : modified shell corrections.
- Ec::Float64 : cutoff energy between the analytical and calculated shell corrections
  formulas, in unit of the electron rest mass energy.
- p_n::Vector{Float64} (size 6) : fit parameters for the analytical formula.
- S_0::Float64 : integrated longitudinal generalized oscillator strenghts (GOS).
-------------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------------
    FILE: mott_data_boschini_2013.jld2
-------------------------------------------------------------------------------------------
REFERENCE: 
    - Boschini, M. J., Consolandi, C., Gervasi, M., Giani, S., Grandi, D., Ivanchenko, V.,
      ... & Tacconi, M. (2013). An expression for the Mott cross section of electrons and 
      positrons on nuclei with Z up to 118. Radiation Physics and Chemistry, 90, 39-66.
-------------------------------------------------------------------------------------------
DATA STRUCTURE:
├── [particle]
│   ├── [Z] -> bjk
where:
- particle::String = "electrons" or "positrons" : type of particle. 
- Z::Int64 = 1 or 2 or ... or 118 : atomic number.
- bjk::Array{Float64} (size 5 x 6) : Parameter b_{j,k} from Boschini paper, with j between
  0 and 4 and k between 1 and 6.
-------------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------------
    FILE: relaxation_JENDL5.jld2
-------------------------------------------------------------------------------------------
REFERENCE: 
    - JENDL-5 library on https://www-nds.iaea.org/exfor/endf.htm.
-------------------------------------------------------------------------------------------
DATA STRUCTURE:
├── ["Fluorescence"]
│   ├── [Z]
│   |   ├── [subshell]
│   |   |   ├── ["ΔE"] -> ΔE
│   |   |   ├── ["P"] -> P
│   |   |   ├── ["Secondary_shells"] -> SS
├── ["Auger"]
│   ├── [Z]
│   |   ├── [subshell]
│   |   |   ├── ["Tertiary_shells"] -> TS
│   |   |   ├── (same field than with "Fluorescence")
where:
- relaxation_type::String = "Fluorescence" or "Auger" : type de relaxation.
- Z::Int64 = 1 or 2 or ... or 99 : atomic number.
- subshell::String = "K", "L1", "L2", ... : subshell.
- ΔE::Vector{Float64} : vector of the transition energy from the subshell and the others.
- P:::Vector{Float64} : vector of the transition probability from the subshell and the
  others.
- SS::Vector{String} : secondary shell implied in the transition.
- TS::Vector{String} : tertiary shell implied in the transition (Auger only).
-------------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------------
    FILE: rayleigh_factors_JENDL5.jld2
-------------------------------------------------------------------------------------------
REFERENCE: 
    - JENDL-5 library on https://www-nds.iaea.org/exfor/endf.htm.
    - Trkov, A., Herman, M., & Brown, D. A. (2012). ENDF-6 formats manual. Brookhaven
      National Laboratory, 80.
-------------------------------------------------------------------------------------------
DATA STRUCTURE:
├── ["x"]
|   ├── [Z] -> x
├── ["F"]
|   ├── [Z] -> F
├── ["E_real"]
|   ├── [Z] -> E_real
├── ["f_real"]
|   ├── [Z] -> f_real
├── ["E_imag"]
|   ├── [Z] -> E_imag
├── ["f_imag"]
|   ├── [Z] -> f_imag
where:
- Z::Int64 = 1 or 2 or ... or 99 : atomic number.
- x::Vector{Float64} : momentum transfer [x = 2*E/(hc)*sqrt((1-mu)/2)] associated to F.
- F::Vector{Float64} : Atomic form factor.
- E_real::Vector{Float64} : energy coordinates associated with f_real.
- f_real::Vector{Float64} : real part of the anomalous scattering factor.
- E_imag::Vector{Float64} : energy coordinates associated with f_imag.
- f_imag::Vector{Float64} : imaginary part of the anomalous scattering factor.
-------------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------------
    FILE: photoelectric_JENDL5.jld2
-------------------------------------------------------------------------------------------
REFERENCE: 
    - JENDL-5 library on https://www-nds.iaea.org/exfor/endf.htm.
-------------------------------------------------------------------------------------------
DATA STRUCTURE:
├── ["E"]
|   ├── [Z]
|   |   ├── [subshell] -> E
├── ["σ"]
|   ├── [Z]
|   |   ├── [subshell] -> σ
where:
- Z::Int64 = 1 or 2 or ... or 99 : atomic number.
- E::Vector{Float64} : energy associated with σ for the given atom and subshell.
- σ::Vector{Float64} : photoelectric cross-section for the given atom and subshell.
-------------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------------
    FILE: photoelectric_biggs_lighthill_1988.jld2
-------------------------------------------------------------------------------------------
REFERENCE: 
    - Biggs, F., & Lighthill, R. (1988). Analytical approximations for X-ray cross sections
      III (No. SAND-87-0070). Sandia National Lab.(SNL-NM), Albuquerque, NM
      (United States).
-------------------------------------------------------------------------------------------
DATA STRUCTURE:
├── ["E"]
|   ├── [Z] -> E
├── ["M"]
|   ├── [Z] -> M
where:
- Z::Int64 = 1 or 2 or ... or 99 : atomic number.
- E::Vector{Float64} : energy associated with M for the given atom.
- M::Array{Float64} : cross-section interpolation parameters for the given atom.
-------------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------------
    FILE: pair_production_JENDL5.jld2
-------------------------------------------------------------------------------------------
REFERENCE: 
    - JENDL-5 library on https://www-nds.iaea.org/exfor/endf.htm.
-------------------------------------------------------------------------------------------
DATA STRUCTURE:
├── ["E"]
|   ├── [Z] -> E
├── ["σ"]
|   ├── [Z] -> σ
where:
- Z::Int64 = 1 or 2 or ... or 99 : atomic number.
- E::Vector{Float64} : energy associated with σ for the given atom.
- σ::Vector{Float64} : total cross-section for the given atom.
-------------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------------
    FILE: compton_factors_JENDL5.jld2
-------------------------------------------------------------------------------------------
REFERENCE: 
    - JENDL-5 library on https://www-nds.iaea.org/exfor/endf.htm.
    - Trkov, A., Herman, M., & Brown, D. A. (2012). ENDF-6 formats manual. Brookhaven
      National Laboratory, 80.
-------------------------------------------------------------------------------------------
DATA STRUCTURE:
├── ["x"]
|   ├── [Z] -> x
├── ["F"]
|   ├── [Z] -> F
where:
- Z::Int64 = 1 or 2 or ... or 99 : atomic number.
- x::Vector{Float64} : momentum transfer [x = 2*E/(hc)*sqrt((1-mu)/2)] associated to F.
- F::Vector{Float64} : incoherent scattering function.
-------------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------------
    FILE: bremsstrahlung_photons_distribution_poskus_2019.jld2
-------------------------------------------------------------------------------------------
REFERENCE: 
    - Poškus, A. (2019). Shape functions and singly differential cross sections of
      bremsstrahlung at electron energies from 10 eV to 3 MeV for Z= 1–100. Atomic Data and
      Nuclear Data Tables, 129, 101277.
    - Bienvenue, C., Naceur, A., Hébert, A., & Carrier, J. F. (2025). Toward highly
      accurate multigroup coupled photon-electron-positron cross-sections for the Boltzmann
      Fokker-Planck equation. Journal of Computational Physics, 524, 113740.
-------------------------------------------------------------------------------------------
DATA STRUCTURE:
├── ["A"]
|   ├── [Z] -> A
├── ["B"]
|   ├── [Z] -> B
├── ["C"]
|   ├── [Z] -> C
├── ["E"]
|   ├── [Z] -> E
├── ["r"]
|   ├── [Z] -> r
where:
- Z::Int64 = 1 or 2 or ... or 99 : atomic number.
- A::Array{Float64} : Interpolation parameter.
- B::Array{Float64} : Interpolation parameter.
- C::Array{Float64} : Interpolation parameter.
- E::Vector{Float64} : Energies.
- r::Array{Float64} : Indexes.
-------------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------------
    FILE: bremsstrahlung_data_seltzer_berger_1986.jld2
-------------------------------------------------------------------------------------------
REFERENCE: 
    - Seltzer, S. M., & Berger, M. J. (1986). Bremsstrahlung energy spectra from electrons
      with kinetic energy 1 keV–10 GeV incident on screened nuclei and orbital electrons of
      neutral atoms with Z= 1–100. Atomic data and nuclear data tables, 35(3), 345-418.

-------------------------------------------------------------------------------------------
DATA STRUCTURE:
├── ["incident_electron_energy"] -> Ei
├── ["radiative_energy_fraction"] -> r
├── ["scaled_stopping_powers"]
|   ├── [Z] -> scaled_sp
├── ["scaled_cross_sections"]
|   ├── [Z] -> scaled_cs

where:
- Z::Int64 = 1 or 2 or ... or 99 : atomic number.
- Ei::Vector{Float64} : Incoming electron energy.
- r::Vector{Float64} : radiative energy fraction.
- scaled_sp::Vector{Float64} : scaled stopping power.
- scaled_cs::Array{Float64} : matrix of scaled cross-sections.
-------------------------------------------------------------------------------------------

