# Reference(s)
# - Varelas, C. and Biersack, J. (1970), Reflection of energetic particles from
#   atomic or ionic chains in single crystals, Nuclear Instruments and Methods
#   79(2), 213-218. doi:10.1016/0029-554X(70)90141-2
#   Reciprocal interpolation formula: Sec. 4.2, Eq. (6), p. 216.
# - Andersen, H. H. and Ziegler, J. F. (1977), Hydrogen: Stopping Powers and
#   Ranges in All Elements, The Stopping and Ranges of Ions in Matter vol. 3,
#   Pergamon Press, New York.
# - Ziegler, J. F. (1977), Helium: Stopping Powers and Ranges in All Elemental
#   Matter, The Stopping and Ranges of Ions in Matter vol. 4, Pergamon Press,
#   New York.
# - ICRU Report 49 (1993), Stopping Powers and Ranges for Protons and Alpha
#   Particles, ICRU, Bethesda, MD.
# - Salvat, F. (2022), Bethe stopping-power formula and its corrections,
#   Physical Review A 106, 032809. doi:10.1103/PhysRevA.106.032809
#   Free-exponent low-energy parametrization: Eqs. (134)-(135); fitting window
#   and objective: Eq. (136) and the text above it. Source of the fitting
#   prescription, not of the functional form -- see icru49_stopping.

"""
    low_energy_stopping(Z::Vector{Int64},ωz::Vector{Float64},atz::Vector{Float64},
    ρ::Float64,N_density::Float64,Ei::Float64,particle::Particle,
    density_correction::String,state_of_matter::String,I_eff::Float64,Ecut::Float64,
    Emax::Float64)

Gives the collisional stopping power below `Ecut`, where the corrected Bethe formula is no
longer valid.

Everything Radiant computes below `Ecut` lives in this file. `bethe`'s low-energy branch is
a single call into it and holds no low-energy physics of its own. There is one dispatch,
on `low_energy_model`:

  - `"analytic"` : the analytical extension anchored on the Bethe branch, below.
  - `"icru49"` : the fitted ICRU 49 / Andersen-Ziegler form, `icru49_stopping`.

The analytical extension solves both coefficients of

    S = exp(A - B*(log t)^1.5)        t = Ei/Emax ≥ 1
    S = 1.5*exp(A)*sqrt(t)*(1 - t/3)  t < 1

from the value and the slope of the Bethe branch at `Ecut`, so the join is smooth. Both
branches have zero slope at `t = 1`, which puts the stopping-power maximum at `Emax`. It
extrapolates from that single anchor, so its error grows with distance from `Ecut`.

The fitted form covers Z = 1-92, protons and alphas only. For any other particle, or any
element outside that range, this falls back to the analytical extension silently: a muon
transport run must keep working when a user selects `"icru49"` for their protons.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `ωz::Vector{Float64}` : weight fraction of the elements composing the material.
- `atz::Vector{Float64}` : atomic fraction of the elements composing the material.
- `ρ::Float64` : material density.
- `N_density::Float64` : atom density of the material [in cm⁻³].
- `Ei::Float64` : incoming particle energy [in mₑc²].
- `particle::Particle` : incoming particle.
- `density_correction::String` : density-effect correction type.
- `state_of_matter::String` : material state.
- `I_eff::Float64` : user-defined mean excitation energy override [in mₑc²]; NaN ⟹ table.
- `Ecut::Float64` : energy below which this branch is used [in mₑc²].
- `Emax::Float64` : energy of the stopping-power maximum of the extension [in mₑc²].
- `atomic_weights::Union{Nothing,Vector{Float64}}` : element atomic weights [u].
- `is_bethe_correction::Bool` : boolean to enable or not the Bethe correction.
- `is_density_correction::Bool` : boolean to enable or not the density correction.
- `is_shell_correction::Bool` : boolean to enable or not the shell correction.
- `is_lindhard_sorensen::Bool` : boolean to enable or not the Lindhard-Sorensen correction.
- `is_barkas::Bool` : boolean to enable or not the Barkas correction.
- `user_ecut::Float64` : the cutoff `bethe` was given [in MeV], forwarded unchanged to the
  anchor calls so they resolve the same `Ecut` and take the corrected-Bethe branch.
- `low_energy_model::String` : low-energy model, `"analytic"` or `"icru49"`.
- `low_energy_source::String` : coefficient source of the fitted model.

# Output Argument(s)
- `S::Float64` : stopping power.

# Reference(s)
- Salvat (2022), Bethe stopping-power formula and its corrections.
- Salvat and Andreo (2023), SBETHE: Stopping powers of materials for swift charged
  particles from the corrected Bethe formula.

"""
function low_energy_stopping(Z::Vector{Int64},ωz::Vector{Float64},atz::Vector{Float64},ρ::Float64,N_density::Float64,Ei::Float64,particle::Particle,density_correction::String,state_of_matter::String,I_eff::Float64,Ecut::Float64,Emax::Float64; atomic_weights::Union{Nothing,Vector{Float64}}=nothing,is_bethe_correction::Bool=true,is_density_correction::Bool=true,is_shell_correction::Bool=true,is_lindhard_sorensen::Bool=true,is_barkas::Bool=true,user_ecut::Float64=NaN,low_energy_model::String="analytic",low_energy_source::String="srim")

    mₑ = 0.51099895069

    #----
    # Fitted ICRU 49 / Andersen-Ziegler branch
    #----
    # The unit change sits on this side of the dispatch: the fitted form works in MeV and
    # returns a mass stopping power, while the analytic extension already works in bethe's
    # units and returns a linear one.
    if low_energy_model != "analytic" && has_low_energy_data(Z,particle)
        return icru49_stopping(Z,ωz,Ei*mₑ,particle,state_of_matter,low_energy_source) * ρ/mₑ
    end

    #----
    # Estimate parameters A and B
    #----
    # The anchor energies are both ≥ Ecut, so these two calls take the corrected-Bethe
    # branch and never come back here: the recursion is one level deep by construction.
    # They pass low_energy_model = "analytic" whatever the caller selected, because under
    # the handover window the fitted form is what bethe returns at Ecut -- the blend weight
    # is zero there -- so an anchor taken with the caller's setting would anchor this
    # extension on the fitted form rather than on the Bethe branch. The correction switches
    # are forwarded unchanged, or the extension would be anchored on a different formula
    # from the one it is joined to.
    h = 0.00005
    Sc = bethe(Z,ωz,atz,ρ,N_density,Ecut,particle,density_correction,state_of_matter,I_eff;
        atomic_weights=atomic_weights,
        is_bethe_correction=is_bethe_correction,
        is_density_correction=is_density_correction,
        is_shell_correction=is_shell_correction,
        is_lindhard_sorensen=is_lindhard_sorensen,
        is_barkas=is_barkas,
        user_ecut=user_ecut,
        low_energy_model="analytic")
    Sc_plus_h = bethe(Z,ωz,atz,ρ,N_density,Ecut+h,particle,density_correction,state_of_matter,I_eff;
        atomic_weights=atomic_weights,
        is_bethe_correction=is_bethe_correction,
        is_density_correction=is_density_correction,
        is_shell_correction=is_shell_correction,
        is_lindhard_sorensen=is_lindhard_sorensen,
        is_barkas=is_barkas,
        user_ecut=user_ecut,
        low_energy_model="analytic")
    dSdE = (Sc_plus_h - Sc)/h
    A = (3*log(Sc)*Sc - 2*dSdE*Ecut*log(Ecut/Emax))/(3*Sc)
    B = -2*Ecut*dSdE/(3*Sc*sqrt(log(Ecut/Emax)))

    #----
    # Compute total cross-sections using analytical formula
    #----
    t = Ei/Emax
    if t ≥ 1
        S = exp(A-B*(log(t))^1.5)
    else
        S = 1.5*exp(A)*sqrt(t)*(1-t/3)
    end
    return S
end

"""
    emax_stopping_maximum(Z::Vector{Int64},atz::Vector{Float64},particle::Particle)

Least-squares power-law fit to the measured stopping-power maximum of all 92 elements,
against SRIM-2013. Returns the energy [in MeV].

Both branches of the analytic low-energy extension have zero slope at `t = Ei/Emax = 1`,
so `Emax` is where that extension puts the stopping-power maximum. A single constant
asserts the same peak energy in hydrogen as in uranium; measured against SRIM-2013 the
real peak runs 41-164 keV for protons and 0.42-1.19 MeV for alphas, rising with the
target's atomic number, and the proton constant is low by about a factor of two on
average.

`Zbar = sum(atz .* Z)` is the mean atomic number per atom, the same definition
`barkas_correction` uses for `CB`. The fit was derived on pure elements, where the atom-
and mass-weighted means coincide, so the compound case is unconstrained either way and
consistency decides it.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `atz::Vector{Float64}` : atomic fraction of the elements composing the material.
- `particle::Particle` : incoming particle.

# Output Argument(s)
- `Emax::Float64` : energy of the stopping-power maximum [in MeV].

"""
function emax_stopping_maximum(Z::Vector{Int64},atz::Vector{Float64},particle::Particle)
    Zbar = sum(atz .* Z)
    if is_proton(particle) || is_antiproton(particle)
        return 0.0481 * Zbar^0.190
    elseif is_alpha(particle)
        return 0.4829 * Zbar^0.114
    else
        error("No stopping-maximum fit for the given particle.")
    end
end

"""
    low_energy_tables(particle::Particle,source::String)

Gives the pair of coefficient tables an ICRU 49 lookup resolves against, and the list of
elements that are gases at STP.

The base table is the SRIM-fitted one and covers Z = 1-92 for both particles. The override
is consulted first, per element and per phase slot, so a per-element alternative source can
replace part of the base table without any other rule changing. With `"srim"` there is no
override and every element comes from the base table.

`"iaea"` supplies coefficients for the same form, in the same units, with the same `a2`
bound, fitted to the IAEA electronic stopping power database instead. Coverage is partial —
41 elements for protons and 30 for alphas — so it is an override on top of full SRIM
coverage rather than a table in its own right. Every element it omits is served by the SRIM
entry through the same lookup, so `"iaea"` is well defined for all of Z = 1-92 and both
particles, never fails and never warns. A material made entirely of elements outside it
produces results bit-identical to `"srim"`.

The two tables carry different phase conventions, each stated in its own provenance string
in the data file, so the slot is resolved once per table by `low_energy_slot`.

# Input Argument(s)
- `particle::Particle` : incoming particle.
- `source::String` : coefficient source, `"srim"` or `"iaea"`.

# Output Argument(s)
- `base::Dict` : the SRIM-fitted per-particle table, covering all 92 elements.
- `override::Dict` : the per-element override, empty for `"srim"`.
- `gas_elements::Vector{Int64}` : elements that are gases at STP.

"""
function low_energy_tables(particle::Particle,source::String)
    data = fast_load("low_energy_ion.jld2")
    key  = is_alpha(particle) ? "alpha" : "proton"
    base = data["fitted_exponent"][key]
    if source == "srim"
        override = Dict{String,Any}()
    elseif source == "iaea"
        override = data["fitted_exponent_iaea"][key]
    else
        error("Unknown '$source' low-energy coefficient source.")
    end
    # gas_elements is a fact about the elements rather than a borrowing from either table,
    # and is the same eleven either way.
    return base, override, data["fitted_exponent"]["gas_elements"]
end

"""
    low_energy_slot(table::Dict,Z::Int64,is_gas::Bool,gas_elements::Vector{Int64})

Gives the phase slot of `table` an element should be read from.

The phase is resolved *before* the source, and the coverage question is then asked of the
slot rather than of the element. The tempting shortcut — is this element in the override
anywhere, and if so use it — is a different question and gives the wrong answer for every
gas-flagged element in a condensed material: the IAEA `"condensed"` slot is empty, because
no condensed-phase measurement set for a gas-flagged element could be determined, so that
shortcut would hand a condensed target oxygen's *gas* coefficients, about 8% high, and
never reach the perfectly good SRIM condensed entry sitting beside it. Nothing in the file
fails loudly if this is got wrong.

The two tables do not share a convention. The SRIM table carries a `"gas"` slot covering
all 92 elements and uses it for any target declared gaseous. The IAEA table has no `"gas"`
slot at all: its `"coefficients"` are the measured fits, which for a gas-flagged element
are the gas-phase ones, and its `"condensed"` slot is what a condensed target of such an
element would need. So each table is asked its own question, and a single shared expression
would silently drop the SRIM gas slot.

# Input Argument(s)
- `table::Dict` : coefficient table to resolve the slot in.
- `Z::Int64` : atomic number of the element.
- `is_gas::Bool` : whether the material is gaseous.
- `gas_elements::Vector{Int64}` : elements that are gases at STP.

# Output Argument(s)
- `slot::String` : phase slot, `"gas"`, `"condensed"` or `"coefficients"`.

"""
function low_energy_slot(table,Z::Int64,is_gas::Bool,gas_elements::Vector{Int64})
    if is_gas
        # The haskey on "gas" is what makes this reader work against a table written before
        # that slot existed, and the slot is not restricted to the elements that are gases
        # at STP.
        return (haskey(table,"gas") && haskey(table["gas"],Z)) ? "gas" : "coefficients"
    else
        return (Z in gas_elements && haskey(table,"condensed")) ? "condensed" : "coefficients"
    end
end

"""
    has_low_energy_data(Z::Vector{Int64},particle::Particle)

Tells whether the fitted ICRU 49 coefficient table covers every element of the material for
the given particle. The table covers Z = 1-92, protons and alphas only.

This predicate guards both the branch below `Ecut` and the handover blend above it. The two
must agree about a material: if they did not, a muon would get the analytic extension below
`Ecut` and a fitted-form blend above it, which is worse than either choice alone.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `particle::Particle` : incoming particle.

# Output Argument(s)
- `is_covered::Bool` : whether the fitted table covers the material.

"""
function has_low_energy_data(Z::Vector{Int64},particle::Particle)
    if !(is_proton(particle) || is_alpha(particle)) return false end
    return all(1 .≤ Z .≤ 92)
end

"""
    icru49_stopping(Z::Vector{Int64},ωz::Vector{Float64},E_MeV::Float64,particle::Particle,
    state_of_matter::String,source::String)

Gives the mass stopping power [in MeV cm²/g] from the fitted ICRU 49 / Andersen-Ziegler
form, at an energy `E_MeV` [in MeV].

    S = s_l*s_h/(s_l + s_h),   s_l = a1*T^a2,   s_h = (a3/T)*log(1 + a4/T + a5*T)

Reciprocal interpolation form: Varelas & Biersack (1970), Sec. 4.2, Eq. (6). There it
interpolates the *per-collision* electronic energy loss dEe at fixed impact parameter b,
between the Firsov low-velocity limit and the Lindhard-modified Bethe-Bloch high-velocity
limit, which leave a gap near the stopping maximum. The high-energy branch already has the
three-term logarithm B*ln(eps + 1 + 5/eps)/eps with eps proportional to E; Andersen &
Ziegler (1977) recast this as a macroscopic stopping-power fit and promoted the two
constants inside the logarithm to free parameters A4, A5. The E^0.45 low branch is an
Andersen-Ziegler empirical adjustment, not from the 1970 paper. The exponent of that low
branch is fitted per element here rather than held at 0.45; the form with it fitted, and
the window and objective used, are Salvat (2022), Eqs. (134)-(136), where it is introduced
as ICRU Report 49's own.

`T` is the kinetic energy in the unit the coefficients were fitted on, keV for protons and
MeV for alphas. The factor is read from the table rather than implied by the particle: a
missing factor of 1000 gives a smooth, plausible, wrong curve rather than an error.

The fitted exponent is what removes the need for a second branch below 10 keV, so there is
no `E_LSS` test in front of the form: `s_l -> 0` as `T^a2` and it runs to zero energy on
its own. The result is smooth through 10 keV and through 40 keV, where the older
fixed-exponent tables have their handover to a `sqrt(E)` branch.

Compounds are always built from the elements with mass fractions, `S = Σⱼ ωⱼ Sⱼ`. This
branch is unconditionally Bragg-additive: there is no option and no compound-level fit,
which is why `set_is_bragg_additivity` is scoped to the Bethe branch.

Phase matching is automatic and has no switch. It substitutes fitted coefficients for the
11 elements that are gases at STP; it is not a percentage applied to the stopping power.
SRIM's electronic stopping is about 8% lower in the condensed phase than the gas phase for
the same element, which is what motivated the swap, but nothing multiplies by it. The gas
slot covers all 92 elements, not only those 11: SRIM's gas flag selects a stopping
parametrisation and applies to whatever target is declared gaseous, so a gaseous compound
takes every constituent in its own phase. Mixing acetylene from condensed carbon instead of
gas carbon costs 5.9% at 0.75 MeV.

Below about 5 keV this branch is badly wrong for some elements, and no choice of exponent
fixes it: helium is worst, at 153% median below 10 keV. Treat 10 keV, the bottom of the
fitting window, as the floor of validity.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `ωz::Vector{Float64}` : weight fraction of the elements composing the material.
- `E_MeV::Float64` : incoming particle energy [in MeV].
- `particle::Particle` : incoming particle.
- `state_of_matter::String` : material state.
- `source::String` : coefficient source.

# Output Argument(s)
- `S::Float64` : mass stopping power [in MeV cm²/g].

# Reference(s)
- Varelas and Biersack (1970), Reflection of energetic particles from atomic or ionic
  chains in single crystals. Reciprocal combination rule only, Sec. 4.2, Eq. (6).
- Andersen and Ziegler (1977), Hydrogen: Stopping Powers and Ranges in All Elements.
- Ziegler (1977), Helium: Stopping Powers and Ranges in All Elemental Matter.
- ICRU Report 49 (1993), Stopping Powers and Ranges for Protons and Alpha Particles.
- Salvat (2022), Bethe stopping-power formula and its corrections. Eqs. (134)-(136), for
  the free-exponent parametrization and the fitting prescription.

"""
function icru49_stopping(Z::Vector{Int64},ωz::Vector{Float64},E_MeV::Float64,particle::Particle,state_of_matter::String,source::String)

    base, override, gas_elements = low_energy_tables(particle,source)
    T = E_MeV/base["T_unit_MeV"]
    is_gas = state_of_matter == "gaz"

    S = 0.0
    for i in eachindex(Z)

        # Resolve the phase first, then the source, and ask the coverage question of the
        # slot rather than of the element -- see low_energy_slot. Each table is asked with
        # its own phase convention. Selecting a source does not switch the whole material
        # over: it is resolved per element, at the point of use, and a compound routinely
        # ends up drawing from both tables at once.
        base_slot = low_energy_slot(base,Z[i],is_gas,gas_elements)
        over_slot = low_energy_slot(override,Z[i],is_gas,gas_elements)
        coef = haskey(override,over_slot) ? get(override[over_slot],Z[i],base[base_slot][Z[i]]) :
                                            base[base_slot][Z[i]]

        a1,a2,a3,a4,a5 = coef
        s_l = a1 * T^a2
        s_h = (a3/T) * log(1 + a4/T + a5*T)
        S += ωz[i] * s_l*s_h/(s_l + s_h)
    end
    return S
end

"""
    blend_weight(E_MeV::Float64,Ecut_MeV::Float64,E_top_MeV::Float64)

Weight on the Bethe branch: 0 at `Ecut_MeV`, 1 at `E_top_MeV`, with zero derivative at both
ends (3t² - 2t³). That both derivatives vanish is what makes the blended result C¹ — at
`Ecut` the slope is the fitted form's own, at `E_top_MeV` it is the Bethe branch's.

The window sits *above* `Ecut`, not below it, so the fitted form is the answer at `Ecut`
and for an octave above it and the Bethe branch is never evaluated below its own cutoff.
That is the opposite of the obvious arrangement and it is deliberate: blending below `Ecut`
needs the Bethe branch at energies where it refuses to evaluate, which means recursing with
an artificially low cutoff and swallowing the exception when the shell-correction table
runs out. Blending above needs none of it, and measured against NIST over the 25 elements
it tabulates the window above scores better on every column.

# Input Argument(s)
- `E_MeV::Float64` : incoming particle energy [in MeV].
- `Ecut_MeV::Float64` : bottom of the handover window [in MeV].
- `E_top_MeV::Float64` : top of the handover window [in MeV].

# Output Argument(s)
- `w::Float64` : weight on the Bethe branch.

"""
function blend_weight(E_MeV::Float64,Ecut_MeV::Float64,E_top_MeV::Float64)
    t = (E_MeV - Ecut_MeV)/(E_top_MeV - Ecut_MeV)
    t ≤ 0 && return 0.0
    t ≥ 1 && return 1.0
    return t^2*(3 - 2*t)
end
