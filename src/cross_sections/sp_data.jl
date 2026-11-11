# Data-driven stopping powers: SRIM-2013 tabulated values, and user-supplied tables.
#
# Both replace the computed electronic stopping power with numbers read from a file rather
# than evaluated, and both interpolate log-log with the same monotone cubic. They differ in
# what may be trusted: the SRIM tables ship with Radiant, are keyed on (Z, ion, phase) and
# cannot be mistyped, while a user table arrives at runtime, is keyed on a path, and must be
# validated line by line before it is believed. Keep that split visible when editing --
# nothing in the SRIM section needs a parse error, and nothing in the user section may
# assume a well-formed file.
#
# Reference(s)
# - Ziegler, J. F., Ziegler, M. D. and Biersack, J. P. (2010), SRIM - The
#   Stopping and Range of Ions in Matter, Nuclear Instruments and Methods B
#   268(11-12), 1818-1823. doi:10.1016/j.nimb.2010.02.091
#   Tables here were produced with SRIM-2013 (SR Module reporting
#   SRIM-2012.01), output units 6 = MeV/(g/cm2).
# - Ziegler, J. F. and Manoyan, J. M. (1988), The stopping of ions in compounds,
#   Nuclear Instruments and Methods B 35(3-4), 215-228.
#   doi:10.1016/0168-583X(88)90273-X
#   The "cores and bonds" model behind SRIM's compound correction. SRIM stores
#   only the per-molecule bond inventory (Data/Compound.dat); the bond stopping
#   coefficients are not in any shipped data file.
# - Fritsch, F. N. and Carlson, R. E. (1980), Monotone piecewise cubic
#   interpolation, SIAM Journal on Numerical Analysis 17(2), 238-246.
#   doi:10.1137/0717021
#   Already implemented as cubic_hermite_spline in src/tools/spline.jl; applied
#   here to log(E), log(S), on both paths.
#
# The stopping-power values in a user-supplied table have no reference here, by
# construction: they are supplied by the user at runtime. Provenance for a given
# table belongs in that table's own comment lines, not in this file. That absence
# is deliberate and is not an incomplete reference block.

#----
# SRIM-2013 tabulated stopping powers
#----
"""
    srim_phase(state_of_matter::String)

Maps Radiant's state of matter onto SRIM's binary phase flag. `"gaz"` selects `"gas"`;
`"solid"` and `"liquid"` both select `"solid"`, because SRIM has no liquid.

The phase follows the material the user declared and is never decided from the element.
All 92 elements have a different gas table, not only the ones that are gases at room
temperature: 87 of 92 differ by more than 5% for protons, and the largest differences sit
on lead, titanium, gadolinium and barium at about 8.8%. SRIM applies a near-uniform
gas-phase enhancement to whatever is declared a gas, so a solid containing oxygen must not
pick up oxygen's gas table.

# Input Argument(s)
- `state_of_matter::String` : material state.

# Output Argument(s)
- `phase::String` : `"solid"` or `"gas"`.

"""
function srim_phase(state_of_matter::String)
    if state_of_matter == "gaz"
        return "gas"
    elseif state_of_matter ∈ ["solid","liquid"]
        return "solid"
    else
        error("The state of matter is either solid, liquid or gaz.")
    end
end

"""
    srim_stopping_power(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,Ei::Float64,
    particle::Particle,state_of_matter::String)

Gives the electronic stopping power from SRIM-2013's own tabulated values, mixed from the
elements the way SRIM does it internally. Returns a *linear* stopping power [in mₑc² per
cm], as `bethe` does.

**The point is fidelity to SRIM, not accuracy.** A user turns this on to reproduce SRIM
numbers — to compare against a SRIM-based workflow, or to feed a code that was tuned on
them. Everything in the corrected Bethe path is about getting the physics right; this is
about getting *SRIM* right, and the two are not the same thing.

SRIM's compound calculation is, in full,

    S_elec = ( Σⱼ ωⱼ S_elecⱼ(E,phase) ) * ( 1 + (CC - 1)*g(E) )

and this implements it with `CC` pinned to 1, so the second factor is identically 1 and
never appears. Every part of the first factor was measured rather than assumed: Bragg
mixing is exact at CC = 1 (H2O from separate elemental runs against SRIM's own two-element
compound run, 114 points, max 0.058%), and the phase flag is nothing but a choice of
elemental table (C2H2 and CH4 in both phases, max 0.076%). `ωⱼ` are mass fractions, the
same rule Bragg additivity uses and for the same reason; this is not the Eq. (53) electron
weighting that applies inside the Bethe bracket.

**The compound correction is deliberately not implemented, and there is no way to set it.**
CC is a property of the molecule, not of the interaction, so one setting on
`Inelastic_Collision` would apply one molecule's bond correction to every material in the
geometry. The numbers are not available either: SRIM computes CC inside SR.exe from the
bond inventory in Data/Compound.dat by the Ziegler-Manoyan "cores and bonds" model, and the
bond stopping coefficients are in no shipped SRIM data file — the only way to obtain one is
to read `Bragg Correction = X%` off SRIM's Compound Dictionary by hand. And what SRIM would
supply is questionable anyway: its own help states the correction is computed for He ions
and that *"It is not known what is the compound correction for heavy ions (Z > 3)"*, yet it
applies that He-derived number to protons unchanged.

**What that costs, plainly.** For a molecular target the correction SRIM's dictionary would
have applied reaches about 6% below 1 MeV — liquid water is -6.00% flat from 10 keV to
400 keV — fading to zero above about 2.5 MeV/u. So below a few MeV this reproduces *SRIM
run on a composition*, which is not the same thing as *SRIM's compound library entry*. A
user comparing against a number read off the SRIM GUI for a compound in the dictionary will
see that difference.

**SRIM has no Fermi density effect.** Not a crude one — none. Above roughly 500 MeV/u these
tables are therefore wrong in a way that grows with energy: for liquid water the density
effect removes about 1.8% of the stopping number at 2 GeV and about 8.5% at 10 GeV, and
none of it is here. It is also why SRIM cannot tell liquid water from water vapour except
through the gas flag. The tables still run to 1 GeV/u because a user comparing against SRIM
wants SRIM's numbers over SRIM's whole range, including where they are wrong — but this is
the wrong switch for high-energy work, and the corrected Bethe path is the right one.

Interpolation is monotone cubic Hermite (Fritsch-Carlson) on `log(E)` against `log(S)`,
built once per (element, particle, phase) and cached. It is chosen over a natural cubic
spline for robustness rather than accuracy: it is local, so one bad point cannot ring
through the table, it is non-oscillatory by construction, and `natural_cubic_spline` would
assemble and solve a dense 4n x 4n system — 3120 x 3120 for 780 points — per element per
phase for no measurable gain. The shipped grid is never subsampled or resampled: it carries
exact nodes at SRIM's internal handovers, 1 MeV for protons and 8 MeV for alphas, with
helium having a third of its own at 2 keV.

Outside the tabulated range this errors rather than extrapolating or falling back to
`bethe`, because a quiet fallback puts a discontinuity in the middle of a stopping-power
curve where the user cannot see it.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `ωz::Vector{Float64}` : weight fraction of the elements composing the material.
- `ρ::Float64` : material density.
- `Ei::Float64` : incoming particle energy [in mₑc²].
- `particle::Particle` : incoming particle.
- `state_of_matter::String` : material state.

# Output Argument(s)
- `S::Float64` : linear stopping power [in mₑc² per cm].

# Reference(s)
- Ziegler, Ziegler and Biersack (2010), SRIM - The Stopping and Range of Ions in Matter.
- Ziegler and Manoyan (1988), The stopping of ions in compounds.
- Fritsch and Carlson (1980), Monotone piecewise cubic interpolation.

"""
function srim_stopping_power(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,Ei::Float64,particle::Particle,state_of_matter::String)

    mₑ = 0.51099895069

    #----
    # Resolve the ion and the phase
    #----
    if is_proton(particle)
        ion = "proton"
    elseif is_alpha(particle)
        ion = "alpha"
    else
        error("There is no SRIM stopping-power table for the '$(get_type(particle))' particle. " *
              "set_use_srim_sp is only supported for Proton and Alpha; turn it off, or restrict " *
              "the calculation to those particles.")
    end
    phase = srim_phase(state_of_matter)

    data = fast_load("srim_stopping_power.jld2")
    E_grid = data[ion]["E"]
    S_elec = data[ion][phase]["S_elec"]

    #----
    # Range check, on the energy and on the elements
    #----
    E_MeV = Ei*mₑ
    if E_MeV < E_grid[1] || E_MeV > E_grid[end]
        error("Energy $E_MeV MeV is outside the SRIM table for $ion, which covers " *
              "$(E_grid[1]) to $(E_grid[end]) MeV. Narrow the group structure, or turn " *
              "set_use_srim_sp off.")
    end
    for z in Z
        if z < 1 || z > 92
            error("SRIM stopping-power tables cover Z = 1 to 92 and the material contains " *
                  "Z = $z. Turn set_use_srim_sp off for this material.")
        end
    end

    #----
    # Bragg mix over the elements, in mass fractions
    #----
    # The interpolant is built once per (element, particle, phase) and cached: building it
    # inside sp would rebuild a 780-point derivative array at every quadrature point of
    # every energy group.
    logE = log(E_MeV)
    S = 0.0
    for i in eachindex(Z)
        spline = fast_spline(log.(E_grid),log.(S_elec[Z[i],:]),Z[i],string("srim_sp_",ion,"_",phase))
        S += ωz[i] * exp(spline(logE))
    end

    # The tables hold a mass stopping power and sp returns a linear one. The density here is
    # the material's, which is the correct one; the 1.0 g/cm³ of the generating runs never
    # enters, because SRIM's mass stopping power is exactly density-independent.
    return S * ρ/mₑ
end


#----
# User-supplied stopping-power tables
#----
"""
    resolve_user_stopping_power_file(file::String)

Resolves the path of a user-supplied stopping-power table. An absolute path is used as
given; anything else is taken relative to `data/` in the package root.

# Input Argument(s)
- `file::String` : file name or path.

# Output Argument(s)
- `path::String` : resolved absolute path.

"""
function resolve_user_stopping_power_file(file::String)
    dir  = joinpath(find_package_root(),"data")
    path = isabspath(file) ? file : joinpath(dir,file)
    if !isfile(path)
        error("The stopping-power file '$file' was not found. Tried '$path'. A relative name " *
              "is taken relative to '$dir', which is where such a table is expected; give an " *
              "absolute path to read it from anywhere else.")
    end
    return abspath(path)
end

"""
    parse_user_stopping_power(path::String)

Parses and validates a two-column stopping-power table.

Column 1 is the projectile's **total kinetic energy in MeV** — not energy per nucleon, and
not the `mₑc²` that `bethe` works in. Column 2 is the **mass** stopping power in
**MeV cm²/g**. The separator is whitespace or a comma, blank lines are skipped, and a line
whose first non-blank character is `#` or `%` is a comment.

Nothing else may be in the file. A third column is an error rather than something to
ignore: a user who pasted in SRIM's output needs to be told, because the extra column there
is *nuclear* stopping and adding it in silently would be wrong by percent at the low end.

Energies must be strictly increasing and are never sorted silently — an out-of-order row is
far more often a corrupt paste than a descending table, and sorting turns it into a smooth
wrong curve.

Only the *file* is validated, never the numbers. A user reaches for this because their
numbers are not Radiant's, so there is no comparison against `bethe`, no warning when the
two disagree, and no check that the curve has a Bragg peak or falls like 1/β² at the top.

# Input Argument(s)
- `path::String` : resolved absolute path.

# Output Argument(s)
- `E_MeV::Vector{Float64}` : energies [in MeV], strictly increasing.
- `S_mass::Vector{Float64}` : mass stopping powers [in MeV cm²/g].

"""
function parse_user_stopping_power(path::String)

    # readlines rather than a streaming eachline: every check below can throw, and a stream
    # left open by an error keeps a lock on the file on Windows.
    E_MeV = Float64[]; S_mass = Float64[]; lines_of = Int64[]
    for (n,line) in enumerate(readlines(path))
        str = strip(line)
        isempty(str) && continue
        (str[1] == '#' || str[1] == '%') && continue
        fields = split(replace(str,","=>" "))
        if length(fields) != 2
            error("Line $n of '$path' has $(length(fields)) fields; a stopping-power table " *
                  "has exactly two, energy [MeV] and mass stopping power [MeV cm²/g]. If this " *
                  "came from SRIM, the extra column is nuclear stopping and must be removed " *
                  "rather than added in.")
        end
        e = tryparse(Float64,fields[1]); s = tryparse(Float64,fields[2])
        if isnothing(e) || isnothing(s)
            error("Line $n of '$path' does not parse as two numbers: '$str'.")
        end
        if !isfinite(e) || !isfinite(s) || e ≤ 0 || s ≤ 0
            error("Line $n of '$path' has a non-positive or non-finite value ($e, $s). The " *
                  "interpolation is in log-log and neither logarithm exists there.")
        end
        if !isempty(E_MeV) && e ≤ E_MeV[end]
            error("Line $n of '$path' has energy $e MeV, which does not exceed the $(E_MeV[end]) " *
                  "MeV on line $(lines_of[end]). Energies must be strictly increasing; they are " *
                  "not sorted for you, because an out-of-order row is more often a corrupt paste " *
                  "than a descending table.")
        end
        push!(E_MeV,e); push!(S_mass,s); push!(lines_of,n)
    end
    if length(E_MeV) < 3
        error("'$path' holds $(length(E_MeV)) data rows; at least 3 are needed, because a " *
              "monotone cubic needs a slope on both sides of a node.")
    end
    return E_MeV, S_mass
end

"""
    user_stopping_power_table(file::String)

Gives the parsed columns of a user-supplied stopping-power table, reading and validating the
file on first use and caching it by resolved path afterwards, so a file registered on ten
materials is read once.

# Input Argument(s)
- `file::String` : file name or path.

# Output Argument(s)
- `E_MeV::Vector{Float64}` : energies [in MeV].
- `S_mass::Vector{Float64}` : mass stopping powers [in MeV cm²/g].

"""
function user_stopping_power_table(file::String)
    path = resolve_user_stopping_power_file(file)
    id = string("user_sp_",path)
    if ~haskey(cache_radiant[],id)
        cache_radiant[][id] = parse_user_stopping_power(path)
    end
    return cache_radiant[][id]
end

"""
    user_stopping_power(file::String,Ei::Float64,ρ::Float64)

Gives the electronic stopping power read from a user-supplied table. Returns a *linear*
stopping power [in mₑc² per cm], as `bethe` does.

**This is the one path with no validation behind it.** Everything else in Radiant's
stopping-power chain was measured against NIST or SRIM; a user table is measured against
nothing, because its correctness is the user's to establish. Radiant validates the *file*;
it cannot validate the *numbers*.

The table must be a **total electronic** stopping power. `sp` computes `Sr = Stot - Sc`
where `Sc` is the catastrophic part it derives from its own differential cross-section, so
feeding in a restricted table subtracts the catastrophic part twice, and a nuclear-stopping
column — or the sum of the two — is the most likely thing to be pasted in by mistake, small
enough at first glance to survive a sanity check while being percent-level and rising below
about 100 keV/u. The soft/catastrophic split is therefore approximate on this path in
exactly the way it is for the `"icru49"` low-energy branch: `Stot` comes from the user's
source and `Sc` from Radiant's analytic cross-section, and the two were never reconciled.

**The density that scales the table is the material's, not the one the table came from.**
For everything except the density effect that is correct, because mass stopping power is
density-independent. But the density effect *is* baked into whatever numbers were supplied,
at whatever density their source assumed, and rescaling cannot touch it. A user who takes a
water-vapour table and hands it to a liquid-water material gets vapour physics scaled to
liquid density, silently. Radiant cannot detect this and does not pretend to.

Interpolation is monotone cubic Hermite (Fritsch-Carlson) on `log(E)` against `log(S)`,
built once and cached by file path. It passes through every tabulated point to round-off,
which is the whole reason a user hands in a table; it cannot overshoot, where a natural
cubic spline rings near the Bragg peak and produces a stopping power that rises where the
data falls; it is local, so one mistyped row corrupts two intervals rather than the whole
curve; and it is two passes over the array rather than a dense 4n x 4n solve. Log-log
because `log S` is close to linear in `log E` above the Bragg peak, so the curve the
interpolant has to follow is nearly straight, and because it guarantees `S > 0` between
nodes.

Outside the tabulated range this errors. It does not extrapolate, does not fall back to
`bethe` and does not clamp to the endpoint: the fix is always to extend the table or to
narrow the group structure, and a silent fallback would put a discontinuity in the middle of
a stopping-power curve where the user cannot see it.

# Input Argument(s)
- `file::String` : file name or path.
- `Ei::Float64` : incoming particle energy [in mₑc²].
- `ρ::Float64` : material density.

# Output Argument(s)
- `S::Float64` : linear stopping power [in mₑc² per cm].

# Reference(s)
- Fritsch and Carlson (1980), Monotone piecewise cubic interpolation.

"""
function user_stopping_power(file::String,Ei::Float64,ρ::Float64)

    mₑ = 0.51099895069
    E_MeV, S_mass = user_stopping_power_table(file)
    E = Ei*mₑ

    if E < E_MeV[1] || E > E_MeV[end]
        error("Energy $E MeV is outside the stopping-power table '$file', which covers " *
              "$(E_MeV[1]) to $(E_MeV[end]) MeV. The table registered with set_stopping_power " *
              "must span the whole energy range of the group structure for that particle: " *
              "extend the table, or narrow the group structure.")
    end

    spline = fast_spline(log.(E_MeV),log.(S_mass),1,string("user_sp_",resolve_user_stopping_power_file(file)))
    return exp(spline(log(E))) * ρ/mₑ
end
