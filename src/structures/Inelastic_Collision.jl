"""
Inelastic_Collision

Structure used to define parameters for production of multigroup inelastic collisional cross-sections.

# Mandatory field(s)
- N/A

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}} = Dict((Positron,Positron) => ["S"],(Positron,Electron) => ["P"],(Electron,Electron) => ["S","P"])` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which values correspond:
    - `(Positron,Positron) => ["S"]` : scattering of incident positrons
    - `(Electron,Electron) => ["S"]` : scattering of incident electrons
    - `(Positron,Electron) => ["P"]` : production of electrons by incident positrons
    - `(Electron,Electron) => ["P"]` : production of electrons by incident electrons

"""
mutable struct Inelastic_Collision <: Interaction

    # Variable(s)
    name::String
    incoming_particle::Vector{Type}
    interaction_particles::Vector{Type}
    interaction_types::Dict{Tuple{Type,Type},Vector{String}}
    is_CSD::Bool
    is_AFP::Bool
    is_AFP_decomposition::Bool
    is_elastic::Bool
    is_subshells_dependant::Bool
    is_shell_correction::Bool
    is_focusing_møller::Bool
    is_hydrogenic_distribution_term::Bool
    is_distant_collision::Bool
    density_correction::String
    scattering_model::String
    user_ecut::Dict{Type,Float64}
    user_emax::Dict{Type,Float64}
    user_emax_fit::Set{Type}
    low_energy_model::String
    low_energy_source::String
    is_bragg_additivity::Bool

    # Constructor(s)
    function Inelastic_Collision()
        this = new()
        this.name = "Inelastic_Collision"
        this.interaction_types = Dict((Positron,Positron) => ["S"],(Positron,Electron) => ["P"],(Electron,Electron) => ["S","P"],(Proton,Proton) => ["S"],(Proton,Electron) => ["P"],(Alpha,Alpha) => ["S"],(Alpha,Electron) => ["P"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.set_density_correction("fano")
        this.is_CSD = true
        this.is_AFP = true
        this.is_AFP_decomposition = false
        this.is_elastic = false
        this.is_focusing_møller = false
        this.is_hydrogenic_distribution_term = false
        this.set_is_subshells_dependant(true)
        this.set_is_shell_correction(true)
        this.set_scattering_model("BFP")
        this.user_ecut = Dict{Type,Float64}()
        this.user_emax = Dict{Type,Float64}()
        this.user_emax_fit = Set{Type}()
        this.set_low_energy_model("analytic")
        this.set_low_energy_source("srim")
        this.set_is_bragg_additivity(false)

        return this
    end
end

# Method(s)
"""
    set_interaction_types(this::Inelastic_Collision,interaction_types::Dict{Tuple{DataType,DataType},Vector{String}})

To define the interaction types for inelastic collision processes.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}}` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `(Positron,Positron) => ["S"]` : scattering of incident positrons
    - `(Electron,Electron) => ["S"]` : scattering of incident electrons
    - `(Positron,Electron) => ["P"]` : production of electrons by incident positrons
    - `(Electron,Electron) => ["P"]` : production of electrons by incident electrons

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_collision = Inelastic_Collision()
julia> elastic_collision.set_interaction_types( Dict((Electron,Electron) => ["S"]) ) # No knock-on electrons.
```
"""
function set_interaction_types(this::Inelastic_Collision,interaction_types)
    this.interaction_types = interaction_types
end

"""
    set_density_correction(this::Inelastic_Collision,density_correction::String)

Set the Fermi density correction.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `density_correction::String` : type of density effect:
    - `fano` : Fano density effect.
    - `sternheimer` : Sternheimer semi-empirical density effect.
    - `none` : no density effect.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_collision = Inelastic_Collision()
julia> elastic_collision.set_density_correction("sternheimer")
```
"""
function set_density_correction(this::Inelastic_Collision,density_correction::String)
    if lowercase(density_correction) ∉ ["fano","sternheimer","none"] error("Unknown '$density_correction' density correction.") end
    this.density_correction = lowercase(density_correction)
end

"""
    set_is_shell_correction(this::Inelastic_Collision,is_shell_correction::Bool)

Activate or desactivate shell correction for stopping powers.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `is_shell_correction::Bool` : activate (true) or desactivate (false) the shell correction.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_collision = Inelastic_Collision()
julia> elastic_collision.set_is_shell_correction(false)
```
"""
function set_is_shell_correction(this::Inelastic_Collision,is_shell_correction::Bool)
    this.is_shell_correction = is_shell_correction
end

"""
    set_is_subshells_dependant(this::Inelastic_Collision,is_subshells_dependant::Bool)

Compute the inelastic cross-sections assuming bounded or unbounded electrons.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `is_subshells_dependant::Bool` : bounded (true) or unbounded (false) electrons.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_collision = Inelastic_Collision()
julia> elastic_collision.set_is_subshells_dependant(false)
```
"""
function set_is_subshells_dependant(this::Inelastic_Collision,is_subshells_dependant::Bool)
    this.is_subshells_dependant = is_subshells_dependant
end

"""
    scattering_model(this::Inelastic_Collision,scattering_model::String)

To define the solver for inelastic scattering.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic structure.
- `solver::String` : solver for inelastic scattering, which can be:
    - `BFP` : Boltzmann Fokker-Planck solver.
    - `FP` : Fokker-Planck solver.
    - `BTE` : Boltzmann solver.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> elastic_collision = Inelastic_Collision()
julia> elastic_collision.set_scattering_model("FP")
```
"""
function set_scattering_model(this::Inelastic_Collision,scattering_model::String)
    if uppercase(scattering_model) ∉ ["BFP","FP","BTE"] error("Unknown scattering model (should be BFP, FP or BTE).") end
    this.scattering_model = uppercase(scattering_model)
end

"""
    set_ecut(this::Inelastic_Collision,particle::Type,energy::Real)

To set the energy above which the corrected Bethe formula is used for a given particle,
and below which the stopping power is taken from the low-energy branch.

Lowering it extends the Bethe formula into a region where its corrections are less well
determined; the shell-correction tables stop at a tabulated γ and raise an error naming
this setter when they run out.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `particle::Type` : particle type, one of `Electron`, `Positron`, `Muon`, `Antimuon`,
  `Proton`, `Antiproton` or `Alpha`.
- `energy::Real` : cutoff energy [in MeV].

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> inelastic_collision = Inelastic_Collision()
julia> inelastic_collision.set_ecut(Proton,1.0)
```
"""
function set_ecut(this::Inelastic_Collision,particle::Type,energy::Real)
    if !haskey(ECUT_DEFAULTS,particle) error("The cutoff energy is undefined for the '$particle' particle. It is defined for $(collect(keys(ECUT_DEFAULTS))).") end
    if energy ≤ 0 error("The cutoff energy should be positive.") end
    this.user_ecut[particle] = float(energy)
end

"""
    get_ecut(this::Inelastic_Collision,particle::Type)

To get the energy above which the corrected Bethe formula is used for a given particle:
the value `set_ecut` stored, or `ECUT_DEFAULTS[particle]` if none was set.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `particle::Type` : particle type.

# Output Argument(s)
- `energy::Float64` : cutoff energy [in MeV].

# Examples
```jldoctest
julia> inelastic_collision = Inelastic_Collision()
julia> inelastic_collision.get_ecut(Proton)
0.75
```
"""
function get_ecut(this::Inelastic_Collision,particle::Type)
    if !haskey(ECUT_DEFAULTS,particle) error("The cutoff energy is undefined for the '$particle' particle. It is defined for $(collect(keys(ECUT_DEFAULTS))).") end
    return get(this.user_ecut,particle,ECUT_DEFAULTS[particle])
end

"""
    set_emax(this::Inelastic_Collision,particle::Type,energy::Real)
    set_emax(this::Inelastic_Collision,particle::Type,mode::String)

To set the energy of the stopping-power maximum of the analytic low-energy extension.

`Emax` is read only inside that extension and is meaningful only there: both of its
branches have zero slope at `Ei/Emax = 1`, so `Emax` is where the extension puts the
stopping-power maximum. Nothing at or above `Ecut` changes, and neither does the
`"icru49"` low-energy branch, which does not use it.

Passing `"f"` selects `emax_stopping_maximum`, a power-law fit in the mean atomic number
to the measured maximum of all 92 elements against SRIM-2013. A single constant asserts
the same peak energy in hydrogen as in uranium, where the real peak runs 41-164 keV for
protons and 0.42-1.19 MeV for alphas. `"f"` is accepted for `Proton` and `Alpha` only;
there is no Z-dependent fit for the others.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `particle::Type` : particle type.
- `energy::Real` : energy of the stopping-power maximum [in MeV].
- `mode::String` : `"f"`, to select the Z-dependent fit.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> inelastic_collision = Inelastic_Collision()
julia> inelastic_collision.set_emax(Proton,"f")
julia> inelastic_collision.set_emax(Alpha,0.9)
```
"""
function set_emax(this::Inelastic_Collision,particle::Type,energy::Real)
    if !haskey(EMAX_DEFAULTS,particle) error("The stopping-power maximum is undefined for the '$particle' particle. It is defined for $(collect(keys(EMAX_DEFAULTS))).") end
    if energy ≤ 0 error("The energy of the stopping-power maximum should be positive.") end
    delete!(this.user_emax_fit,particle)
    this.user_emax[particle] = float(energy)
end

function set_emax(this::Inelastic_Collision,particle::Type,mode::String)
    if lowercase(mode) != "f" error("Unknown '$mode' mode for the stopping-power maximum. Give an energy [in MeV], or \"f\" for the Z-dependent fit.") end
    if particle ∉ [Proton,Alpha] error("The \"f\" mode is only supported for Proton and Alpha; there is no Z-dependent fit for the '$particle' particle.") end
    delete!(this.user_emax,particle)
    push!(this.user_emax_fit,particle)
end

"""
    get_emax(this::Inelastic_Collision,particle::Type)

To get the energy of the stopping-power maximum of the analytic low-energy extension.

Gives the value `set_emax` stored, or `EMAX_DEFAULTS[particle]` if none was set. Returns
`"f"` unchanged when the Z-dependent fit is selected: the resolved value depends on the
material and cannot be known here.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `particle::Type` : particle type.

# Output Argument(s)
- `energy::Union{Float64,String}` : energy of the stopping-power maximum [in MeV], or `"f"`.

# Examples
```jldoctest
julia> inelastic_collision = Inelastic_Collision()
julia> inelastic_collision.get_emax(Proton)
0.05
```
"""
function get_emax(this::Inelastic_Collision,particle::Type)
    if !haskey(EMAX_DEFAULTS,particle) error("The stopping-power maximum is undefined for the '$particle' particle. It is defined for $(collect(keys(EMAX_DEFAULTS))).") end
    if particle in this.user_emax_fit return "f" end
    return get(this.user_emax,particle,EMAX_DEFAULTS[particle])
end

"""
    set_low_energy_model(this::Inelastic_Collision,model::String)

To set the model used for the stopping power below `Ecut`.

- `analytic` : analytical extension anchored on the value and the slope of the Bethe
  branch at `Ecut`. Smooth at the join, but it extrapolates from that single anchor and
  its error grows with distance from it.
- `icru49` : the empirical form of ICRU Report 49, attributed there to Andersen and
  Ziegler, with per-element coefficients fitted to SRIM-2013. The low branch exponent is
  fitted per element rather than pinned at 0.45, which is what removes the separate
  `sqrt(E)` branch the pinned version carries below 10 keV. It is confined to [0.30, 0.80],
  a departure from the parametrization as published: unconstrained, seven of the 92 proton
  fits ran to 0.8-1.31 with 22-43% error extrapolating below the fitting window.

The default is `analytic`, and it stays the default in this change: that the evidence
favours the fitted branch does not make a default swap part of it.

Measured against NIST PSTAR/ASTAR below `Ecut` over the 74 ICRU 49 materials, median
deviation goes from 18.72% to 2.67% for protons and 26.40% to 3.60% for alphas. Against
the same form with the exponent pinned at 0.45 the improvement is decisive for alphas and
close to a coin flip for protons -- 3.33% to 2.67%, with the whole proton gain coming from
below 10 keV, where the fitted exponent replaces the pinned form's bolted-on `sqrt(E)`.
That region is where PSTAR and ASTAR are themselves fits to experiment in this same
empirical form, so the comparison is between two empirical fits of one functional form and
not between a model and theory.

**10 keV is the floor of validity.** Below about 5 keV this branch is badly wrong for some
elements -- helium worst, 153% median below 10 keV -- and no choice of exponent fixes it.

Above `Ecut` the fitted form is blended onto the Bethe branch over one octave, so that the
handover is C¹ rather than leaving a step of up to 13.6% in the CSD operator.

Below `Ecut` with `icru49` selected, `sp` computes `Sr = Stot - Sc` where `Stot` is a fit
to SRIM's *total electronic* stopping while `Sc` comes from the analytic heavy-particle
differential cross-section, so the soft/catastrophic split is approximate in that region.
In practice it rarely bites -- a 750 keV proton can transfer at most about 1.6 keV to an
electron, so `Sc` is zero unless the electron group structure reaches below that -- and it
is not new, since the analytic extension is not a Bethe evaluation either.

The fitted table covers Z = 1-92, protons and alphas only. Any other particle, or any
element outside that range, falls back to the analytic extension silently.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `model::String` : low-energy model, `"analytic"` or `"icru49"`.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> inelastic_collision = Inelastic_Collision()
julia> inelastic_collision.set_low_energy_model("icru49")
```

# Reference(s)
- ICRU Report 49 (1993), Stopping Powers and Ranges for Protons and Alpha Particles.
- Andersen and Ziegler (1977), Hydrogen: Stopping Powers and Ranges in All Elements.
- Salvat (2022), Bethe stopping-power formula and its corrections. Eqs. (134)-(136), for
  the free-exponent parametrization and the fitting prescription.
"""
function set_low_energy_model(this::Inelastic_Collision,model::String)
    if lowercase(model) ∉ ["analytic","icru49"] error("Unknown '$model' low-energy model.") end
    this.low_energy_model = lowercase(model)
end

"""
    get_low_energy_model(this::Inelastic_Collision)

To get the model used for the stopping power below `Ecut`.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.

# Output Argument(s)
- `model::String` : low-energy model.

# Examples
```jldoctest
julia> inelastic_collision = Inelastic_Collision()
julia> inelastic_collision.get_low_energy_model()
"analytic"
```
"""
function get_low_energy_model(this::Inelastic_Collision)
    return this.low_energy_model
end

"""
    set_low_energy_source(this::Inelastic_Collision,source::String)

To choose which coefficients the fitted low-energy branch reads.

- `srim` : fitted to SRIM-2013, all 92 elements from one fit.
- `iaea` : fitted to the IAEA electronic stopping power database — the compiled
  experimental literature, 720 references, publications from 1928 to 2025. Same form, same
  units, same phase convention, same `a2` bound; only the numbers differ.

**It is inert unless `set_low_energy_model("icru49")` is also selected.** There are no
coefficients in the analytic branch to choose between, and results are bit-identical
between the two sources there.

`iaea` is a **partial override on top of the SRIM table**, not a table in its own right:
it covers 41 elements for protons and 30 for alphas, and it is resolved per element at the
point of use, so a compound routinely draws from both at once. An element is missing only
because the measurements cannot determine the form — no points on one side of the stopping
maximum, too few points, or nothing reaching `Ecut` — and every element it omits is served
by the SRIM entry through the same lookup. So `iaea` never fails, never warns, and leaves
most materials untouched.

Each entry carries `n_references` and a leave-one-reference-out cross-validated error
under `["quality"]`, and those are the numbers to read before relying on one. The
cross-validated median is 4.38% for protons and 3.38% for alphas, against in-sample
medians of 2.83% and 2.03%; the tail is real and is not hidden — terbium fits its own alpha
points well and misses by 44% when a reference is held out. Radiant does not gate on these
numbers: which entries to trust is the user's judgement, and a filter applied here would
be invisible to them.

**This is a choice of reference, not an upgrade.** Scored against NIST PSTAR/ASTAR, which
neither table was fitted to, `iaea` is closer for 11 of 23 elements for protons and 6 of
19 for alphas, and the median is *worse* for both — 1.72% to 2.38% and 1.85% to 3.13%. The
gains and the losses are both large where they occur: Pt 10.68% to 3.67% against Pb 2.96%
to 7.19% for protons, Ar 1.42% to 0.42% against Ti 1.04% to 11.27% for alphas. ICRU 49 is
itself an evaluation of experiment in this same functional form, so this is two empirical
fits weighing the same disagreeing literature differently, not a model scored against
theory: `srim` tracks ICRU 49 more closely because SRIM and ICRU 49 are two evaluations by
overlapping communities of the same older measurements, while `iaea` tracks the current
compiled database, which includes measurements neither of them saw. The measurements
disagree among themselves by 14-17%.

The honest reason to select `iaea` is wanting elemental stopping tied to the measurement
compilation rather than to SRIM's evaluation of it, or working on an element where the two
disagree and having a reason to trust the measurements. On the evidence above, for alphas
it is usually worse.

The default stays `srim`: it is the only source covering all 92 elements from one fit, and
defaulting to a partial one would put a provenance seam in the middle of the table where a
trend across Z could no longer be told from a change of data source.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `source::String` : coefficient source, `"srim"` or `"iaea"`.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> inelastic_collision = Inelastic_Collision()
julia> inelastic_collision.set_low_energy_model("icru49")
julia> inelastic_collision.set_low_energy_source("iaea")
```

# Reference(s)
- IAEA, Electronic Stopping Power of Matter for Ions (database),
  https://www-nds.iaea.org/stopping/. Sole source of the "iaea" coefficients.
"""
function set_low_energy_source(this::Inelastic_Collision,source::String)
    if lowercase(source) ∉ ["srim","iaea"] error("Unknown '$source' low-energy coefficient source.") end
    this.low_energy_source = lowercase(source)
end

"""
    get_low_energy_source(this::Inelastic_Collision)

To get the coefficient source of the fitted low-energy branch. Inert unless the `icru49`
low-energy model is selected.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.

# Output Argument(s)
- `source::String` : coefficient source.

# Examples
```jldoctest
julia> inelastic_collision = Inelastic_Collision()
julia> inelastic_collision.get_low_energy_source()
"srim"
```

# Reference(s)
- IAEA, Electronic Stopping Power of Matter for Ions (database),
  https://www-nds.iaea.org/stopping/. Sole source of the "iaea" coefficients.
"""
function get_low_energy_source(this::Inelastic_Collision)
    return this.low_energy_source
end

"""
    set_is_bragg_additivity(this::Inelastic_Collision,is_bragg_additivity::Bool)

To build the stopping power of a compound from its elements rather than from its own
composition. Each element is evaluated separately and the results are combined with mass
fractions, `S = Σⱼ ωⱼ Sⱼ`, which is what mass stopping power is additive in.

Every element is evaluated at the *compound* density, because the density effect in a
compound is set by the compound's electron density and not by what the pure element would
do on its own.

**Scope: the Bethe branch only.** The fitted `icru49` low-energy branch is already
unconditionally additive whatever this flag says, so scoping the flag here is what stops
it meaning different things in different energy regions. With the `analytic` model each
element's extension anchors to that element's own Bethe value and slope at `Ecut`.

**Enabling this overrides the tabulated compound mean excitation energy.**
`effective_mean_excitation_energy` returns the ICRU value for canonical compounds — 75 eV
for liquid water — before falling back to additivity. Building water from H and O
elementwise means each element carries its own I and that 75 eV never enters. For
electrons this is a visible shift, and in the opposite direction from what a user who
trusts ICRU would want, so it is a choice to be made knowingly.

Note also that any additivity-based kernel misses NIST's liquid water by about 28% at
1 keV and 18% at 100 keV. Applying additivity to NIST's own hydrogen and oxygen tables
reproduces the same gap, so that is the assumption failing rather than the implementation.

Applies to all particles.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `is_bragg_additivity::Bool` : activate (true) or deactivate (false) Bragg additivity.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> inelastic_collision = Inelastic_Collision()
julia> inelastic_collision.set_is_bragg_additivity(true)
```
"""
function set_is_bragg_additivity(this::Inelastic_Collision,is_bragg_additivity::Bool)
    this.is_bragg_additivity = is_bragg_additivity
end

"""
    get_is_bragg_additivity(this::Inelastic_Collision)

To know whether compounds are built from their elements by Bragg additivity.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.

# Output Argument(s)
- `is_bragg_additivity::Bool` : whether Bragg additivity is enabled.

# Examples
```jldoctest
julia> inelastic_collision = Inelastic_Collision()
julia> inelastic_collision.get_is_bragg_additivity()
false
```
"""
function get_is_bragg_additivity(this::Inelastic_Collision)
    return this.is_bragg_additivity
end

"""
    in_distribution(this::Inelastic_Collision)

Describe the energy discretization method for the incoming particle in the inelastic collision
interaction.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function in_distribution(this::Inelastic_Collision)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    out_distribution(this::Inelastic_Collision)

Describe the energy discretization method for the outgoing particle in the inelastic collision
interaction.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.

"""
function out_distribution(this::Inelastic_Collision)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    bounds(this::Inelastic_Collision,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String,
    Ec::Float64,Ui::Float64,particle::Particle)

Gives the integration energy bounds for the outgoing particle for inelastic collision
interaction. 

# Input Argument(s)
- `this::Annihilation` : annihilation structure.
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `Ei::Float64` : energy of the incoming particle.
- `type::String` : type of interaction.
- `Ec::Float64` : energy cutoff between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.

# Output Argument(s)
- `Ef⁻::Float64` : upper bound.
- `Ef⁺::Float64` : lower bound.
- `isSkip::Bool` : define if the integration is skipped or not.

"""
function bounds(this::Inelastic_Collision,Ef⁻::Float64,Ef⁺::Float64,Ei::Float64,type::String,Ec::Float64,Ui::Float64,particle::Particle)
    if is_electron(particle)
        # Scattered electron
        if type == "S"
            Ef⁻ = min(Ef⁻,Ec,Ei-Ui)
            Ef⁺ = max(Ef⁺,(Ei-Ui)/2)
        # Knock-on electron
        elseif type == "P"
            Ef⁻ = min(Ef⁻,(Ei-Ui)/2)
            Ef⁺ = max(Ef⁺,0.0)
        else
            error("Unknown type.")
        end
    elseif is_positron(particle)
        # Scattered positron/proton/alpha
        if type == "S"
            Ef⁻ = min(Ef⁻,Ec,Ei-Ui)
            Ef⁺ = max(Ef⁺,0.0)
        # Knock-on electron
        elseif type == "P"
            Ef⁻ = min(Ef⁻,Ei-Ui)
            Ef⁺ = max(Ef⁺,0.0)
        else
            error("Unknown type.")
        end
    elseif is_proton(particle) || is_alpha(particle)
        mₑ = 0.51099895069
        ratio_mass = get_mass(particle)/mₑ
        γ = (Ei + ratio_mass)/ratio_mass
        β² = (γ^2-1)/γ^2
        R = 1/(1+(1/ratio_mass)^2+2*γ/ratio_mass)
        W_ridge = 2*β²*γ^2*R
        # Scattered heavy particle
        if type == "S"
            Ef⁻ = min(Ef⁻,Ec,Ei-Ui)
            Ef⁺ = max(Ef⁺,0.0,Ei-W_ridge)
        # Knock-on electron
        elseif type == "P"
            Ef⁻ = min(Ef⁻,Ei-Ui,W_ridge-Ui)
            Ef⁺ = max(Ef⁺,0.0)
        else
            error("Unknown type.")
        end
    else
        error("Unknown particle")
    end
    if (Ef⁻-Ef⁺ < 0) isSkip = true else isSkip = false end
    return Ef⁻,Ef⁺,isSkip
end

"""
    dcs(this::Inelastic_Collision,L::Int64,Ei::Float64,Ef::Float64,type::String,
    particle::Particle,Ui::Float64,Zi::Real,Ti::Float64)

Gives the Legendre moments of the scattering cross-sections for inelastic collision
interaction. 

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `Z::Int64` : atomic number.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : incoming particle energy.
- `Ef::Float64` : outgoing particle energy.
- `type::String` : type of interaction.
- `particle::Particle` : incoming particle.
- `δi::Int64` : subshell index.
- `Ui::Float64` : binding energy of the subshell.
- `Zi::Int64` : number of electrons in the subshell.
- `Ti::Float64` : mean kinetic energy of electrons in the subshell.
- `ri::Float64` : average radius per subshell.

# Output Argument(s)
- `σl::Vector{Float64}` : Legendre moments of the scattering cross-sections.

"""
function dcs(this::Inelastic_Collision,Z::Int64,L::Int64,Ei::Float64,Ef::Float64,type::String,particle::Particle,δi::Int64,Ui::Float64,Zi::Real,Ti::Float64,ri::Float64)

    #----
    # Initialization
    #----
    σs = 0.0
    σl = zeros(L+1)
    if type == "S"
        E⁺ = Ef
        E⁻ = Ei-Ui-E⁺
        W = E⁻+Ui
    elseif type == "P"
        E⁻ = Ef
        E⁺ = Ei-Ui-E⁻
        W = E⁻+Ui
    else
        error("Unknown type")
    end
    if W < 0 return σl end

    #----
    # Close collisions
    #----
    is_heavy = false
    sigma_reg = 0.0
    A_over_W2 = 0.0
    cache = nothing
    if is_electron(particle)
        σs += moller(Zi,Ei,W,Ui,Ti,this.is_focusing_møller,this.is_hydrogenic_distribution_term)
    elseif is_positron(particle)
        σs += bhabha(Zi,Ei,W)
    elseif is_proton(particle) || is_alpha(particle)
        is_heavy = true
        # Create cache once for all Legendre orders
        cache = HeavyInelasticCache(Zi, Ei, particle)
        sigma_reg, A_over_W2 = heavy_sigma_components(cache, W)
    else
        error("Unknown particle")
    end

    #----
    # Compute the angular distribution
    #----
    if is_heavy
        Wl = angular_inelastic_collision_heavy_particle(Ei,Ef,L,particle)
    else
        Wl = angular_moller(Ei,Ef,L)
    end

    #----
    # Compute the Legendre moments of the cross-section
    #----
    if is_heavy
        for l in range(0,L)
            # Subtract analytic leading 1/W piece to regularize the numeric integrand
            if W > 0
                coef = heavy_leading_1overW_coeff(cache, l)
                leading_term = coef / W
            else
                leading_term = 0.0
            end
            σl[l+1] += Wl[l+1] * sigma_reg + (Wl[l+1] - 1.0) * A_over_W2 - leading_term
        end
    else
        for l in range(0,L) σl[l+1] += Wl[l+1] * σs end
    end
    return σl
end

"""
    tcs(this::Inelastic_Collision,Ei::Float64,Ec::Float64,particle::Particle,Z::Int64)

Gives the total cross-section for inelastic collision interaction. 

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure. 
- `Ei::Float64` : incoming particle energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.
- `Z::Int64` : atomic number.
- `Eout::Vector{Float64}` : energy boundaries associated with outgoing particles.

# Output Argument(s)
- `σt::Float64` : total cross-section.

"""
function tcs(this::Inelastic_Collision,Ei::Float64,Ec::Float64,particle::Particle,Z::Int64)

    #----
    # Initialization
    #----
    σt = 0

    #----
    # Close collisions
    #----
    if is_electron(particle)
        σt += integrate_moller(Z,Ei,0,Ei-Ec,this.is_focusing_møller,this.is_hydrogenic_distribution_term,this.is_subshells_dependant)
    elseif is_positron(particle)
        σt += integrate_bhabha(Z,Ei,0,Ei-Ec,this.is_subshells_dependant)
    elseif is_proton(particle) || is_alpha(particle)
        σt += integrate_inelastic_collision_heavy_particle(Z,Ei,0,particle,Ei-Ec,this.is_subshells_dependant)
    else
        error("Unknown particle")
    end

    return σt
end

"""
    acs(this::Inelastic_Collision,Ei::Float64,Ec::Float64,particle::Particle,Z::Int64)

Gives the absorption cross-section for inelastic collision interaction. 

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure. 
- `Ei::Float64` : incoming particle energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.
- `Z::Int64` : atomic number.
- `Ecutoff::Float64` : cutoff energy.

# Output Argument(s)
- `σa::Float64` : absorption cross-section.

"""
function acs(this::Inelastic_Collision,Ei::Float64,Ec::Float64,particle::Particle,Z::Int64,Ecutoff::Float64)

    #----
    # Close collisions
    #----
    if is_electron(particle)
        σa = integrate_moller(Z,Ei,0,Ei-min(Ec,Ecutoff),this.is_focusing_møller,this.is_hydrogenic_distribution_term,this.is_subshells_dependant)
    elseif is_positron(particle)
        σa = integrate_bhabha(Z,Ei,0,Ei-min(Ec,Ecutoff),this.is_subshells_dependant)
    else
        error("Unknown particle")
    end
    return σa
    
end

"""
    sp(this::Inelastic_Collision,Z::Vector{Int64},ωz::Vector{Float64},atz::Vector{Float64},
    ρ::Float64,N_density::Float64,state_of_matter::String,Ei::Float64,Ec::Float64,
    particle::Particle,I_eff::Float64,atomic_weights::Union{Nothing,Vector{Float64}})

Gives the stopping power for inelastic collision interaction.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure.
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `ωz::Vector{Float64}` : weight fraction of the elements composing the material.
- `atz::Vector{Float64}` : atomic fraction of the elements composing the material.
- `ρ::Float64` : material density.
- `N_density::Float64` : nuclei density of the material [in cm⁻³].
- `state_of_matter::String` : state of matter.
- `Ei::Float64` : incoming particle energy.
- `Ec::Float64` : cutoff energy between soft and catastrophic interactions.
- `particle::Particle` : incoming particle.
- `I_eff::Float64` : effective mean excitation energy override [in mₑc²]; NaN to use tables.
- `atomic_weights::Union{Nothing,Vector{Float64}}` : element atomic weights [u].

# Output Argument(s)
- `S::Float64` : stopping power.

"""
function sp(this::Inelastic_Collision,Z::Vector{Int64},ωz::Vector{Float64},atz::Vector{Float64},ρ::Float64,N_density::Float64,state_of_matter::String,Ei::Float64,Ec::Float64,particle::Particle,I_eff::Float64=NaN,atomic_weights::Union{Nothing,Vector{Float64}}=nothing)

    # Compute the total stopping power
    Stot = bethe(Z,ωz,atz,ρ,N_density,Ei,particle,this.density_correction,state_of_matter,I_eff;
        atomic_weights=atomic_weights,
        is_shell_correction=this.is_shell_correction,
        user_ecut=get(this.user_ecut,get_type(particle),NaN),
        user_emax=get(this.user_emax,get_type(particle),NaN),
        is_emax_fit=(get_type(particle) in this.user_emax_fit),
        low_energy_model=this.low_energy_model,
        low_energy_source=this.low_energy_source,
        is_bragg_additivity=this.is_bragg_additivity)

    # Compute the catastrophic Møller- or Bhabha- derived stopping power
    Sc = 0
    Nz = length(Z)
    for i in range(1,Nz)

        # Close collision
        if is_electron(particle)
            Sc += atz[i] * N_density * integrate_moller(Z[i],Ei,1,Ei-Ec,this.is_focusing_møller,this.is_hydrogenic_distribution_term,this.is_subshells_dependant)
        elseif is_positron(particle)
            Sc += atz[i] * N_density * integrate_bhabha(Z[i],Ei,1,Ei-Ec,this.is_subshells_dependant)
        elseif is_proton(particle) || is_alpha(particle)
            Sc += atz[i] * N_density * integrate_inelastic_collision_heavy_particle(Z[i],Ei,1,particle,Ei-Ec,this.is_subshells_dependant)
        end
    end

    # Compute the soft stopping power
    Sr = Stot - Sc

    return Sr
end

"""
    mt(this::Inelastic_Collision)

Gives the momentum transfer for inelastic collision interaction.

# Input Argument(s)
- `this::Inelastic_Collision` : inelastic collision structure. 

# Output Argument(s)
- `T::Float64` : momentum transfer.

"""
function mt(this::Inelastic_Collision)
    T = 0.0
    return T
end

