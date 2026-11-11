"""
    Material

Structure used to define a material and its properties.

# Mandatory field(s)
- `name::String` : name (or identifier) of the Material structure.
- `density::Float64` : density [in g/cm³].
- `elements::Vector{String}` : vector of the element in the composition of the material.
- `weight_fractions::Vector{Float64}` : vector of the weight fraction for each element in the composition of the material.
- `mass_numbers::Vector{Vector{Int64}}` : isotope mass numbers for each element, if specified.
- `isotope_atomic_fractions::Vector{Vector{Float64}}` : isotope atomic fractions for each element, if specified.

# Optional field(s) - with default values
- `state_of_matter::String = "solid"` : state of the matter.

"""
mutable struct Material

    # Variable(s)
    tag                 ::String
    density             ::Union{Missing,Float64}
    state_of_matter     ::String
    number_of_elements  ::Int64
    elements            ::Vector{String}
    atomic_numbers      ::Vector{Int64}
    weight_fractions    ::Vector{Float64}
    mean_excitation_energy   ::Union{Missing,Float64}
    stopping_power_files     ::Dict{Type,String}
    mass_numbers             ::Vector{Vector{Int64}}
    isotope_atomic_fractions ::Vector{Vector{Float64}}

    # Constructor(s)
    function Material(tag::String)
        this = new()
        this.tag                = tag
        this.density            = missing
        this.state_of_matter    = "solid"
        this.number_of_elements = 0
        this.elements           = Vector{String}()
        this.atomic_numbers     = Vector{Int64}()
        this.weight_fractions   = Vector{Float64}()
        this.mean_excitation_energy      = missing
        this.mass_numbers                = Vector{Vector{Int64}}()
        this.isotope_atomic_fractions    = Vector{Vector{Float64}}()
        this.stopping_power_files        = Dict{Type,String}()
        return this
    end
end

# Method(s)
"""
    println(this::Material)

To print the material properties.

# Input Argument(s)
- `this::Material` : material.

# Output Argument(s)
N/A

"""
function println(this::Material)
    println("Material:")
    println("   Density (g/cm³):              $(this.density)")
    println("   Elements in the compound:     $(this.elements)")
    println("   Weight fractions:             $(this.weight_fractions)")
    if any(!isempty.(this.mass_numbers))
        println("   Isotope mass numbers:          $(this.mass_numbers)")
        println("   Isotope atomic fractions:      $(this.isotope_atomic_fractions)")
    end
    if !isempty(this.stopping_power_files)
        println("   Stopping-power files:         $(this.stopping_power_files)")
    end
    println("   State of matter:              $(this.state_of_matter)")
end

"""
    set_density(this::Material,density::Real)

To set the density of the material.

# Input Argument(s)
- `this::Material` : material.
- `density::Real` : density [in g/cm³].

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> mat.set_density(19.3)
```
"""
function set_density(this::Material,density::Real)
    if density < 0 error("The density should be positive.") end
    this.density = density
end

"""
    set_mean_excitation_energy(this::Material,I::Real)

To override the mean excitation energy I of the material (in eV). When set, this value
takes precedence over the internally-computed one (canonical compound table or Bragg
additivity rule) for the collisional stopping power and density effect. Use it to match
the ICRU-37/ESTAR/TOPAS value of a given compound (e.g. striated muscle ≈ 74.7 eV).

# Input Argument(s)
- `this::Material` : material.
- `I::Real` : mean excitation energy [in eV].

# Examples
```jldoctest
julia> muscle.set_mean_excitation_energy(74.7)
```
"""
function set_mean_excitation_energy(this::Material,I::Real)
    if I ≤ 0 error("The mean excitation energy should be positive.") end
    this.mean_excitation_energy = float(I)
end

"""
    get_mean_excitation_energy(this::Material)

To get the user-defined mean excitation energy [in eV], or `missing` if it was not set.

# Input Argument(s)
- `this::Material` : material.

# Output Argument(s)
- `I::Union{Missing,Float64}` : mean excitation energy [in eV] or `missing`.
"""
function get_mean_excitation_energy(this::Material)
    return this.mean_excitation_energy
end

"""
    set_state_of_matter(this::Material,state_of_matter::String)

To set the state of the matter (solid or liquid or gaz).

# Input Argument(s)
- `this::Material` : material.
- `state_of_matter::String` : state of matter, which value is given by "solid", "liquid" or "gaz".

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> mat.set_state_of_matter("liquid")
```
"""
function set_state_of_matter(this::Material,state_of_matter::String)
    if state_of_matter ∉ ["solid","liquid","gaz"] error("The state of matter is either solid, liquid or gaz.") end
    this.state_of_matter = state_of_matter
end

"""
    add_element(this::Material,symbol::String,weight_fraction::Real=1)

To add an element which is part of the composition of the material.

# Input Argument(s)
- `this::Material` : material.
- `symbol::String` : element symbol, with value such as "H", "He", etc.
- `weight_fraction::Real` : weight fraction of the element in the material (between 0 and 1).

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> water.add_element("H",0.1111)
julia> water.add_element("O",0.8889)
```
"""
function add_element(this::Material,symbol::String,weight_fraction::Real=1)
    if ~(0 ≤ weight_fraction ≤ 1) error("Weight fraction should have values between 0 and 1.") end
    push!(this.elements,lowercase(symbol))
    push!(this.atomic_numbers,atomic_number(lowercase(symbol)))
    push!(this.weight_fractions,weight_fraction)
    push!(this.mass_numbers, Vector{Int64}())
    push!(this.isotope_atomic_fractions, Vector{Float64}())
    this.number_of_elements += 1
    if sum(this.weight_fractions) > 1 error("Weight fraction exceed 1.") end
    if weight_fraction == 1 && ismissing(this.density)
        this.set_density(density(atomic_number(lowercase(symbol))))
    end
end

"""
    add_element_isotopes(this::Material,symbol::String,weight_fraction::Real,
    mass_numbers::Vector{Int64},isotope_atomic_fractions::Vector{Float64})

To add an element and explicitly specify its isotope mass numbers and isotope atomic fractions.

# Input Argument(s)
- `this::Material` : material.
- `symbol::String` : element symbol, with value such as "H", "He", etc.
- `weight_fraction::Real` : weight fraction of the element in the material (between 0 and 1).
- `mass_numbers::Vector{Int64}` : isotope mass numbers for the element.
- `isotope_atomic_fractions::Vector{Float64}` : atomic fractions for the isotopes (should sum to 1).

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> water.add_element_isotopes("H",0.1111,[1,2],[0.999885,0.000115])
```
"""
function add_element_isotopes(this::Material,symbol::String,weight_fraction::Real,mass_numbers::Vector{Int64},isotope_atomic_fractions::Vector{Float64})
    if ~(0 ≤ weight_fraction ≤ 1) error("Weight fraction should have values between 0 and 1.") end
    if length(mass_numbers) != length(isotope_atomic_fractions)
        error("Mass numbers and isotope atomic fractions must have the same length.")
    end
    if sum(isotope_atomic_fractions) > 1 + 1e-12 || sum(isotope_atomic_fractions) < 1 - 1e-12
        error("Isotope atomic fractions must sum to 1.")
    end
    push!(this.elements,lowercase(symbol))
    push!(this.atomic_numbers,atomic_number(lowercase(symbol)))
    push!(this.weight_fractions,weight_fraction)
    push!(this.mass_numbers, copy(mass_numbers))
    push!(this.isotope_atomic_fractions, copy(isotope_atomic_fractions))
    this.number_of_elements += 1
    if sum(this.weight_fractions) > 1 error("Weight fraction exceed 1.") end
    if weight_fraction == 1 && ismissing(this.density)
        this.set_density(density(atomic_number(lowercase(symbol))))
    end
end

"""
    set_element_isotopes(this::Material,symbol::String,
    mass_numbers::Vector{Int64},isotope_atomic_fractions::Vector{Float64})

To set isotope mass numbers and isotope atomic fractions for an existing element in the material.

# Input Argument(s)
- `this::Material` : material.
- `symbol::String` : element symbol, with value such as "H", "He", etc.
- `mass_numbers::Vector{Int64}` : isotope mass numbers for the element.
- `isotope_atomic_fractions::Vector{Float64}` : atomic fractions for the isotopes (should sum to 1).

# Output Argument(s)
N/A
"""
function set_element_isotopes(this::Material,symbol::String,mass_numbers::Vector{Int64},isotope_atomic_fractions::Vector{Float64})
    if length(mass_numbers) != length(isotope_atomic_fractions)
        error("Mass numbers and isotope atomic fractions must have the same length.")
    end
    if sum(isotope_atomic_fractions) > 1 + 1e-12 || sum(isotope_atomic_fractions) < 1 - 1e-12
        error("Isotope atomic fractions must sum to 1.")
    end
    idx = findfirst(==(lowercase(symbol)), this.elements)
    if idx === nothing
        error("Element $(symbol) is not defined in the material.")
    end
    this.mass_numbers[idx] = copy(mass_numbers)
    this.isotope_atomic_fractions[idx] = copy(isotope_atomic_fractions)
end

"""
    get_density(this::Material)

To get the density of the material.

# Input Argument(s)
- `this::Material` : material.

# Output Argument(s)
- `density::Real` : density [in g/cm³].

# Examples
```jldoctest
julia> density = mat.get_density()
```
"""
function get_density(this::Material)
    if ismissing(this.density) error("Density is missing.") end
    return this.density
end

"""
    get_number_of_elements(this::Material)

To get the number of elements in the material composition.

# Input Argument(s)
- `this::Material` : material.

# Output Argument(s)
- `number_of_elements::Int64` : number of elements.

# Examples
```jldoctest
julia> N = mat.get_number_of_elements()
```
"""
function get_number_of_elements(this::Material)
    return this.number_of_elements
end

"""
    get_atomic_numbers(this::Material)

To get the atomic numbers of the elements in the material composition.

# Input Argument(s)
- `this::Material` : material.

# Output Argument(s)
- `atomic_numbers::Vector{Int64}` : atomic numbers.

# Examples
```jldoctest
julia> atomic_number = mat.get_atomic_numbers()
```
"""
function get_atomic_numbers(this::Material)
    return this.atomic_numbers
end

"""
    get_weight_fractions(this::Material)

To get the weight fractions associated with the elements in the material composition.

# Input Argument(s)
- `this::Material` : material.

# Output Argument(s)
- `weight_fractions::Vector{Float64}` : weight fractions.

# Examples
```jldoctest
julia> weight_fractions = mat.weight_fractions()
```
"""
function get_weight_fractions(this::Material)
    return this.weight_fractions
end

"""
    get_mass_numbers(this::Material)

To get the isotope mass numbers for each element in the material composition.

# Input Argument(s)
- `this::Material` : material.

# Output Argument(s)
- `mass_numbers::Vector{Vector{Int64}}` : isotope mass numbers for each element.

"""
function get_mass_numbers(this::Material)
    return this.mass_numbers
end

"""
    get_isotope_atomic_fractions(this::Material)

To get the isotope atomic fractions for each element in the material composition.

# Input Argument(s)
- `this::Material` : material.

# Output Argument(s)
- `isotope_atomic_fractions::Vector{Vector{Float64}}` : isotope atomic fractions for each element.

"""
function get_isotope_atomic_fractions(this::Material)
    return this.isotope_atomic_fractions
end

"""
    get_tag(this::Material)

Get the tag (string identifier) of the material.

# Input Argument(s)
- `this::Material` : material.

# Output Argument(s)
- `tag::String` : tag of the material.

# Examples
```jldoctest
julia> tag = mat.get_tag()
```
"""
function get_tag(this::Material)
    return this.tag
end

"""
    set_tag(this::Material,tag::String)

Set the tag (string identifier) of the material.

# Input Argument(s)
- `this::Material` : material.
- `tag::String` : tag.

# Examples
```jldoctest
julia> mat = Material()
julia> mat.set_tag("water")
```
"""
function set_tag(this::Material,tag::String)
    if isempty(tag) error("Material tag cannot be empty.") end
    this.tag = tag
end

"""
    get_state_of_matter(this::Material)

To get the state of matter of the material.

# Input Argument(s)
- `this::Material` : material.

# Output Argument(s)
- `state_of_matter::String` : state of matter.

# Examples
```jldoctest
julia> state_of_matter = mat.get_state_of_matter()
```
"""
function get_state_of_matter(this::Material)
    return this.state_of_matter
end

"""
    set_stopping_power(this::Material,particle::Type,file::String)

To read the electronic stopping power of this material, for one particle, out of a
two-column text file instead of computing it. On this path `bethe` is not called for the
electronic term and neither is the SRIM lookup: the file *is* the stopping power.

This is the escape hatch. A user with a measured table, a table from a code Radiant does
not implement, or a table their group has agreed on can use it without the answer passing
through any model. It lives on `Material` rather than on `Inelastic_Collision` because the
table describes one material and nothing about it generalises to the next one in the
geometry.

**It is also the one path in Radiant's stopping-power chain with no validation behind it.**
Everything else was measured against NIST or SRIM; a user table is measured against
nothing, because its correctness is the user's to establish. Radiant checks that the file
is well-formed and then does exactly what it says.

# File format

Plain text, two columns, conventionally placed in `data/`:

```
# proton in liquid water, PSTAR, retrieved 2026-01-14
# E [MeV]      S/rho [MeV cm^2/g]
1.000e-03      1.769e+02
1.500e-03      1.984e+02
...
1.000e+04      1.848e+00
```

Column 1 is the projectile's **total kinetic energy in MeV** - not energy per nucleon.
Column 2 is the **mass** stopping power in **MeV cm²/g**. The separator is whitespace or a
comma, blank lines are skipped, and a line starting with `#` or `%` is a comment. Energies
must be strictly increasing and are never sorted for you. A third column is an error, not
something ignored: in SRIM's output that column is *nuclear* stopping, and adding it in
silently would be wrong by percent at the low end.

The file is parsed and validated here, at the call site where the user is standing, rather
than at the first call to `sp` three layers away. An absolute path is used as given;
anything else is taken relative to `data/`.

# What the table must be

A **total electronic** stopping power: not restricted, not soft-only, and with no nuclear
component. `sp` computes `Sr = Stot - Sc` and subtracts its own catastrophic part from
whatever it is given, so a restricted table subtracts that part twice. The soft/catastrophic
split is therefore approximate on this path, in the same way it is for the `icru49`
low-energy branch: `Stot` comes from the user's source and `Sc` from Radiant's analytic
cross-section, and the two were never reconciled. Radiant cannot check any of this.

# What it applies to, and what it overrides

**One table, one particle.** A table registered for `Proton` changes nothing for any other
particle: electrons in the same material still go through `bethe`, and alphas still go
through `bethe` or through SRIM. Radiant runs coupled multi-particle problems, and a single
un-keyed table would quietly apply a proton table to electrons the moment someone added an
electron group. A user with tables for two particles calls this twice.

It **supersedes everything** for that one particle in that one material -
`set_use_srim_sp`, `set_ecut`, `set_emax`, `set_low_energy_model`,
`set_low_energy_source`, `set_is_bragg_additivity`, the Sternheimer, atom-fraction and
Barkas corrections, and `set_mean_excitation_energy` - all of which stay set and stop
having an effect. None of that is an error: handing in a table is switching sources.

It does **not** replace the radiative or nuclear stopping power. `Bremsstrahlung.sp` and
`Elastic_Scattering.sp` still run and still contribute, so a Radiant total will differ from
the table at high energy.

# Density, which is where this will go wrong

The table is scaled by **this material's** density. For everything except the density effect
that is correct, because mass stopping power is density-independent - but the density effect
is baked into whatever numbers were supplied, at whatever density their source assumed, and
rescaling cannot touch it. A water-vapour table handed to a liquid-water material gives
vapour physics scaled to liquid density, silently.

# Range

The table must span the whole energy range of the group structure for that particle. An
energy outside it is an error, not an extrapolation and not a fallback to `bethe`.

# Input Argument(s)
- `this::Material` : material.
- `particle::Type` : particle type the table applies to.
- `file::String` : file name, taken relative to `data/`, or an absolute path.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> mat = Material("water")
julia> mat.set_stopping_power(Proton,"pstar_water_proton.txt")
```
"""
function set_stopping_power(this::Material,particle::Type,file::String)
    if particle ∉ [Photon,Electron,Positron,Proton,Antiproton,Alpha,Muon,Antimuon]
        error("Unknown '$particle' particle. A stopping-power table can be registered for " *
              "Photon, Electron, Positron, Proton, Antiproton, Alpha, Muon or Antimuon.")
    end
    path = resolve_user_stopping_power_file(file)
    user_stopping_power_table(path)   # parse and validate now, not at the first call to sp
    this.stopping_power_files[particle] = path
end

"""
    get_stopping_power(this::Material,particle::Type)
    get_stopping_power(this::Material)

To get the resolved path of the stopping-power table registered for a particle, or
`missing` if none is; or, with no particle, the whole dictionary.

# Input Argument(s)
- `this::Material` : material.
- `particle::Type` : particle type.

# Output Argument(s)
- `file::Union{Missing,String}` : resolved path, or `missing`.

# Examples
```jldoctest
julia> mat.get_stopping_power(Proton)
```
"""
function get_stopping_power(this::Material,particle::Type)
    return get(this.stopping_power_files,particle,missing)
end

function get_stopping_power(this::Material)
    return this.stopping_power_files
end