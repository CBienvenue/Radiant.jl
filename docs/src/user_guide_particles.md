# 3 Particles

In Radiant, particles are represented by the `Radiant.Particle` object. There are two primary ways to instantiate such an object:

**Using Reserved Particle Constructors:** Reserved constructors are designed for predefined particles with known mass and charge. Each of these generates a `Radiant.Particle` object preloaded with the physical properties specific to that particle type.

```julia
photon     = Photon()
electron   = Electron()
positron   = Positron()
proton     = Proton()
antiproton = Antiproton()
alpha      = Alpha()
muon       = Muon()
antimuon   = Antimuon()
```

Optionally, a string tag can be provided to the constructor:

```julia
e_beam = Electron("primary_electron")
```

Setting an explicit tag is recommended when several `Particle` objects of the same kind are used in a single calculation, since the tag is used to distinguish them.

**Using the General Constructor:** The general constructor allows the creation of user-defined particles. The tag is mandatory; the mass (in MeV) and the charge are optional.

```julia
heavy_ion = Particle("heavy_ion",1.0e3,2.0)   # tag, mass [MeV], charge
```

The default reserved particles use the following properties:

| Particle      | Tag           | Mass [MeV]    | Charge |
|---------------|---------------|---------------|--------|
| `Photon`      | `"photon"`    | 0             | 0      |
| `Electron`    | `"electron"`  | 0.511         | -1     |
| `Positron`    | `"positron"`  | 0.511         | +1     |
| `Proton`      | `"proton"`    | 938.272       | +1     |
| `Antiproton`  | `"antiproton"`| 938.272       | -1     |
| `Alpha`       | `"alpha"`     | 3727.379      | +2     |
| `Muon`        | `"muon"`      | 105.658       | -1     |
| `Antimuon`    | `"antimuon"`  | 105.658       | +1     |

## 3.1 Inspecting a Particle

Getter methods give access to a particle's properties and can be used to test its type:

```julia
get_tag(electron)         # "electron"
get_mass(electron)        # 0.51099895069
get_charge(electron)      # -1
get_type(electron)        # Electron (abstract type)

is_photon(electron)       # false
is_electron(electron)     # true
is_positron(electron)     # false
# Similarly: is_proton, is_antiproton, is_alpha, is_muon, is_antimuon
```

The tag of a particle can be reassigned at any time using `set_tag(particle,new_tag)`.

## 3.2 Reserved vs Custom Particles

The reserved particle types (`Photon`, `Electron`, `Positron`, ...) can only be instantiated through their dedicated constructors. This guarantees that interactions and cross-section models (e.g. Compton scattering for photons, Møller for electrons) recognize the particle and apply the appropriate physics. Particles created with the general `Particle(tag, ...)` constructor have type `Custom_Particle` and cannot trigger the built-in physics models.

## 3.3 Summary of the Particle API

| Method                            | Description                                              |
|-----------------------------------|----------------------------------------------------------|
| `Photon()`, `Electron()`, ...     | Reserved-particle constructors.                          |
| `Photon(tag)`, `Electron(tag)`, ...| Same as above with a custom tag.                        |
| `Particle(tag,mass=missing,charge=missing)` | Custom particle constructor.                   |
| `set_tag(p,tag)`                  | Set the particle tag.                                    |
| `get_tag(p)`                      | Get the particle tag.                                    |
| `get_mass(p)`                     | Get the mass [MeV].                                      |
| `get_charge(p)`                   | Get the charge [in units of the elementary charge].      |
| `get_type(p)`                     | Get the particle type (abstract type).                   |
| `is_photon(p)`, `is_electron(p)`, ... | Test particle identity.                              |
