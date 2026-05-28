# The MATXS File Format

MATXS is a generalized, flexible multigroup interface file based on the CCCC (Committee on
Computer Code Coordination) standard conventions and produced by NJOY's `matxsr` module. It
stores cross-section vectors and group-to-group transfer matrices for an arbitrary set of
particles, materials and reactions.

Radiant reads and writes a MATXS file through the `"matxs"` cross-section source and the
`cs.write(path,"matxs")` method (see Section 4.5 of the User's Guide). This page describes the
record structure, the reaction-name keywords and the particle codes used by Radiant, so that
the files can be inspected or exchanged with other tools.

!!! note "Interoperability and round-trip"
    Radiant writes two layers of reactions. A *standard* interop layer uses recognized MATXS
    reaction names (`*tot0`, `*heat`, `*scat`) so that NJOY/TRANSX-style tooling can read the
    main quantities. An *extension* layer (`*edep`, `*abs`, `*momt`, `*stpw`, `*cdep`) carries
    the Radiant-specific quantities that have no standard MATXS slot; tools that do not
    recognize these names simply skip them. Reading a Radiant-written MATXS file back therefore
    reproduces the original `Cross_Sections` exactly.

## File structure

A MATXS file is a sequence of records. The first four records describe the whole file; the
group structures are then given once per particle; and the data of each material follows as a
material-control record and, per submaterial, a set of vector and matrix records.

```
file identification        (always)
file control               (always)
set hollerith id           (always)
file data                  (always)

  (repeat per particle)
    group structure         (always)

  (repeat per material)
    material control        (always)
      (repeat per submaterial)
        vector control      (if the submaterial holds vector reactions)
        vector block(s)
        matrix control      (if the submaterial holds a matrix reaction)
        matrix sub-block(s)
```

In the BCD (ASCII) encoding each record begins with a 4-character tag. The records written by
Radiant are:

| Tag    | Record                       | Content                                                        |
|--------|------------------------------|----------------------------------------------------------------|
| ` 0v ` | File identification          | File name, user id, version.                                   |
| ` 1d ` | File control                 | `npart, ntype, nholl, nmat, maxw, length`.                     |
| ` 2d ` | Set hollerith identification | Free-text description of the set.                              |
| ` 3d ` | File data                    | Particle / data-type / material names, group counts, `jinp`, `joutp`, submaterial counts. |
| ` 4d ` | Group structure              | Group upper bounds in descending order followed by the minimum energy (`Ng+1` values), one record per particle. |
| ` 5d ` | Material control             | Material name and, per submaterial, `temp, sigz, itype, n1d, n2d, locs`. |
| ` 6d ` | Vector control               | Vector reaction names and their first/last group band (`nfg`, `nlg`). |
| ` 7d ` | Vector block                 | The vector reaction data (`e12.5` reals).                      |
| ` 8d ` | Matrix control               | Matrix reaction name, Legendre count `lord`, `jconst`, and the band arrays `jband`, `ijj`. |
| ` 9d ` | Matrix sub-block             | The banded transfer-matrix data (`e12.5` reals).               |

`npart` is the number of particles, `nmat` the number of materials, and `ntype` the number of
data types. Radiant defines one data type per ordered `(incident, outgoing)` particle pair, so
`ntype = npart²`; the type of submaterial `s` carries the matrix for its pair and, on the
"self" pairs (incident = outgoing), the vector reactions of the incident particle.

## Particle codes

The hollerith particle codes (`hprt`, and the prefix of every reaction name) are:

| Code  | Particle              |
|-------|-----------------------|
| `g`   | photon (gamma)        |
| `b`   | electron (beta)       |
| `p`   | positron              |

Data-type names (`htype`) are the concatenation of the incident and outgoing codes, e.g. `gg`
(photon→photon), `gb` (photon→electron), `bb` (electron→electron).

## Reaction-name keywords

Each reaction name is the incident-particle code followed by a reaction suffix. The vector
reactions (record ` 6d `/` 7d `) are written only on the self pairs; the matrix reaction
(record ` 8d `/` 9d `) is written for every pair. Using the photon prefix `g` as an example:

| Keyword   | Layer     | Length  | Radiant quantity                                       |
|-----------|-----------|---------|--------------------------------------------------------|
| `gtot0`   | standard  | `Ng`    | `get_total` — total cross-section (P₀).                |
| `gheat`   | standard  | `Ng`    | Heating / energy deposition (the per-group part of `get_energy_deposition`). |
| `gscat`   | standard  | banded  | `get_scattering` — group-to-group transfer matrix, Legendre orders `0…L`. |
| `gedep`   | extension | `Ng+1`  | `get_energy_deposition` (authoritative; supersedes `gheat` on read). |
| `gabs`    | extension | `Ng`    | `get_absorption`.                                      |
| `gmomt`   | extension | `Ng`    | `get_momentum_transfer`.                               |
| `gstpw`   | extension | `Ng+1`  | `get_boundary_stopping_powers`.                        |
| `gcdep`   | extension | `Ng+1`  | `get_charge_deposition`.                               |

The same suffixes apply with the `b` (electron) and `p` (positron) prefixes. On reading,
`*heat` is ignored because the full `Ng+1` energy-deposition vector is recovered from `*edep`.

### Scattering-matrix banding

The transfer matrix is stored in the compact banded form used by MATXS. For each outgoing
group, `jband` gives the number of contributing incident groups and `ijj` the lowest incident
group number of the band. Within a band the data are ordered by Legendre moment (`P₀` band,
then `P₁` band, …), and the incident groups within each moment are listed in descending order.

## Encodings

Radiant supports two encodings, selected by the `binary` argument of `cs.write`. The reader
auto-detects the encoding from the first bytes of the file.

- **BCD / ASCII** (`binary = false`, default) — tagged text records. Reals are written as
  fixed-width `e12.5` fields (6 per line), integers as `i6` (12 per line) and hollerith
  identifiers as `a8` (8 per line).
- **Binary** (`binary = true`) — each record is written as a single FORTRAN-style unformatted
  record: a 4-byte `Int32` length marker, a payload of native-endian `Int32` integers,
  `Float64` reals and 8-byte hollerith fields, and a trailing length marker.

!!! warning "Caveats"
    - The binary encoding is a **Radiant variant** (`Int32`/`Float64` payloads); it is not
      guaranteed byte-compatible with NJOY's own binary MATXS files. Validate against a
      reference file if true binary interoperability with NJOY is required.
    - The exact standard reaction-name strings and the positron code are provisional. Confirm
      them against the specific TRANSX/MATXS naming of your toolchain when exchanging files.

## Example

```julia
# Write a built library, then read it back into a fresh Cross_Sections
cs.write("./library.matxs","matxs")        # or "matxs", true for the binary encoding

cs2 = Cross_Sections()
cs2.set_source("matxs")
cs2.set_file("./library.matxs")
cs2.set_materials([water])                  # same materials, in the same order
cs2.set_particles([electron,photon])        # same particles
cs2.build()
```
