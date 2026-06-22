"""
    create_radiant_endf_library.jl

Organizes a directory of raw ENDF-6 files into a Radiant ENDF database under
`data/`. Every file found under `source_path` is probed for an ENDF-6 MF=1,
MT=451 header (regardless of its extension); files without one are skipped.
The target isotope and incident particle are read directly from that header
(no reliance on the source filename), using Radiant's own ENDF line parser.
Recognized files are copied (not modified) into
`data/<library_name>/<particle_tag>/`, renamed to match Radiant's ENDF naming
convention (`Radiant.endf_path_for_isotope`).

Run this script once to prepare an ENDF library before using it in a
simulation (e.g. via `init_elastic_scattering_endf_db`).

# Usage
    julia --project=. create_radiant_endf_library.jl <library_name> <source_path>

# Argument(s)
- `library_name` : name of the ENDF database directory to create under `data/`
  (e.g. `"TENDL2023"`).
- `source_path` : directory to search (recursively) for ENDF-structured files.
"""

using Radiant

const DATA_ROOT = @__DIR__
const HEADER_LINES = 8 # TPID + the four MF=1/MT=451 HEAD/CONT records.

# Pads a line to the 80-column width assumed by Radiant's ENDF parser; some
# raw ENDF files are missing the trailing sequence-number field.
pad80(line::AbstractString) = length(line) >= 80 ? line : line * " "^(80 - length(line))

# Reads just the first few lines of `path` (cheap and safe even for large or
# binary files) and returns its target `(Z, A)` and incident-particle `NSUB`
# code, or `nothing` if `path` has no recognizable ENDF MF=1 section.
function endf_header(path::String)
    lines = String[]
    open(path) do io
        for _ in 1:HEADER_LINES
            eof(io) && break
            push!(lines, pad80(readline(io)))
        end
    end
    i1 = Radiant.find_first_mf(lines, 1)
    i1 === nothing && return nothing
    (ZA, _, _, _, _, _, _, _, _, _, i2, _) = Radiant.read_CONT(lines, i1)
    (_, _, _, _, _, _, _, _, _, _, i3, _) = Radiant.read_CONT(lines, i2)
    (_, _, _, _, nsub, _, _, _, _, _, _, _) = Radiant.read_CONT(lines, i3)
    za_int = round(Int, ZA)
    return (Z = za_int ÷ 1000, A = za_int % 1000, NSUB = round(Int, nsub))
end

# Maps an ENDF MF=1/MT=451 NSUB (sub-library) code to the matching Radiant
# Particle, reusing the projectile ZAP codes already defined in `Radiant.ZAP`.
function particle_from_nsub(NSUB::Int)
    projectile_zap = NSUB ÷ 10
    for particle in (Proton(), Alpha())
        if Radiant.ZAP(particle) == projectile_zap
            return particle
        end
    end
    error("no Radiant Particle matches incident-particle ZAP=$(projectile_zap) (NSUB=$(NSUB))")
end

function organize_endf_library(library_name::String, source_path::String)
    isdir(source_path) || error("Source path not found: $(source_path)")

    n_copied = 0
    n_skipped = 0
    for (root, _, files) in walkdir(source_path)
        for fname in files
            src = joinpath(root, fname)

            header = try
                endf_header(src)
            catch
                nothing
            end
            if header === nothing
                n_skipped += 1
                continue
            end

            particle = try
                particle_from_nsub(header.NSUB)
            catch err
                @warn "Skipping $(fname): $(err)"
                n_skipped += 1
                continue
            end

            dest = Radiant.endf_path_for_isotope(library_name, header.Z, header.A, particle;
                                                  data_root = DATA_ROOT)
            mkpath(dirname(dest))
            cp(src, dest; force = true)
            println("$(fname) -> $(relpath(dest, DATA_ROOT))")
            n_copied += 1
        end
    end
    println("Done: $(n_copied) file(s) copied, $(n_skipped) file(s) skipped.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 2 || error("Usage: julia create_radiant_endf_library.jl <library_name> <source_path>")
    organize_endf_library(ARGS[1], ARGS[2])
end
