"""
    write_matxs(cross_sections::Cross_Sections,file::String,binary::Bool=false)

Write a MATXS-format (NJOY/CCCC) cross-sections file from a Cross_Sections structure.

The file follows the CCCC MATXS record structure produced by NJOY's `matxsr` module
(file identification → file control → set hollerith id → file data → group structures →
material control → vector control/blocks → matrix control/sub-blocks). To round-trip the
full Radiant data set without loss, standard and Radiant-specific reactions are written:

- standard vectors use the recognized total (`*tot0`), heating (`*heat`), charge (`*char`),
  and scattering (`*scat`) names;
- Radiant-specific vectors use `*abs` and `*edep`, plus `EMOMTR`/`BSTC`/`PSTC` for charged
  particles. These vectors carry the data required to reconstruct the exact
  `Cross_Sections`; tools that do not recognize them may consume only the standard subset.

The emitted vector names are `gtot0`, `gheat`, `gchar`, `gabs`, `gedep` for photons;
`btot0`, `bheat`, `bchar`, `babs`, `EMOMTR`, `BSTC`, `bedep` for electrons; and
`ptot0`, `pheat`, `pchar`, `pabs`, `PSTC`, `pedep` for positrons. Legacy Radiant-private
names are accepted by `read_matxs` but are not emitted.

Radiant stores energy boundaries internally in MeV and writes them in MATXS eV. Vector and
matrix values are divided by material density when written; `read_matxs` multiplies them by
the supplied material density when rebuilding the library.

Two encodings are supported and share the same record content: the BCD (ASCII) encoding
(`binary=false`, tagged text records with `e12.5`/`i6`/`a8` fields) and a binary encoding
(`binary=true`) writing each record as one native FORTRAN-style unformatted record (4-byte
length markers around `Int32`/`Float64`/8-byte hollerith payloads). The binary encoding is a
Radiant variant and is not guaranteed byte-compatible with NJOY's own binary MATXS files.

# Input Argument(s)
- `cross_sections::Cross_Sections` : cross-sections library.
- `file::String` : file name and directory.
- `binary::Bool` : write the binary encoding instead of BCD/ASCII.

# Output Argument(s)
N/A

# Reference(s)
N/A

"""
function write_matxs(cross_sections::Cross_Sections,file::String,binary::Bool=false)

#----
# Extract information
#----
particles = cross_sections.get_particles()
materials = cross_sections.get_materials()
Npart = cross_sections.get_number_of_particles()
Nmat = cross_sections.get_number_of_materials()
L = cross_sections.get_legendre_order()
lord = L + 1
densities = [materials[m].density for m in range(1,Nmat)]
if any(ismissing.(densities)) || any(densities .<= 0.0)
    error("MATXS output requires strictly positive material densities.")
end
densities = Float64.(densities)

Ng = [cross_sections.get_number_of_groups(particles[i]) for i in range(1,Npart)]
codes = [matxs_particle_code(particles[i]) for i in range(1,Npart)]

# Vector reactions written for each incident particle (fixed order, shared by read_matxs).
vector_names_by_particle = [matxs_vector_names(codes[i]) for i in range(1,Npart)]

#----
# File-level data type table : one type per ordered (incident,outgoing) particle pair.
# Self pairs (jinp==joutp) additionally carry the vector reactions of the incident particle.
#----
ntype,htype,jinp,joutp = matxs_type_table(Npart,codes)
n1d_by_type = [(jinp[t] == joutp[t]) ? length(vector_names_by_particle[jinp[t]]) : 0 for t in range(1,ntype)]
n2d_by_type = ones(Int64,ntype)
type_record_counts = [(n1d_by_type[t] > 0 ? 2 : 0) + 2*n2d_by_type[t] for t in range(1,ntype)]
submaterial_locs = vcat([0],cumsum(type_record_counts)[1:end-1])
material_record_count = 1 + sum(type_record_counts)

#----
# Write the MATXS file
#----
io = open(file, binary ? "w" : "w")
try

    # Record 0 : file identification
    matxs_emit(io,binary," 0v ",[(:holl,["matxs","Radiant"]),(:int,[5])])

    # Record 1 : file control
    nholl = 9
    maxw  = 99999     # fits in i6; blocks are written un-split and read by computed count
    flength = 0
    matxs_emit(io,binary," 1d ",[(:int,[Npart,ntype,nholl,Nmat,maxw,flength])])

    # Record 2 : set hollerith identification (nholl a8 words)
    idwords = matxs_holl_words("Coupled multigroup cross-sections produced with Radiant.jl",nholl)
    matxs_emit(io,binary," 2d ",[(:holl,idwords)])

    # Record 3 : file data
    hprt  = [codes[i] for i in range(1,Npart)]
    hmatn = [string(first(materials[m].tag,8)) for m in range(1,Nmat)]
    nsubm = fill(ntype,Nmat)
    locm  = [(m-1)*material_record_count for m in range(1,Nmat)]
    matxs_emit(io,binary," 3d ",[(:holl,vcat(hprt,htype,hmatn)),
                                 (:int,vcat(Ng,jinp,joutp,nsubm,locm))])

    # Record 4 : group structures (per particle, descending bounds + emin).
    # MATXS stores energy boundaries in eV; Radiant keeps its internal energy grid in MeV.
    for i in range(1,Npart)
        eb = 1.0e6 .* cross_sections.get_energy_boundaries(particles[i])
        matxs_emit(io,binary," 4d ",[(:real,collect(Float64,eb))])
    end

    # Pre-fetch per-particle / per-pair data
    total = [cross_sections.get_total(particles[i]) for i in range(1,Npart)]
    edep  = [cross_sections.get_energy_deposition(particles[i]) for i in range(1,Npart)]
    absor = [cross_sections.get_absorption(particles[i]) for i in range(1,Npart)]
    momt  = [cross_sections.get_momentum_transfer(particles[i]) for i in range(1,Npart)]
    stpw  = [cross_sections.get_boundary_stopping_powers(particles[i]) for i in range(1,Npart)]
    cdep  = [cross_sections.get_charge_deposition(particles[i]) for i in range(1,Npart)]
    scat  = Array{Array{Float64,4}}(undef,Npart,Npart)
    for i in range(1,Npart), o in range(1,Npart)
        scat[i,o] = cross_sections.get_scattering(particles[i],particles[o],L) # (Nmat,Ngi,Ngo,L+1)
    end

    # Records 5-9 : per material
    for m in range(1,Nmat)

        # Record 5 : material control (hmat + amass + per-submaterial temp,sigz,itype,n1d,n2d,locs)
        reals5 = Float64[0.0]   # amass
        ints5  = Int64[]
        for t in range(1,ntype)
            n1d = n1d_by_type[t]
            n2d = n2d_by_type[t]
            append!(reals5,[0.0,0.0])              # temp, sigz
            append!(ints5,[t,n1d,n2d,submaterial_locs[t]]) # itype, n1d, n2d, locs
        end
        matxs_emit(io,binary," 5d ",[(:holl,[string(first(materials[m].tag,8))]),
                                     (:real,reals5),(:int,ints5)])

        # Per submaterial (== per data type)
        for t in range(1,ntype)
            i = jinp[t]; o = joutp[t]
            Ngi = Ng[i]; Ngo = Ng[o]

            # Vector control + block (self pairs only)
            if i == o
                hvps = vector_names_by_particle[i]
                nfg  = ones(Int64,length(hvps))
                nlg  = [name in ("BSTC","CSTC","PSTC") || endswith(name,"edep") ? Ngi+1 : Ngi for name in hvps]
                matxs_emit(io,binary," 6d ",[(:holl,hvps),(:int,vcat(nfg,nlg))])

                vps = Float64[]
                for name in hvps
                    if endswith(name,"tot0")
                        append!(vps,total[i][:,m] ./ densities[m])
                    elseif endswith(name,"heat")
                        append!(vps,edep[i][1:Ngi,m] ./ densities[m])
                    elseif endswith(name,"char")
                        append!(vps,cdep[i][1:Ngi,m] ./ densities[m])
                    elseif endswith(name,"abs")
                        append!(vps,absor[i][:,m] ./ densities[m])
                    elseif name == "EMOMTR"
                        append!(vps,momt[i][:,m] ./ densities[m])
                    elseif name in ("BSTC","CSTC","PSTC")
                        append!(vps,stpw[i][:,m] ./ densities[m])
                    elseif endswith(name,"edep")
                        append!(vps,edep[i][:,m] ./ densities[m])
                    else
                        error("Unsupported MATXS vector name '$name'.")
                    end
                end
                matxs_emit(io,binary," 7d ",[(:real,vps)])
            end

            # Matrix control + sub-block (every pair)
            matview = @view scat[i,o][m,:,:,:]   # (Ngi,Ngo,L+1)
            jband,ijj = matxs_band(matview,Ngi,Ngo)
            matxs_emit(io,binary," 8d ",[(:holl,[codes[i]*"scat"]),
                                         (:int,vcat([lord,0],jband,ijj))])

            mdat = Float64[]
            for gf in range(1,Ngo)
                jb = jband[gf]; ij = ijj[gf]
                if jb == 0 continue end
                for l in range(1,lord), gi in ij:-1:(ij-jb+1)
                    push!(mdat,matview[gi,gf,l] / densities[m])
                end
            end
            matxs_emit(io,binary," 9d ",[(:real,mdat)])
        end
    end

finally
    close(io)
end

end

"""
    matxs_type_table(Npart::Int64,codes::Vector{String})

Build the MATXS data-type table : one ordered (incident,outgoing) particle pair per type.
Returns `(ntype,htype,jinp,joutp)`.

"""
function matxs_type_table(Npart::Int64,codes::Vector{String})
    ntype = Npart*Npart
    htype = Vector{String}(undef,ntype)
    jinp  = Vector{Int64}(undef,ntype)
    joutp = Vector{Int64}(undef,ntype)
    for i in range(1,Npart), o in range(1,Npart)
        t = (i-1)*Npart + o
        htype[t] = codes[i]*codes[o]
        jinp[t]  = i
        joutp[t] = o
    end
    return ntype,htype,jinp,joutp
end

"""
    matxs_band(matview,Ngi::Int64,Ngo::Int64)

Compute the MATXS banding arrays `(jband,ijj)` of a `(Ngi,Ngo,L+1)` scattering matrix view :
for each final group, `jband` is the number of contributing initial groups and `ijj` is
the highest initial-group number in the band.

"""
function matxs_band(matview,Ngi::Int64,Ngo::Int64)
    jband = zeros(Int64,Ngo)
    ijj   = ones(Int64,Ngo)
    for gf in range(1,Ngo)
        gmin = 0; gmax = 0
        for gi in range(1,Ngi)
            if any(matview[gi,gf,:] .!= 0.0)
                if gmin == 0 gmin = gi end
                gmax = gi
            end
        end
        if gmin != 0
            jband[gf] = gmax - gmin + 1
            ijj[gf]   = gmax
        end
    end
    return jband,ijj
end

"""
    matxs_particle_code(particle::Particle)

Return the single-character MATXS hollerith particle code used by the Radiant MATXS
reader/writer (`"g"` photon, `"b"` electron/beta, `"p"` positron).

"""
function matxs_particle_code(particle::Particle)
    if is_photon(particle)   return "g" end
    if is_electron(particle) return "b" end
    if is_positron(particle) return "p" end
    error("Unsupported particle for MATXS format (expected photon, electron or positron).")
end

"""
    matxs_vector_names(code::String)

Return Dragon-facing MATXS vector names for an incident particle code. Records 35 and 36
of FMAC-M are charged-particle quantities; keep them out of photon data.

"""
function matxs_vector_names(code::String)
    if code == "b"
        return ["btot0","bheat","bchar","babs","EMOMTR","BSTC","bedep"]
    elseif code == "g"
        return ["gtot0","gheat","gchar","gabs","gedep"]
    elseif code == "p"
        return ["ptot0","pheat","pchar","pabs","PSTC","pedep"]
    end
    error("Unsupported MATXS particle code '$code'.")
end

"""
    matxs_holl_words(str::String,nwords::Int64)

Split a string into `nwords` fixed `a8` hollerith words (right-padded with spaces).

"""
function matxs_holl_words(str::String,nwords::Int64)
    s = rpad(str,nwords*8)
    return [s[(i-1)*8+1:i*8] for i in 1:nwords]
end

"""
    matxs_emit(io::IO,binary::Bool,tag::String,segments)

Write one MATXS record. `segments` is an ordered list of `(:int|:real|:holl, values)` pairs.
In BCD mode the record tag and first payload fields are written on the same
fixed-format line, following the CCCC/Dragon text reader conventions; in binary mode it
is written as a single FORTRAN-style unformatted record (`Int32` length marker, payload of `Int32`/`Float64`/8-byte hollerith, length marker).

"""
function matxs_emit(io::IO,binary::Bool,tag::String,segments)
    if binary
        # Note : Radiant defines its own `write` (Cross_Sections.write), so byte output must
        # explicitly use `Base.write`.
        buf = IOBuffer()
        for (kind,vals) in segments
            if kind === :int
                for v in vals Base.write(buf,Int32(v)) end
            elseif kind === :real
                for v in vals Base.write(buf,Float64(v)) end
            elseif kind === :holl
                for s in vals Base.write(buf,rpad(first(s,8),8)) end
            end
        end
        bytes = take!(buf)
        marker = Int32(length(bytes))
        Base.write(io,marker); Base.write(io,bytes); Base.write(io,marker)
    else
        _matxs_emit_bcd(io,tag,segments)
    end
end

"""
    _matxs_emit_bcd(io::IO,tag::String,segments)

Write one BCD MATXS record using the fixed-column layout expected by the CCCC/Dragon
reader.  The record tag occupies the first columns of the first payload line; continuation
lines are untagged.

"""
function _matxs_emit_bcd(io::IO,tag::String,segments)
    stag = rpad(tag,4)
    if tag == " 0v "
        holl = String.(segments[1][2])
        ints = Int.(segments[2][2])
        println(io,stag,rpad(first(holl[1],8),8)," ",rpad(first(holl[2],8),8),
                rpad("",8)," ",@sprintf("%6d",ints[1]))
    elseif tag == " 1d "
        print(io,stag,"  ")
        _matxs_write_ints(io,segments[1][2]; first_line=6, cont_line=12)
    elseif tag == " 2d "
        print(io,stag)
        _matxs_write_holl(io,segments[1][2]; first_prefix="    ", first_line=8,
                          cont_line=9)
    elseif tag == " 3d "
        print(io,stag)
        _matxs_write_holl(io,segments[1][2]; first_prefix="    ", first_line=8,
                          cont_line=9)
        _matxs_write_ints(io,segments[2][2])
    elseif tag == " 4d " || tag == " 7d " || tag == " 9d " || tag == "10d "
        print(io,stag)
        _matxs_write_reals(io,segments[1][2]; first_prefix="        ", first_line=5,
                           cont_line=6)
    elseif tag == " 5d "
        holl = String.(segments[1][2])
        reals = Float64.(segments[2][2])
        ints = Int.(segments[3][2])
        println(io,stag,rpad(first(holl[1],8),8),@sprintf("%12.5E",reals[1]))
        nsub = length(ints) ÷ 4
        for i in 1:nsub
            rbase = 2 + 2*(i-1)
            ibase = 1 + 4*(i-1)
            println(io,@sprintf("%12.5E",reals[rbase]),@sprintf("%12.5E",reals[rbase+1]),
                    @sprintf("%6d",ints[ibase]),@sprintf("%6d",ints[ibase+1]),
                    @sprintf("%6d",ints[ibase+2]),@sprintf("%6d",ints[ibase+3]))
        end
    elseif tag == " 6d "
        print(io,stag)
        _matxs_write_holl(io,segments[1][2]; first_prefix="    ", first_line=8,
                          cont_line=9)
        _matxs_write_ints(io,segments[2][2])
    elseif tag == " 8d "
        holl = String.(segments[1][2])
        ints = Int.(segments[2][2])
        println(io,stag,"    ",rpad(first(holl[1],8),8))
        _matxs_write_ints(io,ints)
    else
        error("Unsupported MATXS BCD record tag '$tag'.")
    end
end

"""
    _matxs_write_reals(io::IO,vals)

Write a sequence of reals as fixed-width `e12.5` fields, 6 per line.

"""
function _matxs_write_reals(io::IO,vals; first_prefix::String="", first_line::Int64=6,
                            cont_line::Int64=6)
    fields = [@sprintf("%12.5E",Float64(v)) for v in vals]
    i = 1
    nfirst = min(first_line,length(fields))
    println(io,first_prefix*join(fields[i:i+nfirst-1]))
    i += nfirst
    while i <= length(fields)
        n = min(cont_line,length(fields)-i+1)
        println(io,join(fields[i:i+n-1]))
        i += n
    end
end

"""
    _matxs_write_ints(io::IO,vals)

Write a sequence of integers as fixed-width `i6` fields, 12 per line.

"""
function _matxs_write_ints(io::IO,vals; first_line::Int64=12, cont_line::Int64=12)
    fields = [@sprintf("%6d",Int(v)) for v in vals]
    i = 1
    nfirst = min(first_line,length(fields))
    println(io,join(fields[i:i+nfirst-1]))
    i += nfirst
    while i <= length(fields)
        n = min(cont_line,length(fields)-i+1)
        println(io,join(fields[i:i+n-1]))
        i += n
    end
end

"""
    _matxs_write_holl(io::IO,names)

Write a sequence of hollerith identifiers as fixed-width `a8` fields, 8 per line.

"""
function _matxs_write_holl(io::IO,names; first_prefix::String="", first_line::Int64=8,
                           cont_line::Int64=8)
    fields = [rpad(first(n,8),8) for n in names]
    i = 1
    nfirst = min(first_line,length(fields))
    println(io,first_prefix*join(fields[i:i+nfirst-1]))
    i += nfirst
    while i <= length(fields)
        n = min(cont_line,length(fields)-i+1)
        println(io,join(fields[i:i+n-1]))
        i += n
    end
end
