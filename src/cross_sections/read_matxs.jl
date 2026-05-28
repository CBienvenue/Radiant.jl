"""
    read_matxs(cross_sections::Cross_Sections)

Read a MATXS-format (NJOY/CCCC) cross-sections file and store its content in a Cross_Sections
structure. This is the counterpart of [`write_matxs`](@ref) and reverses both the standard
interop reactions (`*tot0`, `*heat`, `*scat`) and the Radiant-private extension reactions
(`*edep`, `*abs`, `*momt`, `*stpw`, `*cdep`) to rebuild the full multigroup data set.

The BCD (ASCII) and binary encodings produced by `write_matxs` are both supported; the
encoding is auto-detected from the first bytes of the file.

# Input Argument(s)
- `cross_sections::Cross_Sections` : cross-sections library.

# Output Argument(s)
N/A

# Reference(s)
N/A

"""
function read_matxs(cross_sections::Cross_Sections)

#----
# Extract information
#----
file = cross_sections.get_file()
particles = cross_sections.get_particles()
Npart_cs = cross_sections.get_number_of_particles()

# Auto-detect the encoding and build a reader
raw = read(file)
is_bcd = length(raw) >= 4 && raw[1:4] == codeunits(" 0v ")
rd = _Matxs_Reader(is_bcd,raw)

#----
# Record 0 : file identification  +  Record 1 : file control
#----
_matxs_next!(rd," 0v "); _matxs_holl!(rd,2); _matxs_ints!(rd,1)          # id (ignored)
_matxs_next!(rd," 1d ")
ctrl = _matxs_ints!(rd,6)
Npart,ntype,nholl,Nmat = ctrl[1],ctrl[2],ctrl[3],ctrl[4]
if Npart != Npart_cs error("MATXS file holds $Npart particles but $Npart_cs were provided.") end
if Nmat  != cross_sections.get_number_of_materials() error("MATXS file material count mismatch.") end

#----
# Record 2 : set hollerith identification (ignored)
#----
_matxs_next!(rd," 2d "); _matxs_holl!(rd,nholl)

#----
# Record 3 : file data
#----
_matxs_next!(rd," 3d ")
holl = _matxs_holl!(rd,Npart+ntype+Nmat)
hprt = holl[1:Npart]
ints = _matxs_ints!(rd,Npart+2*ntype+2*Nmat)
ngrp_f = ints[1:Npart]
jinp   = ints[Npart+1:Npart+ntype]
joutp  = ints[Npart+ntype+1:Npart+2*ntype]

# Map each file particle slot to the index of the matching provided particle (by type)
slot2cs = Vector{Int64}(undef,Npart)
for j in range(1,Npart)
    code = strip(hprt[j])
    idx = code == "g" ? findfirst(is_photon,particles)   :
          code == "b" ? findfirst(is_electron,particles) :
          code == "p" ? findfirst(is_positron,particles) : nothing
    if isnothing(idx) error("MATXS particle code '$code' has no matching provided particle.") end
    slot2cs[j] = idx
end

#----
# Record 4 : group structures (per particle, file order)
#----
Ng = zeros(Int64,Npart)                          # indexed by provided-particle order
energy_boundaries = Vector{Vector{Float64}}(undef,Npart)
for j in range(1,Npart)
    _matxs_next!(rd," 4d ")
    eb = _matxs_reals!(rd,ngrp_f[j]+1)
    n = slot2cs[j]
    Ng[n] = ngrp_f[j]
    energy_boundaries[n] = eb
end

#----
# Records 5-9 : per material
#----
multigroup_cross_sections = Array{Multigroup_Cross_Sections}(undef,Npart,Nmat)
legendre_order = 0

for mfile in range(1,Nmat)

    # Record 5 : material control
    _matxs_next!(rd," 5d ")
    _matxs_holl!(rd,1)                            # hmat (ignored)
    _matxs_reals!(rd,1+2*ntype)                   # amass + temp/sigz (ignored)
    ints5 = _matxs_ints!(rd,4*ntype)              # itype,n1d,n2d,locs per submaterial
    sub_itype = [ints5[(s-1)*4+1] for s in 1:ntype]
    sub_n1d   = [ints5[(s-1)*4+2] for s in 1:ntype]
    sub_n2d   = [ints5[(s-1)*4+3] for s in 1:ntype]

    vecs = Dict{Int64,Dict{String,Vector{Float64}}}()   # provided incident index -> suffix -> data
    mats = Dict{Tuple{Int64,Int64},Array{Float64,3}}()  # (incident,outgoing) -> scattering

    for s in range(1,ntype)
        t = sub_itype[s]
        i = slot2cs[jinp[t]]; o = slot2cs[joutp[t]]
        Ngi = Ng[i]; Ngo = Ng[o]

        # Vector control + block
        if sub_n1d[s] > 0
            _matxs_next!(rd," 6d ")
            hvps = _matxs_holl!(rd,sub_n1d[s])
            bands = _matxs_ints!(rd,2*sub_n1d[s])
            nfg = bands[1:sub_n1d[s]]
            nlg = bands[sub_n1d[s]+1:end]
            lengths = nlg .- nfg .+ 1
            _matxs_next!(rd," 7d ")
            vps = _matxs_reals!(rd,sum(lengths))
            d = get!(vecs,i,Dict{String,Vector{Float64}}())
            pos = 0
            for k in range(1,sub_n1d[s])
                sfx = strip(hvps[k])[2:end]       # drop the particle-code prefix
                d[sfx] = vps[pos+1:pos+lengths[k]]
                pos += lengths[k]
            end
        end

        # Matrix control + sub-block
        if sub_n2d[s] > 0
            _matxs_next!(rd," 8d ")
            _matxs_holl!(rd,1)                     # hmtx (ignored)
            hdr = _matxs_ints!(rd,2+2*Ngo)
            lord = hdr[1]
            legendre_order = max(legendre_order,lord-1)
            jband = hdr[3:2+Ngo]
            ijj   = hdr[3+Ngo:2+2*Ngo]
            kmax = lord*sum(jband)
            _matxs_next!(rd," 9d ")
            mdat = _matxs_reals!(rd,kmax)
            scat = zeros(Ngi,Ngo,lord)
            pos = 0
            for gf in range(1,Ngo)
                jb = jband[gf]; ij = ijj[gf]
                if jb == 0 continue end
                for l in range(1,lord), gi in (ij+jb-1):-1:ij
                    pos += 1
                    scat[gi,gf,l] = mdat[pos]
                end
            end
            mats[(i,o)] = scat
        end
    end

    # Assemble Multigroup_Cross_Sections for each incident particle of this material
    for n in range(1,Npart)
        mcs = Multigroup_Cross_Sections(Ng[n])
        d = vecs[n]
        mcs.set_total(d["tot0"])
        mcs.set_absorption(d["abs"])
        mcs.set_momentum_transfer(d["momt"])
        mcs.set_energy_deposition(d["edep"])      # authoritative Ng+1 (heat is interop-only)
        mcs.set_boundary_stopping_powers(d["stpw"])
        mcs.set_charge_deposition(d["cdep"])
        for nout in range(1,Npart)                # push scattering in outgoing order
            mcs.set_scattering(mats[(n,nout)])
        end
        multigroup_cross_sections[n,mfile] = mcs
    end
end

#----
# Store in the Cross_Sections structure
#----
energy = [(energy_boundaries[n][1] + energy_boundaries[n][2])/2 for n in range(1,Npart)]
cutoff = [energy_boundaries[n][end] for n in range(1,Npart)]

cross_sections.set_energy(energy)
cross_sections.set_cutoff(cutoff)
cross_sections.set_number_of_groups(Ng)
cross_sections.set_legendre_order(legendre_order)
cross_sections.set_energy_boundaries(energy_boundaries)
cross_sections.set_multigroup_cross_sections(multigroup_cross_sections)

end

#----
# Unified MATXS record reader (BCD/ASCII or binary)
#----

mutable struct _Matxs_Reader
    binary ::Bool
    lines  ::Vector{String}   # BCD : record lines (columns preserved)
    li     ::Int64
    io     ::IO               # binary : record stream
    buf    ::Vector{UInt8}    # binary : current record payload
    pos    ::Int64
end

"""
    _Matxs_Reader(is_bcd::Bool,raw::Vector{UInt8})

Build a MATXS record reader over the raw file bytes, for BCD/ASCII (`is_bcd=true`) or binary.

"""
function _Matxs_Reader(is_bcd::Bool,raw::Vector{UInt8})
    if is_bcd
        # Keep every line (do NOT drop blanks : an all-space line is a legitimate hollerith
        # continuation word). Strip only a trailing carriage return.
        lines = String[]
        for line in split(String(raw),'\n')
            push!(lines,endswith(line,"\r") ? String(line[1:end-1]) : String(line))
        end
        return _Matxs_Reader(false,lines,1,IOBuffer(),UInt8[],1)
    else
        return _Matxs_Reader(true,String[],1,IOBuffer(raw),UInt8[],1)
    end
end

"""
    _matxs_next!(rd::_Matxs_Reader,tag::String)

Advance to the next record. In BCD mode the current line must be the record `tag`; in binary
mode the next FORTRAN-style record is loaded (the `tag` is ignored).

"""
function _matxs_next!(rd::_Matxs_Reader,tag::String)
    if rd.binary
        m = read(rd.io,Int32)
        rd.buf = read(rd.io,Int(m))
        read(rd.io,Int32)                          # trailing length marker
        rd.pos = 1
    else
        if strip(rd.lines[rd.li]) != strip(tag) error("Expected MATXS record '$tag' at line $(rd.li): $(rd.lines[rd.li])") end
        rd.li += 1
    end
end

"""
    _matxs_reals!(rd::_Matxs_Reader,n::Int64)

Read `n` reals from the current record (BCD : 12-char `e12.5`, 6 per line; binary : `Float64`).

"""
function _matxs_reals!(rd::_Matxs_Reader,n::Int64)
    out = Vector{Float64}(undef,n)
    if rd.binary
        for k in 1:n
            out[k] = reinterpret(Float64,rd.buf[rd.pos:rd.pos+7])[1]; rd.pos += 8
        end
    else
        got = 0
        while got < n
            line = rd.lines[rd.li]; rd.li += 1
            k = min(6,n-got)
            for j in 1:k
                got += 1
                out[got] = parse(Float64,strip(line[(j-1)*12+1:j*12]))
            end
        end
    end
    return out
end

"""
    _matxs_ints!(rd::_Matxs_Reader,n::Int64)

Read `n` integers from the current record (BCD : 6-char `i6`, 12 per line; binary : `Int32`).

"""
function _matxs_ints!(rd::_Matxs_Reader,n::Int64)
    out = Vector{Int64}(undef,n)
    if rd.binary
        for k in 1:n
            out[k] = Int(reinterpret(Int32,rd.buf[rd.pos:rd.pos+3])[1]); rd.pos += 4
        end
    else
        got = 0
        while got < n
            line = rd.lines[rd.li]; rd.li += 1
            k = min(12,n-got)
            for j in 1:k
                got += 1
                out[got] = parse(Int,strip(line[(j-1)*6+1:j*6]))
            end
        end
    end
    return out
end

"""
    _matxs_holl!(rd::_Matxs_Reader,n::Int64)

Read `n` hollerith identifiers from the current record (BCD : 8-char `a8`, 8 per line; binary :
8-byte fields).

"""
function _matxs_holl!(rd::_Matxs_Reader,n::Int64)
    out = Vector{String}(undef,n)
    if rd.binary
        for k in 1:n
            out[k] = String(rd.buf[rd.pos:rd.pos+7]); rd.pos += 8
        end
    else
        got = 0
        while got < n
            k = min(8,n-got)
            line = rpad(rd.lines[rd.li],k*8); rd.li += 1   # tolerate trailing-space trimming
            for j in 1:k
                got += 1
                out[got] = line[(j-1)*8+1:j*8]
            end
        end
    end
    return out
end
