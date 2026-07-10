"""
    read_matxs(cross_sections::Cross_Sections)

Read a MATXS-format (NJOY/CCCC) cross-sections file and store its content in a Cross_Sections
structure. This is the counterpart of [`write_matxs`](@ref) and reverses both the standard
interop reactions (`*tot0`, `*heat`, `*char`, `*scat`) and the Dragon-facing charged-particle
reactions (`EMOMTR`, `BSTC`, `PSTC`). Older Radiant-private extension reactions
(`*momt`, `*stpw`, `*cdep`) are still accepted when present.

Radiant requires the `*tot0`, `*abs`, `*edep`, and either `*cdep` or `*char` vector for
each incident particle in order to reconstruct a complete cross-section library. The
`*momt` and `*stpw` vectors are optional and default to zero when absent. Standard-only
MATXS files that omit Radiant's required extension vectors are therefore rejected with an
explicit error rather than being reconstructed with guessed values.

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
densities = cross_sections.get_densities()
if any(densities .<= 0.0)
    error("MATXS input requires strictly positive material densities.")
end

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
    eb = _matxs_reals!(rd,ngrp_f[j]+1) ./ 1.0e6
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
                sfx = matxs_vector_key(strip(hvps[k]))
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
                for l in range(1,lord), gi in ij:-1:(ij-jb+1)
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
        d = get(vecs,n,Dict{String,Vector{Float64}}())
        rho = densities[mfile]
        code = matxs_particle_code(particles[n])
        total = _matxs_required_vector(d,"tot0",Ng[n],code)
        absorption = _matxs_required_vector(d,"abs",Ng[n],code)
        energy_deposition = _matxs_required_vector(d,"edep",Ng[n]+1,code)
        momentum_transfer = _matxs_optional_vector(d,"momt",Ng[n],code)
        stopping_powers = _matxs_optional_vector(d,"stpw",Ng[n]+1,code)
        if haskey(d,"cdep")
            charge = _matxs_vector(d,"cdep",Ng[n]+1,code)
        elseif haskey(d,"char")
            charge = vcat(_matxs_vector(d,"char",Ng[n],code),0.0)
        else
            error("MATXS file is missing required charge-deposition vector 'cdep' or 'char' for incident particle '$code'.")
        end
        mcs.set_total(total .* rho)
        mcs.set_absorption(absorption .* rho)
        mcs.set_momentum_transfer(momentum_transfer .* rho)
        mcs.set_energy_deposition(energy_deposition .* rho) # authoritative Ng+1 (heat is interop-only)
        mcs.set_boundary_stopping_powers(stopping_powers .* rho)
        charge = charge .* rho
        mcs.set_charge_deposition(charge)
        for nout in range(1,Npart)                # push scattering in outgoing order
            mcs.set_scattering(mats[(n,nout)] .* rho)
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

function _matxs_required_vector(data::Dict{String,Vector{Float64}},name::String,
                                length_expected::Int64,code::String)
    if !haskey(data,name)
        error("MATXS file is missing required vector '*$name' for incident particle '$code'.")
    end
    return _matxs_vector(data,name,length_expected,code)
end

function _matxs_optional_vector(data::Dict{String,Vector{Float64}},name::String,
                                length_expected::Int64,code::String)
    if !haskey(data,name)
        return zeros(length_expected)
    end
    return _matxs_vector(data,name,length_expected,code)
end

function _matxs_vector(data::Dict{String,Vector{Float64}},name::String,
                       length_expected::Int64,code::String)
    vector = data[name]
    if length(vector) != length_expected
        error("MATXS vector '*$name' for incident particle '$code' has length $(length(vector)); expected $length_expected.")
    end
    return vector
end

#----
# Unified MATXS record reader (BCD/ASCII or binary)
#----

mutable struct _Matxs_Reader
    binary ::Bool
    lines  ::Vector{String}   # BCD : record lines (columns preserved)
    li     ::Int64
    rec    ::Vector{String}   # BCD : current record payload lines
    ri     ::Int64
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
        return _Matxs_Reader(false,lines,1,String[],1,IOBuffer(),UInt8[],1)
    else
        return _Matxs_Reader(true,String[],1,String[],1,IOBuffer(raw),UInt8[],1)
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
        if rd.li > length(rd.lines) error("Expected MATXS record '$tag' at end of file.") end
        if !_matxs_is_tag_line(rd.lines[rd.li],tag)
            error("Expected MATXS record '$tag' at line $(rd.li): $(rd.lines[rd.li])")
        end
        first_line = rd.lines[rd.li]
        rd.li += 1
        payload = String[]
        while rd.li <= length(rd.lines) && !_matxs_starts_record(rd.lines[rd.li])
            push!(payload,rd.lines[rd.li])
            rd.li += 1
        end
        rd.rec = _matxs_bcd_payload(tag,first_line,payload)
        rd.ri = 1
    end
end

function _matxs_is_tag_line(line::String,tag::String)
    return length(line) >= 4 && strip(line[1:4]) == strip(tag)
end

function _matxs_starts_record(line::String)
    tags = ("0v","1d","2d","3d","4d","5d","6d","7d","8d","9d","10d")
    return length(line) >= 4 && strip(line[1:4]) in tags
end

function _matxs_bcd_payload(tag::String,first_line::String,continuation::Vector{String})
    line = rpad(first_line,4)
    rest = line[5:end]
    if strip(rest) == ""
        return continuation
    end
    if tag == " 0v "
        rest = rpad(rest,35)
        return String[rest[1:8] * rest[10:17], rest[end-5:end]]
    elseif tag == " 1d "
        return vcat(String[rest[3:end]],continuation)
    elseif tag == " 2d " || tag == " 3d " || tag == " 6d "
        return vcat(String[rest[5:end]],continuation)
    elseif tag == " 4d " || tag == " 7d " || tag == " 9d " || tag == "10d "
        return vcat(String[rest[9:end]],continuation)
    elseif tag == " 5d "
        rest = rpad(rest,20)
        real_lines = String[rest[9:end]]
        int_lines = String[]
        for cont in continuation
            cont = rpad(cont,48)
            push!(real_lines,cont[1:24])
            push!(int_lines,cont[25:48])
        end
        return vcat(String[rest[1:8]],real_lines,int_lines)
    elseif tag == " 8d "
        return vcat(String[rpad(rest,12)[5:12]],continuation)
    else
        error("Unsupported MATXS BCD record tag '$tag'.")
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
            line = rd.rec[rd.ri]; rd.ri += 1
            k = min(length(line) ÷ 12,n-got)
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
            line = rd.rec[rd.ri]; rd.ri += 1
            k = min(length(line) ÷ 6,n-got)
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
            line = rd.rec[rd.ri]; rd.ri += 1
            k = min(max(length(line) ÷ 8,1),n-got)
            line = rpad(line,k*8)   # tolerate trailing-space trimming
            for j in 1:k
                got += 1
                out[got] = line[(j-1)*8+1:j*8]
            end
        end
    end
    return out
end

function matxs_vector_key(name::AbstractString)
    if name == "EMOMTR"
        return "momt"
    elseif name in ("BSTC","CSTC","PSTC")
        return "stpw"
    elseif length(name) >= 2
        sfx = name[2:end]
        return sfx == "char" ? "char" : sfx
    end
    return name
end
