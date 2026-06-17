"""
    endf_reading

Module for reading and processing ENDF (Evaluated Nuclear Data File) format data.

"""

"""
    ZAP(particle::Particle)

Gives the ENDF product identifier associated with a Radiant particle.

# Input Argument(s)
- `particle::Particle` : Radiant particle.

# Output Argument(s)
- `zap::Union{Nothing,Int}` : ENDF ZAP identifier, or `nothing` when the particle is not
  mapped.

# Reference(s)
- ENDF-6 Formats Manual, product identifier conventions.
"""
function ZAP(particle::Particle)::Union{Nothing,Int}
    if is_proton(particle)
        return 1001
    elseif is_alpha(particle)
        return 2004
    end

    particle_name = lowercase(string(get_type(particle)))
    if particle_name == "deuteron"
        return 1002
    elseif particle_name == "triton"
        return 1003
    elseif particle_name == "alpha"
        return 2004
    end
    return nothing
end

"""
    interpret_law5_list_record(E::Float64, LTP::Int, NW::Int, NL::Int,
    A::Vector{Float64}, LIDP::Int=-1)

Interprets one ENDF LAW=5 LIST record at a single incident energy and validates the
reported LIST length against the expected LAW=5 payload size.

# Input Argument(s)
- `E::Float64` : incident energy associated with the LIST record.
- `LTP::Int` : ENDF LAW=5 representation flag.
- `NW::Int` : number of values reported in the ENDF LIST record; compared with the
  expected LAW=5 value and warned on mismatch.
- `NL::Int` : highest nuclear partial-wave order for `LTP <= 2`, or number of tabulated cosines for `LTP > 2`.
- `A::Vector{Float64}` : raw LIST payload values; its length must match the expected
  LAW=5 payload size.
- `LIDP::Int` : ENDF identical-particle flag.

# Output Argument(s)
- `record::NamedTuple` : interpreted LAW=5 data with fields `E`, `LTP`, `NL`, `mu`,
  `pni`, `b`, `a`, and `c`.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering.
"""
function interpret_law5_list_record(E::Float64, LTP::Int, NW::Int, NL::Int,
                                    A::Vector{Float64}, LIDP::Int = -1)
    local n_expected

    if LTP > 2
        n_expected = 2 * NL
        if NW != n_expected
            @warn "LAW=5 LTP=$(LTP) reports NW=$(NW), expected $(n_expected) for NL=$(NL)."
        end
        if length(A) != n_expected
            error("LAW=5 LTP=$(LTP) expected $(n_expected) values in LIST array, got $(length(A)).")
        end

        mu = A[1:2:end]
        pni_raw = A[2:2:end]
        pni = LTP == 14 ? exp.(pni_raw) : pni_raw
        return (E=E, LTP=LTP, NL=NL, mu=mu, pni=pni, b=Float64[], a=ComplexF64[], c=Float64[])
    elseif LTP == 1
        if LIDP == 0
            nb = 2 * NL + 1
            n_expected = 4 * NL + 3
        elseif LIDP == 1
            nb = NL + 1
            n_expected = 3 * NL + 3
        else
            error("LAW=5 LTP=1 requires LIDP=0 or LIDP=1, got LIDP=$(LIDP).")
        end
        if NW != n_expected
            @warn "LAW=5 LTP=1 reports NW=$(NW), expected $(n_expected) for NL=$(NL), LIDP=$(LIDP)."
        end
        if length(A) != n_expected
            error("LAW=5 LTP=1 expected $(n_expected) values in LIST array, got $(length(A)).")
        end

        b = A[1:nb]
        a = ComplexF64[]
        coeffs = A[(nb + 1):n_expected]
        for k in 1:2:length(coeffs)-1
            push!(a, complex(coeffs[k], coeffs[k + 1]))
        end
        return (E=E, LTP=LTP, NL=NL, mu=Float64[], pni=Float64[], b=b, a=a, c=Float64[])
    elseif LTP == 2
        n_expected = NL + 1
        if NW != n_expected
            @warn "LAW=5 LTP=2 reports NW=$(NW), expected $(n_expected) for NL=$(NL)."
        end
        if length(A) != n_expected
            error("LAW=5 LTP=2 expected $(n_expected) values in LIST array, got $(length(A)).")
        end

        c = A[1:n_expected]
        return (E=E, LTP=LTP, NL=NL, mu=Float64[], pni=Float64[], b=Float64[], a=ComplexF64[], c=c)
    else
        error("Unsupported LAW=5 LTP=$(LTP).")
    end
end

"""
    validate_section_ids(mat::Int, mf::Int, mt::Int,
    expected_mat::Int, expected_mf::Int, expected_mt::Int,
    record_name::AbstractString)

Checks that a parsed ENDF record still belongs to the expected MAT/MF/MT section.

# Input Argument(s)
- `mat::Int` : parsed ENDF MAT identifier.
- `mf::Int` : parsed ENDF MF identifier.
- `mt::Int` : parsed ENDF MT identifier.
- `expected_mat::Int` : expected ENDF MAT identifier.
- `expected_mf::Int` : expected ENDF MF identifier.
- `expected_mt::Int` : expected ENDF MT identifier.
- `record_name::AbstractString` : record label used in error messages.

# Output Argument(s)
- `nothing` : returns `nothing` when the identifiers match; otherwise throws an error.
"""
function validate_section_ids(mat::Int, mf::Int, mt::Int,
                              expected_mat::Int, expected_mf::Int, expected_mt::Int,
                              record_name::AbstractString)
    if mat != expected_mat || mf != expected_mf || mt != expected_mt
        error("$(record_name) crossed ENDF section boundary: expected MAT/MF/MT " *
              "$(expected_mat)/$(expected_mf)/$(expected_mt), got $(mat)/$(mf)/$(mt).")
    end
end

"""
    validate_interp_metadata(NBT::Vector{Int}, INT::Vector{Int}, n_points::Int,
    record_name::AbstractString; supported_ints=nothing)

Checks ENDF interpolation-region metadata for consistent lengths and boundaries.

# Input Argument(s)
- `NBT::Vector{Int}` : interpolation-region upper boundaries.
- `INT::Vector{Int}` : interpolation laws for each region.
- `n_points::Int` : number of tabulated points or subrecords governed by the metadata.
- `record_name::AbstractString` : record label used in error messages.
- `supported_ints` : optional collection of interpolation laws allowed for this record.

# Output Argument(s)
- `nothing` : returns `nothing` when the metadata is consistent; otherwise throws an error.
"""
function validate_interp_metadata(NBT::Vector{Int}, INT::Vector{Int}, n_points::Int,
                                  record_name::AbstractString; supported_ints=nothing)
    length(NBT) == length(INT) ||
        error("$(record_name) has inconsistent interpolation metadata: " *
              "length(NBT)=$(length(NBT)), length(INT)=$(length(INT)).")

    isempty(NBT) && return nothing

    last_nbt = 0
    for nbt in NBT
        nbt > last_nbt ||
            error("$(record_name) interpolation boundaries must be strictly increasing.")
        last_nbt = nbt
    end

    NBT[end] == n_points ||
        error("$(record_name) final interpolation boundary NBT[end]=$(NBT[end]) " *
              "does not match number of points $(n_points).")

    if supported_ints !== nothing
        for itp in INT
            itp in supported_ints ||
                error("$(record_name) has unsupported ENDF interpolation INT=$(itp).")
        end
    end

    return nothing
end

"""
    parse_endf_float(field::AbstractString)

Parses an ENDF floating-point field written in the fixed-width ENDF notation.

# Input Argument(s)
- `field::AbstractString` : eleven-character ENDF floating-point field.

# Output Argument(s)
- `value::Float64` : parsed floating-point value.

# Reference(s)
- ENDF-6 Formats Manual, numeric field representation.
"""
function parse_endf_float(field::AbstractString)::Float64
    s = strip(field)
    isempty(s) && return 0.0

    # Find exponent sign near the end (not the leading sign)
    pos = 0
    for i in lastindex(s):-1:2
        if s[i] == '+' || s[i] == '-'
            pos = i
            break
        end
    end
    if pos == 0
        return parse(Float64, s)
    end
    mant = strip(s[1:pos-1])
    expo = strip(s[pos:end])
    m = isempty(mant) ? 0.0 : parse(Float64, mant)
    e = parse(Int, expo)
    return m * 10.0^e
end

"""
    parse_endf_line(line::AbstractString)

Parses one ENDF-6 line into six data fields and the MAT/MF/MT/NS identifiers.

# Input Argument(s)
- `line::AbstractString` : ENDF-6 formatted line.

# Output Argument(s)
- `fields::Tuple` : six raw eleven-character fields followed by `mat`, `mf`, `mt`, and
  `ns` identifiers.

# Reference(s)
- ENDF-6 Formats Manual, record layout.
"""
function parse_endf_line(line::AbstractString)
    if length(line) < 80
        error("Line too short to be ENDF-6: '$line'")
    end

    f1 = line[1:11];  f2 = line[12:22]; f3 = line[23:33]
    f4 = line[34:44]; f5 = line[45:55]; f6 = line[56:66]

    mat = parse(Int, strip(line[67:70]))
    mf = parse(Int, strip(line[71:72]))
    mt = parse(Int, strip(line[73:75]))
    ns = parse(Int, strip(line[76:80]))

    return (f1,f2,f3,f4,f5,f6, mat,mf,mt,ns)
end

"""
    read_CONT(lines::Vector{String}, i::Int)

Reads one ENDF CONT or HEAD-like record.

# Input Argument(s)
- `lines::Vector{String}` : ENDF file lines.
- `i::Int` : index of the record to read.

# Output Argument(s)
- `record::Tuple` : parsed `C1`, `C2`, `L1`, `L2`, `N1`, `N2`, identifiers, next line
  index, and raw line.

# Reference(s)
- ENDF-6 Formats Manual, CONT and HEAD records.
"""
function read_CONT(lines::Vector{String}, i::Int)
    ln = lines[i]
    (f1,f2,f3,f4,f5,f6,mat,mf,mt,ns) = parse_endf_line(ln)
    C1 = parse_endf_float(f1)
    C2 = parse_endf_float(f2)
    L1 = parse(Int, strip(f3))
    L2 = parse(Int, strip(f4))
    N1 = parse(Int, strip(f5))
    N2 = parse(Int, strip(f6))
    return (C1,C2,L1,L2,N1,N2, mat,mf,mt, ns, i+1, ln)
end

"""
    read_interp_pairs(lines::Vector{String}, i::Int, NR::Int)

Reads ENDF interpolation boundary pairs from TAB1 or TAB2 records.

# Input Argument(s)
- `lines::Vector{String}` : ENDF file lines.
- `i::Int` : index of the first interpolation-pair line.
- `NR::Int` : number of interpolation regions.

# Output Argument(s)
- `(NBT, INT, i_next)::Tuple` : interpolation boundaries, interpolation laws, and the
  next unread line index.

# Reference(s)
- ENDF-6 Formats Manual, TAB1 and TAB2 interpolation tables.
"""
function read_interp_pairs(lines::Vector{String}, i::Int, NR::Int)
    NBT = Int[]
    INT = Int[]
    j = i
    while length(NBT) < NR
        (f1,f2,f3,f4,f5,f6,mat,mf,mt,ns) = parse_endf_line(lines[j])
        raw = (strip(f1),strip(f2),strip(f3),strip(f4),strip(f5),strip(f6))
        vals = Int[]
        for s in raw
            push!(vals, isempty(s) ? 0 : parse(Int, s))
        end
        for p in 1:2:6
            if length(NBT) >= NR
                break
            end
            push!(NBT, vals[p])
            push!(INT, vals[p+1])
        end
        j += 1
    end
    return (NBT, INT, j)
end

"""
    read_TAB1(lines::Vector{String}, i::Int)

Reads an ENDF TAB1 record, including interpolation metadata and tabulated `(X, Y)` pairs.

# Input Argument(s)
- `lines::Vector{String}` : ENDF file lines.
- `i::Int` : index of the TAB1 header record.

# Output Argument(s)
- `record::Tuple` : TAB1 header values, interpolation arrays, tabulated values,
  identifiers, next line index, and raw header line.

# Reference(s)
- ENDF-6 Formats Manual, TAB1 records.
"""
function read_TAB1(lines::Vector{String}, i::Int)
    (C1,C2,L1,L2,NR,NP, mat,mf,mt,ns, j, ln0) = read_CONT(lines, i)

    (NBT, INT, j2) = read_interp_pairs(lines, j, NR)
    validate_interp_metadata(NBT, INT, NP, "TAB1"; supported_ints=1:6)
    j = j2

    X = Float64[]
    Y = Float64[]
    while length(X) < NP
        (f1,f2,f3,f4,f5,f6,mat2,mf2,mt2,ns2) = parse_endf_line(lines[j])
        validate_section_ids(mat2, mf2, mt2, mat, mf, mt, "TAB1 data")
        vals = (parse_endf_float(f1), parse_endf_float(f2),
                parse_endf_float(f3), parse_endf_float(f4),
                parse_endf_float(f5), parse_endf_float(f6))
        for k in 1:2:6
            if length(X) >= NP
                break
            end
            push!(X, vals[k])
            push!(Y, vals[k+1])
        end
        j += 1
    end

    return (C1,C2,L1,L2,NR,NP,NBT,INT,X,Y, mat,mf,mt,ns, j, ln0)
end

"""
    read_TAB2(lines::Vector{String}, i::Int)

Reads an ENDF TAB2 record, including interpolation metadata for a following sequence of
subrecords.

# Input Argument(s)
- `lines::Vector{String}` : ENDF file lines.
- `i::Int` : index of the TAB2 header record.

# Output Argument(s)
- `record::Tuple` : TAB2 header values, interpolation arrays, identifiers, next line
  index, and raw header line.

# Reference(s)
- ENDF-6 Formats Manual, TAB2 records.
"""
function read_TAB2(lines::Vector{String}, i::Int)
    (C1,C2,L1,L2,NR,NP, mat,mf,mt,ns, j, ln0) = read_CONT(lines, i)
    (NBT, INT, j2) = read_interp_pairs(lines, j, NR)
    validate_interp_metadata(NBT, INT, NP, "TAB2")
    j = j2
    return (C1,C2,L1,L2,NR,NP,NBT,INT, mat,mf,mt,ns, j, ln0)
end

"""
    read_LIST(lines::Vector{String}, i::Int)

Reads an ENDF LIST record and its floating-point payload.

# Input Argument(s)
- `lines::Vector{String}` : ENDF file lines.
- `i::Int` : index of the LIST header record.

# Output Argument(s)
- `record::Tuple` : LIST header values, payload array, identifiers, next line index, and
  raw header line.

# Reference(s)
- ENDF-6 Formats Manual, LIST records.
"""
function read_LIST(lines::Vector{String}, i::Int)
    (C1,C2,L1,L2,NPL,N2, mat,mf,mt,ns, j, ln0) = read_CONT(lines, i)
    B = Float64[]
    while length(B) < NPL
        (f1,f2,f3,f4,f5,f6,mat2,mf2,mt2,ns2) = parse_endf_line(lines[j])
        validate_section_ids(mat2, mf2, mt2, mat, mf, mt, "LIST data")
        append!(B, (parse_endf_float(f1), parse_endf_float(f2), parse_endf_float(f3),
                    parse_endf_float(f4), parse_endf_float(f5), parse_endf_float(f6)))
        j += 1
    end
    resize!(B, NPL)
    return (C1,C2,L1,L2,NPL,N2,B, mat,mf,mt,ns, j, ln0)
end

# =============================================================================
# Navigation utilities
# =============================================================================

"""
    find_section_start(lines::Vector{String}, mf_target::Int, mt_target::Int; start_i::Int=1)

Finds the first line of an ENDF section with the requested MF and MT identifiers.

# Input Argument(s)
- `lines::Vector{String}` : ENDF file lines.
- `mf_target::Int` : requested MF number.
- `mt_target::Int` : requested MT number.
- `start_i::Int` : line index at which the search starts.

# Output Argument(s)
- `i::Union{Nothing,Int}` : section start line, or `nothing` when no section is found.
"""
function find_section_start(lines::Vector{String}, mf_target::Int, mt_target::Int; start_i::Int=1)
    for i in start_i:length(lines)
        ln = lines[i]
        if length(ln) < 80; continue; end
        (_,_,_,_,_,_, _, mf, mt, _) = parse_endf_line(ln)
        if mf == mf_target && mt == mt_target
            return i
        end
    end
    return nothing
end

"""
    find_first_mf(lines::Vector{String}, mf_target::Int)

Finds the first line in an ENDF file with the requested MF identifier.

# Input Argument(s)
- `lines::Vector{String}` : ENDF file lines.
- `mf_target::Int` : requested MF number.

# Output Argument(s)
- `i::Union{Nothing,Int}` : first matching line, or `nothing` when no MF section is found.
"""
function find_first_mf(lines::Vector{String}, mf_target::Int)
    for (i,ln) in enumerate(lines)
        if length(ln) < 80; continue; end
        (_,_,_,_,_,_, _, mf, _, _) = parse_endf_line(ln)
        if mf == mf_target
            return i
        end
    end
    return nothing
end

# =============================================================================
# Clean extraction functions (logic without verbose output)
# =============================================================================

"""
    parse_law5_section(lines::Vector{String}, i::Int)

Parses an ENDF MF=6 LAW=5 angular-distribution section.

# Input Argument(s)
- `lines::Vector{String}` : ENDF file lines.
- `i::Int` : index of the LAW=5 TAB2 header.

# Output Argument(s)
- `section::Tuple` : spin, identical-particle flag, LTP value, LTP>2 tables, LTP=1
  coefficient tables, LTP=2 coefficient tables, interpolation metadata, and next line
  index.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=5 charged-particle elastic scattering.
"""
function parse_law5_section(lines::Vector{String}, i::Int)
    (SPI,_,LIDP,_,NRt,NE,NBTt,INTt, mat,mf,mt,ns, i, _) = read_TAB2(lines, i)

    mupni_tables = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
    ba_tables_ltp1 = Dict{Float64, Tuple{Int, Vector{Float64}, Vector{ComplexF64}}}()
    c_tables_ltp2 = Dict{Float64, Tuple{Int, Vector{Float64}}}()
    LTP_seen = nothing

    for ib in 1:NE
        # Read LIST: [0.0, Ein, LTP, 0, NW, NL / A(E)]
        (C1l,Ein,LTP,_,NW,NL,Aarr, _, _, _, _, i, _) = read_LIST(lines, i)

        if LTP_seen === nothing
            LTP_seen = LTP
        elseif LTP != LTP_seen
            error("Mixed LAW=5 LTP values in one section: first LTP=$(LTP_seen), got LTP=$(LTP) at E=$(Ein).")
        end

        # Interpret LIST array
        rec = interpret_law5_list_record(Ein, LTP, NW, NL, Aarr, LIDP)
        if !isempty(rec.mu)
            mupni_tables[Ein] = (rec.mu, rec.pni)
        end
        if !isempty(rec.b) || !isempty(rec.a)
            ba_tables_ltp1[Ein] = (rec.NL, rec.b, rec.a)
        end
        if !isempty(rec.c)
            c_tables_ltp2[Ein] = (rec.NL, rec.c)
        end
    end

    return (SPI, LIDP, LTP_seen, mupni_tables, ba_tables_ltp1, c_tables_ltp2, NBTt, INTt, i)
end

"""
    parse_product_subsection(lines::Vector{String}, i::Int, ZA::Float64, AWR::Float64)

Parses one ENDF MF=6 product subsection, including its product TAB1 record. LAW=5
subsections are fully parsed; unsupported LAW values return only the product header and
yield/interpolation metadata, leaving any LAW-specific records unread.

# Input Argument(s)
- `lines::Vector{String}` : ENDF file lines.
- `i::Int` : index of the product TAB1 record.
- `ZA::Float64` : target ZA value from the MF=6 section header.
- `AWR::Float64` : target atomic-weight ratio from the MF=6 section header.

# Output Argument(s)
- `(product, i_next)::Tuple` : product NamedTuple and the next unread line index. For
  unsupported LAW values, `i_next` points immediately after the product TAB1 record.

# Reference(s)
- ENDF-6 Formats Manual, MF=6 product subsections.
"""
function parse_product_subsection(lines::Vector{String}, i::Int, ZA::Float64, AWR::Float64)
    # Read product TAB1
    (ZAP,AWP,LIP,LAW,NR,NP,NBT,INT,Eint,yi, _, _, _, _, i, _) = read_TAB1(lines, i)

    mupni_tables = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
    ba_tables_ltp1 = Dict{Float64, Tuple{Int, Vector{Float64}, Vector{ComplexF64}}}()
    c_tables_ltp2 = Dict{Float64, Tuple{Int, Vector{Float64}}}()
    SPI_seen = nothing
    LIDP_seen = nothing
    LTP_seen = nothing
    NBT_seen = Int[]
    INT_seen = Int[]

    # Handle LAW-specific structure
    if LAW == 5
        (SPI_seen, LIDP_seen, LTP_seen, mupni_tables, ba_tables_ltp1, c_tables_ltp2, NBT_seen, INT_seen, i) =
            parse_law5_section(lines, i)
    end

    return (ZAP=ZAP, AWP=AWP, LIP=LIP, LAW=LAW, Eint=Eint, yi=yi,
            mupni_tables=mupni_tables, ba_tables_ltp1=ba_tables_ltp1, c_tables_ltp2=c_tables_ltp2,
            SPI=SPI_seen, LIDP=LIDP_seen, LTP=LTP_seen,
            mf6_nbt=NBT_seen, mf6_int=INT_seen,
            ZA=ZA, AWR=AWR), i
end

"""
    read_mf3_mt(lines::Vector{String}, mt::Int)

Reads an ENDF MF=3 cross-section table for a requested MT reaction.

# Input Argument(s)
- `lines::Vector{String}` : ENDF file lines.
- `mt::Int` : requested MT reaction number.

# Output Argument(s)
- `data::Tuple` : target `ZA`, `AWR`, energy grid, cross-section values, and TAB1
  interpolation metadata.

# Reference(s)
- ENDF-6 Formats Manual, MF=3 reaction cross sections.
"""
function read_mf3_mt(lines::Vector{String}, mt::Int)
    i_mf3 = find_first_mf(lines, 3)
    i_mf3 === nothing && error("MF=3 not found.")

    i0 = find_section_start(lines, 3, mt; start_i=i_mf3)
    i0 === nothing && error("MF=3 MT=$mt not found.")

    (ZA,AWR,JP,LCT,NK,zero, mat,mf,mt_check,ns, i, _) = read_CONT(lines, i0)
    validate_section_ids(mat, mf, mt_check, mat, 3, mt, "MF=3 HEAD")
    (C1,C2,L1,L2,NR,NP,NBT,INT,E,S, mat2,mf2,mt2,ns2, j, _) = read_TAB1(lines, i)
    validate_section_ids(mat2, mf2, mt2, mat, 3, mt, "MF=3 TAB1")

    return (ZA, AWR, E, S, NBT, INT)
end

"""
    read_mf6_mt(lines::Vector{String}, mt::Int; source::AbstractString="ENDF data")

Reads ENDF MF=6 product energy-angle distributions for a requested MT reaction. LAW=5
products are parsed fully; a final unsupported product is allowed as a header-only recoil
or residual product, while unsupported non-final products throw an error because their
record length is not known here.

# Input Argument(s)
- `lines::Vector{String}` : ENDF file lines.
- `mt::Int` : requested MT reaction number.
- `source::AbstractString="ENDF data"` : source label used in error messages.

# Output Argument(s)
- `products::Vector{NamedTuple}` : parsed product subsections for the requested MF=6/MT
  section.

# Reference(s)
- ENDF-6 Formats Manual, MF=6 product distributions.
"""
function read_mf6_mt(lines::Vector{String}, mt::Int; source::AbstractString="ENDF data")
    i_mf6 = find_first_mf(lines, 6)
    i_mf6 === nothing && error("MF=6 not found in $(source).")

    i0 = find_section_start(lines, 6, mt; start_i=i_mf6)
    i0 === nothing && error("MF=6 MT=$mt not found in $(source).")

    # Read MF=6 HEAD (A-section)
    (ZA,AWR,JP,LCT,NK,zero, mat,mf,mt_check,ns, i, _) = read_CONT(lines, i0)
    validate_section_ids(mat, mf, mt_check, mat, 6, mt, "MF=6 HEAD")

    products = NamedTuple[]
    for kprod in 1:NK
        prod, i = parse_product_subsection(lines, i, ZA, AWR)
        if prod.LAW != 5 && kprod != NK
            error("Unsupported MF=6 product LAW=$(prod.LAW) before the final product; cannot skip LAW-specific data safely.")
        end
        push!(products, prod)
    end

    return products
end

# =============================================================================
# Interpolation helpers
# =============================================================================

# ENDF interpolation types commonly used in TAB1:
# 1 = histogram, 2 = lin-lin, 3 = lin-log, 4 = log-lin, 5 = log-log,
# 6 = charged-particle low-energy interpolation.
"""
    interp_pair(x::Float64, x1::Float64, y1::Float64, x2::Float64, y2::Float64,
    itp::Int; T::Float64=0.0)

Interpolates between two ENDF tabulated points using an ENDF interpolation law, including
the charged-particle `INT=6` law in its logarithm-free form.

# Input Argument(s)
- `x::Float64` : interpolation point.
- `x1::Float64` : lower tabulated abscissa.
- `y1::Float64` : value at `x1`.
- `x2::Float64` : upper tabulated abscissa.
- `y2::Float64` : value at `x2`.
- `itp::Int` : ENDF interpolation law (`1` histogram, `2` lin-lin, `3` lin-log,
  `4` log-lin, `5` log-log, `6` charged-particle low-energy interpolation).
- `T::Float64` : threshold energy used by `INT=6`, defaulting to zero for elastic
  charged-particle data; endothermic reactions should pass their threshold explicitly.

# Output Argument(s)
- `y::Float64` : interpolated value.

# Reference(s)
- ENDF-6 Formats Manual, interpolation laws including the charged-particle `INT=6` law.
"""
function interp_pair(x::Float64, x1::Float64, y1::Float64, x2::Float64, y2::Float64,
                     itp::Int; T::Float64=0.0)
    if x2 == x1
        return y1
    end

    if itp == 1
        return y1
    elseif itp == 2
        t = (x - x1) / (x2 - x1)
        return y1 + t * (y2 - y1)
    elseif itp == 3
        x > 0 && x1 > 0 && x2 > 0 || error("linear-log interpolation requires positive x.")
        t = (log(x) - log(x1)) / (log(x2) - log(x1))
        return y1 + t * (y2 - y1)
    elseif itp == 4
        y1 > 0 && y2 > 0 || error("log-linear interpolation requires positive y.")
        t = (x - x1) / (x2 - x1)
        return exp(log(y1) + t * (log(y2) - log(y1)))
    elseif itp == 5
        x > 0 && x1 > 0 && x2 > 0 && y1 > 0 && y2 > 0 ||
            error("log-log interpolation requires positive x and y.")
        t = (log(x) - log(x1)) / (log(x2) - log(x1))
        return exp(log(y1) + t * (log(y2) - log(y1)))
    elseif itp == 6
        x > T && x1 > T && x2 > T || error("INT=6 interpolation requires x, x1, and x2 greater than T.")
        x > 0 && x1 > 0 && x2 > 0 || error("INT=6 interpolation requires positive energies.")
        y1 >= 0 && y2 >= 0 || error("INT=6 interpolation requires non-negative values.")
        if y1 == 0.0 || y2 == 0.0
            return 0.0
        end
        u = inv(sqrt(x - T))
        u1 = inv(sqrt(x1 - T))
        u2 = inv(sqrt(x2 - T))
        alpha = (u - u1) / (u2 - u1)
        return ((y2 * x2)^alpha * (y1 * x1)^(1.0 - alpha)) / x
    else
        error("Unsupported ENDF interpolation INT=$itp")
    end
end

"""
    select_interpolation(nbt::Vector{Int}, int::Vector{Int}, idx::Int)

Selects the ENDF interpolation law for a TAB1/TAB2 interval index.

# Input Argument(s)
- `nbt::Vector{Int}` : ENDF interpolation-region upper boundaries.
- `int::Vector{Int}` : ENDF interpolation laws for each region.
- `idx::Int` : interval index used by the caller.

# Output Argument(s)
- `itp::Int` : selected ENDF interpolation law, defaulting to lin-lin (`2`) when
  metadata is missing.

# Reference(s)
- ENDF-6 Formats Manual, TAB1 interpolation metadata.
"""
function select_interpolation(nbt::Vector{Int}, int::Vector{Int}, idx::Int)::Int
    if isempty(nbt) || isempty(int)
        return 2
    end
    itp = int[end]
    for r in 1:length(nbt)
        if (idx + 1) <= nbt[r]
            itp = int[r]
            break
        end
    end
    return itp
end

# Interpolate a TAB1 curve (E[], Y[]) with region boundaries NBT[] and INT[].
"""
    interp_TAB1(x::Float64, X::Vector{Float64}, Y::Vector{Float64}, NBT::Vector{Int}, INT::Vector{Int})

Interpolates an ENDF TAB1 curve using its region boundaries and interpolation laws.
For `INT=6`, this helper uses the elastic charged-particle default `T=0.0`.

# Input Argument(s)
- `x::Float64` : interpolation point.
- `X::Vector{Float64}` : sorted TAB1 abscissa grid.
- `Y::Vector{Float64}` : TAB1 values.
- `NBT::Vector{Int}` : interpolation-region upper boundaries.
- `INT::Vector{Int}` : interpolation laws for each region; `INT=6` is evaluated with
  `T=0.0` by this helper.

# Output Argument(s)
- `y::Float64` : interpolated value, clamped to endpoint values outside the tabulated
  range.

# Reference(s)
- ENDF-6 Formats Manual, interpolation laws including the charged-particle `INT=6` law.
"""
function interp_TAB1(x::Float64, X::Vector{Float64}, Y::Vector{Float64}, NBT::Vector{Int}, INT::Vector{Int})
    length(X) == length(Y) || error("TAB1 interpolation requires length(X) == length(Y).")
    !isempty(X) || error("TAB1 interpolation requires at least one point.")

    if x <= X[1]
        return Y[1]
    elseif x >= X[end]
        return Y[end]
    end

    i = searchsortedlast(X, x)
    i = clamp(i, 1, length(X) - 1)

    itp = select_interpolation(NBT, INT, i)

    return interp_pair(x, X[i], Y[i], X[i + 1], Y[i + 1], itp)
end

"""
    interp_linear(x::Vector{Float64}, y::Vector{Float64}, xq::Float64)

Interpolates a tabulated function linearly on a sorted grid.

# Input Argument(s)
- `x::Vector{Float64}` : sorted abscissa grid.
- `y::Vector{Float64}` : tabulated values.
- `xq::Float64` : interpolation point.

# Output Argument(s)
- `yq::Float64` : linearly interpolated value, clamped to endpoint values outside the
  tabulated range.
"""
function interp_linear(x::Vector{Float64}, y::Vector{Float64}, xq::Float64)::Float64
    length(x) == length(y) || error("Linear interpolation requires length(x) == length(y).")
    n = length(x)
    n >= 2 || error("Linear interpolation requires at least two points.")

    if xq <= x[1]
        return y[1]
    elseif xq >= x[end]
        return y[end]
    end

    j = searchsortedfirst(x, xq)
    j = clamp(j, 2, n)
    x1 = x[j - 1]
    x2 = x[j]
    y1 = y[j - 1]
    y2 = y[j]
    t = (xq - x1) / (x2 - x1)
    return y1 + t * (y2 - y1)
end
