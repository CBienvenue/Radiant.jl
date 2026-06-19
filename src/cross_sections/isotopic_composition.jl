"""
    isotopic_composition(Z::Int64)

Return the natural isotopic composition for atomic number Z.

# Input Argument(s)
- `Z::Int64`: atomic number of the element.

# Output Argument(s)
- `comp::Vector{Tuple{Int, Float64}}`: vector of (A, abundance) pairs, where A is the
  mass number and abundance is the natural isotopic fraction.

# Notes
- Abundances should sum to 1.0 for each Z.
- Populate `isotopic_composition_dict` with data for Z = 1..100.

# Reference(s)
- National Institute of Standards and Technology (NIST), Atomic Weights and
  Isotopic Compositions for All Elements,
  https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl

"""
function isotopic_composition(Z::Int64)
    if haskey(isotopic_composition_dict, Z)
        return isotopic_composition_dict[Z]
    else
        error("No natural isotopic composition data for Z=$(Z). Please define them.")
    end
end

"""
    material_isotopic_composition(materials, Z::Vector{Vector{Int64}})

Builds the isotope mass-number and atomic-fraction arrays for all materials. Missing
isotope entries for an element are filled with natural isotopic composition.

# Input Argument(s)
- `materials` : material objects used in the cross-section generation.
- `Z::Vector{Vector{Int64}}` : atomic numbers for each material.

# Output Argument(s)
- `A::Vector{Vector{Vector{Int64}}}` : isotope mass numbers for each material element.
- `atpercentA::Vector{Vector{Vector{Float64}}}` : isotope atomic fractions matching `A`.
"""
function material_isotopic_composition(materials, Z::Vector{Vector{Int64}})
    Nmat = length(materials)
    A = Vector{Vector{Vector{Int64}}}(undef,Nmat)
    atpercentA = Vector{Vector{Vector{Float64}}}(undef,Nmat)

    for n in range(1,Nmat)
        A[n] = [copy(v) for v in materials[n].get_mass_numbers()]
        atpercentA[n] = [copy(v) for v in materials[n].get_isotope_atomic_fractions()]
        Nz = length(Z[n])

        if length(A[n]) < Nz
            append!(A[n], [Int64[] for _ in 1:(Nz - length(A[n]))])
        end
        if length(atpercentA[n]) < Nz
            append!(atpercentA[n], [Float64[] for _ in 1:(Nz - length(atpercentA[n]))])
        end

        for i in range(1,Nz)
            if isempty(A[n][i]) || isempty(atpercentA[n][i])
                println("Isotopes not provided for Z=$(Z[n][i]) in material $(n); using natural composition.")
                pairs = isotopic_composition(Z[n][i])
                A[n][i] = [p[1] for p in pairs]
                atpercentA[n][i] = [p[2] for p in pairs]
            end
        end
    end

    return A, atpercentA
end

"""
    unique_isotopes(Z::Vector{Vector{Int64}}, A::Vector{Vector{Vector{Int64}}})

Returns the sorted unique `(Z, A)` isotope pairs used by the material isotope arrays.

# Input Argument(s)
- `Z::Vector{Vector{Int64}}` : atomic numbers for each material.
- `A::Vector{Vector{Vector{Int64}}}` : isotope mass numbers for each material element.

# Output Argument(s)
- `isotopes::Vector{Tuple{Int,Int}}` : sorted unique isotope pairs.
"""
function unique_isotopes(Z::Vector{Vector{Int64}}, A::Vector{Vector{Vector{Int64}}})
    unique_isotopes_set = Set{Tuple{Int,Int}}()
    for n in range(1,length(Z))
        for i in range(1,length(Z[n]))
            for Ai in A[n][i]
                push!(unique_isotopes_set, (Z[n][i], Ai))
            end
        end
    end
    return sort!(collect(unique_isotopes_set), by=x -> (x[1], x[2]))
end
# Dictionary mapping atomic numbers to natural isotopic compositions.
# Each value is a vector of (A, abundance) pairs.
const isotopic_composition_dict = Dict{Int, Vector{Tuple{Int, Float64}}}(
    1 => [(1, 0.999885), (2, 0.000115)],
    2 => [(3, 1.34e-06), (4, 0.99999866)],
    3 => [(6, 0.0759), (7, 0.9241)],
    4 => [(9, 1)],
    5 => [(10, 0.199), (11, 0.801)],
    6 => [(12, 0.9893), (13, 0.0107)],
    7 => [(14, 0.99636), (15, 0.00364)],
    8 => [(16, 0.99757), (17, 0.00038), (18, 0.00205)],
    9 => [(19, 1)],
    10 => [(20, 0.9048), (21, 0.0027), (22, 0.0925)],
    11 => [(23, 1)],
    12 => [(24, 0.7899), (25, 0.1), (26, 0.1101)],
    13 => [(27, 1)],
    14 => [(28, 0.92223), (29, 0.04685), (30, 0.03092)],
    15 => [(31, 1)],
    16 => [(32, 0.9499), (33, 0.0075), (34, 0.0425), (36, 0.0001)],
    17 => [(35, 0.7576), (37, 0.2424)],
    18 => [(36, 0.003336), (38, 0.000629), (40, 0.996035)],
    19 => [(39, 0.932581), (40, 0.000117), (41, 0.067302)],
    20 => [(40, 0.96941), (42, 0.00647), (43, 0.00135), (44, 0.02086), (46, 4e-05), (48, 0.00187)],
    21 => [(45, 1)],
    22 => [(46, 0.0825), (47, 0.0744), (48, 0.7372), (49, 0.0541), (50, 0.0518)],
    23 => [(50, 0.0025), (51, 0.9975)],
    24 => [(50, 0.04345), (52, 0.83789), (53, 0.09501), (54, 0.02365)],
    25 => [(55, 1)],
    26 => [(54, 0.05845), (56, 0.91754), (57, 0.02119), (58, 0.00282)],
    27 => [(59, 1)],
    28 => [(58, 0.68077), (60, 0.26223), (61, 0.011399), (62, 0.036346), (64, 0.009255)],
    29 => [(63, 0.6915), (65, 0.3085)],
    30 => [(64, 0.4917), (66, 0.2773), (67, 0.0404), (68, 0.1845), (70, 0.0061)],
    31 => [(69, 0.60108), (71, 0.39892)],
    32 => [(70, 0.2057), (72, 0.2745), (73, 0.0775), (74, 0.365), (76, 0.0773)],
    33 => [(75, 1)],
    34 => [(74, 0.0089), (76, 0.0937), (77, 0.0763), (78, 0.2377), (80, 0.4961), (82, 0.0873)],
    35 => [(79, 0.5069), (81, 0.4931)],
    36 => [(78, 0.00355), (80, 0.02286), (82, 0.11593), (83, 0.115), (84, 0.56987), (86, 0.17279)],
    37 => [(85, 0.7217), (87, 0.2783)],
    38 => [(84, 0.0056), (86, 0.0986), (87, 0.07), (88, 0.8258)],
    39 => [(89, 1)],
    40 => [(90, 0.5145), (91, 0.1122), (92, 0.1715), (94, 0.1738), (96, 0.028)],
    41 => [(93, 1)],
    42 => [(92, 0.1453), (94, 0.0915), (95, 0.1584), (96, 0.1667), (97, 0.096), (98, 0.2439), (100, 0.0982)],
    44 => [(96, 0.0554), (98, 0.0187), (99, 0.1276), (100, 0.126), (101, 0.1706), (102, 0.3155), (104, 0.1862)],
    45 => [(103, 1)],
    46 => [(102, 0.0102), (104, 0.1114), (105, 0.2233), (106, 0.2733), (108, 0.2646), (110, 0.1172)],
    47 => [(107, 0.51839), (109, 0.48161)],
    48 => [(106, 0.0125), (108, 0.0089), (110, 0.1249), (111, 0.128), (112, 0.2413), (113, 0.1222), (114, 0.2873), (116, 0.0749)],
    49 => [(113, 0.0429), (115, 0.9571)],
    50 => [(112, 0.0097), (114, 0.0066), (115, 0.0034), (116, 0.1454), (117, 0.0768), (118, 0.2422), (119, 0.0859), (120, 0.3258), (122, 0.0463), (124, 0.0579)],
    51 => [(121, 0.5721), (123, 0.4279)],
    52 => [(120, 0.0009), (122, 0.0255), (123, 0.0089), (124, 0.0474), (125, 0.0707), (126, 0.1884), (128, 0.3174), (130, 0.3408)],
    53 => [(127, 1)],
    54 => [(124, 0.000952), (126, 0.00089), (128, 0.019102), (129, 0.264006), (130, 0.04071), (131, 0.212324), (132, 0.269086), (134, 0.104357), (136, 0.088573)],
    55 => [(133, 1)],
    56 => [(130, 0.00106), (132, 0.00101), (134, 0.02417), (135, 0.06592), (136, 0.07854), (137, 0.11232), (138, 0.71698)],
    57 => [(138, 0.0008881), (139, 0.9991119)],
    58 => [(136, 0.00185), (138, 0.00251), (140, 0.8845), (142, 0.11114)],
    59 => [(141, 1)],
    60 => [(142, 0.27152), (143, 0.12174), (144, 0.23798), (145, 0.08293), (146, 0.17189), (148, 0.05756), (150, 0.05638)],
    62 => [(144, 0.0307), (147, 0.1499), (148, 0.1124), (149, 0.1382), (150, 0.0738), (152, 0.2675), (154, 0.2275)],
    63 => [(151, 0.4781), (153, 0.5219)],
    64 => [(152, 0.002), (154, 0.0218), (155, 0.148), (156, 0.2047), (157, 0.1565), (158, 0.2484), (160, 0.2186)],
    65 => [(159, 1)],
    66 => [(156, 0.00056), (158, 0.00095), (160, 0.02329), (161, 0.18889), (162, 0.25475), (163, 0.24896), (164, 0.2826)],
    67 => [(165, 1)],
    68 => [(162, 0.00139), (164, 0.01601), (166, 0.33503), (167, 0.22869), (168, 0.26978), (170, 0.1491)],
    69 => [(169, 1)],
    70 => [(168, 0.00123), (170, 0.02982), (171, 0.1409), (172, 0.2168), (173, 0.16103), (174, 0.32026), (176, 0.12996)],
    71 => [(175, 0.97401), (176, 0.02599)],
    72 => [(174, 0.0016), (176, 0.0526), (177, 0.186), (178, 0.2728), (179, 0.1362), (180, 0.3508)],
    73 => [(180, 0.0001201), (181, 0.9998799)],
    74 => [(180, 0.0012), (182, 0.265), (183, 0.1431), (184, 0.3064), (186, 0.2843)],
    75 => [(185, 0.374), (187, 0.626)],
    76 => [(184, 0.0002), (186, 0.0159), (187, 0.0196), (188, 0.1324), (189, 0.1615), (190, 0.2626), (192, 0.4078)],
    77 => [(191, 0.373), (193, 0.627)],
    78 => [(190, 0.00012), (192, 0.00782), (194, 0.3286), (195, 0.3378), (196, 0.2521), (198, 0.07356)],
    79 => [(197, 1)],
    80 => [(196, 0.0015), (198, 0.0997), (199, 0.1687), (200, 0.231), (201, 0.1318), (202, 0.2986), (204, 0.0687)],
    81 => [(203, 0.2952), (205, 0.7048)],
    82 => [(204, 0.014), (206, 0.241), (207, 0.221), (208, 0.524)],
    83 => [(209, 1)],
    90 => [(232, 1)],
    91 => [(231, 1)],
    92 => [(234, 5.4e-05), (235, 0.007204), (238, 0.992742)]
)