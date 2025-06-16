"""
    baro(Z::Int64,Ei::Float64,Ef::Float64,F₀_type::String="baro",
    is_triplet_production::Bool=true)

Gives the differential pair production cross-section based on the model of Baró et al.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `Ei::Float64` : incoming photon energy.
- `Ef::Float64` : outgoing lepton energy.
- `F₀_type::String` : type of correction factor, which can takes the following values:
    - `"baro"` : correction factor approximation from Baró et al.
    - `"epdl97"` : interpolation based on the total cross-section from EPDL97 library.
    - `"none"` : no correction factor. 
- `is_triplet_production::Bool` : boolean to consider or not triplet production.

# Output Argument(s)
- `σs::Float64` : pair production differential cross-section.

# Reference(s)
- Tsai (1974), Pair production and bremsstrahlung of charged leptons.
- Baró et al. (1994), Analytical cross sections for Monte Carlo simulation of photon
  transport.

"""
function baro(Z::Int64,Ei::Float64,Ef::Float64,F₀_type::String="epdl97",is_triplet_production::Bool=true)

    #----
    # Initialization
    #----
    if Ei-2 < 0 return 0 end
    rₑ = 2.81794092e-13 # (in cm)
    α = 1/137
    a = α*Z
    rs, n∞ = baro_coefficient(Z)
    fc = high_energy_Coulomb_correction(Z)
    Cr = 1.0093 # High-energy limit of Mork and Olsen's radiative correction

    #----
    # Compute the triplet production factor
    #----
    if is_triplet_production && Ei > 4
        ν = (0.2840-0.1909*a)*log(4/Ei) + (0.1095+0.2206*a)*log(4/Ei)^2 + (0.02888-0.04269*a)*log(4/Ei)^3 + (0.002527+0.002623*a)*log(4/Ei)^4
        η = (1-exp(-ν))*n∞
    else
        η = 0
    end

    #----
    # Correcting factor to fit experimental values
    #----
    if F₀_type == "epdl97"
        F₀ = baro_correction_factor(Z,Ei)
    elseif F₀_type == "baro"
        F₀ = (-0.1774-12.10*a+11.18*a^2)*sqrt(2/Ei) + (8.523+73.26*a-44.41*a^2)*(2/Ei) - (13.52+121.1*a-96.41*a^2)*sqrt(2/Ei)^3 + (8.946+62.05*a-63.41*a^2)*(2/Ei)^2
    elseif F₀_type == "none"
        F₀ = 0.0
    else
        error("Unknow type of correction factor.")
    end

    #----
    # Compute the differential cross-section
    #----
    ϵ = (Ef+1)/Ei
    ϵ₀ = 1/Ei
    b = rs/2 * ϵ₀/(ϵ*(1-ϵ))
    g1 = 7/3 - 2*log(1+b^2) - 6*b*atan(1/b) - b^2*(4 - 4*b*atan(1/b) - 3*log(1+1/b^2))
    g2 = 11/6 - 2*log(1+b^2) - 3*b*atan(1/b) + b^2/2*(4 - 4*b*atan(1/b) - 3*log(1+1/b^2))
    g0 = 4*log(rs) - 4*fc + F₀
    ϕ₁ = g1 + g0
    ϕ₂ = g2 + g0
    if F₀_type != "none" ϕ₁ = max(ϕ₁,0); ϕ₂ = max(ϕ₂,0) end
    σs = 2/3 * α * rₑ^2 * Z*(Z+η) * Cr * (2*(1/2-ϵ)^2*ϕ₁+ϕ₂)/Ei

    return σs
end

"""
    baro_correction_factor(Z::Int64,Ei::Float64)

Compute the correction factor from Baró et al. using experimental data from EPDL97 library.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `Ei::Float64` : incoming photon energy.

# Output Argument(s)
- `F₀::Float64` : correction factor.

# Reference(s)
- Baró et al. (1994), Analytical cross sections for Monte Carlo simulation of photon
  transport.

"""
function baro_correction_factor(Z::Int64,Ei::Float64)

    # Extract data
    data = fast_load("pair_production_EPDL97.jld2")
    E_exp = data["E"][Z]
    σt_exp = data["σ"][Z]

    # Compute F₀ parameters from experimental data following Baró methodology
    if ~haskey(data,"F0") || ~haskey(data["F0"],Z)
        rₑ = 2.81794092e-13 # (in cm)
        α = 1/137
        a = α*Z
        rs, n∞ = baro_coefficient(Z)
        Cr = 1.0093 # High-energy limit of Mork and Olsen's radiative correction
        Np = 64
        N_exp = length(E_exp)
        F₀_exp = zeros(N_exp)
        for i in range(1,N_exp)
            Ei_exp = E_exp[i]
            σti_exp = σt_exp[i]
            σti_theory = 0

            if Ei_exp > 4
                ν = (0.2840-0.1909*a)*log(4/Ei_exp) + (0.1095+0.2206*a)*log(4/Ei_exp)^2 + (0.02888-0.04269*a)*log(4/Ei_exp)^3 + (0.002527+0.002623*a)*log(4/Ei_exp)^4
                η = (1-exp(-ν))*n∞
            else
                η = 0
            end

            # Compute total cross-section
            Eout = log_energy_group_structure(80,Ei_exp,2)
            Ngf = length(Eout)-1
            Np = 8
            u,w = quadrature(Np,"gauss-legendre")
            for gf in range(1,Ngf+1)
                Ef⁻ = Eout[gf]
                if (gf != Ngf+1) Ef⁺ = Eout[gf+1] else Ef⁺ = 0.0 end
                Ef⁻ = min(Ef⁻,Ei_exp-2)
                ΔEf = Ef⁻ - Ef⁺
                if ΔEf > 0
                    for n in range(1,Np)
                        Ef = (u[n]*ΔEf + (Ef⁻+Ef⁺))/2
                        σs = baro(Z,Ei_exp,Ef,"none",true)
                        σti_theory += ΔEf/2 * w[n] * σs
                    end
                end
            end

            # Compute F₀ that fit the experimental data
            A = 2/3 * α * rₑ^2 * Z*(Z+η) * Cr
            F₀_exp[i] = 6*Ei_exp^3/(7*Ei_exp^3-18*Ei_exp^2+12*Ei_exp-8) * (σti_exp - σti_theory)/A

        end
        if ~haskey(data,"F0")
            data["F0"] = Dict()
        end
        data["F0"][Z] = cubic_hermite_spline(E_exp,F₀_exp)
        cache_radiant[]["pair_production_EPDL97"] = data
    end

    # Interpolate F₀ parameters
    F₀_spline = data["F0"][Z]
    F₀ = F₀_spline(Ei)
    return F₀
end

"""
    baro_coefficient(Z::Int64)

Extract the screening parameter (rs) and the asymptotic triplet contribution factor (n∞) 
for Baró cross-sections calculations.

# Input Argument(s)
- `Z::Int64`: atomic number of the element.

# Output Argument(s)
- `rs::Float64`: screening parameter.
- `n∞::Float64`: asymptotic triplet contribution factor.

# Reference(s)
- Baró et al. (1994), Analytical cross sections for Monte Carlo simulation of photon
  transport.
- Salvat (2019), PENELOPE-2018: A Code System for Monte Carlo Simulation of Electron and
  Photon Transport.

"""
function baro_coefficient(Z::Int64)
  if haskey(baro_coefficient_dict, Z)
      return baro_coefficient_dict[Z]
  else
      error("The input atomic number Z is invalid (Z ∈ {1,99}).")
  end
end

const baro_coefficient_dict = Dict(
    1 => (122.81, 1.157),
    2 => (73.167, 1.169),
    3 => (69.228, 1.219),
    4 => (67.301, 1.201),
    5 => (64.696, 1.189),
    6 => (61.228, 1.174),
    7 => (57.524, 1.176),
    8 => (54.033, 1.169),
    9 => (50.787, 1.163),
    10 => (47.851, 1.157),
    11 => (46.373, 1.174),
    12 => (45.401, 1.183),
    13 => (44.503, 1.186),
    14 => (43.815, 1.184),
    15 => (43.074, 1.180),
    16 => (42.321, 1.178),
    17 => (41.586, 1.175),
    18 => (40.953, 1.170),
    19 => (40.524, 1.180),
    20 => (40.256, 1.187),
    21 => (39.756, 1.184),
    22 => (39.144, 1.180),
    23 => (38.462, 1.177),
    24 => (37.778, 1.166),
    25 => (37.174, 1.169),
    26 => (36.663, 1.166),
    27 => (35.986, 1.164),
    28 => (35.317, 1.162),
    29 => (34.688, 1.154),
    30 => (34.197, 1.156),
    31 => (33.786, 1.157),
    32 => (33.422, 1.158),
    33 => (33.068, 1.157),
    34 => (32.740, 1.158),
    35 => (32.438, 1.158),
    36 => (32.143, 1.158),
    37 => (31.884, 1.166),
    38 => (31.622, 1.173),
    39 => (31.438, 1.174),
    40 => (31.142, 1.175),
    41 => (30.950, 1.170),
    42 => (30.758, 1.169),
    43 => (30.561, 1.172),
    44 => (30.285, 1.169),
    45 => (30.097, 1.168),
    46 => (29.832, 1.164),
    47 => (29.581, 1.167),
    48 => (29.411, 1.170),
    49 => (29.247, 1.172),
    50 => (29.085, 1.174),
    51 => (28.930, 1.175),
    52 => (28.721, 1.178),
    53 => (28.580, 1.179),
    54 => (28.442, 1.180),
    55 => (28.312, 1.187),
    56 => (28.139, 1.194),
    57 => (27.973, 1.197),
    58 => (27.819, 1.196),
    59 => (27.675, 1.194),
    60 => (27.496, 1.194),
    61 => (27.285, 1.194),
    62 => (27.093, 1.194),
    63 => (26.911, 1.194),
    64 => (26.705, 1.196),
    65 => (26.516, 1.197),
    66 => (26.304, 1.196),
    67 => (26.108, 1.197),
    68 => (25.929, 1.197),
    69 => (25.730, 1.198),
    70 => (25.577, 1.198),
    71 => (25.403, 1.200),
    72 => (25.245, 1.201),
    73 => (25.100, 1.202),
    74 => (24.941, 1.204),
    75 => (24.790, 1.205),
    76 => (24.655, 1.206),
    77 => (24.506, 1.208),
    78 => (24.391, 1.207),
    79 => (24.262, 1.208),
    80 => (24.145, 1.212),
    81 => (24.039, 1.215),
    82 => (23.922, 1.218),
    83 => (23.813, 1.221),
    84 => (23.712, 1.224),
    85 => (23.621, 1.227),
    86 => (23.523, 1.230),
    87 => (23.430, 1.237),
    88 => (23.331, 1.243),
    89 => (23.238, 1.247),
    90 => (23.139, 1.250),
    91 => (23.048, 1.251),
    92 => (22.967, 1.252),
    93 => (22.833, 1.255),
    94 => (22.694, 1.256),
    95 => (22.624, 1.257),
    96 => (22.545, 1.259),
    97 => (22.446, 1.262),
    98 => (22.358, 1.262),
    99 => (22.264, 1.265)
)