"""
    mean_excitation_energy(Z::Int64)

Extract the mean excitation energy for atomic number Z ∈ {1,100}.

# Input Argument(s)
- `Z::Int64`: atomic number of the element.

# Output Argument(s)
- `I::Float64`: mean excitation energy [in mₑc²].

# Reference(s)
- Seltzer and Berger (1981) : Evaluation of the Collision Stopping Power of Elements and
  Compounds for Electrons and Positrons.

"""
function mean_excitation_energy(Z::Int64)

    # Define the mean excitation energies (in eV)
    mean_excitation_energies = Dict(
        1 => 19.2, 2 => 41.8, 3 => 40.0, 4 => 63.7, 5 => 76.0,
        6 => 78.0, 7 => 82.0, 8 => 95.0, 9 => 115.0, 10 => 137.0,
        11 => 149.0, 12 => 156.0, 13 => 166.0, 14 => 173.0, 15 => 173.0,
        16 => 180.0, 17 => 174.0, 18 => 188.0, 19 => 190.0, 20 => 191.0,
        21 => 216.0, 22 => 233.0, 23 => 245.0, 24 => 257.0, 25 => 272.0,
        26 => 286.0, 27 => 297.0, 28 => 311.0, 29 => 322.0, 30 => 330.0,
        31 => 334.0, 32 => 350.0, 33 => 347.0, 34 => 348.0, 35 => 343.0,
        36 => 352.0, 37 => 363.0, 38 => 366.0, 39 => 379.0, 40 => 393.0,
        41 => 417.0, 42 => 424.0, 43 => 428.0, 44 => 441.0, 45 => 449.0,
        46 => 470.0, 47 => 470.0, 48 => 469.0, 49 => 488.0, 50 => 488.0,
        51 => 487.0, 52 => 485.0, 53 => 474.0, 54 => 482.0, 55 => 488.0,
        56 => 491.0, 57 => 501.0, 58 => 523.0, 59 => 535.0, 60 => 546.0,
        61 => 560.0, 62 => 574.0, 63 => 580.0, 64 => 591.0, 65 => 614.0,
        66 => 628.0, 67 => 650.0, 68 => 658.0, 69 => 674.0, 70 => 684.0,
        71 => 694.0, 72 => 705.0, 73 => 718.0, 74 => 727.0, 75 => 736.0,
        76 => 746.0, 77 => 757.0, 78 => 790.0, 79 => 790.0, 80 => 800.0,
        81 => 810.0, 82 => 823.0, 83 => 823.0, 84 => 830.0, 85 => 825.0,
        86 => 794.0, 87 => 827.0, 88 => 826.0, 89 => 841.0, 90 => 847.0,
        91 => 878.0, 92 => 890.0, 93 => 902.0, 94 => 921.0, 95 => 934.0,
        96 => 939.0, 97 => 952.0, 98 => 966.0, 99 => 980.0, 100 => 994.0
    )

    # Check if Z exists in the dictionary
    if haskey(mean_excitation_energies, Z)
        mₑc² = 0.510999
        return mean_excitation_energies[Z] / (1e6 * mₑc²)
    else
        error("The input atomic number Z is invalid (Z ∉ {1,100}).")
    end
end

function effective_mean_excitation_energy(Z::Vector{Int64},ωz::Vector{Float64},I_override::Float64=NaN)

    # 1) Explicit per-material override (in mₑc²), set via Material.set_mean_excitation_energy
    if !isnan(I_override)
        Ieff = I_override

    else
        # 2) Canonical compound from the tabulated ICRU values
        Itab = tabulated_mean_excitation_energy(Z,ωz)
        if !isnan(Itab)
            Ieff = Itab

        # 3) General case : Bragg additivity rule
        else
            ZoverA = ωz .* Z ./ standard_atomic_weight.(Z)
            Ieff = exp(  sum(ZoverA .* log.(mean_excitation_energy.(Z)))  /  sum(ZoverA)  )
        end
    end

    return Ieff
end

"""
    tabulated_mean_excitation_energy(Z::Vector{Int64},ωz::Vector{Float64})

Return the ICRU-37/44 recommended mean excitation energy (in mₑc²) for a few canonical
compositions (liquid water, dry air, ...), or `NaN` if the composition does not match a
tabulated one. For any other material, the Bragg additivity rule applies, and a precise
value can be imposed with `Material.set_mean_excitation_energy` (recommended for tissues,
e.g. striated muscle ≈ 74.7 eV, whose CT-derived composition does not match a template).
"""
function tabulated_mean_excitation_energy(Z::Vector{Int64},ωz::Vector{Float64})
    mₑc² = 0.510999
    # Water keeps its original (lenient) detection and value (78 eV)
    if is_water(Z,ωz) return 78.0/(1e6*mₑc²) end
    # (atomic number => mass fraction) templates, with mean excitation energy in eV
    templates = (
        (Dict(6=>0.000124, 7=>0.755267, 8=>0.231781, 18=>0.012827), 85.7),  # dry air (sea level)
    )
    for (comp,I_eV) in templates
        if matches_composition(Z,ωz,comp) return I_eV/(1e6*mₑc²) end
    end
    return NaN
end

"""
    matches_composition(Z,ωz,comp;rtol,atol)

Test whether the material `(Z,ωz)` matches a tabulated composition `comp` (atomic number
=> mass fraction): identical element set and weight fractions within tolerance.
"""
function matches_composition(Z::Vector{Int64},ωz::Vector{Float64},comp::Dict{Int64,Float64};rtol::Float64=0.02,atol::Float64=0.005)
    if Set(Z) != Set(keys(comp)) return false end
    for (z,w_ref) in comp
        w = ωz[findfirst(==(z),Z)]
        if abs(w - w_ref) > atol + rtol*w_ref return false end
    end
    return true
end

function is_water(Z::Vector{Int64},ωz::Vector{Float64})
    # Check if Z contain H and O and if weight fraction correspond
    return any(x -> x == 1, Z) && any(x -> x == 8, Z) && (abs(ωz[Z .== 1][1] - 0.1119) < 0.01 && abs(ωz[Z .== 8][1] - 0.8881) < 0.01)
end