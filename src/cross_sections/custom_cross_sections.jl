function custom_cross_sections(cross_sections::Cross_Sections)

    Nmat = cross_sections.get_number_of_materials()

    multigroup_cross_sections = Array{Multigroup_Cross_Sections}(undef,1,Nmat)
    for i in range(1,1), n in range(1,Nmat)
        mcs = Multigroup_Cross_Sections(1)
        Σt = zeros(1); Σa = zeros(1); Σs = zeros(1); Σe = zeros(1+1); Σc = zeros(1+1); S = zeros(1+1); α = zeros(1);
        Σsℓ = zeros(1,1,0+1)

        Σt[1] = cross_sections.custom_absorption[n] + cross_sections.custom_scattering[n]
        Σsℓ[1,1,1] = cross_sections.custom_scattering[n]

        mcs.set_scattering(Σsℓ)
        mcs.set_total(Σt)
        mcs.set_absorption(Σa)
        mcs.set_stopping_powers(S)
        mcs.set_momentum_transfer(α)
        mcs.set_energy_deposition(Σe)
        mcs.set_charge_deposition(Σc)
        multigroup_cross_sections[i,n] = mcs
    end

    cross_sections.set_number_of_groups(1)
    cross_sections.set_multigroup_cross_sections(multigroup_cross_sections)

end