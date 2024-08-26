"""
    generate_cross_sections(cross_sections::Cross_Sections)

Generate cross sections and save their content in a Cross_Sections structure.

# Input Argument(s)
- 'cross_sections::Cross_Sections': cross sections informations.

# Output Argument(s)
- 'cross_sections::Cross_Sections': cross sections with actualized informations.

# Reference(s)
N/A

"""
function generate_cross_sections(cross_sections::Cross_Sections)

#----
# Initialization
#----
Ng = cross_sections.get_number_of_groups()
E₁ = cross_sections.get_energy()
Ec = cross_sections.get_cutoff()
L = cross_sections.get_legendre_order()
group_structure = cross_sections.get_group_structure()
Nmat = cross_sections.get_number_of_materials()
Npart = cross_sections.get_number_of_particles()
materials = cross_sections.get_materials()
particles = cross_sections.get_particles()
interactions = cross_sections.get_interactions()
material_names = Vector{String}()
ρ = zeros(Nmat)
state_of_matter = Vector{String}()
Z = Vector{Vector{Int64}}(undef,Nmat)
ωz = Vector{Vector{Float64}}(undef,Nmat)
for n in range(1,Nmat)
    ρ[n] = materials[n].get_density()
    Z[n] = materials[n].get_atomic_numbers()
    ωz[n] = materials[n].get_weight_fractions()
    push!(material_names,materials[n].get_name())
    push!(state_of_matter,materials[n].get_state_of_matter())
end

#----
# Generate the group structure
#----
Eᵇ = Vector{Vector{Float64}}(undef,Npart)
for n in range(1,Npart)
    Eᵇ[n] = energy_group_structure(Ng[n],E₁,Ec,group_structure[n])
end

#----
# Multigroup cross-sections production
#----
multigroup_cross_sections = Array{Multigroup_Cross_Sections}(undef,Npart,Nmat)
for i in range(1,Npart), n in range(1,Nmat)
    mcs = Multigroup_Cross_Sections(Ng[i])
    Σt = zeros(Ng[i]); Σa = zeros(Ng[i]); Σs = zeros(Ng[i]); Σe = zeros(Ng[i]); Σc = zeros(Ng[i]); S = zeros(Ng[i]+1); α = zeros(Ng[i]);
    for j in range(1,Npart)
        Σsℓ = zeros(Ng[i],Ng[j],L+1)
        for interaction in interactions
            for pin in interaction.get_in_particles(), pout in interaction.get_out_particles()
                if pin != particles[i] || pout != particles[j] continue end
                for type in interaction.get_types(pin,pout)
                    println([type,particles[i],particles[j]])
                    Σsℓi, Σti, Σai, Σsi, Σei, Σci, Si, αi = multigroup(Z[n],ωz[n],ρ[n],state_of_matter[n],Eᵇ[i],Eᵇ[j],L,interaction,type,pin,pout,particles,8,true,interactions)
                    Σsℓ .+= Σsℓi; Σt .+= Σti; Σa .+= Σai; Σs .+= Σsi; Σe .+= Σei; Σc .+= Σci; S .+= Si; α .+= αi;
                end
            end
        end
        mcs.set_scattering(Σsℓ)
    end
    mcs.set_total(Σt)
    mcs.set_absorption(Σa)
    mcs.set_stopping_powers(S)
    mcs.set_momentum_transfer(α)
    mcs.set_energy_deposition(Σe)
    mcs.set_charge_deposition(Σc)
    multigroup_cross_sections[i,n] = mcs
end

cross_sections.set_energy(E₁)
cross_sections.set_cutoff(Ec)
cross_sections.set_number_of_groups(Ng)
cross_sections.set_legendre_order(L)
cross_sections.set_multigroup_cross_sections(multigroup_cross_sections)
cross_sections.set_energy_boundaries(Eᵇ)

end