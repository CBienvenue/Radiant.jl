"""
    generate_cross_sections(cross_sections::Cross_Sections)

Generate cross sections and save their content in a Cross_Sections structure.

# Input Argument(s)
- `cross_sections::Cross_Sections`: cross sections informations.

# Output Argument(s)
N/A

# Reference(s)
N/A

"""
function generate_cross_sections(cross_sections::Cross_Sections)

#----
# Initialization
#----
L = cross_sections.get_legendre_order()
group_structure = cross_sections.get_group_structure()
Nmat = cross_sections.get_number_of_materials()
Npart = cross_sections.get_number_of_particles()
materials = cross_sections.get_materials()
particles = cross_sections.get_particles()
interactions = cross_sections.get_interactions()
ρ = zeros(Nmat)
state_of_matter = Vector{String}()
Z = Vector{Vector{Int64}}(undef,Nmat)
ωz = Vector{Vector{Float64}}(undef,Nmat)
for n in range(1,Nmat)
    ρ[n] = materials[n].get_density()
    Z[n] = materials[n].get_atomic_numbers()
    ωz[n] = materials[n].get_weight_fractions()
    push!(state_of_matter,materials[n].get_state_of_matter())
end

#----
# Generate the group structure
#----
E₁ = Vector{Float64}(undef,Npart)
Ec = Vector{Float64}(undef,Npart)
Ng = Vector{Int64}(undef,Npart)
Eᵇ = Vector{Vector{Float64}}(undef,Npart)
for n in range(1,Npart)
    if haskey(group_structure,particles[n])
        Eᵇ[n] = group_structure[particles[n]]
    elseif haskey(group_structure,"default")
        Eᵇ[n] = group_structure["default"]
    else
        error("Group structure is not defined for $(particles[n].type).")
    end
    E₁[n] = (Eᵇ[n][1]+Eᵇ[n][2])/2
    Ec[n] = Eᵇ[n][end]
    Ng[n] = length(Eᵇ[n]) - 1
end

#----
# Deal with interaction interdependances
#----
interaction_interdependances(interactions,particles)

#----
# Multigroup cross-sections production
#----
multigroup_cross_sections = Array{Multigroup_Cross_Sections}(undef,Npart,Nmat)
for i in range(1,Npart), n in range(1,Nmat)
    mcs = Multigroup_Cross_Sections(Ng[i])
    Σt = zeros(Ng[i]); Σa = zeros(Ng[i]); Σe = zeros(Ng[i]+1); Σc = zeros(Ng[i]+1); Sb = zeros(Ng[i]+1); S = zeros(Ng[i]); T = zeros(Ng[i]);
    for j in range(1,Npart)
        Σsl = zeros(Ng[i],Ng[j],L+1)
        for interaction in interactions
            for pin in interaction.get_in_particles(), pout in interaction.get_out_particles()
                if pin != get_type(particles[i]) || pout != get_type(particles[j]) continue end
                for type in interaction.get_types(pin,pout)
                    println("\n Interaction: $(typeof(interaction)) | Type: $type | Incoming particle: $(get_type(particles[i])) | Outgoing particle: $(get_type(particles[j]))")
                    @time Σsli, Σti, Σai, Σei, Σci, Sbi, Si, Ti = multigroup(Z[n],ωz[n],ρ[n],state_of_matter[n],Eᵇ[i],Eᵇ[j],L,interaction,type,particles[i],particles[j],particles,interactions)
                    Σsl .+= Σsli; Σt .+= Σti; Σa .+= Σai; Σe .+= Σei; Σc .+= Σci; Sb .+= Sbi; S .+= Si; T .+= Ti;
                end
            end
        end
        mcs.set_scattering(Σsl)
    end
    mcs.set_total(Σt)
    mcs.set_absorption(Σa)
    mcs.set_boundary_stopping_powers(Sb)
    mcs.set_stopping_powers(S)
    mcs.set_momentum_transfer(T)
    mcs.set_energy_deposition(Σe)
    mcs.set_charge_deposition(Σc)
    multigroup_cross_sections[i,n] = mcs
end

cross_sections.set_energy(E₁)
cross_sections.set_cutoff(Ec)
cross_sections.set_number_of_groups(Ng)
cross_sections.set_multigroup_cross_sections(multigroup_cross_sections)
cross_sections.set_energy_boundaries(Eᵇ)

end