"""
    electron_cascades(ΔE_auger::Vector{Vector{Float64}},
    ΔE_fluorescence::Vector{Vector{Float64}},η_auger::Vector{Vector{Float64}},
    η_fluorescence::Vector{Vector{Float64}},ΔE⁻::Vector{Float64},η⁻::Vector{Float64},
    primary_shell⁻::Vector{String},secondary_shell⁻::Vector{String},
    tertiary_shell⁻::Vector{String},ΔE::Vector{Float64},η::Vector{Float64},
    PS::Vector{String},SS::Vector{String},TS::Vector{String},Nshells::Int64,
    subshells::Vector{String},PS₀::Vector{String},ηmin::Float64)

Identify atomic transition and accumulate the probability and energy of the produced
radiation they generate.

# Input Argument(s)
- 'ΔE_auger::Vector{Vector{Float64}}' : energies of produced Auger electrons associated to
  a vacancy in each energy shell.
- 'ΔE_fluorescence::Vector{Vector{Float64}}' : energies of produced fluorescence photons
  associated to a vacancy in each energy shell.
- 'η_auger::Vector{Vector{Float64}}' : probability of produced Auger electrons associated
  to a vacancy in each energy shell.
- 'η_fluorescence::Vector{Vector{Float64}}' : probability of produced fluorescence photons
  associated to a vacancy in each energy shell.
- 'ΔE⁻::Vector{Float64}' : energy transition of the cascade.
- 'η⁻::Vector{Float64}' : transition probability of the cascade.
- 'primary_shell⁻::Vector{String}' : subshell which contain a vacancy.
- 'secondary_shell⁻::Vector{String}' : subshell which electron is filling the vacancy.
- 'tertiary_shell⁻::Vector{String}' : subshell which electron is ejected (Auger electron) 
  following the filling of the vacancy.
- 'ΔE::Vector{Float64}' : energy transition of the following cascade.
- 'η::Vector{Float64}' : transition probability of the following cascade.
- 'PS::Vector{String}' : subshell which contain a vacancy.
- 'SS::Vector{String}' : subshell which electron is filling the vacancy.
- 'TS::Vector{String}' : subshell which electron is ejected (Auger electron) 
  following the filling of the vacancy.
- 'Nshells::Int64' : Number of subshells.
- 'subshells::Vector{String}' : List of subshells in the atom.
- 'PS₀::Vector{String}' : subshell which contain the first vacancy (before cascades).
- 'ηmin::Float64' : minimum probability of the particle production to be considered.

# Output Argument(s)
- 'ΔE_auger::Vector{Vector{Float64}}' : energies of produced Auger electrons associated to
  a vacancy in each energy shell.
- 'ΔE_fluorescence::Vector{Vector{Float64}}' : energies of produced fluorescence photons
  associated to a vacancy in each energy shell.
- 'η_auger::Vector{Vector{Float64}}' : probability of produced Auger electrons associated
  to a vacancy in each energy shell.
- 'η_fluorescence::Vector{Vector{Float64}}' : probability of produced fluorescence photons
  associated to a vacancy in each energy shell.

# Reference(s)
N/A

"""
function electron_cascades(ΔE_auger::Vector{Vector{Float64}},ΔE_fluorescence::Vector{Vector{Float64}},η_auger::Vector{Vector{Float64}},η_fluorescence::Vector{Vector{Float64}},ΔE⁻::Vector{Float64},η⁻::Vector{Float64},primary_shell⁻::Vector{String},secondary_shell⁻::Vector{String},tertiary_shell⁻::Vector{String},ΔE::Vector{Float64},η::Vector{Float64},PS::Vector{String},SS::Vector{String},TS::Vector{String},Nshells::Int64,subshells::Vector{String},PS₀::Vector{String},ηmin::Float64)

    Nt⁻ = length(ΔE⁻)
    Nt = length(ΔE)
    @inbounds for i in range(1,Nt⁻)

        # Accumulate the probability that an initial transition i results in the production of a specific transition in an outer shell
        for δi in range(1,Nshells)
            if PS₀[i] == subshells[δi] && ηmin < η⁻[i]
                # Accumulate fluorescence production
                if TS[i] == ""
                    push!(ΔE_fluorescence[δi],ΔE⁻[i])
                    push!(η_fluorescence[δi],η⁻[i])
                # Accumulate Auger electron production
                else
                    push!(ΔE_auger[δi],ΔE⁻[i])
                    push!(η_auger[δi],η⁻[i])
                end
            end
        end

        # Compute the next generation of cascades
        Nt⁺ = 0
        ΔE⁺ = Vector{Float64}()
        η⁺ = Vector{Float64}()
        PS⁺ = Vector{String}()
        SS⁺ = Vector{String}()
        TS⁺ = Vector{String}()
        PS₀⁺ = Vector{String}()
        for j in range(1,Nt)
            if secondary_shell⁻[i] == PS[j] && ηmin < η⁻[i]
                push!(ΔE⁺,ΔE[j])
                push!(η⁺,η[j]*η⁻[i])
                push!(PS⁺,PS[j])
                push!(SS⁺,SS[j])
                push!(TS⁺,TS[j])
                push!(PS₀⁺,PS₀[i])
                Nt⁺ += 1
            end
            if tertiary_shell⁻[i] == PS[j] && ηmin < η⁻[i]
                push!(ΔE⁺,ΔE[j])
                push!(η⁺,η[j]*η⁻[i])
                push!(PS⁺,PS[j])
                push!(SS⁺,SS[j])
                push!(TS⁺,TS[j])
                push!(PS₀⁺,PS₀[i])
                Nt⁺ += 1
            end
        end

        # If there is additionnal relaxation available with outer-shells electron filling the remaining vacancies, then continue, otherwise, it is the end
        if Nt⁺ != 0
            ΔE_auger,ΔE_fluorescence,η_auger,η_fluorescence = electron_cascades(ΔE_auger,ΔE_fluorescence,η_auger,η_fluorescence,ΔE⁺,η⁺,PS⁺,SS⁺,TS⁺,ΔE,η,PS,SS,TS,Nshells,subshells,PS₀⁺,ηmin)
        end
    end
    return ΔE_auger,ΔE_fluorescence,η_auger,η_fluorescence
end

"""
    atomic_electron_cascades(type::String,Z::Vector{Int64},Ecutoff::Float64,
    ηmin::Float64=0.001)

Accumulate the probability and energy of the radiation produced by vacancies in the atomic
structure for each atom in the material.

# Input Argument(s)
- 'type::String' : type of the radiation (auger or fluorescence).
- 'Z::Vector{Int64}' : atomic number of the elements in the material.
- 'Ecutoff::Float64' : cutoff energy.
- 'ηmin::Float64' : minimum probability of the production of specific Auger electrons
  following electron cascades.

# Output Argument(s)
- 'vec_ΔE::Vector{Vector{Vector{Float64}}}' : energy of the particles produced, by material
  and by a vacancy in each subshell. 
- 'vec_η::Vector{Vector{Vector{Float64}}}' : probability of the particles production, by
  material and by a vacancy in each subshell. 

# Reference(s)
N/A

"""
function atomic_electron_cascades(type::String,Z::Vector{Int64},Ecutoff::Float64,ηmin::Float64=0.001)

    # Initialization
    if type ∉ ["auger","fluorescence"] error("Type of radiation following electron cascades is either auger of fluorescence.") end
    Nz = length(Z)
    vec_ΔE = Vector{Vector{Vector{Float64}}}(undef,Nz)
    vec_η = Vector{Vector{Vector{Float64}}}(undef,Nz)
    path = joinpath(find_package_root(), "data", "relaxation_JENDL5.jld2")
    data = load(path)
    
    for iz in range(1,Nz)

        fluorescence = data["Fluorescence"][Z[iz]]
        auger = data["Auger"][Z[iz]]

        # Produce lists
        Nt = 0
        ΔE = Vector{Float64}()
        η = Vector{Float64}()
        PS = Vector{String}()
        SS = Vector{String}()
        TS = Vector{String}()
        Nshells,Zi,Ui,Ti,ri,subshells = electron_subshells(Z[iz])
        for δi in range(1,Nshells)
            fi = fluorescence[subshells[δi]]
            ai = auger[subshells[δi]]
            Nti = length(fi["ΔE"])
            for δj in range(1,Nti)
                Nt += 1
                push!(ΔE,fi["ΔE"][δj])
                push!(η,fi["P"][δj])
                push!(PS,subshells[δi])
                push!(SS,fi["Secondary_shells"][δj])
                push!(TS,"")
            end
            Nti = length(ai["ΔE"])
            for δj in range(1,Nti)
                Nt += 1
                push!(ΔE,ai["ΔE"][δj])
                push!(η,ai["P"][δj])
                push!(PS,subshells[δi])
                push!(SS,ai["Secondary_shells"][δj])
                push!(TS,ai["Tertiary_shells"][δj])
            end
        end

        # Extract electron cascades
        ΔE_auger_temp = Vector{Vector{Float64}}(undef,Nshells)
        η_auger_temp = Vector{Vector{Float64}}(undef,Nshells)
        ΔE_fluorescence_temp = Vector{Vector{Float64}}(undef,Nshells)
        η_fluorescence_temp = Vector{Vector{Float64}}(undef,Nshells)
        for δi in range(1,Nshells)
            ΔE_auger_temp[δi] = Vector{Float64}()
            ΔE_fluorescence_temp[δi] = Vector{Float64}()
            η_auger_temp[δi] = Vector{Float64}()
            η_fluorescence_temp[δi] = Vector{Float64}()
        end
        ΔE_auger_temp,ΔE_fluorescence_temp,η_auger_temp,η_fluorescence_temp = electron_cascades(ΔE_auger_temp,ΔE_fluorescence_temp,η_auger_temp,η_fluorescence_temp,ΔE,η,PS,SS,TS,ΔE,η,PS,SS,TS,Nshells,subshells,PS,ηmin)

        if type == "auger"
            ΔE_auger = Vector{Vector{Float64}}(undef,Nshells)
            η_auger = Vector{Vector{Float64}}(undef,Nshells)
            for δi in range(1,Nshells)
                ΔE_auger[δi] = Vector{Float64}()
                η_auger[δi] = Vector{Float64}()
                Nt = length(ΔE_auger_temp[δi])
                for i in range(1,Nt)
                    if ΔE_auger_temp[δi][i] > Ecutoff
                        push!(ΔE_auger[δi],ΔE_auger_temp[δi][i])
                        push!(η_auger[δi],η_auger_temp[δi][i])
                    end
                end
            end
            vec_ΔE[iz] = ΔE_auger
            vec_η[iz] = η_auger
        elseif type == "fluorescence"
            ΔE_fluorescence = Vector{Vector{Float64}}(undef,Nshells)
            η_fluorescence = Vector{Vector{Float64}}(undef,Nshells)
            for δi in range(1,Nshells)
                ΔE_fluorescence[δi] = Vector{Float64}()
                η_fluorescence[δi] = Vector{Float64}()
                Nt = length(ΔE_fluorescence_temp[δi])
                for i in range(1,Nt)
                    if ΔE_fluorescence_temp[δi][i] > Ecutoff
                        push!(ΔE_fluorescence[δi],ΔE_fluorescence_temp[δi][i])
                        push!(η_fluorescence[δi],η_fluorescence_temp[δi][i])
                    end
                end
            end
            vec_ΔE[iz] = ΔE_fluorescence
            vec_η[iz] = η_fluorescence
        end
    end

    return vec_ΔE, vec_η
end