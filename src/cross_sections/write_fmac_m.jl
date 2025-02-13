"""
    write_fmac_m(cross_sections::Cross_Sections)

Write FMAC-M cross sections file.

# Input Argument(s)
- 'cross_sections::Cross_Sections': cross sections informations.

# Output Argument(s)
N/A

# Reference(s)
N/A

"""
function write_fmac_m(cross_sections::Cross_Sections,fmac_file::String)
    

# Write FMAC-M line per line
open(fmac_file,"w") do file

    # Record 1. File identification
    println(file,"FMAC-M                       5")

    # Record 2. File label
    println(file,"\n")

    # Record 3. Integer parameters
    v0 = zeros(Int64,10)
    v0[1] = 1
    print_int(file,v0)

    # Record 4. Character parameters
    println(file,"NGCAP")

    #----
    # Record II + 1 : Integer parameters
    #----
    v1 = zeros(Int64,40)
    v1[1] = sum(cross_sections.get_number_of_groups())  # Cumulative number of groups
    v1[2] = cross_sections.get_number_of_particles()    # Number of particles
    v1[15] = cross_sections.get_number_of_materials()   # Number of materials
    v1[17] = 38                                         # Number of parameters in Record II + 2
    v1[20] = cross_sections.get_legendre_order() + 1    # Legendre order + 1
    print_int(file,v1)

    #----
    # Record II + 2 : Control integer that define the length of each Record III + i
    #----
    v2 = zeros(Int64,38)
    v2[1] = 12              # Information text
    v2[2] = v1[2]           # Number of energy groups per particle type
    v2[3] = v1[2]           # Particle names
    v2[5] = v1[1] + v1[2]   # Energy mesh boundaries
    v2[15] = v1[1]          # Total cross-sections
    v2[35] = v1[1] + v1[2]  # Stopping powers
    v2[36] = v1[1]          # Momentum transfers
    v2[37] = v1[1]          # Energy deposition cross-sections
    v2[38] = v1[1]          # Charge deposition cross-sections
    print_int(file,v2)

    #----
    # Record III + 1 : Information text
    #----
    if v2[1] != 0
        println(file,"Coupled multigroup cross-sections produced with Radiant.jl")
    end

    #----
    # Record III + 2 : Number of groups of energy mesh by particle type
    #----
    if v2[2] != 0
        number_of_groups = cross_sections.get_number_of_groups()
        print_int(file,number_of_groups)
    end

    #----
    # Record III + 3 : Particle names by particle types
    #----
    if v2[3] != 0
        Npart = v1[2]
        particle_names = cross_sections.get_particles()
        for i in range(1,Npart)
            if is_electron(particle_names[i])
                fmac_particles_names = "BETA"
            elseif is_photon(particle_names[i])
                fmac_particles_names = "GAMA"
            elseif is_positron(particle_names[i])
                fmac_particles_names = "POSITR"
            else
                error("Unkown particle.")
            end
            print(file,string(fmac_particles_names,"  "))
        end
        println(file)
    end

    #----
    # Record III + 5 : Energy mesh boundaries
    #----
    if v2[5] != 0
        Npart = v1[2]
        particle_names = cross_sections.get_particles()
        mesh_boundaries = Vector{Float64}()
        for i in range(1,Npart)
            append!(mesh_boundaries,cross_sections.get_energy_boundaries(particle_names[i])[1:end-1])
        end
        for i in range(1,Npart)
            append!(mesh_boundaries,cross_sections.get_energy_boundaries(particle_names[Npart-i+1])[end])
        end
        print_float(file,mesh_boundaries)
    end

    #----
    # Record III + 15 : Total cross-sections
    #----
    if v2[15] != 0
        Npart = v1[2]
        Nmat = v1[15]
        particle_names = cross_sections.get_particles()
        total_cross_sections = Vector{Float64}()
        for i in range(1,Npart)
            tcs = cross_sections.get_total(particle_names[i])
            for n in range(1,Nmat)
                append!(total_cross_sections,tcs[:,n])
            end
        end
        print_float(file,total_cross_sections)
    end

    #----
    # Record III + 35 : Stopping powers
    #----
    if v2[35] != 0
        Npart = v1[2]
        Nmat = v1[15]
        particle_names = cross_sections.get_particles()
        stopping_powers = Vector{Float64}()
        for i in range(1,Npart)
            sp = cross_sections.get_boundary_stopping_powers(particle_names[i])
            for n in range(1,Nmat)
                append!(stopping_powers,sp[1:end-1,n])
            end
        end
        for i in range(1,Npart)
            sp = cross_sections.get_boundary_stopping_powers(particle_names[Npart-i+1])
            for n in range(1,Nmat)
                append!(stopping_powers,sp[end,n])
            end
        end
        print_float(file,stopping_powers)
    end

    #----
    # Record III + 36 : Momentum transfers
    #----
    if v2[36] != 0
        Npart = v1[2]
        Nmat = v1[15]
        particle_names = cross_sections.get_particles()
        momentum_transfer = Vector{Float64}()
        for i in range(1,Npart)
            mt = cross_sections.get_momentum_transfer(particle_names[i])
            for n in range(1,Nmat)
                append!(momentum_transfer,mt[:,n])
            end
        end
        print_float(file,momentum_transfer)
    end

    #----
    # Record III + 37 : Energy deposition cross-sections
    #----
    if v2[37] != 0
        Npart = v1[2]
        Nmat = v1[15]
        particle_names = cross_sections.get_particles()
        energy_deposition = Vector{Float64}()
        for i in range(1,Npart)
            ed = cross_sections.get_energy_deposition(particle_names[i])
            for n in range(1,Nmat)
                append!(energy_deposition,ed[1:end-1,n])
            end
        end
        print_float(file,energy_deposition)
    end

    #----
    # Record III + 38 : Charge deposition cross-sections
    #----
    if v2[38] != 0
        Npart = v1[2]
        Nmat = v1[15]
        particle_names = cross_sections.get_particles()
        charge_deposition = Vector{Float64}()
        for i in range(1,Npart)
            cd = cross_sections.get_charge_deposition(particle_names[i])
            for n in range(1,Nmat)
                append!(charge_deposition,cd[1:end-1,n])
            end
        end
        print_float(file,charge_deposition)
    end

    #----
    # Record III + x > 39 : Legendre moments of scattering cross sections
    #----
    Npart = v1[2]
    Nmat = v1[15]
    particle_names = cross_sections.get_particles()
    Ng = cross_sections.get_number_of_groups()
    L = cross_sections.get_legendre_order()
    for j in range(1,Npart), gf in range(1,Ng[j]), i in range(1,Npart), gi in range(1,Ng[i])

        dcs = cross_sections.get_scattering(particle_names[i],particle_names[j],L)
        dcs_ij = dcs[:,gi,gf,:]
        
        # Integer parameters
        p = zeros(Int64,3+Nmat)
        p[1] = sum(Ng[1:i-1]) + gi
        p[2] = sum(Ng[1:j-1]) + gf
        for n in range(1,Nmat)
            if all(x -> x == 0.0,dcs_ij[n,:])
                p[2+n] = 0
                p[end] += 1
            else
                L_eff = L+1
                for i in range(0,L)
                    if (dcs_ij[n,L+1-i] != 0.0) L_eff = L-i; break; end
                end
                p[2+n] = -(L_eff+1)
                p[end] += L_eff+1
            end
        end
        print_int(file,p)

        # Legendre moments of the differential cross sections
        scat = Vector{Float64}()
        for n in range(1,Nmat)
            if p[2+n] == 0
                append!(scat,0.0)
            else
                L_eff = -p[2+n]-1
                append!(scat,dcs_ij[n,1:L_eff+1])
            end
        end
        print_float(file,scat)

    end

end 
end

"""
    print_int(file::IOStream,i::Int64)

Write a single integer to the FMAC-M file.

# Input Argument(s)
- 'file::IOStream' : file.
- 'i::Int64' : integer to write in the file.

# Output Argument(s)
N/A

# Reference(s)
N/A

"""
function print_int(file::IOStream,i::Int64)
    space_length = 12 - length(string(i))
    print(file,string(join(fill(" ",space_length)),string(i)))
end

"""
    print_float(file::IOStream,i::Float64)

Write a single real to the FMAC-M file.

# Input Argument(s)
- 'file::IOStream' : file.
- 'i::Float64' : real to write in the file.

# Output Argument(s)
N/A

# Reference(s)
N/A

"""
function print_float(file::IOStream,i::Float64)

    if (sign(i) == -1) sng = "-" else sng = " " end
    if i != 0 exponent = floor(Int64,log10(abs(i))) else exponent = 0 end
    if (sign(exponent) == -1) sng_exp = "-" else sng_exp = "+" end
    mantissa = abs(i)/10.0^exponent
    exponent = abs(exponent)

    if exponent < 10
        exponent_string = string("0",exponent)
    elseif exponent < 99
        exponent_string = string(exponent)
    else
        if sng_exp == 1
            error("Value is too big.")
        else
            exponent_string = string("00")
            mantissa = 0.0
            sng = " "
        end
    end
    mantissa_string = string(round.(mantissa; digits=5))
    mantissa_length = 7 - length(mantissa_string)
    if mantissa_length > 0
        print(file,string(sng,mantissa_string,join(fill("0",mantissa_length)),"E",sng_exp,exponent_string))
    else
        print(file,string(sng,mantissa_string,"E",sng_exp,exponent_string))
    end

end

"""
    print_int(file::IOStream,v::Vector{Int64})

Write a vector of integer to the FMAC-M file.

# Input Argument(s)
- 'file::IOStream' : file.
- 'v::Vector{Int64}' : integer vector to write in the file.

# Output Argument(s)
N/A

# Reference(s)
N/A

"""
function print_int(file::IOStream,v::Vector{Int64})
    Nint = 0
    for i in range(1,length(v))
        print_int(file,v[i])
        Nint += 1
        if mod(Nint,6) == 0 println(file) end
    end
    if mod(Nint,6) != 0 println(file) end
end

"""
    print_float(file::IOStream,v::Vector{Float64})

Write a vector of real to the FMAC-M file.

# Input Argument(s)
- 'file::IOStream' : file.
- 'v::Vector{Float64}' : real vector to write in the file.

# Output Argument(s)
N/A

# Reference(s)
N/A

"""
function print_float(file::IOStream,v::Vector{Float64})
    Nint = 0
    for i in range(1,length(v))
        print_float(file,v[i])
        Nint += 1
        if mod(Nint,6) == 0 println(file) end
    end
    if mod(Nint,6) != 0 println(file) end
end