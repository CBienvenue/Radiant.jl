"""
    read_fmac_m(cross_sections::Cross_Sections)

Read FMAC-M cross sections file and save their content in a Cross_Sections structure.

# Input Argument(s)
- 'cross_sections::Cross_Sections': cross sections informations.

# Output Argument(s)
N/A

# Reference(s)
N/A

"""
function read_fmac_m(cross_sections::Cross_Sections)

#----
# Extract information
#----
file = cross_sections.get_file()
materials = cross_sections.get_materials()
particles = cross_sections.get_particles()
Nmat = cross_sections.get_number_of_materials()
densities = zeros(Nmat)
for n in range(1,Nmat)
    densities[n] = materials[n].get_density()
end

#----
# Initialization
#----

numberOfGroups = Vector{Int64}()
Nscat = Vector{Int64}()
Tmin = Vector{Int64}()
Tmax = Vector{Int64}()
legendreOrderScattering = Vector{Int64}()
restEnergies = Vector{Float64}()
energyBoundaries = Vector{Float64}()
energyBoundaries_temp = Vector{Float64}()
groupVelocities = Vector{Float64}()
totalCrossSections_temp = Vector{Float64}()
totalCrossSections = Vector{Float64}()
absorptionCrossSections_temp = Vector{Float64}()
absorptionCrossSections = Vector{Float64}()
stoppingPowers = Vector{Float64}()
stoppingPowers_temp = Vector{Float64}()
continuousScatteringCrossSections_temp = Vector{Float64}()
continuousScatteringCrossSections = Vector{Float64}()
energyDepositionCrossSections_temp = Vector{Float64}()
energyDepositionCrossSections = Vector{Float64}()
chargeDepositionCrossSections_temp = Vector{Float64}()
chargeDepositionCrossSections = Vector{Float64}()
scatteringCrossSections = Vector{Float64}()

P = []
particleNames = []
cumulativeNumberOfGroups = []
numberOfParticles = []
numberOfMaterials = []
legendreOrder = []
title = []
scatDict = Dict{Int64,Dict{Int64,Dict{Int64,Vector{Float64}}}}()
scatNDict = Dict{Int64,Dict{Int64,Dict{Int64}}}()
p = 1


# Positional variables
line_index = 1
L = []
L_II2 = nothing
index_III = 1
index_III_line = 1
index_scat = 1
record_III = []

#----
# Read the FMAC-M file line per line
#----

open(file) do f
for line in readlines(f)
line = one_space(line)
if line != ""

# FMAC-M format validation
if line_index == 1
    if line[1:6] != "FMAC-M" error("The file is not a FMAC-M format file.") end

# N/A
elseif line_index >= 2 && line_index <= 6

# Record II + 1 : Integer parameters
elseif line_index >=7 && line_index <= 13
    append!(P,map(x->parse(Int,x),split(line," ")))
    # Store informations
    if line_index == 13
    cumulativeNumberOfGroups = P[1]
    numberOfParticles = P[2]
    numberOfMaterials = P[15]
    L_II2 = 13+ceil(Int64,P[17]/6)
    legendreOrder = P[20] - 1
    end

# Record II + 2 : Lengths of records III
elseif line_index >= 14 && line_index <= L_II2
    append!(L,map(x->parse(Int,x),split(line," ")))
    if line_index == L_II2
    for i in L 
        if i != 0
            append!(record_III,index_III)
        end
        index_III += 1
    end
    index_III = 1
    end

# Record III
elseif line_index > L_II2 && index_III <= length(record_III)

# Initalization of vectors
if line_index == L_II2 + 1
    energyBoundaries = zeros(cumulativeNumberOfGroups+numberOfParticles)
    totalCrossSections = zeros(cumulativeNumberOfGroups,numberOfMaterials)
    absorptionCrossSections = zeros(cumulativeNumberOfGroups,numberOfMaterials)
    stoppingPowers = zeros(cumulativeNumberOfGroups+numberOfParticles,numberOfMaterials)
    continuousScatteringCrossSections = zeros(cumulativeNumberOfGroups,numberOfMaterials)
    energyDepositionCrossSections = zeros(cumulativeNumberOfGroups,numberOfMaterials)
    chargeDepositionCrossSections = zeros(cumulativeNumberOfGroups,numberOfMaterials)

end

# Record III + 1 : Title
if record_III[index_III] == 1
    title = line
    index_III += 1

# Record III + 2 : Number of groups of energy mesh by particle types
elseif record_III[index_III] == 2
    index_III, index_III_line,numberOfGroups = read_Int!(index_III,index_III_line,L[2],line,numberOfGroups)

# Record III + 3 : Particle names
elseif record_III[index_III] == 3
    particleNames = split(line," ")
    particles_temp = Vector{Particle}()
    for i in range(1,numberOfParticles)
        if particleNames[i] == "GAMA" && any(is_photon.(particles)) push!(particles_temp,particles[findfirst(is_photon.(particles))]) end
        if particleNames[i] == "BETA" && any(is_electron.(particles)) push!(particles_temp,particles[findfirst(is_electron.(particles))]) end
        if particleNames[i] == "POSITR" && any(is_positron.(particles)) push!(particles_temp,particles[findfirst(is_positron.(particles))]) end
    end
    particles = particles_temp
    index_III += 1

# Record III + 4 : Particle rest energies by particle types
elseif record_III[index_III] == 4
    index_III, index_III_line,restEnergies = read_Float!(index_III,index_III_line,L[4],line,restEnergies)

# Record III + 5 : Energy mesh boundaries
elseif record_III[index_III] == 5
    index_III, index_III_line,energyBoundaries_temp = read_Float!(index_III,index_III_line,L[5],line,energyBoundaries_temp)
    # To put the values in order
    if index_III_line == 1
    for n in range(1,numberOfParticles), ig in range(1,numberOfGroups[n]+1)
        if ig == numberOfGroups[n] + 1
        energyBoundaries[sum(numberOfGroups[1:n].+1)] = energyBoundaries_temp[end-n+1]    
        elseif n != 1
        energyBoundaries[sum(numberOfGroups[1:n-1].+1)+ig] = energyBoundaries_temp[sum(numberOfGroups[1:n-1])+ig]
        else
        energyBoundaries[ig] = energyBoundaries_temp[ig]
        end
    end
    end

# Record III + 6 : Group velocities
elseif record_III[index_III] == 6
    index_III, index_III_line,groupVelocities = read_Float!(index_III,index_III_line,L[6],line,groupVelocities)

#Record III + 11 : Array that define the minimal group number from which a transition to group g is available
elseif record_III[index_III] == 11
    index_III, index_III_line,Tmin = read_Int!(index_III,index_III_line,L[11],line,Tmin)

#Record III + 12 : Array that define the maximal group number from which a transition to group g is available
elseif record_III[index_III] == 12
    index_III, index_III_line,Tmax = read_Int!(index_III,index_III_line,L[12],line,Tmax)

#Record III + 14 : Number of scattering cross-section moments used to calculate scattering from the g-th group
elseif record_III[index_III] == 14
    index_III, index_III_line,Nscat = read_Int!(index_III,index_III_line,L[14],line,Nscat)

# Record III + 15 : Total cross sections
elseif record_III[index_III] == 15
    index_III, index_III_line,totalCrossSections_temp = read_Float!(index_III,index_III_line,L[15],line,totalCrossSections_temp)
    if index_III_line == 1
        for imat in range(1,numberOfMaterials)
            totalCrossSections[:,imat] = totalCrossSections_temp[imat:numberOfMaterials:end]
        end
    end

# Record III + 16 : Absorption cross sections
elseif record_III[index_III] == 16
    index_III, index_III_line,absorptionCrossSections_temp = read_Float!(index_III,index_III_line,L[16],line,absorptionCrossSections_temp)
    if index_III_line == 1
        for imat in range(1,numberOfMaterials)
            absorptionCrossSections[:,imat] = absorptionCrossSections_temp[imat:numberOfMaterials:end]
        end
    end

# Record III + 35 : Restricted stopping power 
elseif record_III[index_III] == 35
    index_III, index_III_line,stoppingPowers_temp = read_Float!(index_III,index_III_line,L[35],line,stoppingPowers_temp)
    # To put the values in order
    if index_III_line == 1
    for imat in range(1,numberOfMaterials)
    for n in range(1,numberOfParticles), ig in range(1,numberOfGroups[n]+1)
        if ig == numberOfGroups[n] + 1
        stoppingPowers[sum(numberOfGroups[1:n].+1),imat] = stoppingPowers_temp[end-n*numberOfMaterials+imat]    
        elseif n != 1
        stoppingPowers[sum(numberOfGroups[1:n-1].+1)+ig,imat] = stoppingPowers_temp[((ig-1)*numberOfMaterials+imat)+sum(numberOfGroups[1:n-1])*numberOfMaterials]
        else
        stoppingPowers[ig,imat] = stoppingPowers_temp[((ig-1)*numberOfMaterials+imat)]
        end
    end
    end
    end

# Record III + 36 : Restricted continuous scattering cross-sections
elseif record_III[index_III] == 36
    index_III, index_III_line,continuousScatteringCrossSections_temp = read_Float!(index_III,index_III_line,L[36],line,continuousScatteringCrossSections_temp)
    if index_III_line == 1
        for imat in range(1,numberOfMaterials)
            continuousScatteringCrossSections[:,imat] = continuousScatteringCrossSections_temp[imat:numberOfMaterials:end]
        end
    end

# Record III + 37 : Energy deposition cross-sections
elseif record_III[index_III] == 37
    index_III, index_III_line,energyDepositionCrossSections_temp = read_Float!(index_III,index_III_line,L[37],line,energyDepositionCrossSections_temp)
    if index_III_line == 1
        for imat in range(1,numberOfMaterials)
            energyDepositionCrossSections[:,imat] = energyDepositionCrossSections_temp[imat:numberOfMaterials:end]
        end
    end

# Record III + 38 : Charge deposition cross-sections
elseif record_III[index_III] == 38
    index_III, index_III_line,chargeDepositionCrossSections_temp = read_Float!(index_III,index_III_line,L[38],line,chargeDepositionCrossSections_temp)
    if index_III_line == 1
        for imat in range(1,numberOfMaterials)
            chargeDepositionCrossSections[:,imat] = chargeDepositionCrossSections_temp[imat:numberOfMaterials:end]
        end
    end
end

# Record III + 39 : Legendre moments of scattering cross sections
else
    # Legendre orders
    if index_scat % 2 == 1
    index_III, index_III_line,legendreOrderScattering = read_Int!(index_III,index_III_line,3+numberOfMaterials,line,legendreOrderScattering)
    if index_III_line == 1
        for imat in range(1,numberOfMaterials)
            if haskey(scatNDict,imat)
                if haskey(scatNDict[imat],legendreOrderScattering[2])
                    scatNDict[imat][legendreOrderScattering[2]][legendreOrderScattering[1]] =abs.(legendreOrderScattering[2+imat])
                else
                    scatNDict[imat][legendreOrderScattering[2]] =  Dict(legendreOrderScattering[1] => abs.(legendreOrderScattering[2+imat]))
                end
            else
                scatNDict[imat] = Dict(legendreOrderScattering[2] => Dict(legendreOrderScattering[1] => abs.(legendreOrderScattering[2+imat])))
            end
        end
        scatteringCrossSections=Vector{Float64}()
        index_scat += 1
    end
        
    # Scattering cross sections
    else
    index_III, index_III_line,scatteringCrossSections = read_Float!(index_III,index_III_line,legendreOrderScattering[3+numberOfMaterials],line,scatteringCrossSections)
    if index_III_line == 1
        for imat in range(1,numberOfMaterials)
            pos1 = scatNDict[imat][legendreOrderScattering[2]][legendreOrderScattering[1]]
            pos2 = 0
            if imat != 1
                for imat2 in range(1,imat-1)
                    pos2 = pos2 + scatNDict[imat2][legendreOrderScattering[2]][legendreOrderScattering[1]]
                end
            end
            if haskey(scatDict,imat)
                if haskey(scatDict[imat],legendreOrderScattering[2])
                    scatDict[imat][legendreOrderScattering[2]][legendreOrderScattering[1]] = scatteringCrossSections[pos2+1:pos2+pos1]
                else
                    scatDict[imat][legendreOrderScattering[2]] = Dict(legendreOrderScattering[1] => scatteringCrossSections[pos2+1:pos2+pos1])
                end
            else
                scatDict[imat] = Dict(legendreOrderScattering[2] => Dict(legendreOrderScattering[1] => scatteringCrossSections[pos2+1:pos2+pos1]))
            end
        end
        legendreOrderScattering=Vector{Int64}()
        index_scat += 1
    end
    end
end
end
line_index += 1
end
end

#----
# Storing in data structure
#----

# Initialize new instance
multigroup_cross_sections = Array{Multigroup_Cross_Sections}(undef,numberOfParticles,numberOfMaterials)
energy_boundaries = Vector{Vector{Float64}}(undef,numberOfParticles)

# Per particle
for n in range(1,numberOfParticles)
    index = 0
    for i in range(1,n-1) index = index + numberOfGroups[i] end

    energy_boundaries[n] = energyBoundaries[index+1+(n-1):index+(numberOfGroups[n]+1)+(n-1)]

    # Per material
    for imat in range(1,numberOfMaterials) 

        mcs = Multigroup_Cross_Sections(numberOfGroups[n])
        mcs.set_total(totalCrossSections[index+1:index+numberOfGroups[n],imat])
        mcs.set_absorption(absorptionCrossSections[index+1:index+numberOfGroups[n],imat])
        mcs.set_stopping_powers(stoppingPowers[index+1+(n-1):index+(numberOfGroups[n]+1)+(n-1),imat])
        mcs.set_momentum_transfer(continuousScatteringCrossSections[index+1:index+numberOfGroups[n],imat])

        sp_cutoff = stoppingPowers[index+1+(n-1):index+(numberOfGroups[n]+1)+(n-1),imat][end]
        charge = get_charge(particles[n])
        mcs.set_energy_deposition(push!(energyDepositionCrossSections[index+1:index+numberOfGroups[n],imat],sp_cutoff*energy_boundaries[n][end]/(energy_boundaries[n][end-1]-energy_boundaries[n][end])))
        mcs.set_charge_deposition(push!(chargeDepositionCrossSections[index+1:index+numberOfGroups[n],imat],sp_cutoff*(-charge)))

        # Per scattered particle
        for m in range(1,numberOfParticles)

            index2 = 0
            for i in range(1,m-1)
                index2 = index2 + numberOfGroups[i]
            end
            
            scat = zeros(numberOfGroups[n],numberOfGroups[m],legendreOrder+1)
            for ig_i in range(1,numberOfGroups[n]), ig_f in range(1,numberOfGroups[m])
                if haskey(scatDict[imat][ig_f+index2],ig_i+index)
                    Ls = length(scatDict[imat][ig_f+index2][ig_i+index])
                    for ℓ in range(1,Ls)
                        scat[ig_i,ig_f,ℓ] = scatDict[imat][ig_f+index2][ig_i+index][ℓ]
                    end
                end
            end

            mcs.set_scattering(scat)

        end

        multigroup_cross_sections[n,imat] = mcs

    end

end

cross_sections.set_energy((energyBoundaries[1] + energyBoundaries[2])/2)
cross_sections.set_cutoff(energyBoundaries[end])
cross_sections.set_number_of_groups(numberOfGroups)
cross_sections.set_legendre_order(legendreOrder)
cross_sections.set_multigroup_cross_sections(multigroup_cross_sections)
cross_sections.set_energy_boundaries(energy_boundaries)

end

"""
    read_Int!(index_III::Int64,index_III_line::Int64,Li::Int64,line::String,
    vector::Vector{Int64})

Read sequential integer values in FMAC-M file.

# Input Argument(s)
- 'index_III::Int64' : index of the line in the FMAC-M file.
- 'index_III_line::Int64' : index in the data to be collected.
- 'Li::Int64' : length of the data to be collected.
- 'line::String' : line content.
- 'vector::vector{Int64}' : vector that contain the integer values.

# Output Argument(s)
- 'index_III::Int64' : index of the line in the FMAC-M file.
- 'index_III_line::Int64' : index in the data to be collected.
- 'vector::vector{Int64}' : vector that contain the integer values.

# Reference(s)
N/A

"""
function read_Int!(index_III::Int64,index_III_line::Int64,Li::Int64,line::String,vector::Vector{Int64})
    append!(vector,map(x->parse(Int,x),split(line," ")))
    vector=Int64.(vector)
    # Is it the last line of the record ?
    if ceil(Int,Li/6) == index_III_line
        index_III += 1 
        index_III_line = 1
    else
        index_III_line += 1
    end
    return index_III,index_III_line,vector
end

"""
    read_Float!(index_III::Int64,index_III_line::Int64,Li::Int64,line::String,
    vector::Vector{Float64})

Read sequential real values in FMAC-M file.

# Input Argument(s)
- 'index_III::Int64' : index of the line in the FMAC-M file.
- 'index_III_line::Int64' : index in the data to be collected.
- 'Li::Int64' : length of the data to be collected.
- 'line::String' : line content.
- 'vector::vector{Float64}' : vector that contain the real values.

# Output Argument(s)
- 'index_III::Int64' : index of the line in the FMAC-M file.
- 'index_III_line::Int64' : index in the data to be collected.
- 'vector::vector{Float64}' : vector that contain the real values.

# Reference(s)
N/A

"""
function read_Float!(index_III::Int64,index_III_line::Int64,Li::Int64,line::String,vector::Vector{Float64})
    if string(line[1]) != "-" line = "+" * line end
    # Is it the last line of the record ?
    if ceil(Int,Li/6) == index_III_line
        index_III += 1
        index_III_line = 1
    else
        index_III_line += 1
    end
    if index_III_line != 1
        nlin = 6
    else
        nlin = Li%6
        if nlin==0
            nlin=6
        end
    end
    for i in 1:nlin
        string(line[(i-1)*12+1])=="-" ? sign=-1.0 : sign=1.0
        mantissa = sign*parse(Float64,line[(i-1)*12+2:(i-1)*12+8])
        string(line[(i-1)*12+10])=="-" ? sign=-1.0 : sign=1.0
        expo = sign*parse(Float64,line[(i-1)*12+11:(i-1)*12+12])
        append!(vector,mantissa * 10^expo)
        vector=Float64.(vector)
    end
    return index_III,index_III_line,vector
end