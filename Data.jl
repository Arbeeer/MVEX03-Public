using XLSX, DataFrames, LinearAlgebra, MathOptInterface

# Constants
noOfDays = 20
noOfTimeSteps = 32
lambda = 1
competencies = ["BB", "Faktura_Chatt", "Mob", "Reskontra", "Reskontra_Chatt"]

#Agents
# The Array will consist of 1s (yes) and 0s (no). So every row i will tell us which competencies agent i has and/or doesn't have
function Agents(competencies)
    #load xlsx file
    file_path = "C:\\Users\\Asus\\OneDrive\\Skrivbord\\AgentSkills.xlsx"
    xlsx_file = XLSX.readxlsx(file_path)
    sheet = XLSX.getdata(xlsx_file[1])

    agentComp = Dict() #agent competence
    #agentWork = Dict() #employment type (percentage hours)
    for i in 1:size(sheet,1)
        for j in competencies
            if occursin(j, String(sheet[i]))
                if haskey(agentComp, i)
                    agentComp[i] = agentComp[i] * "," * j
                else
                    agentComp[i] = j
                end
            end
        end
        #agentWork[i] = sheet[i,2]
    end
        return agentComp#, agentWork
end

#dictionary with agents as keys and their competencies as output in the order of the "competencies" - array
comps = Agents(competencies)
#workh = Agents(competencies)[2]
AGENT_I = zeros(Int, length(comps), length(competencies))
 for i in 1:length(comps)
    for j in 1:length(competencies)
        if occursin(competencies[j], comps[i])
            AGENT_I[i,j] = 1
        end
    end
end

noOfAgentsPerComp = Dict()
noOfAgentsPerComp[1] = sum(AGENT_I[:,1])
noOfAgentsPerComp[2] = sum(AGENT_I[:,2])
noOfAgentsPerComp[3] = sum(AGENT_I[:,3])
noOfAgentsPerComp[4] = sum(AGENT_I[:,4])
noOfAgentsPerComp[5] = sum(AGENT_I[:,5])


#SHIFT_J = [1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]'
SHIFT_J = [0 1 1 1 0 1 1 1 1 1 1 1 0 0 0 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 0]'
#SHIFT_J = [1 1 1 1 0 1 1 1 1 1 1 1 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 0 0]'
#SHIFT_J = [1 1 1 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]'

## Days in planning period
DAY_L = 1:noOfDays

## Time steps in planning period
TIMESTEP_H = 1:noOfTimeSteps

COMPETENCE_K = 1:length(competencies)

t_h = 0.25



function FindDemand(competence::Int, F_hk::Int) #finds demand per time step over the year 2022 for a given competence.
    #load xlsx file
    file_path = "C:\\Users\\Asus\\OneDrive\\Skrivbord\\Competence.xlsx"
    #file_path = "C:\\Users\\Asus\\OneDrive\\Skrivbord\\CompetenceForecast.xlsx"
    xlsx_file = XLSX.readxlsx(file_path)
    sheet = XLSX.getdata(xlsx_file[competence])
    
    daysInYear = Dict() #key=competence, val=number of days in full year in excel sheet
    daysInYear[1] = 320 #1280/4
    daysInYear[2] = 261 #1044/4
    daysInYear[3] = 273 #1092/4
    daysInYear[4] = 256 #1024/4
    daysInYear[5] = 259 #1036/4


    summerTime = Dict()
    summerTime[1] = [73, 265] #[292/4, 1060/4] #start and end of summer time, e.g. 292/4=73 --> summer time starts on day 74
    summerTime[2] = [62, 217] #[248/4, 868/4]
    summerTime[3] = [64, 228] #[256/4, 912/4]
    summerTime[4] = [61, 212] #[244/4, 848/4]
    summerTime[5] = [60, 215] #[240/4, 860/4]


    summation=0
    counter = 0
    noOfCalls = zeros(daysInYear[competence], noOfTimeSteps) # = d_lhk
    for i in 1:daysInYear[competence]
        if (competence < 5)
            if (i <= summerTime[competence][1] || i > summerTime[competence][2])
                startRow = 37 #winter time
                endRow = 68
                offSet = 36
            else
                startRow = 33 #summer time
                endRow = 64
                offSet = 32
            end
        elseif (competence == 5)
            if (i <= summerTime[competence][1] || i > summerTime[competence][2])
                startRow = 17 #winter time
                endRow = 48
                offSet = 16
            else
                startRow = 13 #summer time
                endRow = 44
                offSet = 12
            end
        end


        for j in startRow:endRow #avg AHT
            if (!ismissing(sheet[j,3+4*(i-1)]) && sheet[j,3+4*(i-1)] != 0)
                summation = summation + (sheet[j,4+4*(i-1)])/sheet[j,3+4*(i-1)] #+ sheet[j,5+4*(i-1)]
                counter = counter+1
            end

            acceptedCalls = zeros(Int, 1, noOfTimeSteps)
            for k in startRow:endRow
                if (ismissing(sheet[k,3+4*(i-1)]))
                    acceptedCalls[k-offSet] = 0
                else
                    acceptedCalls[k-offSet] = sheet[k,3+4*(i-1)]
                end
            end

                if (sum(acceptedCalls) != 0)
                    if (ismissing(sheet[j, 2+4*(i-1)]))
                    noOfCalls[i, j-offSet] = 0
                else
                    noOfCalls[i, j-offSet] = sheet[j, 2+4*(i-1)]
                end
            else
                noOfCalls[i, j-offSet] = Inf
            end
        end
    end

    rowsWithInf = []
    for i in 1:daysInYear[competence]
        if (sum(noOfCalls[i,:]) == Inf)
            push!(rowsWithInf, i)
        end
    end

    rowsToKeep = filter(x -> !(x in rowsWithInf), [i for i in 1:daysInYear[competence]])

    noOfCalls = noOfCalls[rowsToKeep, :] #remove rows from noOfCalls

    #demand = zeros(Int, daysInYear[competence], noOfTimeSteps)
    demand = (noOfCalls .* (round(Int, summation/counter)+180))./F_hk # noOfCalls * AHT / F_hk, should maybe subtract no of agents by 1 to get maxwait- rather than maxhandling time, in that case do so in the model file.
    #println("AHT: ", round(Int, summation/counter), "s")

    return demand, round(Int, summation/counter)+180, noOfCalls #+180 = avg. after call work.
end

#D_lhk = FindDemand(1, 900)[1]
#sum(D_lhk[1,:])

# a = [i for i in 38:75]

function TeliaDemand(competence::Int) #finds demand per time step over the year 2022 for a given competence.
    #load xlsx file
    file_path = "C:\\Users\\Asus\\OneDrive\\Skrivbord\\ErlangTelia.xlsx"
    xlsx_file = XLSX.readxlsx(file_path)
    sheet = XLSX.getdata(xlsx_file[competence])
    
    days = Dict() #key=competence, val=number of days in excel sheet
    days[1] = 21
    days[2] = 21
    days[3] = 21
    days[4] = 21
    days[5] = 21

    noOfCalls = zeros(days[competence], noOfTimeSteps) # = d_lhk forecast
    noOfStaff = zeros(days[competence], noOfTimeSteps) # basically sum_j w_jk
    startRow = 1
    endRow = 335
    startCol = 1
    if (competence == 2 || competence == 5)
        endCol = 28
    else 
        endCol = 32
    end
    for i in 1:days[competence]
        for j in startRow:endRow
            for k in startCol:endCol
                noOfCalls[i,k] = sheet[1+16*(i-1),k]
                noOfStaff[i,k] = sheet[14+16*(i-1),k]
            end
        end
    end
    return noOfCalls, noOfStaff
end


function FindDemand1(competence::Int, F_hk::Int)
    noOfCalls = TeliaDemand(competence)[1]
    AHT = FindDemand(competence, F_hk)[2]
    demand = ((noOfCalls .* AHT)./F_hk)
    return demand
end

 #a = [i for i in 25:62]
function FindDemand2(competence::Int, F_hk::Int) #finds demand per time step over the year 2022 for a given competence.
    #load xlsx file
    file_path = "C:\\Users\\Asus\\OneDrive\\Skrivbord\\Competence_2023.xlsx"
    #file_path = "C:\\Users\\Asus\\OneDrive\\Skrivbord\\CompetenceForecast.xlsx"
    xlsx_file = XLSX.readxlsx(file_path)
    sheet = XLSX.getdata(xlsx_file[competence])
    
    daysInYear = Dict() #key=competence, val=number of days in 0.5 year in excel sheet
    daysInYear[1] = 118 #472/4
    daysInYear[2] = 103 #412/4
    daysInYear[3] = 116 #464/4
    daysInYear[4] = 103 #412/4
    daysInYear[5] = 103 #412/4


    summerTime = Dict()
    summerTime[1] = 73 #[292/4, 1060/4] #start of summer time, e.g. 292/4=73 --> summer time starts on day 74
    summerTime[2] = 62 #[248/4, 868/4]
    summerTime[3] = 74 #[296/4, 912/4]
    summerTime[4] = 61 #[244/4, 848/4]
    summerTime[5] = 64 #[256/4, 860/4]


    summation=0
    counter = 0
    noOfCalls = zeros(daysInYear[competence], noOfTimeSteps) # = d_lhk
    for i in 1:daysInYear[competence]
        if (competence == 1 || competence == 3 || competence == 4)
            if (i <= summerTime[competence])
                startRow = 37 #winter time
                endRow = 68
                offSet = 36
            else
                startRow = 33 #summer time
                endRow = 64
                offSet = 32
            end
        elseif (competence == 5 || competence == 2)
            if (i <= summerTime[competence])
                startRow = 25 #winter time
                endRow = 56
                offSet = 24
            else
                startRow = 21 #summer time
                endRow = 52
                offSet = 20
            end
        end


        if (competence == 1 || competence == 3)
            for j in startRow:endRow #avg AHT
                if (!ismissing(sheet[j,7+4*(i-1)]) && sheet[j,7+4*(i-1)] != 0)
                    summation = summation + (sheet[j,8+4*(i-1)] + sheet[j,9+4*(i-1)])/sheet[j,7+4*(i-1)]
                    counter = counter+1
                end

                acceptedCalls = zeros(Int, 1, noOfTimeSteps)
                for k in startRow:endRow
                    if (ismissing(sheet[k,7+4*(i-1)]))
                        acceptedCalls[k-offSet] = 0
                    else
                        acceptedCalls[k-offSet] = sheet[k,7+4*(i-1)]
                    end
                end

                    if (sum(acceptedCalls) != 0)
                        if (ismissing(sheet[j, 6+4*(i-1)]))
                        noOfCalls[i, j-offSet] = 0
                    else
                        noOfCalls[i, j-offSet] = sheet[j, 6+4*(i-1)]
                    end
                else
                    noOfCalls[i, j-offSet] = Inf
                end
            end
        else
            for j in startRow:endRow #avg AHT
                if (!ismissing(sheet[j,3+4*(i-1)]) && sheet[j,3+4*(i-1)] != 0)
                    summation = summation + (sheet[j,4+4*(i-1)] + sheet[j,5+4*(i-1)])/sheet[j,3+4*(i-1)]
                    counter = counter+1
                end

                acceptedCalls = zeros(Int, 1, noOfTimeSteps)
                for k in startRow:endRow
                    if (ismissing(sheet[k,3+4*(i-1)]))
                        acceptedCalls[k-offSet] = 0
                    else
                        acceptedCalls[k-offSet] = sheet[k,3+4*(i-1)]
                    end
                end

                    if (sum(acceptedCalls) != 0)
                        if (ismissing(sheet[j, 2+4*(i-1)]))
                        noOfCalls[i, j-offSet] = 0
                    else
                        noOfCalls[i, j-offSet] = sheet[j, 2+4*(i-1)]
                    end
                else
                    noOfCalls[i, j-offSet] = Inf
                end
            end
        end
    end

    rowsWithInf = []
    for i in 1:daysInYear[competence]
        if (sum(noOfCalls[i,:]) == Inf)
            push!(rowsWithInf, i)
        end
    end

    rowsToKeep = filter(x -> !(x in rowsWithInf), [i for i in 1:daysInYear[competence]])

    noOfCalls = noOfCalls[rowsToKeep, :] #remove rows from noOfCalls

    #demand = zeros(Int, daysInYear[competence], noOfTimeSteps)
    demand = (noOfCalls .* round(Int, summation/counter))./F_hk # noOfCalls * AHT / F_hk, should maybe subtract no of agents by 1 to get maxwait- rather than maxhandling time, in that case do so in the model file.
    #println("AHT: ", round(Int, summation/counter), "s")

    return demand, round(Int, summation/counter), noOfCalls
end