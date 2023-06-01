# Shifts of length >4h but <6h.
using LinearAlgebra

noOfTimeSteps = 32

function GenerateMatricesMedium(n::Int)
    matrices = []
    for k in 1:n-2, l in k+1:n-1, m in l+1:n
        matrix = fill(1, 1, n)
        matrix[1, [k,l,m]] .= 0
        if (IsValidMedium(vec(matrix)))
            push!(matrices, matrix)
        end
    end
    return matrices
end

#show(GenerateMatrices(23))

# define the restrictions
function IsValidMedium(arr)
    indices = findall(x -> x == 0, arr)

    sum(arr) == length(arr)-3&&

    # check that the first three time-steps sum to 3 (earliest B1 start 45min)
    sum(arr[1:3]) == 3 &&

    #  check that the first 11 time-steps sum to 10 (latest B1 start 2h 45min)
    sum(arr[1:11]) == 10 &&

    # check that the time-steps between the first and second break sum up to at least five, but no more than 13 (at least 1h 15min between breaks, but no more than 3h)
    5 <= sum(arr[indices[1]:indices[2]]) <= 12 &&

    #check if second and third 0 are consecutive (30min lunch break)
    ContainsConsecutiveZerosMedium(arr) &&

    #last 3 elements must be 1
    sum(arr[length(arr)-2:length(arr)]) == 3
end

# helper function to check for consecutive sequences of zeros
function ContainsConsecutiveZerosMedium(arr)
    indices = findall(x -> x == 0, arr)
    if  (indices[3] - indices[2] == 1)
        return true
    else
        return false
    end
end

function ShiftRight(arr::Vector{T}, k::Int) where T #shift each element k steps to the right
    n = length(arr)
    new_arr = similar(arr)
    for i in 1:n
        new_i = i + k
        if new_i > n
            new_i = new_i - n
        end
        new_arr[new_i] = arr[i]
    end
    return new_arr
end

function GenerateMediumShifts(H) # H = length of shift (4.25h - 5.75h)
    matrices = GenerateMatricesMedium(H)
    validShifts = zeros(Int, 1, H) #init
    for i in 1:length(matrices)
        validShifts = vcat(validShifts, matrices[i])
    end
    validShifts = validShifts[2:end, :] #rm first (init) row
    return validShifts
end

function ReducedCostsMedium(shift::Vector, u::Vector) #computes the reduced costs of all shifts that belong to the same type as `shift', with dual variable `u'
    shiftCopy = copy(shift)
    costVec = []
    while (shiftCopy[length(shiftCopy)] != 1)
        costVec = push!(costVec, dot(shiftCopy,u))
        shiftCopy = ShiftRight(shiftCopy,1)
    end
    costVec = push!(costVec, dot(shiftCopy,u))
    #costVec = costVec[2:end] #rm first elem
    return costVec
end


function BestColumnMedium(u::Vector{Float64}, H::Int)
    shiftColl = GenerateMediumShifts(H) #generate shifts
    fillTheDay = zeros(Int, size(shiftColl,1), noOfTimeSteps-H) #38-H zeroes i.e. as many zeroes as are left in the opening hours
    allShifts = hcat(shiftColl, fillTheDay) #all (unique) shifts
    costDict = Dict() #dictionary mapping biggest reduced cost to its corresponding shift type
    minCostVec = [] #init
    for n in 1:size(allShifts,1)
        minRedCosts = maximum(ReducedCostsMedium(allShifts[n,:], u))
        minCostVec = push!(minCostVec, minRedCosts)
        costDict[minRedCosts] = n #store the best shift type in a dictionary
    end
    #minCostVec = minCostVec[2:end, :]
    bestCost = maximum(minCostVec)

    #if (bestCost < 0) #if the min-cost is <0 we find the shift variation that gave it
        bestShiftType = costDict[bestCost]
        k = 0
        while (dot(ShiftRight(allShifts[bestShiftType,:], k), u) != bestCost && k<size(allShifts,2)) #add termination criteria k<shift length for safety (the most amount of shift variations per type is 34, which is for the shortest shift), also the ShiftRight function can only handle the case when k <= arr_length. 
            k = k+1
        end
        return ShiftRight(allShifts[bestShiftType,:], k)
   # else
        #return nothing
    #end
end

function BestColMedium(u::Vector{Float64})
    minCostOfShifts = [] #8 so that all i below fit
    colDict = Dict()
    for i in 0:6 #17+6=23
        result = BestColumnMedium(u, 17+i) #shortest shift is 17 time steps
        if result === nothing
            minCostOfShifts[i+1] = Inf
        else
            minCostOfShifts = push!(minCostOfShifts, dot(result, u))
            colDict[dot(result, u)] = result
        end
    end
    minCost = maximum(minCostOfShifts)
    #if minCost < 0
        return colDict[minCost]
    #else return nothing
    #end
end

#u = vec(2*rand(1, noOfTimeSteps).-1)
#u = vec(zeros(1,noOfTimeSteps))
#result = BestColMedium(u)'
#sum(result)
#result = BestColumnLong(u, 36)'

#@elapsed BestColumnLong(u,36)
#@elapsed BestColumnMedium(u)
#@elapsed GenerateMatrices(36)

#if(result === nothing)
    #println("No better column found")
#else
    #println("Best column: ", result)
    #println("Smallest reduced cost: ", dot(result, u))
#end