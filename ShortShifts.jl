# Shifts of length >=3h but <=4h.
using LinearAlgebra

noOfTimeSteps = 32

function GenerateMatricesShort(n::Int)
    matrices = []
    for k in 1:n-1
        matrix = fill(1, 1, n)
        matrix[1, [k]] .= 0
        if (IsValidShort(vec(matrix)))
            push!(matrices, matrix)
        end
    end
    return matrices
end

#GenerateMatricesShort(16)
p=16/32
H=64
q=6/32
#(p*H-q*H)*(H-p*H)
(H^2)*(p-q)*(1-p)
(p-q)*(1-p)

# define the restrictions
function IsValidShort(arr)

    sum(arr) == length(arr)-1&&

    # check that the first three time-steps sum to 3 (earliest B1 start 45min)
    sum(arr[1:3]) == 3 &&

    #last 3 elements must be 1
    sum(arr[length(arr)-2:length(arr)]) == 3
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

function GenerateShortShifts(H) # H = length of shift (4.25h - 5.75h)
    matrices = GenerateMatricesShort(H)
    validShifts = zeros(Int, 1, H) #init
    for i in 1:length(matrices)
        validShifts = vcat(validShifts, matrices[i])
    end
    validShifts = validShifts[2:end, :] #rm first (init) row
    return validShifts
end

function ReducedCostsShort(shift::Vector, u::Vector) #computes the reduced costs of all shifts that belong to the same type as `shift', with dual variable `u'
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


function BestColumnShort(u::Vector{Float64}, H::Int)
    shiftColl = GenerateShortShifts(H) #generate shifts
    fillTheDay = zeros(Int, size(shiftColl,1), noOfTimeSteps-H) #32-H zeroes i.e. as many zeroes as are left in the opening hours
    allShifts = hcat(shiftColl, fillTheDay) #all (unique) shifts
    costDict = Dict() #dictionary mapping smallest reduced cost to its corresponding shift type
    minCostVec = [] #init
    for n in 1:size(allShifts,1)
        maxRedCosts = maximum(ReducedCostsShort(allShifts[n,:], u))
        minCostVec = push!(minCostVec, maxRedCosts)
        costDict[maxRedCosts] = n #store the best shift type in a dictionary
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
    #else
        #return nothing
    #end
end

function BestColShort(u::Vector{Float64})
    minCostOfShifts = []
    colDict = Dict()
    for i in 0:4 #12+4=16
        result = BestColumnShort(u, 12+i) #shortest shift is 12 time steps, including break
        if result === nothing
            minCostOfShifts[i+1] = Inf
        else
            minCostOfShifts = push!(minCostOfShifts, dot(result, u))
            colDict[dot(result, u)] = result
        end
    end
    maxCost = maximum(minCostOfShifts)
    #if minCost < 0
        return colDict[maxCost]
    #else return nothing
    #end
end

#u = vec(2*rand(1, noOfTimeSteps).-1)
#u = vec(zeros(1,45))
#result = BestColShort(u)'
#sum(result)


#if(result === nothing)
    #println("No better column found")
#else
    #println("Best column: ", result)
    #println("Smallest reduced cost: ", dot(result, u))
#end