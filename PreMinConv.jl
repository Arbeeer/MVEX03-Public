using JuMP, Gurobi, Clp, Pkg, Cbc
using LinearAlgebra, Printf
include("ShortShifts.jl")
include("MediumShifts.jl")
include("LongShifts.jl")
include("Data.jl")
#include("SmallData.jl")

#CopySHIFT_J = SHIFT_J

function MinConvComb(SHIFT_Jc, D_hk, lambda)
    a_hj = SHIFT_Jc
    m = Model() #Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
    @variable(m, w_kj[COMPETENCE_K, 1:size(SHIFT_Jc,2)] >= 0)
    @variable(m, yu_hk[TIMESTEP_H, COMPETENCE_K] >= 0)

    
    @objective(m, Min, 
    lambda*sum(sum(yu_hk[h,k] for h in TIMESTEP_H) for k in COMPETENCE_K) +
    (1-lambda)*sum(sum(sum((a_hj[h,j])*w_kj[k,j] for h in TIMESTEP_H) for j in 1:size(SHIFT_Jc,2)) for k in COMPETENCE_K) +
    sum(sum(((findlast(x -> x == 1, BitArray(SHIFT_Jc[:,j])) - findnext(BitArray(SHIFT_Jc[:,j]), 1)) >= 24 ? 2*(1-lambda)*w_kj[k,j] : (1-lambda)*w_kj[k,j]) for k in COMPETENCE_K) for j in 1:size(SHIFT_Jc,2))
    )
    
    #@objective(m, Min, sum(sum(q*yu_hk[h,k] for h in TIMESTEP_H) for k in COMPETENCE_K) 
     #   + sum(sum(sum((a_hj[h,j]*c)*w_kj[k,j] for h in TIMESTEP_H) for j in 1:size(SHIFT_Jc,2)) for k in COMPETENCE_K)
      #  )

    @constraint(m, availability[h in TIMESTEP_H, k in COMPETENCE_K],
        sum(a_hj[h,j]*w_kj[k,j] for j in 1:size(SHIFT_Jc,2)) + yu_hk[h,k] >= D_hk[h,k])

    @constraint(m, limitedAgentsPerCompetence[k in COMPETENCE_K],
        sum(w_kj[k,j] for j in 1:size(SHIFT_Jc,2)) <= noOfAgentsPerComp[k])
    
    @constraint(m, limitedAgentsPerDay,
        sum(sum(w_kj[k,j] for j in 1:size(SHIFT_Jc,2)) for k in COMPETENCE_K) <= size(AGENT_I,1))
    

    set_optimizer(m, Gurobi.Optimizer)
    println("solving BuildModelAndSolveIt...")
    optimize!(m)
    println("Termination status, problem: ", termination_status(m))
    #print(value.(w_kjl))

    dualVal = value.(dual.(availability))
    unregister(m, :w_kj)
    unregister(m, :availability)

    return dualVal, objective_value(m), value.(w_kj), value.(yu_hk)
end

function ColumnGeneration(SHIFT_J, D_hk, lambda)
    CopySHIFT_J = SHIFT_J
    k=1
    q=1
    while(k==1 && q<=300) #introduce more termination criterias
        termCrit = size(CopySHIFT_J, 2)
        #solution_summary(m)
        #has_duals(m)
        #value.(x_ijl)
    
        dualVal = MinConvComb(CopySHIFT_J, D_hk, lambda)[1]
        for i in 1:length(competencies) #get optimal schedule for each competence in that day
            dualVec = vec(dualVal[:,i])
            u = vec(zeros(length(dualVec),1))
            for j in 1:length(dualVec)
                u[j] = dualVec[j]
            end
            bestShort = BestColShort(u) # might have to include SHIFT_J as an input in order to not allow for duplicate shifts, but I don't think that will be necessary since if an existing column in SHIFT_J has a negative reduced cost it means it improves the objective, but since it already is in SHIFT_J it should have improved the objective already, thus existing columns can never have negative reduced costs.
            #println("Elapsed time short: ", @elapsed BestShortShift(u))
            bestLong = BestColLong(u)
            #println("Elapsed time long: ", @elapsed BestLongShift(u))
            bestMedium = BestColMedium(u)

            if (bestShort === nothing)
                redShort = Inf
            else
                redShort = (1-lambda)*(sum(bestShort)+1) - dot(bestShort, u) #paid break (lunch is not paid)
            end
            
            if (bestMedium === nothing)
                redMedium = Inf
            else
                redMedium = (1-lambda)*(sum(bestMedium)+1) - dot(bestMedium, u)
            end

            if (bestLong === nothing)
                redLong = Inf
            else
                redLong = (1-lambda)*(sum(bestLong)+2) - dot(bestLong, u)
            end

            redCosts = [redShort, redMedium, redLong]
            #println("redCosts: ", redCosts)
            #println("udotshort: ", dot(bestShort, u))
            #println("dual value: ", u)
            bestCost = minimum(redCosts)

            if (bestCost <= 0)
                if (bestCost == redShort) #prioritise adding shorter shifts by having them in this order in the if-statement
                    a = 0
                    for j in 1:size(CopySHIFT_J,2) #add only unique columns
                        if (bestShort == CopySHIFT_J[:,j])
                            a = a + 1
                        end
                    end
                    if (a==0)
                        CopySHIFT_J = hcat(CopySHIFT_J, bestShort)
                        println("Added short shift: ", bestShort', "competence: ", i)
                    end
                elseif (bestCost == redMedium)
                    a = 0
                    for j in 1:size(CopySHIFT_J,2) #add only unique columns
                        if (bestMedium == CopySHIFT_J[:,j])
                            a = a + 1
                        end
                    end
                    if (a==0)
                        CopySHIFT_J = hcat(CopySHIFT_J, bestMedium)
                        println("Added medium shift: ", bestMedium', "competence: ", i)
                    end
                else
                    a = 0
                    for j in 1:size(CopySHIFT_J,2) #add only unique columns
                        if (bestLong == CopySHIFT_J[:,j])
                            a = a + 1
                        end
                    end
                    if (a==0)
                        CopySHIFT_J = hcat(CopySHIFT_J, bestLong) #prio long shifts
                        println("Added long shift: ", bestLong', "competence: ", i)
                    end
                end
            end
        end

        q=q+1
        println("Amount of shifts: ", size(CopySHIFT_J, 2))
        println("Iteration: ", q-1)
        if (termCrit == size(CopySHIFT_J,2)) #check if the size of SHIFT_J has changed - if it hasn't we are at optimum
            k=0
        end

    end
    println("Added ", size(CopySHIFT_J,2)-size(SHIFT_J, 2), " total shifts.")
    println("Final objective value: ", MinConvComb(CopySHIFT_J, D_hk, lambda)[2])
    return CopySHIFT_J
end
#println("Total time for column generation: ", @elapsed ColumnGeneration(SHIFT_J))

function IntMinConvComb(SHIFT_Jc, D_hk, lambda)
    a_hj = SHIFT_Jc
    m = Model() #Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
    @variable(m, w_kj[COMPETENCE_K, 1:size(SHIFT_Jc,2)] >= 0, Int)
    @variable(m, yu_hk[TIMESTEP_H, COMPETENCE_K] >= 0, Int)

    @objective(m, Min, 
    lambda*sum(sum(yu_hk[h,k] for h in TIMESTEP_H) for k in COMPETENCE_K) +
    (1-lambda)*sum(sum(sum((a_hj[h,j])*w_kj[k,j] for h in TIMESTEP_H) for j in 1:size(SHIFT_Jc,2)) for k in COMPETENCE_K) +
    sum(sum(((findlast(x -> x == 1, BitArray(SHIFT_Jc[:,j])) - findnext(BitArray(SHIFT_Jc[:,j]), 1)) >= 24 ? 2*(1-lambda)*w_kj[k,j] : (1-lambda)*w_kj[k,j]) for k in COMPETENCE_K) for j in 1:size(SHIFT_Jc,2))
    )
    
    #@objective(m, Min, sum(sum(q*yu_hk[h,k] for h in TIMESTEP_H) for k in COMPETENCE_K) 
     #   + sum(sum(sum((a_hj[h,j]*c)*w_kj[k,j] for h in TIMESTEP_H) for j in 1:size(SHIFT_Jc,2)) for k in COMPETENCE_K)
      #  )

    @constraint(m, availability[h in TIMESTEP_H, k in COMPETENCE_K],
        sum(a_hj[h,j]*w_kj[k,j] for j in 1:size(SHIFT_Jc,2)) + yu_hk[h,k] >= ceil.(Int, D_hk[h,k]))

    @constraint(m, limitedAgentsPerCompetence[k in COMPETENCE_K],
        sum(w_kj[k,j] for j in 1:size(SHIFT_Jc,2)) <= noOfAgentsPerComp[k])

    @constraint(m, limitedAgentsPerDay,
        sum(sum(w_kj[k,j] for j in 1:size(SHIFT_Jc,2)) for k in COMPETENCE_K) <= size(AGENT_I,1))

    set_optimizer(m, Gurobi.Optimizer)
    println("solving MIPModel...")
    optimize!(m)
    println("Termination status, problem: ", termination_status(m))

    return value.(w_kj), value.(yu_hk), objective_value(m)
end


function MILPmodel(w_kjl, SHIFT_Jc, days)
    SHIFTS = [SHIFT_Jc[l] for l in 1:days] 
    m1 = Model() #Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
    @variable(m1, x_ijl[1:size(AGENT_I,1), 1:maximum(size.(SHIFTS,2)), 1:days], Bin)
    @variable(m1, z_i[1:size(AGENT_I,1)], Bin)
    @variable(m1, S_ikjhl[1:size(AGENT_I,1), COMPETENCE_K, 1:maximum(size.(SHIFTS,2)), TIMESTEP_H, 1:days] >= 0)
    a_hj = SHIFTS #\ell dependent

    @objective(m1, Min, sum(z_i[i] for i in 1:size(AGENT_I,1)))

    @constraint(m1, workOrNot[i in 1:size(AGENT_I,1)],
       sum(sum(x_ijl[i,j,l] for j in 1:size(SHIFTS[l],2)) for l in 1:days) <= maximum(size.(SHIFTS,2))*length([1:days])*z_i[i])

    @constraint(m1, matchDemand[l in 1:days, k in COMPETENCE_K, j in 1:size(SHIFTS[l],2), h in 1:noOfTimeSteps],
       sum(S_ikjhl[i,k,j,h,l] for i in findall(x -> x == 1, AGENT_I[:,k])) >= a_hj[l][h,j]*w_kjl[l][k,j])

    @constraint(m1, timedistr[i in 1:size(AGENT_I,1), l in 1:days, j in 1:size(SHIFTS[l],2), h in TIMESTEP_H],
        sum(S_ikjhl[i,k,j,h,l] for k in findall(x -> x == 1, AGENT_I[i,:])) == a_hj[l][h,j]*x_ijl[i,j,l])

    @constraint(m1, maxhours[i in 1:size(AGENT_I,1)], 
        sum(sum(sum(t_h*a_hj[l][h,j]*x_ijl[i,j,l] for h in TIMESTEP_H) for j in 1:size(SHIFTS[l],2)) for l in 1:days) <= 160)

    @constraint(m1, oneshift[i in 1:size(AGENT_I,1), l in 1:days], 
        sum(x_ijl[i,j,l] for j in 1:size(SHIFTS[l],2)) <= 1)


    set_optimizer(m1, Gurobi.Optimizer)
    println("solving MILPmodel...")
    optimize!(m1)
    println("Termination status, problem: ", termination_status(m1))

    return value.(x_ijl), value.(S_ikjhl), objective_value(m1)
end

function wVariable(intSol)
    w = zeros(Int, size(intSol,1), size(intSol,2))
    for i in 1:size(intSol,1)
        for j in 1:size(intSol,2)
            w[i,j] = round(value.(intSol[i,j]))
        end
    end
    return w
end

function obtainResults(days, offset, lambda)
    maxAHT = 900
    D_lhk = zeros(days, length(TIMESTEP_H), length(competencies)) 

    for l in 1:days
        for k in 1:length(competencies)
            D_lhk[l, :, k] = FindDemand(k, maxAHT)[1][l+offset,:] #in order to take 20 days other than the first 20 you offset [l,:] with eg. [l+20,:] - this gives the next 20 days
        end
    end

    results = []
    genShifts = []
    yvals = []
    for i in 1:days
        println("DAY: ", i)
        colGen = ColumnGeneration(SHIFT_J, D_lhk[i,:,:], lambda)
        intSol = IntMinConvComb(colGen, D_lhk[i,:,:], lambda)
        w = wVariable(intSol[1])
        results = push!(results, w)
        genShifts = push!(genShifts, colGen)
        yvals = push!(yvals, intSol[2])
    end
    #println("SHIFTRESULTS: ", genShifts)
    #println("W-RESULTS: ", results)
    return results, genShifts, yvals, D_lhk
end

r = obtainResults(noOfDays, 0, lambda)
#finalSolution = MILPmodel(r[1], r[2], noOfDays)
#f1 = round.(Int, finalSolution[1])
#f2 = finalSolution[2]
#f3 = finalSolution[3]

#r = obtainResults(20, 0, c[1], c[2])
#finalSolution = round.(Int, MILPmodel(r[1], r[2], 20))
#show(r[2][1]')
#r[1][1]
#sum(r)

#sum(r[1][1])
#sum(sum(r[2][1][h,j]*r[1][1][1,j] for h in 1:38) for j in 1:29)
#sum(r[4][1,:,1])


#sum(finalSolution)

#show(w)

#df = DataFrame(finalSolution[:,:,1], :auto)

#CSV.write("result.csv", df, delim=',')

#IntMinConvComb(r[2][1], D_lhk[1,:,:], lambda)[3]