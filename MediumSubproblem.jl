# Generate medium shifts
using JuMP, Gurobi, LinearAlgebra, SparseArrays, AxisArrays
include("Data.jl")
WORKBLOCK_NBM = 1:5
TIMESTEPH = 2:noOfTimeSteps+1
BREAKLM = [1, 2] #two types of breaks

function BuildAndSolveMedium(u_h, lambda)
    m = Model()
    @variable(m, v_nhs[n in WORKBLOCK_NBM, h in 2:noOfTimeSteps+2, s in 1:h-1], Bin)
    @variable(m, a_h[h in 1:noOfTimeSteps] >= 0)
    
    #@objective(m, Min, c[2]*(sum(s*v_nhs[5, noOfTimeSteps+2, s] for s in 1:noOfTimeSteps+1) - sum(h*v_nhs[1,h,1] for h in 2:noOfTimeSteps+2)-2) #-2 for unpaid time-steps
        #- sum(sum(sum(u_h[h-1]*v_nhs[n,h,s] for s in 1:h-1) for h in TIMESTEPH) for n in WORKBLOCK_NBM))

    #@objective(m, Min, (1-lambda)*(sum(s*v_nhs[5, noOfTimeSteps+2, s] for s in 1:noOfTimeSteps+1) - sum(h*v_nhs[1,h,1] for h in 2:noOfTimeSteps+2)-2) #-2 for unpaid time-steps
     #   - sum(sum(sum(u_h[h-1]*v_nhs[n,h,s] for s in 1:h-1) for h in TIMESTEPH) for n in WORKBLOCK_NBM))

    @objective(m, Min, (1-lambda)*(sum(s*v_nhs[5, noOfTimeSteps+2, s] for s in 1:noOfTimeSteps+1) - sum(h*v_nhs[1,h,1] for h in 2:noOfTimeSteps+2)-2) #-2 for unpaid time-steps
        - sum(u_h[h]*a_h[h] for h in 1:noOfTimeSteps))

    @constraint(m, translate[h in TIMESTEPH],
        a_h[h-1] == sum(sum(sum(v_nhs[n,r,s] for r in h+1:noOfTimeSteps+2) for s in 1:h) for n in 2:5) - sum(sum(v_nhs[5,r,s] for r in h:noOfTimeSteps+2) for s in 1:h-1)
        )

    @constraint(m, 
        sum(v_nhs[1,h,1] for h in 2:noOfTimeSteps+2) == 1)

    @constraint(m, oneArcInFirstBlock[r in 2:noOfTimeSteps+1],
        sum(v_nhs[1,h,r] for h in r+1:noOfTimeSteps+2) == 0)
    
    @constraint(m, start[h in TIMESTEPH],
        sum(v_nhs[1,h,s] for s in 1:h-1) == sum(v_nhs[2,r,h] for r in h+1:noOfTimeSteps+2))

    @constraint(m, firstBreak[h in TIMESTEPH],
        sum(v_nhs[2,h,s] for s in 1:h-1) == sum(v_nhs[3,r,h+BREAKLM[1]] for r in h+1+BREAKLM[1]:noOfTimeSteps+2))

    @constraint(m, secondBreak[h in TIMESTEPH],
        sum(v_nhs[3,h,s] for s in 1:h-1) == sum(v_nhs[4,r,h+BREAKLM[2]] for r in h+1+BREAKLM[2]:noOfTimeSteps+2))

    @constraint(m, ending[h in TIMESTEPH],
        sum(v_nhs[4,h,s] for s in 1:h-1) == sum(v_nhs[5,r,h] for r in h+1:noOfTimeSteps+2))
        
    @constraint(m, earliestBreak1[h in 4:noOfTimeSteps+2],
        sum(v_nhs[2,h,s] for s in h-3:h-1) == 0)

    @constraint(m, lastBreakNotTooLate[h in 4:noOfTimeSteps+2],
        sum(v_nhs[4,h,s] for s in h-3:h-1) == 0)

    @constraint(m, earliestBreak2[h in 6:noOfTimeSteps+2],
        sum(v_nhs[3,h,s] for s in h-5:h-1) == 0)
        
    @constraint(m, latestBreak1[h in 12:noOfTimeSteps+2],
        sum(v_nhs[2,h,s] for s in 1:h-11) == 0)

    @constraint(m, latestBreak2[h in 13:noOfTimeSteps+2],
        sum(v_nhs[3,h,s] for s in 1:h-12) == 0)

    @constraint(m,
        v_nhs[2,2,1] + v_nhs[3,2,1] + v_nhs[4,2,1] == 0)

    @constraint(m,
        v_nhs[2,3,1] + v_nhs[3,3,1] + v_nhs[4,3,1] == 0)

    @constraint(m,
        v_nhs[2,3,2] + v_nhs[3,3,2] + v_nhs[4,3,2] == 0)

    @constraint(m,
        v_nhs[3,4,1] + v_nhs[4,4,1] == 0)
        
    @constraint(m,
        v_nhs[3,5,1] + v_nhs[4,5,1] == 0)

    @constraint(m,
        v_nhs[3,5,2] + v_nhs[4,5,2] == 0)

    @constraint(m,
        v_nhs[3,5,3] + v_nhs[4,5,3] == 0)

    @constraint(m,
        v_nhs[3,5,4] + v_nhs[4,5,4] == 0)

    @constraint(m,
        v_nhs[3,4,2] + v_nhs[4,4,2] == 0)

    @constraint(m,
        v_nhs[3,4,3] + v_nhs[4,4,3] == 0)

    @constraint(m,
        17 <= sum(s*v_nhs[5, noOfTimeSteps+2, s] for s in 1:noOfTimeSteps+1) - sum(h*v_nhs[1,h,1] for h in 2:noOfTimeSteps+2) <= 23)

    @constraint(m,
        sum(h*v_nhs[1,h,s] for h in 2:noOfTimeSteps+2, s in 1:h-1) - sum(s*v_nhs[1,h,s] for h in 2:noOfTimeSteps+2, s in 1:h-1)
        + sum(h*v_nhs[2,h,s] for h in 2:noOfTimeSteps+2, s in 1:h-1) - sum(s*v_nhs[2,h,s] for h in 2:noOfTimeSteps+2, s in 1:h-1)
        + sum(h*v_nhs[3,h,s] for h in 2:noOfTimeSteps+2, s in 1:h-1) - sum(s*v_nhs[3,h,s] for h in 2:noOfTimeSteps+2, s in 1:h-1)
        + sum(h*v_nhs[4,h,s] for h in 2:noOfTimeSteps+2, s in 1:h-1) - sum(s*v_nhs[4,h,s] for h in 2:noOfTimeSteps+2, s in 1:h-1)
        + sum(h*v_nhs[5,h,s] for h in 2:noOfTimeSteps+2, s in 1:h-1) - sum(s*v_nhs[5,h,s] for h in 2:noOfTimeSteps+2, s in 1:h-1) == noOfTimeSteps+2 - 3 - 1) #3 time steps are breaks, and we start on 1
    
    
    #@constraint(m,
      #  19 <= sum(s*v_nhs[4, 47, s] for s in 1:46) - sum(h*v_nhs[1,h,1] for h in 2:47) <= 32)  # for some reason the problem is infeasible when including this constr and the one below
    
    #@constraint(m,
     #   10 <= sum(s*v_nhs[3, 47, s] for s in 1:46) - sum(h*v_nhs[1,h,1] for h in 2:47) <= 19)

    
    #@constraint(m, sum(v_nhs[1,h,s] for h in 2:47, s in 1:h-1) == 1)
    @constraint(m, sum(v_nhs[2,h,s] for h in 2:noOfTimeSteps+2, s in 1:h-1) == 1)
    @constraint(m, sum(v_nhs[3,h,s] for h in 2:noOfTimeSteps+2, s in 1:h-1) == 1)
    @constraint(m, sum(v_nhs[4,h,s] for h in 2:noOfTimeSteps+2, s in 1:h-1) == 1)
    @constraint(m, sum(v_nhs[5,h,s] for h in 2:noOfTimeSteps+2, s in 1:h-1) == 1)

    set_optimizer(m, Gurobi.Optimizer)
    optimize!(m)
    println("Termination status, problem: ", termination_status(m))
        


    val = value.(v_nhs)
    obj = objective_value(m)
    a = value.(a_h)

    return val, obj, a
end

function GetMediumShift(mediumArcs)
    indexMatrix = zeros(Int, 3,1)
    for n in WORKBLOCK_NBM
        for h in 2:noOfTimeSteps+2
            for s in 1:h-1
                if (RoundTol(mediumArcs[n,h,s]) == 1) #rounding is needed due to numerical errors
                    indexMatrix = hcat(indexMatrix, [n,h,s])
                end
            end
        end
    end
    return indexMatrix[:, 2:end]
end

function RoundTol(x::Real)
    tol = 1e-3
    if abs(x - round(x)) <= tol
        return round(x)
    else
        return 0
    end
end


function TranslateToMediumShift(indexMatrix)
    shift = zeros(Int, 1,noOfTimeSteps)
    block1 = indexMatrix[:,2]
    block2 = indexMatrix[:,3]
    block3 = indexMatrix[:,4]
    for i in block1[3]:block1[2]-1
        shift[i] = 1
    end
    for i in block2[3]:block2[2]-1
        shift[i] = 1
    end
    for i in block3[3]:block3[2]-1
        shift[i] = 1
    end
    return shift
end

#u = vec(3*rand(1, noOfTimeSteps))
#u = vec([0.038461538461538464 0.0 0.0 0.029461538461538456 0.0 0.0 0.0 0.0 0.0 0.0 0.003538461538461539 0.0 0.0 0.0 0.0 0.009000000000000008 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.032999999999999995 0.009000000000000008 0.0 0.0 0.0 0.0 0.0 0.0 0.0])
function BestMediumShift(u, lambda)
    mediumArcs = BuildAndSolveMedium(u, lambda)
    #mediumShift = GetMediumShift(mediumArcs[1])
    #translatedMediumShift = TranslateToMediumShift(mediumShift)
    return round.(Int, mediumArcs[3]), mediumArcs[2]
end

#B = round.(Int, BuildAndSolveMedium(u, 0.7)[3])'
#mediumShift = GetMediumShift(mediumArcs)
#translatedMediumShift = TranslateToMediumShift(mediumShift)

#A = BestMediumShift(u, 0.7)[1]'
