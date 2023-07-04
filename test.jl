a_hj = r[2]
master = Model() #Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
@variable(m, w_kj[COMPETENCE_K, 1:size(SHIFT_Jc,2)] >= 0, Int)
@variable(m, yu_hk[TIMESTEP_H, COMPETENCE_K] >= 0)


@objective(master, Min, 
lambda*(sum(sum(yu_hk[h,k] for h in TIMESTEP_H) for k in COMPETENCE_K)) +
(1-lambda)*(sum(sum(sum((a_hj[h,j])*w_kj[k,j] for h in TIMESTEP_H) for j in 1:size(SHIFT_Jc,2)) for k in COMPETENCE_K) +
sum(sum(((findlast(x -> x == 1, BitArray(SHIFT_Jc[:,j])) - findnext(BitArray(SHIFT_Jc[:,j]), 1)) >= 24 ? 2*(1-lambda)*w_kj[k,j] : (1-lambda)*w_kj[k,j]) for k in COMPETENCE_K) for j in 1:size(SHIFT_Jc,2)))
)

#@objective(m, Min, sum(sum(q*yu_hk[h,k] for h in TIMESTEP_H) for k in COMPETENCE_K) 
 #   + sum(sum(sum((a_hj[h,j]*c)*w_kj[k,j] for h in TIMESTEP_H) for j in 1:size(SHIFT_Jc,2)) for k in COMPETENCE_K)
  #  )

@constraint(master, availability[h in TIMESTEP_H, k in COMPETENCE_K],
    sum(a_hj[h,j]*w_kj[k,j] for j in 1:size(SHIFT_Jc,2)) + yu_hk[h,k] >= D_hk[h,k])

@constraint(master, limitedAgentsPerCompetence[k in COMPETENCE_K],
    sum(w_kj[k,j] for j in 1:size(SHIFT_Jc,2)) <= noOfAgentsPerComp[k])

@constraint(master, limitedAgentsPerDay,
    sum(sum(w_kj[k,j] for j in 1:size(SHIFT_Jc,2)) for k in COMPETENCE_K) <= size(AGENT_I,1))


set_optimizer(master, Gurobi.Optimizer)

unregister(m, :w_kj)
unregister(m, :availability)


function BnP(fracW, SHIFT_Jc, D_hk, lambda)
while !is_integer_feasible
    a_hj = SHIFT_Jc

    # Select a fractional variable from the master problem solution
    fracVar = findfirst(x -> !isinteger(x), fracW)
    fracVarIndex = findfirst(isequal(fracVar), fracW)
    
    # Create child subproblems by branching on the fractional variable
    subProblem1 = Model()
    subProblem2 = Model()

    # Create variables for subproblems
    @variable(subProblem1, w_kj[COMPETENCE_K, 1:size(SHIFT_Jc,2)] >= 0)
    @variable(subProblem1, yu_hk[TIMESTEP_H, COMPETENCE_K] >= 0)
    @variable(subProblem2, w_kj[COMPETENCE_K, 1:size(SHIFT_Jc,2)] >= 0)
    @variable(subProblem2, yu_hk[TIMESTEP_H, COMPETENCE_K] >= 0)

    # Define the subproblem objectives
    @objective(subProblem1, Min, 
        lambda*(sum(sum(yu_hk[h,k] for h in TIMESTEP_H) for k in COMPETENCE_K)) +
        (1-lambda)*(sum(sum(sum((a_hj[h,j])*w_kj[k,j] for h in TIMESTEP_H) for j in 1:size(SHIFT_Jc,2)) for k in COMPETENCE_K) +
        sum(sum(((findlast(x -> x == 1, BitArray(SHIFT_Jc[:,j])) - findnext(BitArray(SHIFT_Jc[:,j]), 1)) >= 24 ? 2*(1-lambda)*w_kj[k,j] : (1-lambda)*w_kj[k,j]) for k in COMPETENCE_K) for j in 1:size(SHIFT_Jc,2)))
    )

    @objective(subProblem2, Min, 
        lambda*(sum(sum(yu_hk[h,k] for h in TIMESTEP_H) for k in COMPETENCE_K)) +
        (1-lambda)*(sum(sum(sum((a_hj[h,j])*w_kj[k,j] for h in TIMESTEP_H) for j in 1:size(SHIFT_Jc,2)) for k in COMPETENCE_K) +
        sum(sum(((findlast(x -> x == 1, BitArray(SHIFT_Jc[:,j])) - findnext(BitArray(SHIFT_Jc[:,j]), 1)) >= 24 ? 2*(1-lambda)*w_kj[k,j] : (1-lambda)*w_kj[k,j]) for k in COMPETENCE_K) for j in 1:size(SHIFT_Jc,2)))
    )

    # Add constraints to subproblems
    @constraint(subProblem1, availability[h in TIMESTEP_H, k in COMPETENCE_K],
        sum(a_hj[h,j]*w_kj[k,j] for j in 1:size(SHIFT_Jc,2)) + yu_hk[h,k] >= D_hk[h,k])
    @constraint(subProblem1, limitedAgentsPerCompetence[k in COMPETENCE_K],
        sum(w_kj[k,j] for j in 1:size(SHIFT_Jc,2)) <= noOfAgentsPerComp[k])
    @constraint(subProblem1, limitedAgentsPerDay,
        sum(sum(w_kj[k,j] for j in 1:size(SHIFT_Jc,2)) for k in COMPETENCE_K) <= size(AGENT_I,1))
    @constraint(subproblem1, x[fracVarIndex] <= floor(Int, fracVar))

    @constraint(subProblem2, availability[h in TIMESTEP_H, k in COMPETENCE_K],
        sum(a_hj[h,j]*w_kj[k,j] for j in 1:size(SHIFT_Jc,2)) + yu_hk[h,k] >= D_hk[h,k])
    @constraint(subProblem2, limitedAgentsPerCompetence[k in COMPETENCE_K],
        sum(w_kj[k,j] for j in 1:size(SHIFT_Jc,2)) <= noOfAgentsPerComp[k])
    @constraint(subProblem2, limitedAgentsPerDay,
        sum(sum(w_kj[k,j] for j in 1:size(SHIFT_Jc,2)) for k in COMPETENCE_K) <= size(AGENT_I,1))
    @constraint(subproblem2, x[fracVarIndex] >= ceil(Int, fracVar) + 1)


    set_optimizer(subProblem1, Gurobi.Optimizer)
    set_optimizer(subProblem2, Gurobi.Optimizer)
    
    # Solve the subproblems with column generation and update the master problem with the subproblem solutions
    optimize!(subProblem1)
    status1 = termination_status(subProblem1)
    optimize!(subProblem2)
    status2 = termination_status(subProblem2)
    if termination_status(status1) == MOI.OPTIMAL
        sol1 = value.(subProblem1[:w_kj])
        @constraint(master, w_kj[fracVarIndex] == sol1[fracVarIndex])

    elseif termination_status(status2) == MOI.OPTIMAL
        sol2 = value.(subProblem2[:w_kj])
        @constraint(master, w_kj[fracVarIndex] == sol1[fracVarIndex])
    end

    optimize!(master)

      # Retrieve the updated master problem solution
    master_solution = value.(x)
    println("Master Problem Solution: ", master_solution)

    # Check if the updated solution is integer-feasible
    is_integer_feasible = all(x -> isinteger(x), master_solution)

    @constraint(master, w_kj[fracVarIndex] == sol2[fracVarIndex])
end
return objective_value(master), value.(w_kj)
end