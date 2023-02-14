using JuMP
using CPLEX
using PlotlyJS
include("instance_generator.jl")


function generic_solving(n, m, p, q, c_min, c_max)


    #Préparation des données
    alpha, proba_list, c = instance_genarator(n,m,p,q,c_min,c_max)

    model = Model(CPLEX.Optimizer)

    #Variables
    @variable(model, x[i in 1:n, j in 1:m], Bin)
    @variable(model, y[i in 1:n, j in 1:m], Bin)

    #Contraintes espèces protégées
    @constraint(model, [k in 1:q], sum(y[i,j]*log(1-proba_list[k][i,j]) for i in 1:n, j in 1:m) <= log(1-alpha[k]))

    #Contraintes espèces communes
    @constraint(model, [k in q+1:p], sum(x[i,j]*log(1-proba_list[k][i,j]) for i in 1:n, j in 1:m) <= log(1-alpha[k]))

    #Contraintes zones centrales
    
    # Pas de zones centrales sur les bords
    for i in 1:n
        @constraint(model, y[i,1] == 0)
        @constraint(model, y[i,m] == 0)
    end

    for j in 1:m
        @constraint(model, y[1,j] == 0)
        @constraint(model, y[n,j] == 0)
    end

        # Zones centrales encerclées par zones communes
    for i in 2:(n-1)
        for j in 2:(m-1)
            @constraint(model, y[i,j] <= x[i-1,j-1])
            @constraint(model, y[i,j] <= x[i-1,j])
            @constraint(model, y[i,j] <= x[i-1,j+1])
            @constraint(model, y[i,j] <= x[i,j-1])
            @constraint(model, y[i,j] <= x[i,j])
            @constraint(model, y[i,j] <= x[i,j+1])
            @constraint(model, y[i,j] <= x[i+1,j-1])
            @constraint(model, y[i,j] <= x[i+1,j])
            @constraint(model, y[i,j] <= x[i+1,j+1])
        end
    end

    #Objectif
    @objective(model, Min, sum(c[i,j]*x[i,j] for i in 1:n, j in 1:m))

    # Désactive les sorties de CPLEX (optionnel)
    set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0)

    start = time()
    #Résolution
    optimize!(model)

    computation_time = time() - start

    # Récupération du status de la résolution
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL
    if !feasibleSolutionFound 
        println("Infaisable")
        return -1
    end

    return computation_time

end
