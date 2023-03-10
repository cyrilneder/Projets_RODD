using JuMP
using CPLEX
using PlotlyJS

path = "data.txt"

function static_solving(inputFile::String)

    include(inputFile)

    #Préparation des données
    proba_list = [prob1, prob2, prob3, prob4, prob5, prob6]

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

    #Résolution
    optimize!(model)

    # Récupération du status de la résolution
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL
    if feasibleSolutionFound
            # Récupération des valeurs d’une variable
            vx = JuMP.value.(x)
            vy = JuMP.value.(y)
            vOpt = JuMP.objective_value(model)
    end

    #Affichage de solution
    println("Valeur optimale :", vOpt)

    proba_survie = zeros(p)
    for k in 1:q
        proba_survie[k] = 1 - prod([1 - vy[i,j]*proba_list[k][i,j] for i in 1:n, j in 1:m])
    end

    for k in q+1:p
        proba_survie[k] = 1 - prod([1 - vx[i,j]*proba_list[k][i,j] for i in 1:n, j in 1:m])
    end

    println("Proba de survie pour chaque espèce : ",proba_survie)

    plot(heatmap(z=(vx+vy)[end:-1:1,:], x = 1:n, y = 1:m))
end

static_solving(path)