using JuMP
using CPLEX
using PlotlyJS

path = "data.txt"

function static_solving(inputFile::String)

    include(inputFile)

    println(prob1)
    println(prob1[6,2])

    m = Model(CPLEX.Optimizer)

    #Variables
    @variable(m, x[i in 1:n, j in 1:n], Bin)
    @variable(m, y[i in 1:n, j in 1:n], Bin)

    #Contraintes espèces protégées
    #@constraint(m, [k in Pd], sum(y[i,j] for (i,j) in S[k]) >= 1)

    @constraint(m, sum(y[i,j]*log(1-prob1[i,j]) for i in 1:n, j in 1:n) <= log(1-alpha[1]))
    @constraint(m, sum(y[i,j]*log(1-prob2[i,j]) for i in 1:n, j in 1:n) <= log(1-alpha[2]))
    @constraint(m, sum(y[i,j]*log(1-prob3[i,j]) for i in 1:n, j in 1:n) <= log(1-alpha[3]))

    #Contraintes espèces communes
    #@constraint(m, [k in Pc], sum(x[i,j] for (i,j) in S[k]) >= 1)

    @constraint(m, sum(x[i,j]*log(1-prob4[i,j]) for i in 1:n, j in 1:n) <= log(1-alpha[4]))
    @constraint(m, sum(x[i,j]*log(1-prob5[i,j]) for i in 1:n, j in 1:n) <= log(1-alpha[5]))
    @constraint(m, sum(x[i,j]*log(1-prob6[i,j]) for i in 1:n, j in 1:n) <= log(1-alpha[6]))

    #Contraintes zones centrales
    
        # Pas de zones centrales sur les bords

        for i in 1:n
            @constraint(m, y[1,i] == 0)
            @constraint(m, y[10,i] == 0)
            @constraint(m, y[i,1] == 0)
            @constraint(m, y[i,10] == 0)
        end
    
        # Zones centrales encerclées par zones communes
    for i in 2:(n-1)
        for j in 2:(n-1)
            @constraint(m, y[i,j] <= x[i-1,j-1])
            @constraint(m, y[i,j] <= x[i-1,j])
            @constraint(m, y[i,j] <= x[i-1,j+1])
            @constraint(m, y[i,j] <= x[i,j-1])
            @constraint(m, y[i,j] <= x[i,j])
            @constraint(m, y[i,j] <= x[i,j+1])
            @constraint(m, y[i,j] <= x[i+1,j-1])
            @constraint(m, y[i,j] <= x[i+1,j])
            @constraint(m, y[i,j] <= x[i+1,j+1])
        end
    end

    #Objectif
    @objective(m, Min, sum(c[i,j]*x[i,j] for i in 1:n, j in 1:n))

    #Résolution
    optimize!(m)

    # Récupération du status de la résolution
    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL
    if feasibleSolutionFound
            # Récupération des valeurs d’une variable
            vx = JuMP.value.(x)
            vy = JuMP.value.(y)
            vOpt = JuMP.objective_value(m)
    end

    #Affichage de solution
    println("Valeur optimale :", vOpt)
    proba_survie1 = 1 - prod([1 - vy[i,j]*prob2[i,j] for i in 1:n, j in 1:n])
    println("Proba de survie pour k = 1 : ",proba_survie1)
    plot(heatmap(z=vx, x = 1:n, y = 1:n))
    #plot(heatmap(z=prob1, x = 1:p, y = 1:p))
    #lot(heatmap(z=prob5, x = 1:p, y = 1:p))
    #plot(heatmap(z=prob6, x = 1:p, y = 1:p))
end

static_solving(path)