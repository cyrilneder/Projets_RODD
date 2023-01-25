using JuMP
using CPLEX
using PlotlyJS

path = "data.txt"
include(path)

function static_solving(inputFile::String)

    include(inputFile)

    Pd = 1:q
    Pc = (q+1):p

    S1 = Tuple.(findall(!iszero, prob1))
    S2 = Tuple.(findall(!iszero, prob2))
    S3 = Tuple.(findall(!iszero, prob3))
    S4 = Tuple.(findall(!iszero, prob4))
    S5 = Tuple.(findall(!iszero, prob5))
    S6 = Tuple.(findall(!iszero, prob6))

    S = [S1, S2, S3, S4, S5, S6]

    m = Model(CPLEX.Optimizer)

    #Variables
    @variable(m, x[i in 1:n, j in 1:n], Bin)
    @variable(m, y[i in 1:n, j in 1:n], Bin)

    #Contraintes espèces protégées
    @constraint(m, [k in Pd], sum(y[i,j] for (i,j) in S[k]) >= 1)

    @constraint(m, sum(y[i,j]*log(1-proba1[i,j]) for (i,j) in S1) <= log(1-alpha[1]))
    @constraint(m, sum(y[i,j]*log(1-proba2[i,j]) for (i,j) in S2) <= log(1-alpha[2]))
    @constraint(m, sum(y[i,j]*log(1-proba3[i,j]) for (i,j) in S3) <= log(1-alpha[3]))

    #Contraintes espèces communes
    @constraint(m, [k in Pc], sum(x[i,j] for (i,j) in S[k]) >= 1)

    @constraint(m, sum(x[i,j]*log(1-proba4[i,j]) for (i,j) in S4) <= log(1-alpha[4]))
    @constraint(m, sum(x[i,j]*log(1-proba5[i,j]) for (i,j) in S5) <= log(1-alpha[5]))
    @constraint(m, sum(x[i,j]*log(1-proba6[i,j]) for (i,j) in S6) <= log(1-alpha[6]))

    #Contraintes zones centrales
    for i in 2:(n-1)
        for j in 2:(n-1)
            @constraint(m, x[i-1,j-1] >= y[i,j])
            @constraint(m, x[i-1,j] >= y[i,j])
            @constraint(m, x[i-1,j+1] >= y[i,j])
            @constraint(m, x[i,j-1] >= y[i,j])
            @constraint(m, x[i,j] >= y[i,j])
            @constraint(m, x[i,j+1] >= y[i,j])
            @constraint(m, x[i+1,j-1] >= y[i,j])
            @constraint(m, x[i+1,j] >= y[i,j])
            @constraint(m, x[i+1,j+1] >= y[i,j])
        end
    end

    #Objectif
    @objective(m, Min, sum(c[i,j]*x[i,j] for i in 1:n,j in 1:n))

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
    println("x = ", vx)
    println("y = ", vy)
    plot(heatmap(z=vx))
    plot(heatmap(z=vy))
    #plot(heatmap(z=vx+vy))
end

static_solving(path)