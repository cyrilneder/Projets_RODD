using JuMP
using CPLEX

path = "DivGenetique.txt"
include(path)


h = 50
theta1 = 0.001

# Renvoie true si l'individu i possède deux fois le même allèle sur le gène g
function two_a(individu,i,g,a)
    return (individu[i][1][g][1]==a)&&(individu[i][1][g][2]==a)
end

# Renvoie le nombre d'allèle a sur le gène g de l'individu i
function nb_a(individu,i,g,a)
    return length(findall(x->x==a,individu[i][1][g]))
end

function model(inputFile::String, h::Int64, theta1::Float64)

    include(inputFile)

    m = Model(CPLEX.Optimizer)

    #Variables
    @variable(m, 0 <= x[i in 1:N] <= 3, Int)
    @variable(m, p[g in 1:G, a in 1:A] >= 0)
    @variable(m, t[g in 1:G, a in 1:A] >= 0)

    #Contraintes espèces protégées

    @constraint(m, sum(x[i] for i in 1:Nm) == N)
    @constraint(m, sum(x[i] for i in Nm+1:Nm+Nf) == N)


    @constraint(m, [g in 1:G, a in 1:A], p[g,a] >= t[g,a] - sum(x[i] for i in 1:N if two_a(individu,i,g,a)))


    for r in 1:h
        theta = theta1^((h-r)/(h-1))
        println(theta)
        @constraint(m, [g in 1:G, a in 1:A], log(theta) + 1/theta * (t[g,a]-theta)  >= sum(x[i]*log(1-nb_a(individu,i,g,a)/2) for i in 1:N if !(two_a(individu,i,g,a))))
    end
    
    #Objectif
    @objective(m, Min, sum(p[g,a] for g in 1:G, a in 1:A))

    #Résolution
    optimize!(m)

    # Récupération du status de la résolution
    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL
    if feasibleSolutionFound
            # Récupération des valeurs d’une variable
            vx = JuMP.value.(x)
            vOpt = JuMP.objective_value(m)
    end

    #Affichage de solution
    println("Valeur optimale :", vOpt)
    println("Solution x :", vx)
    proba_b = prod([(1 - nb_a(individu,i,2,2)/2)^vx[i] for i in 1:N])
    proba_e = prod([(1 - nb_a(individu,i,5,2)/2)^vx[i] for i in 1:N])
    proba = sum(prod([(1 - nb_a(individu,i,g,a)/2)^vx[i] for i in 1:N]) for g in 1:G, a in 1:A)
    println(proba_b)
    #println("Proba de survie pour k = 1 : ",proba_survie1)
    #plot(heatmap(z=(vx+vy)[end:-1:1,:], x = 1:n, y = 1:n))
    #plot(heatmap(z=prob1, x = 1:p, y = 1:p))
    #lot(heatmap(z=prob5, x = 1:p, y = 1:p))
    #plot(heatmap(z=prob6, x = 1:p, y = 1:p))
end

model(path,h,theta1)