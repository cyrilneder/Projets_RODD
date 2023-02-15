# Durable forest exploitation

using JuMP
using CPLEX
using PlotlyJS

include("ExplForet_ampl.txt")

N, M = size(t)
w1 = 1
w2 = 5
gL = 3 * 1.26157

function get_A_ij(i, j)
    if i == 1
        if j == 1
            return [(i+1, j), (i, j+1)]
        elseif j == N 
            return [(i+1, j), (i, j-1)]
        else
            return [(i+1, j), (i, j+1), (i, j-1)]
        end
    elseif i == M
        if j == 1
            return [(i-1, j), (i, j+1)]
        elseif j == N 
            return [(i-1, j), (i, j-1)]
        else
            return [(i, j-1), (i, j+1), (i-1, j)]
        end
    else
        if j == 1
            return [(i-1, j), (i+1, j), (i, j+1)]
        elseif j == N 
            return [(i-1, j), (i, j-1), (i+1, j)]
        else
            return [(i-1, j), (i, j-1), (i+1, j), (i, j+1)]
        end
    end
end

function first_model(max_coupe::Bool=false)
    m = Model(CPLEX.Optimizer)

    # Variables
    @variable(m, x[i in 1:M, j in 1:N], Bin)
    @variable(m, d[i in 1:M, j in 1:N] >= 0)

    # Contraintes zones coupées adjascentes
    @constraint(m, [i in 1:M, j in 1:N], 
        sum(x[k,l] for (k,l) in get_A_ij(i,j)) - length(get_A_ij(i,j)) * (1 - x[i,j]) <= d[i,j]
    )

    # Contraintes minimum de zones non-coupées
    if max_coupe
        @constraint(m, sum(x) >= 60)
    end

    #Objectif
    @objective(m, Max, 
        w1 * sum(t .* (1 .- x)) + w2 * gL * sum((4 .* x) .- d)
    )

    #Résolution
    optimize!(m)

    # Récupération du status de la résolution
    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL
    if feasibleSolutionFound
            # Récupération des valeurs d’une variable
            vx = JuMP.value.(x)
            vx = round.(Int, vx)
            vd = JuMP.value.(d)
            vOpt = JuMP.objective_value(m)
    end

    #Affichage de solution
    println("Valeur optimale :", round(vOpt, digits=1))
    println("Solution x :")
    display(vx)
    return vx
end

function second_model(max_coupe::Bool=false)
    m = Model(CPLEX.Optimizer)

    # Variables
    @variable(m, x[i in 1:M, j in 1:N], Bin)
    @variable(m, y[i in 1:M, j in 1:N, (k,l) in get_A_ij(i,j)] >= 0)

    # Contraintes linéarisation
    for i in 1:M, j in 1:N
        for (k,l) in get_A_ij(i,j)
            @constraint(m, 1 + y[i, j, (k, l)] >= x[i,j] + x[k,l])
        end
    end

    # Contraintes minimum de zones non-coupées
    if max_coupe
        @constraint(m, sum(x) >= 60)
    end

    #Objectif
    @objective(m, Max, 
        w1 * sum(t .* (1 .- x)) + w2 * gL * (sum(4 .* x) - sum(y))
    )

    write_to_file(m, "model.lp")

    #Résolution
    optimize!(m)

    # Récupération du status de la résolution
    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL
    if feasibleSolutionFound
            # Récupération des valeurs d’une variable
            vx = JuMP.value.(x)
            vx = round.(Int, vx)
            vy = JuMP.value.(y)
            vOpt = JuMP.objective_value(m)
    end

    #Affichage de solution
    println("Valeur optimale :", round(vOpt, digits=1))
    println("Solution x :")
    display(vx)
    return vx
end

# plot(heatmap(z=first_model()[end:-1:1,:]))
plot(heatmap(z=second_model(true)[end:-1:1,:]))
