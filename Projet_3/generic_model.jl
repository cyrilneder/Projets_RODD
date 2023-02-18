using JuMP
using CPLEX
include("instance_generator.jl")

# Renvoie true si l'individu i possède deux fois le même allèle sur le gène g
function two_a(individu,i,g,a)
    return (individu[i][1][g][1]==a)&&(individu[i][1][g][2]==a)
end

# Renvoie le nombre d'allèle a sur le gène g de l'individu i
function nb_a(individu,i,g,a)
    return length(findall(x->x==a,individu[i][1][g]))
end

function generic_model(N, G, A, T, init)

    Nm = trunc(Int,N/2)
    Nf = N - Nm

    individu = instance_genarator(N, G, A)

    m = Model(CPLEX.Optimizer)

    #Variables
    @variable(m, 0 <= x[i in 1:N] <= trunc(Int, 3/8 * N), Int)
    @variable(m, p[g in 1:G, a in 1:A] >= 0)
    @variable(m, t[g in 1:G, a in 1:A] >= 0)

    #Contraintes espèces protégées

    @constraint(m, sum(x[i] for i in 1:Nm) == N)
    @constraint(m, sum(x[i] for i in Nm+1:Nm+Nf) == N)


    @constraint(m, [g in 1:G, a in 1:A], p[g,a] >= t[g,a] - sum(x[i] for i in 1:N if two_a(individu,i,g,a)))


    for r in 1:T
        theta = init^((T-r)/(T-1))

        @constraint(m, [g in 1:G, a in 1:A], log(theta) + 1/theta * (t[g,a]-theta)  >= sum(x[i]*log(1-nb_a(individu,i,g,a)/2) for i in 1:N if !(two_a(individu,i,g,a))))
    end
    
    #Objectif
    @objective(m, Min, sum(p[g,a] for g in 1:G, a in 1:A))

    start = time()

    # Désactive les sorties de CPLEX (optionnel)
    set_optimizer_attribute(m, "CPX_PARAM_SCRIND", 0)

    #Résolution
    optimize!(m)

    computation_time = time()-start

    # Récupération du status de la résolution
    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL
    if feasibleSolutionFound
        return computation_time
    else 
        return -1
    end

end

function factor_analysis()
    
    generic_model(8,5,2,50,0.001)
    
    T_min, T_max = 100,1000
    pas = 100
    
    T_list = T_min:pas:T_max
    
    T_length = length(T_list)
    
    N = 8
    G = 5
    A = 2
    
    nb_iter = 30
    
    avg_time_table = zeros(T_length)
    
    counter = 0
    
    while counter <= nb_iter
        println("Current status counter : ", counter)
        time = @. generic_model(N,G,A,T_list,0.001*50/T_list)
        avg_time_table = avg_time_table + time
        counter += 1
    end
    
    avg_time_table /= nb_iter
    
    println(avg_time_table)
end


factor_analysis()