using JuMP
using CPLEX


function uls_glissant(L, T)
    
    M = 6

    E_max = 3
    
    d = 20 .+ rand(T)*(70 - 20)
    
    f = [10, 30, 60, 75, 90, 105]
    e = [8, 6, 4, 3, 2, 1]
    
    h_t = 1
    p_t = 0

    model = Model(CPLEX.Optimizer)

    #Variables
    @variable(model, x[t in 1:T, m in 1:M] >= 0)
    @variable(model, y[t in 1:T, m in 1:M], Bin)
    @variable(model, s[t in 0:T] >= 0)

    #Contrainte de bilan des stocks
    @constraint(model, s[0] == 0)
    @constraint(model, [t in 1:T], sum(x[t,m] for m in 1:M) - s[t] + s[t-1] == d[t])

    #Contrainte activation des y[t,m]
    @constraint(model, [t in 1:T, m in 1:M], x[t,m] <= sum(d[t1] for t1 in t:T) * y[t,m])
    

    #Contrainte d'émission carbone
    @constraint(model, [t in L:T], sum((E_max - e[m]) * x[t1,m] for m in 1:M, t1 in (t-L+1):t) >= 0)

    #Objectif
    @objective(model, Min, sum(p_t * x[t,m] + f[m] * y[t,m] for t in 1:T, m in 1:M) + sum(h_t * s[t] for t in 1:T))

    
    # Désactive les sorties de CPLEX (optionnel)
    set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0)
    set_time_limit_sec(model, 90)
    
    #Résolution
    optimize!(model)

    # Récupération du status de la résolution
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL
    if feasibleSolutionFound
        # Récupération des valeurs d’une variable
        vx = JuMP.value.(x)
        vOpt = JuMP.objective_value(model)

        emission_moyenne = 1/T * sum(e[m]*vx[t,m] for t in 1:T, m in 1:M)
        return vOpt, emission_moyenne
    else 
        println("Infaisable ou non borné")
    end
end

# nb_iter = 100
# T = 12
# L_list = [L for L in 1:T,i in 1:nb_iter]

# println(L_list)

uls_glissant(4, 12)

# truc = uls_glissant.(L_list,T)
# println(truc)
