using JuMP
using CPLEX
using IterTools
using PlotlyJS
using Random

include("data.txt")

n = 25

println("Instance : ")
cout = rand(1:10, n, n)
display(cout)
println()

parameters = [(3*n,3.5*n,9.2*n), (2*n,2.1*n,5.2*n), (7*n,7.5*n,35*n)]

function get_distance(i::Int, j::Int, k::Int, l::Int)
    return sqrt((i-k)^2 + (j-l)^2)
end

println("Computing distances")
distances = zeros(n, n, n, n)
for i in 1:n, j in 1:n, k in 1:n, l in 1:n 
    distances[i, j, k, l] = get_distance(i, j, k, l)
end

function get_linear_model(lambda::Real, A_min::Real, A_max::Real, B::Real)
    # Model object
    model = Model(CPLEX.Optimizer)
    # Variables
    @variable(model, x[i in 1:n, j in 1:n], Bin)
    @variable(model, y[i in 1:n, j in 1:n, k in 1:n, l in 1:n], Bin)
    # Constraints
    @constraint(model, A_min <= sum(x) <= A_max)
    @constraint(model, sum(cout .* x) <= B)
    @constraint(model, [i in 1:n, j in 1:n], sum(y[i,j,k,l] for k in 1:n for l in 1:n if k != i || l != j) == x[i,j])
    @constraint(model, [i in 1:n, j in 1:n], y[i,j,:,:] .<= x[i,j])
    @constraint(model, [k in 1:n, l in 1:n], y[:,:,k,l] .<= x[k,l])
    # Objective
    @objective(model, Max, sum(distances .* (1 .- y)) - lambda * sum(x))
    write_to_file(model, "model.lp")
    return model
end

function compute_v_lambda(lambda::Real, A_min::Real, A_max::Real, B::Real)
    # Solving model
    model = get_linear_model(lambda, A_min, A_max, B)
    set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0)
    optimize!(model)

    feasible_solution_found = primal_status(model) == MOI.FEASIBLE_POINT
    is_optimal = termination_status(model) == MOI.OPTIMAL
    if !(feasible_solution_found || is_optimal)
        return -1, -1, -1, -1
    end

    # Getting values
    v_lambda = objective_value(model)
    x_lambda = value.(model[:x])
    y_lambda = value.(model[:y])
    f_x_lambda = sum(distances .* (1 .- y_lambda))
    g_x_lambda = sum(x_lambda)
    dmppv = sum(distances .* y_lambda) / g_x_lambda

    if v_lambda < 1e-5
        println("Found DMPPV : $dmppv")
    end

    return v_lambda, x_lambda, f_x_lambda, g_x_lambda, node_count(model)
end

# v_lambda, x_lambda, f_x_lambda, g_x_lambda = compute_v_lambda(1, 30, 35, 920)
# println(v_lambda)
# display(x_lambda)
# println
# println(f_x_lambda)
# println(g_x_lambda)

function dinkelbach(lambda_0::Real, A_min::Real, A_max::Real, B::Real)
    println("\nDinkelbach Result : lambda_0 = $lambda_0, A_min = $A_min, A_max = $A_max, B = $B")
    # Resolution statistics
    nb_iter, nb_nodes, start_time = 0, 0, time()
    # Lambda initialization
    lambda = lambda_0
    # Computing first v_lambda and getting the solution x_lambda
    v_lambda, x_lambda, f_x_lambda, g_x_lambda, nodes = compute_v_lambda(lambda, A_min, A_max, B)
    # Updating statistics
    nb_iter += 1
    nb_nodes += nodes
    # If v_lambda > 0, the solution isn't optimal
    while v_lambda > 1e-5
        # Updating lambda
        lambda = f_x_lambda / g_x_lambda
        # Computing again v_lambda and related values
        v_lambda, x_lambda, f_x_lambda, g_x_lambda, nodes = compute_v_lambda(lambda, A_min, A_max, B)
        # Updating statistics
        nb_iter += 1
        nb_nodes += nodes
    end
    time_spent = round(time() - start_time, digits=2)
    x_lambda = round.(Int, x_lambda)
    println("Number of selected boxes : $(sum(x_lambda))")
    println("Resolution statistics : iterations = $nb_iter, nodes = $nb_nodes, time = $(time_spent)s\n")
    return x_lambda
end

# x_lambda = dinkelbach(1, 70, 75, 350)
# display(x_lambda)

println()
# plot(heatmap(z=x_lambda[end:-1:1,:]))

println("Instances resolution : ")
for (A_min, A_max, B) in parameters
    println("Parameters A_min = $A_min, A_max = $A_max, B = $B")
    dinkelbach(1, A_min, A_max, B)
end