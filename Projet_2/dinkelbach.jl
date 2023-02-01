using JuMP
using CPLEX
using IterTools
using PlotlyJS

include("data.txt")

println("Instance : ")
display(cout)
println()

parameters = [(30,35,92), (20,21,52), (70,75,350)]

function get_distance(i::Int, j::Int, k::Int, l::Int)
    return sqrt((i-k)^2 + (j-l)^2)
end

n = 10

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
    @constraint(model, sum(cout[i,j] * x[i,j] for i in 1:n for j in 1:n) <= B)
    @constraint(model, [i in 1:n, j in 1:n], sum(y[i,j,k,l] for k in 1:n for l in 1:n if k != i || l != j) == x[i,j])
    @constraint(model, [i in 1:n, j in 1:n, k in 1:n, l in 1:n], y[i,j,k,l] <= x[i,j])
    @constraint(model, [i in 1:n, j in 1:n, k in 1:n, l in 1:n], y[i,j,k,l] <= x[k,l])
    # Objective
    @objective(model, Max, sum(distances[i,j,k,l] * y[i,j,k,l] for i in 1:n for j in 1:n for k in 1:n for l in 1:n) - lambda * sum(x))
    # write_to_file(model, "model.lp")
    return model
end

function compute_v_lambda(lambda::Real, A_min::Real, A_max::Real, B::Real)
    # Solving model
    model = get_linear_model(lambda, A_min, A_max, B)
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
    f_x_lambda = sum(distances .* y_lambda)
    g_x_lambda = lambda * sum(x_lambda)

    return v_lambda, x_lambda, f_x_lambda, g_x_lambda
end

# v_lambda, x_lambda, f_x_lambda, g_x_lambda = compute_v_lambda(1, 30, 35, 920)
# println(v_lambda)
# display(x_lambda)
# println
# println(f_x_lambda)
# println(g_x_lambda)

function dinkelbach(lambda_0::Real, A_min::Real, A_max::Real, B::Real)
    # Lambda initialization
    lambda = lambda_0
    # Computing first v_lambda and getting the solution x_lambda
    v_lambda, x_lambda, f_x_lambda, g_x_lambda = compute_v_lambda(lambda, A_min, A_max, B)
    # If v_lambda > 0, the solution isn't optimal
    while v_lambda > 0
        # Updating lambda
        lambda = f_x_lambda / g_x_lambda
        # Computing again v_lambda and related values
        v_lambda, x_lambda, f_x_lambda, g_x_lambda = compute_v_lambda(lambda, A_min, A_max, B)
    end
    plot(heatmap(z=(x_lambda)[end:-1:1,:]))
    return x_lambda
end

display(dinkelbach(1, 30, 35, 92))
println()

# println("Instances resolution : ")
# for (A_min, A_max, B) in parameters
#     println("Parameters A_min = $A_min, A_max = $A_max, B = $B")
#     display(dinkelbach(1, A_min, A_max, B))
# end