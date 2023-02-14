include("generic_solver.jl")

p_min, p_max = 10, 20

n = 10

c_min, c_max = 3, 8

nb_iter = 10

avg_time_table = zeros(p_max-p_min+1)

@simd for p in p_min:p_max

    println("Current status : p = ", p)

    q = trunc(Int, p/2)

    s = 0
    counter = 0
    while counter <= nb_iter
        time = generic_solving(n,n,p,q,c_min,c_max)
        if time != -1
            s += time
            counter += 1
        end
    end
    avg_time_table[p-p_min+1] = s/nb_iter
end

println(avg_time_table)