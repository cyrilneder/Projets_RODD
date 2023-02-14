
function instance_genarator(n,m,p, q, c_min,c_max)
    alpha = rand([0.5,0.6,0.7,0.8,0.9],p)
    
    proba_list = []
    for k in 1:p
        probk = rand([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0.2, 0.3, 0.4, 0.5], (n, m))
        if k <= q
            probk[1,:] = zeros(m)
            probk[n,:] = zeros(m)
            probk[:,1] = zeros(n)
            probk[:,m] = zeros(n)
        end
        push!(proba_list,probk)
    end

    c = rand(c_min:c_max, (n,m))
    return alpha, proba_list, c
end