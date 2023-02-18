
function instance_genarator(N, G, A)
    
    individu = []
    
    for ind in 1:N 
        genotype = []
        
        for g in 1:G
            for a in 1:A
                push!(genotype, rand(1:A,A))        
            end
        end
        push!(individu, [genotype])
    end
    
    return individu
end