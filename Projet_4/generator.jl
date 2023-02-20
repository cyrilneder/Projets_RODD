# Instance generator for forest management

using Random

function generate_instance(M::Int, N::Int)
    t = rand(60:100, M, N)
    open(string(@__DIR__, "\\ExplForet$(M)x$(N).txt"), "w") do file
        println(file, "t = $t")
    end
end

# generate_instance(50, 50)
