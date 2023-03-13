# Figure pour la q2

using Plots

include("uls_intervalle_glissant.jl")

T_vector = [T for T in 10:5:40]
C_vector = Float64[]
E_vector = Float64[]

for T in T_vector
    println("Tests for T = $T")
    C, E = 0, 0
    for i in 1:5
        # 5 tirs pour enlever les effets de l'aléatoire
        c_i, e_i = uls_glissant(4, T)
        print("($c_i, $e_i), ")
        C += c_i / 5
        E += e_i / 5
    end
    println()
    push!(C_vector, C)
    push!(E_vector, E)
end

println("C_vector = $C_vector")
println("E_vector = $E_vector")

# plot(T_vector, C_vector, title="Cout en fonction de T", label="cout")
# savefig("cout_fonction_T_Emax3.png")

# plot(T_vector, E_vector, title="Emission en fonction de T", label="Emission CO2")
# savefig("emission_fonction_T_Emax3.png")