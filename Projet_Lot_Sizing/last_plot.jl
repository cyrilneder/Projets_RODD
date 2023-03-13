using Plots

include("res.txt")

plot(T_vector, [C_vector3, C_vector6m], title="Cout en fonction de T", label=["4 modes" "6 modes"])
savefig("cout_fonction_T_6m.png")

plot(T_vector, [E_vector3, E_vector6m], title="Emission en fonction de T", label=["4 modes" "6 modes"])
savefig("emission_fonction_T_6m.png")