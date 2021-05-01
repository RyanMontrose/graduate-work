using Plots
using GraphRecipes
pyplot()

# 𝒜 = Matrix(adjacency_matrix(random_regular_graph(N_homes, connection_density)))
𝒜 = Matrix(adjacency_matrix(barabasi_albert(N_homes, 3)))
# 𝒜 = Matrix(adjacency_matrix(star_graph(N_homes)))

graphplot(𝒜, nodesize=0.3, nodeshape=:circle, names=1:N_homes, curvature_scalar=0.008, linecolor=:black, fontsize=14, markercolor=:lightblue)