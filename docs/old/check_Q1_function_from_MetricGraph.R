library(Matrix)
library(MetricGraph)


# Function to build a tadpole graph and create a mesh
gets.graph.tadpole <- function(flip_edge = FALSE){
  if(flip_edge) {
    edge1 <- rbind(c(0,0),c(1,0))[c(2,1),]
  } else {
    edge1 <- rbind(c(0,0),c(1,0))}
  theta <- seq(from=-pi,to=pi,length.out = 10000)
  edge2 <- cbind(1+1/pi+cos(theta)/pi,sin(theta)/pi)
  edges <- list(edge1, edge2)
  graph <- metric_graph$new(edges = edges, verbose = 0)
  graph$set_manual_edge_lengths(edge_lengths = c(1,2))
  #graph$build_mesh(h = h)
  return(graph)
}

graph <- gets.graph.tadpole(flip_edge = TRUE)
graph_copy <- graph$clone()

graph$build_mesh(h = 0.02)

mesh_loc <- graph$get_mesh_locations() %>% 
  as.data.frame() %>% 
  mutate(y = 1) %>% 
  rename(edge_number = V1, distance_on_edge = V2)

graph_copy$add_observations(mesh_loc,
                            edge_number = "edge_number",
                            distance_on_edge = "distance_on_edge",
                            data_coords = "PtE",
                            normalized = TRUE, 
                            clear_obs = TRUE)
graph_copy$observation_to_vertex()


Q <- MetricGraph:::Qalpha1(
  theta = c(1, 5), 
  graph = graph_copy, 
  BC = 3000, 
  build = TRUE)

Sigma <- solve(Q, Diagonal(n = nrow(Q)))
graph$plot_function(X = Sigma[,1], type = "plotly", line_color = "red", interpolate_plot = FALSE, name = "True", showlegend = TRUE)


# graph_copy$plot()
# graph_copy$plot_connections()
# graph_copy$nE
# graph$buildDirectionalConstraints(alpha = 1)
# CoB <- graph$buildC(alpha = 2)
# 
# aux <- graph$CoB
# aux
# CoB
