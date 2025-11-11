library(Matrix)
library(MetricGraph)


capture.output(
  knitr::purl(here::here("all_functions.Rmd"), output = here::here("all_functions.R")),
  file = here::here("old/purl_log.txt")
)
source(here::here("all_functions.R"))

graph <- gets.graph.tadpole()
graph_copy <- graph$clone()

graph$build_mesh(h = 0.2)

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


Q <- Qalpha1(theta = c(1,1), graph = graph_copy, BC = 1, build = TRUE)
Q

graph_copy$plot()
graph_copy$plot_connections()
graph_copy$nE

