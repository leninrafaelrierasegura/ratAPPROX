library(Matrix)
library(MetricGraph)


capture.output(
  knitr::purl(here::here("all_functions.Rmd"), output = here::here("all_functions.R")),
  file = here::here("old/purl_log.txt")
)
source(here::here("all_functions.R"))

graph <- gets.graph.tadpole(h = 0.5)
