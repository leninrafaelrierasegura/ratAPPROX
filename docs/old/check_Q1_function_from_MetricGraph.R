library(Matrix)
library(MetricGraph)

Qalpha1 <- function(theta, graph, BC = 1, build = TRUE) {
  
  kappa <- theta[2]
  tau <- theta[1]
  i_ <- j_ <- x_ <- rep(0, dim(graph$V)[1]*4) # it has length 4*nV
  count <- 0
  for(i in 1:graph$nE){ # loop over edges
    l_e <- graph$edge_lengths[i]
    c1 <- exp(-kappa*l_e)
    c2 <- c1^2
    one_m_c2 = 1-c2
    c_1 = 0.5 + c2/one_m_c2
    c_2 = -c1/one_m_c2
    
    if (graph$E[i, 1] != graph$E[i, 2]) { # This is for non-circular edges and codes I(overline(e) != underline(e))
      
      i_[count + 1] <- graph$E[i, 1]
      j_[count + 1] <- graph$E[i, 1]
      x_[count + 1] <- c_1
      
      i_[count + 2] <- graph$E[i, 2]
      j_[count + 2] <- graph$E[i, 2]
      x_[count + 2] <- c_1
      
      
      i_[count + 3] <- graph$E[i, 1]
      j_[count + 3] <- graph$E[i, 2]
      x_[count + 3] <- c_2
      
      i_[count + 4] <- graph$E[i, 2]
      j_[count + 4] <- graph$E[i, 1]
      x_[count + 4] <- c_2
      count <- count + 4
    }else{ # This is for circular edges and codes I(overline(e) = underline(e))
      i_[count + 1] <- graph$E[i, 1]
      j_[count + 1] <- graph$E[i, 1]
      x_[count + 1] <- tanh(0.5 * kappa * l_e)
      count <- count + 1
    }
  }
  if(BC == 1){
    #does this work for circle?
    i.table <- table(i_[1:count])
    index = as.integer(names(which(i.table < 3)))
    i_ <- c(i_[1:count], index)
    j_ <- c(j_[1:count], index)
    x_ <- c(x_[1:count], rep(0.5, length(index)))
    count <- count + length(index)
  }else if(BC==2){
    
    dV <- graph$get_vertices()$degree
    index <- 1:length(dV)
    i_ <- c(i_[1:count], index)
    j_ <- c(j_[1:count], index)
    x_ <- c(x_[1:count], -0.5*dV + 1)
    count <- count + length(index)
    
  }
  if(build){
    Q <- Matrix::sparseMatrix(i = i_[1:count],
                              j = j_[1:count],
                              x = (2 * kappa * tau^2) * x_[1:count], # This is the 2kappa*tau^2 factor
                              dims = c(graph$nV, graph$nV))
    
    
    return(Q)
  } else {
    return(list(i = i_[1:count],
                j = j_[1:count],
                x = (2 * kappa * tau^2) * x_[1:count],
                dims = c(graph$nV, graph$nV)))
  }
}

capture.output(
  knitr::purl(here::here("control_functionality.Rmd"), output = here::here("control_functionality.R")),
  file = here::here("purl_log.txt")
)
source(here::here("control_functionality.R"))

overkill_graph <- gets.graph.tadpole(h = overkill_h)
