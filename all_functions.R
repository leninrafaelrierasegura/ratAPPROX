## ---------------------------------------------------------------------------------------------------------------------------------------------------
# remotes::install_github("davidbolin/rspde", ref = "devel")
# remotes::install_github("davidbolin/metricgraph", ref = "devel")
library(rSPDE)
library(MetricGraph)
library(grateful)

library(ggplot2)
library(reshape2)
library(plotly)


## ---------------------------------------------------------------------------------------------------------------------------------------------------
# Function to build a tadpole graph and create a mesh
gets.graph.tadpole <- function(){
  edge1 <- rbind(c(0,0),c(1,0))#[c(2,1),]
  theta <- seq(from=-pi,to=pi,length.out = 10000)
  edge2 <- cbind(1+1/pi+cos(theta)/pi,sin(theta)/pi)
  edges <- list(edge1, edge2)
  graph <- metric_graph$new(edges = edges, verbose = 0)
  graph$set_manual_edge_lengths(edge_lengths = c(1,2))
  #graph$build_mesh(h = h)
  return(graph)
}


## ---------------------------------------------------------------------------------------------------------------------------------------------------
# Eigenfunctions for the tadpole graph
tadpole.eig <- function(k,graph){
  x1 <- c(0,graph$get_edge_lengths()[1]*graph$mesh$PtE[graph$mesh$PtE[,1]==1,2]) 
  x2 <- c(0,graph$get_edge_lengths()[2]*graph$mesh$PtE[graph$mesh$PtE[,1]==2,2]) 
  
  if(k==0){ 
    f.e1 <- rep(1,length(x1)) 
    f.e2 <- rep(1,length(x2)) 
    f1 = c(f.e1[1],f.e2[1],f.e1[-1], f.e2[-1]) 
    f = list(phi=f1/sqrt(3)) 
    
  } else {
    f.e1 <- -2*sin(pi*k*1/2)*cos(pi*k*x1/2) 
    f.e2 <- sin(pi*k*x2/2)                  
    
    f1 = c(f.e1[1],f.e2[1],f.e1[-1], f.e2[-1]) 
    
    if((k %% 2)==1){ 
      f = list(phi=f1/sqrt(3)) 
    } else { 
      f.e1 <- (-1)^{k/2}*cos(pi*k*x1/2)
      f.e2 <- cos(pi*k*x2/2)
      f2 = c(f.e1[1],f.e2[1],f.e1[-1],f.e2[-1]) 
      f <- list(phi=f1,psi=f2/sqrt(3/2))
    }
  }
  
  return(f)
}


# Function to compute the true covariance matrix
gets_true_cov_mat <- function(graph, kappa, tau, alpha, n.overkill){
  Sigma.kl <- matrix(0,nrow = dim(graph$mesh$V)[1],ncol = dim(graph$mesh$V)[1])
  for(i in 0:n.overkill){
    phi <- tadpole.eig(i,graph)$phi
    Sigma.kl <- Sigma.kl + (1/(kappa^2 + (i*pi/2)^2)^(alpha))*phi%*%t(phi)
    if(i>0 && (i %% 2)==0){ 
      psi <- tadpole.eig(i,graph)$psi
      Sigma.kl <- Sigma.kl + (1/(kappa^2 + (i*pi/2)^2)^(alpha))*psi%*%t(psi)
    }
    
  }
  Sigma.kl <- Sigma.kl/tau^2
  return(Sigma.kl)
}


## ---------------------------------------------------------------------------------------------------------------------------------------------------
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
  print(i_)
  print(j_)
  }
  if(BC == 1){
    #does this work for circle?
    i.table <- table(i_[1:count])
    index = as.integer(names(which(i.table < 3)))
    i_ <- c(i_[1:count], index)
    j_ <- c(j_[1:count], index)
    x_ <- c(x_[1:count], rep(0.5, length(index))) # here is where we add the 0.5 for degree one vertices
    count <- count + length(index)
    # print(i_)
    # print(j_)
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


## ---------------------------------------------------------------------------------------------------------------------------------------------------
gives.indices <- function(graph, factor, constant){
  index.obs1 <- sapply(graph$PtV, 
                       function(i){
                         idx_temp <- i == graph$E[,1]
                         idx_temp <- which(idx_temp)
                         return(idx_temp[1])}
                       )
  index.obs1 <- (index.obs1 - 1) * factor + 1
  index.obs4 <- NULL
  na_obs1 <- is.na(index.obs1)
  if(any(na_obs1)){
    idx_na <- which(na_obs1)
    PtV_NA <- graph$PtV[idx_na]
    index.obs4 <- sapply(PtV_NA, 
                         function(i){
                           idx_temp <- i == graph$E[,2]
                           idx_temp <- which(idx_temp)
                           return(idx_temp[1])}
                         )
    index.obs1[na_obs1] <- (index.obs4 - 1 ) * factor + constant                                                                      
  }
  return(index.obs1)
}

conditioning <- function(graph, alpha = 1){
  i_  =  rep(0, 2 * graph$nE)
  j_  =  rep(0, 2 * graph$nE)
  x_  =  rep(0, 2 * graph$nE)

  count_constraint <- 0
  count <- 0
  for (v in 1:graph$nV) {
    edges_leaving_v  <- which(graph$E[, 1] %in% v) 
    edges_entering_v  <- which(graph$E[, 2] %in% v)
    n_leaving_edges <- length(edges_leaving_v)
    n_entering_edges <- length(edges_entering_v)
    n_e <- n_leaving_edges + n_entering_edges
    if (n_e > 1) { # the alternative is n_e = 1, which means v is a degree one vertex and so no conditioning is needed 
      if (n_entering_edges == 0) {
        edges <- cbind(edges_leaving_v, 1)
      } else if(n_leaving_edges == 0){
        edges <- cbind(edges_entering_v, 2)
      }else{
        edges <- rbind(cbind(edges_leaving_v, 1),
                       cbind(edges_entering_v, 2))
      }
      for (i in 2:n_e) {
        i_[count + 1:2] <- count_constraint + 1
        j_[count + 1:2] <- c(2 * (edges[i-1,1] - 1) + edges[i-1, 2],
                             2 * (edges[i,1]   - 1) + edges[i,   2])
        x_[count + 1:2] <- c(1,-1)
        count <- count + 2
        count_constraint <- count_constraint + 1
      }
    }
  }
  K <- Matrix::sparseMatrix(i = i_[1:count],
                            j = j_[1:count],
                            x = x_[1:count],
                            dims = c(count_constraint, 2*graph$nE))
                         
  CB <- MetricGraph:::c_basis2(K)
  return(CB)
}


## ---------------------------------------------------------------------------------------------------------------------------------------------------
# the one that translates from Vaibhav's
gets_cov_mat_rat_approx_alpha_1_to_2 <- function(graph, kappa, tau, alpha, m){

  # get rational approximation coefficients
  coeff <- rSPDE:::interp_rational_coefficients(
    order = m, 
    type_rational_approx = "chebfun", 
    type_interp = "spline", 
    alpha = alpha)
  
  r <- coeff$r
  p <- coeff$p
  k <- coeff$k
  
  # reparameterization
  nu <- alpha - 1/2
  sigma <- sqrt(gamma(nu) / (tau^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))
  c_alpha <- gamma(alpha)/gamma(alpha - 0.5)
  c_1 <- gamma(floor(alpha))/gamma(floor(alpha) - 0.5)
  
  # get edge lengths
  L_e <- graph$edge_lengths
  
  Qtilde_i <- list() 
  for(order in 0:m){
    if(order == 0){
      P <- order
      ALPHA <- floor(alpha)
      FACTOR <- (2*tau^2)/(k*kappa)
      r00_inverse <- solve(matern.p.joint(s = 0, t = 0, kappa = kappa, p = P, alpha = ALPHA))
      correction_term <- rbind(cbind(r00_inverse, matrix(0, floor(alpha), floor(alpha))),
                               cbind(matrix(0, floor(alpha), floor(alpha)), r00_inverse))
    } else {
      P <- p[order]
      ALPHA <- alpha
      FACTOR <- (2*c_alpha*sqrt(pi)*tau^2)/(r[order] * kappa)
      r00_inverse <- solve(matern.p.joint(s = 0, t = 0, kappa = kappa, p = P, alpha = ALPHA))
      correction_term <- rbind(cbind(r00_inverse, matrix(0, ceiling(alpha), ceiling(alpha))),
                               cbind(matrix(0, ceiling(alpha), ceiling(alpha)), r00_inverse))
    }
    Qtilde_i[[paste0("m=",order)]] <- list()
    for(e in 1:length(L_e)){
      Q_e <- matern.p.precision(loc = c(0, L_e[e]), #if (order == 0) c(0, L_e[e]) else c(L_e[e], 0), 
                                kappa = kappa, 
                                p = P,
                                equally_spaced = TRUE, 
                                alpha = ALPHA)$Q
      Qtilde_i[[paste0("m=",order)]][[e]] <- (Q_e - 0.5 * correction_term)*FACTOR*kappa^(2*alpha)
    }
    Qtilde_i[[paste0("m=",order)]] <- bdiag(Qtilde_i[[paste0("m=",order)]])
  }
  
  Qtilde_0 <- Qtilde_i[[paste0("m=",0)]] # extract Qtilde_0
  Qtilde_i <- Qtilde_i[-1] # remove Qtilde_0
  

  #####################################
  ## CASE m = 0
  #####################################
  COND_0 <- conditioning(graph = graph, alpha = 1)
  index.obs_0 <- gives.indices(graph = graph, factor = 2, constant = 2)
  nc_0 <- 1:length(COND_0$S) # number of constraints
  T_0 <- COND_0$T # change of basis matrix
  W_0 <- Diagonal(2*floor(alpha)*graph$nE)[,-nc_0] # matrix to remove constraints
  Qtilde_0_star_UU <- t(W_0) %*% t(T_0) %*% Qtilde_0 %*% (T_0) %*% W_0 
  A0 <- T_0[index.obs_0, -nc_0] # observation matrix after conditioning
  
  # Qtilde_0_star_UU <- MetricGraph:::Qalpha1(theta = c(tau, kappa), graph = graph, BC = 1, build = TRUE)*c_1*kappa^(2*alpha)/(2 * k * c_alpha * kappa * sigma^2 * tau^2)#(2*tau^2*kappa^(2*alpha))/(k * kappa)
  # A0 <- graph$.__enclos_env__$private$A()
  #####################################
  ## CASE m > 0
  #####################################
  graph$buildC(alpha = 2, edge_constraint = TRUE)
  COND_i <- graph$CoB
  index.obs_i <- gives.indices(graph = graph, factor = 4, constant = 3)
  n_const <- length(COND_i$S)
  ind.const <- c(1:n_const)
  Tc <- COND_i$T[-ind.const, ]
  Qtilde_i_star_UU <- lapply(Qtilde_i, function(Q) Tc %*% Q %*% t(Tc)) 
  Ai <- t(Tc)[index.obs_i, ] # observation matrix after conditioning
  
  # graph$buildC(alpha = 2)
  # COND_i <- graph$CoB
  # index.obs_i <- gives.indices(graph = graph, factor = 4, constant = 3)
  # nc_i <- 1:length(c(1,COND_i$S)) # number of constraints
  # T_i <- COND_i$T # change of basis matrix
  # T_i <- t(T_i)[, c(ncol(T_i), 1:(ncol(T_i) - 1))] # column reordering
  # W_i <- Diagonal(2*ceiling(alpha)*graph$nE)[,-nc_i] # matrix to remove constraints
  # Qtilde_i_star_UU <- lapply(Qtilde_i, function(Q) t(W_i) %*% t(T_i) %*% Q %*% T_i %*% W_i) 
  # Ai <- T_i[index.obs_i, -nc_i] # observation matrix after conditioning
  
  #####################################
  ## Build matrix A and Q_UU
  #####################################
  A <- cbind(A0, do.call(cbind, rep(list(Ai), m)))
  Q_UU <- bdiag(Qtilde_0_star_UU, do.call(bdiag, Qtilde_i_star_UU))
  # Return Sigma
  Sigma <- A %*% solve(Q_UU, t(A)) 
  return(Sigma)
}

