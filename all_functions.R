## --------------------------------------------------------------------------------------------------
# remotes::install_github("davidbolin/rspde", ref = "devel")
# remotes::install_github("davidbolin/metricgraph", ref = "devel")
library(rSPDE)
library(MetricGraph)
library(grateful)

library(ggplot2)
library(reshape2)
library(plotly)


## --------------------------------------------------------------------------------------------------
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


## --------------------------------------------------------------------------------------------------
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


## --------------------------------------------------------------------------------------------------
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


## --------------------------------------------------------------------------------------------------
# Typically, factor = 4, constant = 3
gives.indices <- function(graph, factor, constant){
  # Here, after doing graph$observation_to_vertex() 
  # graph$PtV is just from 1 to graph$nV in some order
  if(is.null(graph$PtV)){
    stop("graph$PtV is NULL, please run graph$observation_to_vertex() first")
  }
  verticesInSomeOrderWeAreRestricted <- graph$PtV
  index.obs1 <- sapply(verticesInSomeOrderWeAreRestricted, 
                       function(i) which(i == graph$E[,1])[1]) 
  # The above just identifies the first index (that is why [1]) in graph$E[,1]
  # where vertex v_i with index graph$PtV[i] appears, it may have NA's
  # as not all vertices are the start of some edge
  index.obs1 <- (index.obs1 - 1) * factor + 1
  #index.obs4 <- NULL
  na_obs1 <- is.na(index.obs1)
  if(any(na_obs1)){
    idx_na <- which(na_obs1)
    PtV_NA <- verticesInSomeOrderWeAreRestricted[idx_na]
    index.obs4 <- sapply(PtV_NA, 
                         function(i) which(i == graph$E[,2])[1])
    index.obs1[na_obs1] <- (index.obs4 - 1 ) * factor + constant                                                                      
  }
  return(index.obs1)
}

indicesForOrderedVertices <- function(graph, alpha){
  if(is.null(graph$PtV)){
    stop("graph$PtV is NULL, please run graph$observation_to_vertex() first")
  }
  # If alpha is not integer, then stop
  if(alpha %% 1 != 0){
    stop("alpha is not an integer, only works for integer alpha")
  }
  initialIndices <- gives.indices(graph, factor = 4, constant = 3)
  adjustedIndices <- (initialIndices - 1) * alpha/2 + 1
  return(adjustedIndices)
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


## --------------------------------------------------------------------------------------------------
computesListOfMatricesQTildeUnconstraint <- function(p,
                                                     kappa, 
                                                     alpha, 
                                                     edgeLengths){
  ca <- ceiling(alpha)
  # initialize Qtilde_i, a list containing block diagonal matrices with blocks Qtilde_{i,e} for each i
  Qtilde_i <- list() 
  for(i in seq_along(p)){
    
    # compute r_(0,0)
    r00 <- matern.p.joint(
      s = 0, 
      t = 0, 
      kappa = kappa, 
      p = p[i], 
      alpha = alpha)
    
    # compute r_(0,0)^(-1)
    r00_inverse <- solve(r00, Diagonal(ca))
    
    # define zero block 
    zero_block <- matrix(0, ca, ca)
    
    # build correction term
    correction_term <- rbind(
      cbind(r00_inverse, zero_block),
      cbind(zero_block, r00_inverse))
    
    # initialize Qtilde_i[[i]], a list containing Qtilde_{i,e} for each edge e
    Qtilde_i[[i]] <- list()
    for(e in seq_along(edgeLengths)){
      
      # compute Q_{i,e}
      Q_e <- matern.p.precision(
        loc = c(0, edgeLengths[e]),
        kappa = kappa, 
        p = p[i],
        equally_spaced = FALSE, 
        alpha = alpha)$Q
      
      # store Qtilde_{i,e}
      Qtilde_i[[i]][[e]] <- Q_e - 0.5 * correction_term
    }
    # build block diagonal matrix Qtilde_i[[i]]
    Qtilde_i[[i]] <- bdiag(Qtilde_i[[i]])
  }
  return(Qtilde_i)
}


## --------------------------------------------------------------------------------------------------
# This is the correct version, it is corrected the constants
gets_cov_mat_rat_approx_alpha_1_to_2 <- function(graph, kappa, tau, alpha, m, build_cov){
  
  if(alpha == 2){
    Q_unconstraint <- MetricGraph:::Qalpha2(theta = c(tau, kappa),
                                         graph = graph,
                                         BC = 3000,
                                         build = TRUE)
    graph$buildC(alpha = alpha, edge_constraint = TRUE) # should always be TRUE
    COND <- graph$CoB
    Tc <- COND$T[-c(1:length(COND$S)), ]

    Q_U <-  Tc %*% Q_unconstraint %*% t(Tc)

    index.obs_i <- gives.indices(graph = graph, factor = 4, constant = 3)
    A <- t(Tc)[index.obs_i, ] 
    if(build_cov){
    Sigma <- A %*% solve(Q_U, t(A))
    return(Sigma)
    }
    return(list(A = A, Q = Q_U))
  }

  # get rational approximation coefficients
  coeff <- rSPDE:::interp_rational_coefficients(
    order = m, 
    type_rational_approx = "chebfun", 
    type_interp = "spline", 
    alpha = alpha)
  
  r <- coeff$r
  p <- coeff$p
  k <- coeff$k
  
  # compute parameters
  fa <- floor(alpha)
  ca <- ceiling(alpha)
  
  nu <- alpha - 1/2
  sigma <- sqrt(gamma(nu) / (tau^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))
  c_alpha <- gamma(alpha)/gamma(alpha - 0.5)
  c_1 <- gamma(fa)/gamma(fa - 0.5)
  
  # get edge lengths
  edgeLengths <- graph$edge_lengths

  Qtilde_i <- computesListOfMatricesQTildeUnconstraint(p, kappa,  alpha,  edgeLengths)
  # --------------------------------------------------
  # CASE i = 0
  # --------------------------------------------------
  
  # When I use MetricGraph:::Qalpha1, I am assuming that tau and sigma are related by tau^2 = gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)) where nu = 1/2
  
  inv_factor_0 <- 1/(k*c_alpha/c_1)
  NU <- fa - 0.5
  TAU <- sqrt(gamma(NU) / (sigma^2 * kappa^(2*NU) * (4*pi)^(1/2) * gamma(NU + 1/2)))
    
  Qtilde_0_star_UU <- MetricGraph:::Qalpha1(
    theta = c(TAU, kappa), 
    graph = graph, 
    BC = 3000, 
    build = TRUE) * inv_factor_0
  
  A_0 <- graph$.__enclos_env__$private$A()

  # --------------------------------------------------
  # CASE i = 1,...,m
  # --------------------------------------------------
  
  # build conditioning matrix
  graph$buildC(alpha = 2, edge_constraint = TRUE) # should always be TRUE
  COND_i <- graph$CoB
  Tc <- COND_i$T[-c(1:length(COND_i$S)), ]
  
  inv_factor_i <- 1/(r*sigma^2)

  Qtilde_i_star_UU <- purrr::map2(
    Qtilde_i, 
    inv_factor_i, 
    function(Q, x) Tc %*% Q %*% t(Tc) * x)

  index.obs_i <- gives.indices(graph = graph, factor = 4, constant = 3)
  
  # Compare the above to 
  # A <- buildMatrixAWhichMapsUToUv(graph, 2)
  # index.obs_i <- which(A == 1, arr.ind = TRUE)[, "col"]
  
  A_i <- t(Tc)[index.obs_i, ] 
  
  # Build matrix A and Q_UU
  A <- cbind(A_0, do.call(cbind, rep(list(A_i), m)))
  Q_UU <- bdiag(Qtilde_0_star_UU, bdiag(Qtilde_i_star_UU))
  if(build_cov){
  # Return Sigma
  Sigma <- A %*% solve(Q_UU, t(A)) 
  return(Sigma)
  }
  return(list(A = A, Q = Q_UU))
}


## --------------------------------------------------------------------------------------------------
gets_cov_mat_rat_approx_alpha_0_to_1 <- function(graph, kappa, tau, alpha, m, build_cov){
  
  if(alpha == 1){
    Q_U <- MetricGraph:::Qalpha1(
      theta = c(tau, kappa),
      graph = graph,
      BC = 3000,
      build = TRUE)
    A <- graph$.__enclos_env__$private$A()
    
    if(build_cov){
    Sigma <- A %*% solve(Q_U, t(A))
    #I <- Matrix::Diagonal(graph$nV)  
    #Sigma <- solve(Q_U, I)
    return(Sigma)
    }
    return(list(A = A, Q = Q_U))
  }

  coeff <- rSPDE:::interp_rational_coefficients(
    order = m, 
    type_rational_approx = "chebfun", 
    type_interp = "spline", 
    alpha = alpha)
  
  r <- coeff$r
  p <- coeff$p
  k <- coeff$k
  
  nu <- alpha - 1/2
  sigma <- sqrt(gamma(nu) / (tau^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))
  
  fa <- floor(alpha)
  ca <- ceiling(alpha)
  
  c_alpha <- gamma(alpha)/gamma(alpha - 0.5)
  
  # --------------------------------------------------
  # CASE i = 0
  # --------------------------------------------------
  
  inv_factor_0 <- 1/(k*c_alpha*sqrt(4*pi)*sigma^2/kappa)
  I <- Matrix::Diagonal(graph$nV)     
  Qtilde_0_star_UU <- I * inv_factor_0
  
  # --------------------------------------------------
  # CASE i = 1,...,m
  # --------------------------------------------------
  
  inv_factor_i <- 1/(r*c_alpha*sqrt(pi)/sqrt(1 - p))
    
  NU <- ca - 0.5
  
  Qtilde_i_star_UU <- list()
  for(i in 1:m){
    
    KAPPA <- kappa*sqrt(1 - p[i])
    TAU <- sqrt(gamma(NU) / (sigma^2 * KAPPA^(2*NU) * (4*pi)^(1/2) * gamma(NU + 1/2)))
    
    Qtilde_i_star_UU[[i]] <- MetricGraph:::Qalpha1(
      theta = c(TAU, KAPPA), 
      graph = graph, 
      BC = 3000, 
      build = TRUE) * inv_factor_i[i]
  
  }
  
  A_0 <- graph$.__enclos_env__$private$A()
  # Build matrix A and Q_UU
  A <- do.call(cbind, rep(list(A_0), m+1))
  Q_UU <- bdiag(Qtilde_0_star_UU, bdiag(Qtilde_i_star_UU))
  # Return Sigma
  
  if(build_cov){
  Sigma <- A %*% solve(Q_UU, t(A))

  # Q_UU <- c(list(Qtilde_0_star_UU), Qtilde_i_star_UU)
  # 
  # Sigma <- Reduce(`+`, lapply(Q_UU, function(Qi) solve(Qi, I)))
  return(Sigma)
  }
  return(list(A = A, Q = Q_UU))
}


## --------------------------------------------------------------------------------------------------
rat_covariance <- function(graph, 
                           kappa, 
                           tau, 
                           alpha,
                           m, 
                           build_cov = TRUE){
  if(alpha <= 0.5){
    stop("alpha = ", alpha, ", alpha should be larger than 0.5")
  }
  else if(alpha > 0.5 && alpha <= 1){
    return(gets_cov_mat_rat_approx_alpha_0_to_1(graph = graph,
                                                kappa = kappa,
                                                tau = tau,
                                                alpha = alpha,
                                                m = m,
                                                build_cov = build_cov))
  }
  else if(alpha > 1 && alpha <= 2){
    return(gets_cov_mat_rat_approx_alpha_1_to_2(graph = graph,
                                                kappa = kappa,
                                                tau = tau,
                                                alpha = alpha,
                                                m = m, 
                                                build_cov = build_cov))
  }
}


## --------------------------------------------------------------------------------------------------

lazy_likelihood_alpha_rat <- function(graph,
                                            kappa,
                                            tau,
                                            sigma_e,
                                            alpha,
                                            m,
                                            Y){
  n <- length(Y)
  Q_and_A <- rat_covariance(
      graph = graph, 
      kappa = kappa, 
      tau = tau, 
      alpha = alpha, 
      m = m, 
      build_cov = FALSE)
  
  AT_Ut <- Q_and_A$A
  
  Q <- Q_and_A$Q
  L <- chol(Q)
  
  Qp <- Q + 1/sigma_e^2 * t(AT_Ut) %*% AT_Ut
  Lp <- chol(Qp)
  
  mu <- 1/sigma_e^2 * solve(Qp, t(AT_Ut) %*% Y)
    
    loglik <- sum(log(diag(L))) -
      n*log(sigma_e) -
      sum(log(diag(Lp))) -
      0.5/sigma_e^2 * t(Y) %*% Y +
      0.5 * t(mu) %*% Qp %*% mu -
          0.5 * n * log(2*pi)
    
    return(loglik)
}


## --------------------------------------------------------------------------------------------------
rat_loglikelihood <- function(graph,
                              theta,
                              alpha,
                              m,
                              BC = 0){
  
  kappa = exp(theta[3])
  sigma_e <- exp(theta[1])
  reciprocal_tau <- exp(theta[2])
  tau <- 1/reciprocal_tau
  
  graph$observation_to_vertex()
  data <- graph$get_data()
  Y <- data$y
  
  if(alpha <= 0.5 || alpha > 2){
    stop("alpha = ", alpha, ", alpha should be in (0.5,2]")
  }
  else if(alpha > 0.5 && alpha <= 2){
    return(lazy_likelihood_alpha_rat(graph = graph, 
                                     kappa = kappa, 
                                     tau = tau, 
                                     sigma_e = sigma_e, 
                                     alpha = alpha, 
                                     m = m, 
                                     Y = Y))
  }
}


## --------------------------------------------------------------------------------------------------
FEM_loglikelihood <- function(object, y, X_cov, repl, A_list, sigma_e, beta_cov) {
  m <- object$m

  Q <- object$Q
  R <- Matrix::Cholesky(Q)

  prior.ld <- c(determinant(R, logarithm = TRUE, sqrt = TRUE)$modulus)

  repl_val <- unique(repl)

  l <- 0

  for (i in repl_val) {
    ind_tmp <- (repl %in% i)
    y_tmp <- y[ind_tmp]

    if (ncol(X_cov) == 0) {
      X_cov_tmp <- 0
    } else {
      X_cov_tmp <- X_cov[ind_tmp, , drop = FALSE]
    }

    na_obs <- is.na(y_tmp)

    y_ <- y_tmp[!na_obs]


    n.o <- length(y_)
    A_tmp <- A_list[[as.character(i)]]
    Q.p <- Q + t(A_tmp) %*% A_tmp / sigma_e^2

    R.p <- Matrix::Cholesky(Q.p)

    posterior.ld <- c(determinant(R.p, logarithm = TRUE, sqrt = TRUE)$modulus)

    l <- l + prior.ld - posterior.ld - n.o * log(sigma_e)

    v <- y_

    if (ncol(X_cov) > 0) {
      X_cov_tmp <- X_cov_tmp[!na_obs, , drop = FALSE]
      v <- v - X_cov_tmp %*% beta_cov
    }

    mu.p <- solve(R.p, as.vector(t(A_tmp) %*% v / sigma_e^2), system = "A")

    v <- v - A_tmp %*% mu.p

    l <- l - 0.5 * (t(mu.p) %*% Q %*% mu.p + t(v) %*% v / sigma_e^2) -
      0.5 * n.o * log(2 * pi)
  }

  return(as.double(l))
}


## --------------------------------------------------------------------------------------------------
gets_De_from_Uv <- function(graph, alpha){
  E  <- graph$E
  nV <- graph$nV
  nE <- graph$nE
  
  De_list <- vector("list", nE)
  
  for (e in 1:nE) {
    e_ini_ter <- E[e,]
    
    De_list[[e]] <- Matrix::sparseMatrix(
      i = c(seq_len(alpha), seq_len(alpha) + alpha),
      j = c(rep(e_ini_ter[1], alpha), rep(e_ini_ter[2], alpha)),
      x = 1,
      dims = c(2 * alpha, nV)
    )
  }
  
  return(De_list)
}


## --------------------------------------------------------------------------------------------------
gets_De_from_U <- function(graph, alpha){
  nE <- graph$nE 
  
  De_list <- vector("list", nE)
  
  for (e in seq_len(nE)) {
    rows <- seq_len(2 * alpha)
    cols <- 2 * alpha * (e - 1) + rows
    
    De_list[[e]] <- Matrix::sparseMatrix(
      i = rows,
      j = cols,
      x = 1,
      dims = c(2 * alpha, 2 * alpha * nE)
    )
  } 
  
  return(De_list)
}


## --------------------------------------------------------------------------------------------------
buildKirchooffConditioningMatrixCaseAlphaEqualOne <- function(graph) {
  edgeMatrix <- graph$E
  degrees <- graph$get_degrees()
  numberOfEdges <- nrow(edgeMatrix)
  
  numberOfConstraints <- sum(degrees[degrees > 1] - 1)
  
  edgeMatrixFlatten <- c(t(edgeMatrix))
  
  numberOfNonZero <- 2 * numberOfConstraints
  i_idx <- integer(numberOfNonZero)
  j_idx <- integer(numberOfNonZero)
  x_val <- numeric(numberOfNonZero)
  
  indexCounter <- 0
  rowCounter <- 0
  
  for (vertex in seq_along(degrees)) {
    degreeOfVertex <- degrees[vertex]
    if (degreeOfVertex < 2) next
    
    indicesOfVertex <- which(edgeMatrixFlatten == vertex)
    
    colIndicesForOne <- indicesOfVertex[-degreeOfVertex]
    colIndicesForMinusOne <- indicesOfVertex[-1]
    
    howManyContinuityConditions <- degreeOfVertex - 1
    
    rowIndicesForBothOneAndMinusOne <- seq_len(howManyContinuityConditions) + rowCounter
    
    # fill +1 entries
    idx_range <- (indexCounter + 1):(indexCounter + howManyContinuityConditions)
    i_idx[idx_range] <- rowIndicesForBothOneAndMinusOne
    j_idx[idx_range] <- colIndicesForOne
    x_val[idx_range] <- 1
    
    # Note that the amount of indices to fill is 2 * howManyContinuityConditions
    # Half for +1 entries, that is, from 1 to howManyContinuityConditions
    # The remaining half for -1 entries, that is, from howManyContinuityConditions + 1 to 2 * howManyContinuityConditions
    
    # fill -1 entries
    idx_range <- (indexCounter + howManyContinuityConditions + 1):
                 (indexCounter + 2 * howManyContinuityConditions)
    i_idx[idx_range] <- rowIndicesForBothOneAndMinusOne
    j_idx[idx_range] <- colIndicesForMinusOne
    x_val[idx_range] <- -1
    
    indexCounter <- indexCounter + 2 * howManyContinuityConditions
    rowCounter <- rowCounter + howManyContinuityConditions
  }
  
  K <- Matrix::sparseMatrix(
    i = i_idx,
    j = j_idx,
    x = x_val,
    dims = c(numberOfConstraints, 2 * numberOfEdges)
  )
  
  return(K)
}


## --------------------------------------------------------------------------------------------------
buildKirchooffConditioningMatrixCaseAlphaEqualThree <- function(graph) {
  alpha <- 2
  n <- 2*alpha*graph$nE
  
  if(is.null(graph$C)) graph$buildC(alpha = alpha)
  K2 <- graph$C
  # indicesOfProcessValues <- seq(1, n - 1, by = alpha)
  # aux <- K2[, indicesOfProcessValues]
  # K1 <- aux[rowSums(aux != 0) > 0, , drop = FALSE]
  
  # This may fail, if so, comment the line below and uncomment the above
  # 
  K1 <- buildKirchooffConditioningMatrixCaseAlphaEqualOne(graph) 
  
  n <- ncol(K1)
  K3 <- Matrix::Matrix(0, nrow(K1), 3*n, sparse = TRUE)
  K3[, seq(3, 3*n, by = 3)] <- K1
  
  n <- ncol(K2)
  K2WithColumnsExpanded <- Matrix::Matrix(0, nrow(K2), n + n/2, sparse = TRUE)
  K2WithColumnsExpanded[, as.vector(rbind(seq(1, n + n/2, by = 3),
                        seq(2, n + n/2, by = 3)))] <- K2

  return(rbind(K2WithColumnsExpanded, K3))
}


## --------------------------------------------------------------------------------------------------
buildMatrixAWhichMapsUToUv <- function(graph, alpha){
  edgeMatrix <- graph$E
  edgeMatrixFlattened <- c(t(edgeMatrix))
  firstVertexOccurrence <- !duplicated(edgeMatrixFlattened)
  verticesOrder <- edgeMatrixFlattened[firstVertexOccurrence]
  positionsOfOriginalVertices <- seq(1, alpha*2*graph$nE - 1, by = alpha)
  colIndicesOfNonZeroEntry <- positionsOfOriginalVertices[firstVertexOccurrence]
  A <- Matrix::sparseMatrix(
    i = seq_len(graph$nV),
    j = colIndicesOfNonZeroEntry,
    x = 1,
    dims = c(graph$nV, alpha*2*graph$nE)
  )
  return(list(A = A, 
              verticesOrder = verticesOrder, 
              colIndicesOfNonZeroEntry = colIndicesOfNonZeroEntry))
}



## --------------------------------------------------------------------------------------------------
getsSmallCovarianceMatrices <- function(D_matrix,
                                        kappa,
                                        tau,
                                        sigma_e,
                                        y_e){
  
  # Pre-compute JointCovarianceMatrix of c(0,ell,\mathbf{t}_e) matrix
  # That is, JointCovarianceMatrix = [Sigma_{X_tX_t} & Sigma_{X_tY} \\ Sigma_{YX_t} & Sigma_{YY}]
  #  where Y = [u(0), u(ell)] and X_t = [u(t_1),...,u(t_{n_e})]

  JointCovarianceMatrix <- MetricGraph:::r_1(D_matrix, kappa = kappa, tau = tau)

  #covariance update see Art p.17
  E.ind <- c(1:2)
  Obs.ind <- -E.ind
  
  Sigma_YY <- JointCovarianceMatrix[E.ind, E.ind, drop = FALSE]
  Sigma_YXt_e <- JointCovarianceMatrix[E.ind, Obs.ind, drop = FALSE]
  Sigma_Xt_eXt_e <- JointCovarianceMatrix[Obs.ind, Obs.ind, drop = FALSE]
  Sigma_Xt_eY <- JointCovarianceMatrix[Obs.ind, E.ind, drop = FALSE]
  
  Se_te_transposed <- solve(Sigma_YY, Sigma_YXt_e)
  
  Sigma_e <- Sigma_Xt_eXt_e - Sigma_Xt_eY  %*% Se_te_transposed
  diag(Sigma_e) <- diag(Sigma_e) + sigma_e^2
  
  chol_Sigma_e <- base::chol(Sigma_e)

  Sigma_e_inverse_Se_te <- backsolve(chol_Sigma_e, forwardsolve(t(chol_Sigma_e), t(Se_te_transposed)))

  Sigma_e_inverse_y_e <- backsolve(chol_Sigma_e, forwardsolve(t(chol_Sigma_e), y_e))
  y_e_transposed_Sigma_e_inverse_y_e <- sum(y_e * Sigma_e_inverse_y_e)
  
  
  Se_te_transposed_Sigma_e_inverse_Se_te <- Se_te_transposed %*% Sigma_e_inverse_Se_te
  Se_te_transposed_Sigma_e_inverse_ye <- t(Sigma_e_inverse_Se_te) %*% y_e
  
  half_log_det_Sigma_e <- sum(log(diag(chol_Sigma_e)))
  
  return(list(Se_te_transposed_Sigma_e_inverse_Se_te = Se_te_transposed_Sigma_e_inverse_Se_te,
              Se_te_transposed_Sigma_e_inverse_ye = Se_te_transposed_Sigma_e_inverse_ye,
              y_e_transposed_Sigma_e_inverse_y_e = y_e_transposed_Sigma_e_inverse_y_e,
              half_log_det_Sigma_e = half_log_det_Sigma_e))
}


## --------------------------------------------------------------------------------------------------
loglikelihoodForAlphaEqualOnePrecompute <- function(theta, 
                                                    graph, 
                                                    precomputeddata,
                                                    data_name = NULL, 
                                                    manual_y = NULL,
                                                    X_cov = NULL, 
                                                    repl, 
                                                    BC, 
                                                    parameterization = "spde") {
  sigma_e <- exp(theta[1])

  repl_vec <- graph$.__enclos_env__$private$data[[".group"]]

  if(is.null(repl)){
    repl <- unique(repl_vec)
  }

  if(parameterization == "matern"){
    kappa = sqrt(8 * 0.5) / exp(theta[3])
  } else{
    kappa = exp(theta[3])
  }

  reciprocal_tau <- exp(theta[2])

  Q.list <- spde_precision(kappa = kappa, tau = 1/reciprocal_tau, alpha = 1,
                           graph = graph, build = FALSE, BC = BC)

  Qv <- Matrix::sparseMatrix(i = Q.list$i, j = Q.list$j, x = Q.list$x, dims = Q.list$dims)
  chol_Qv <- Matrix::Cholesky(Qv, LDL = FALSE, perm = TRUE)
  half_log_det_Qv <- Matrix::determinant(chol_Qv, sqrt=TRUE)$modulus[1]

  #PtE <- graph$get_PtE()
  obs.edges <- precomputeddata$obs.edges

  i_ <- j_ <- x_ <- rep(0, 4 * length(obs.edges))

  if(is.null(repl)){
    u_repl <- unique(graph$.__enclos_env__$private$data[[".group"]])
  } else{
    u_repl <- unique(repl)
  }

  loglik <- 0

  half_log_det_Qp <- NULL
  numberOfObservations <- 0

  # Cache some values used in the loop
  nV <- nrow(graph$V)

  for(j in seq_along(u_repl)){
      curr_repl <- u_repl[j]
      repl_name <- paste0("repl_", curr_repl)

      loglik <- loglik + half_log_det_Qv # we add 0.5 * log|Qv| n_repl times
      count <- 0
      Qp_mu <- numeric(nV)

    for (i in seq_along(obs.edges)) {
      # Use pre-computed replicate indices
      e <- obs.edges[i]
      edge_name <- paste0("edge_", e)

      # Get data for this edge and replicate
      y_e <- precomputeddata$y[[repl_name]][[edge_name]]

      # Skip if no observations
      if(is.null(y_e) || length(y_e) == 0){
        next
      }

      numberOfObservations <- numberOfObservations + length(y_e) # update the counter of number of observations

      if(!is.null(X_cov)){
        n_cov <- ncol(X_cov)
        if(n_cov == 0){
          X_cov_repl <- 0
        } else{
          y_e <- y_e - as.vector(precomputeddata$x[[repl_name]][[edge_name]] %*% theta[4:(3+n_cov)])
        }
      }
      D_matrix <- precomputeddata$D_matrix[[repl_name]][[edge_name]]

      AUX <- getsSmallCovarianceMatrices(D_matrix = D_matrix,
                                         kappa = kappa, 
                                         tau = 1/reciprocal_tau, 
                                         sigma_e = sigma_e,
                                         y_e = y_e)
      # Se_te_transposed_Sigma_e_inverse_Se_te is of size 2\alpha x 2\alpha
      Se_te_transposed_Sigma_e_inverse_Se_te <- AUX$Se_te_transposed_Sigma_e_inverse_Se_te
      # Se_te_transposed_Sigma_e_inverse_ye is of size 2\alpha \times 1, for alpha = 1, is 2x1, one for each vertex connecting e
      Se_te_transposed_Sigma_e_inverse_ye <- AUX$Se_te_transposed_Sigma_e_inverse_ye
      y_e_transposed_Sigma_e_inverse_y_e <- AUX$y_e_transposed_Sigma_e_inverse_y_e
      half_log_det_Sigma_e <- AUX$half_log_det_Sigma_e

      loglik <- loglik - 0.5 * y_e_transposed_Sigma_e_inverse_y_e - half_log_det_Sigma_e
      
      # The above adds - 0.5 * y_e^\top \Sigma_e^{-1}y_e and - 0.5 * log|Sigma_e|
      
      # Note that Qp_mu is of size |V|\times 1, that is why for each edge, we directly update Qp_mu
      # using E, which contains the indices of the vertices for which the edge contributes to
      E <- graph$E[e, ]
      if (E[1] == E[2]) { # for loops
        # Pre-compute matrix product
        Qp_mu[E[1]] <- Qp_mu[E[1]] + sum(as.vector(Se_te_transposed_Sigma_e_inverse_ye))
        # This: Qp_mu[E[1]] <- Qp_mu[E[1]] is the update part, we only 
        i_[count + 1] <- E[1]
        j_[count + 1] <- E[1]
        x_[count + 1] <- sum(Se_te_transposed_Sigma_e_inverse_Se_te)
        count <- count + 1
      } else { # for non loops
        Qp_mu[E] <- Qp_mu[E] + as.vector(Se_te_transposed_Sigma_e_inverse_ye) 
        
        # This: Qp_mu[E] <- Qp_mu[E] is just the update part
        
        idx <- count + 1:4
        i_[idx] <- c(E[1], E[1], E[2], E[2])
        j_[idx] <- c(E[1], E[2], E[1], E[2])
        x_[idx] <- c(Se_te_transposed_Sigma_e_inverse_Se_te[1, 1], Se_te_transposed_Sigma_e_inverse_Se_te[1, 2],
                     Se_te_transposed_Sigma_e_inverse_Se_te[1, 2], Se_te_transposed_Sigma_e_inverse_Se_te[2, 2])
        count <- count + 4
      }

    }

    if(is.null(half_log_det_Qp)){ # This is to compute ONE time only, it assumes all replicates share the same locations
      i_ <- c(Q.list$i, i_[1:count])
      j_ <- c(Q.list$j, j_[1:count])
      x_ <- c(Q.list$x, x_[1:count])

      # Qp = Qv + sum D_e^\top S_e(t_e)^\top \Sigma_e^{-1}S_e(t_e)D_e
      Qp <- Matrix::sparseMatrix(i = i_, j = j_, x = x_, dims = Q.list$dims) 

      chol_Qp <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE) # This does Qp = P^\top L L^\top P

      half_log_det_Qp <- Matrix::determinant(chol_Qp, sqrt=TRUE)$modulus[1]
    }

    loglik <- loglik - half_log_det_Qp # adds 0.5 * log|Qp| n_repl times

    v <- c(as.matrix(Matrix::solve(chol_Qp, Matrix::solve(chol_Qp, Qp_mu, system = "P"), system = "L"))) # This L^{-1}P Qp mu

    loglik <- loglik + 
      0.5  * t(v) %*% v - # This is 0.5 * mu^\top Qp P^\top L^{-\top} L^{-1}P Qp mu = 0.5 * mu^\top Qp Qp^{-1} Qp mu = 0.5 * mu^\top Qp mu
      0.5 * numberOfObservations * log(2*pi)
  }

  return(loglik[1])
}



## --------------------------------------------------------------------------------------------------
loglikelihoodForAlphaEqualTwoPrecompute <- function(theta, 
                                                    precomputed_data, 
                                                    BC = 1, 
                                                    parameterization = "matern") {
  # Extract parameters
  sigma_e <- exp(theta[1])
  reciprocal_tau <- exp(theta[2])
  if(parameterization == "matern"){
    kappa = sqrt(8 * 1.5) / exp(theta[3])
  } else{
    kappa = exp(theta[3])
  }

  # Build Q matrix once
  Q <- spde_precision(kappa = kappa, tau = 1/reciprocal_tau,
                     alpha = 2, graph = precomputed_data$graph, BC=BC)

  R <- Matrix::Cholesky(forceSymmetric(precomputed_data$Tc%*%Q%*%t(precomputed_data$Tc)),
                       LDL = FALSE, perm = TRUE)

  # Get determinant once
  loglik <- 0
  det_R <- Matrix::determinant(R, sqrt=TRUE)$modulus[1]
  det_R_count <- NULL

  # Creating a large pre-allocated array for all potential entries
  total_max_entries <- 16 * length(precomputed_data$obs.edges)
  all_i <- all_j <- all_x <- numeric(total_max_entries)
  all_count <- 0

  n_obs_total <- 0

  # Process each replicate
  for(i in seq_along(precomputed_data$u_repl)) {
    curr_repl <- precomputed_data$u_repl[i]
    repl_name <- paste0("repl_", curr_repl)

    loglik <- loglik + det_R

    # Pre-allocate with exact size needed
    Qpmu <- numeric(4 * precomputed_data$n_edges)

    # Process each edge
    for(j in seq_along(precomputed_data$obs.edges)) {
      e <- precomputed_data$obs.edges[j]
      edge_name <- paste0("edge_", e)

      # Get data for this edge
      y_i <- precomputed_data$y_data[[repl_name]][[edge_name]]

      # Skip if no observations
      if(is.null(y_i) || length(y_i) == 0) {
        next
      }

      # Count observations for log determinant term
      n_obs_total <- n_obs_total + length(y_i)

      # Handle covariates if present
      if(precomputed_data$n_cov > 0) {
        X_cov_e <- precomputed_data$x_data[[repl_name]][[edge_name]]
        n_cov <- ncol(X_cov_e)
        if(n_cov > 0){
          y_i <- y_i - as.vector(X_cov_e %*% theta[4:(3+n_cov)])
        }
      }

      # Get precomputed distance data
      t <- precomputed_data$t_data[[repl_name]][[edge_name]]
      D <- precomputed_data$D_data[[repl_name]][[edge_name]]

      # Pre-allocate matrix
      n_pts <- length(t)
      S <- matrix(0, n_pts + 2, n_pts + 2)

      # Compute all submatrices
      d.index <- c(1,2)
      S[-d.index, -d.index] <- r_2(D, kappa = kappa,
                                  tau = 1/reciprocal_tau, deriv = 0)
      S[d.index, d.index] <- -r_2(as.matrix(dist(c(0,precomputed_data$edge_lengths[e]))),
                                 kappa = kappa, tau = 1/reciprocal_tau,
                                 deriv = 2)
      S[d.index, -d.index] <- -r_2(D[1:2,], kappa = kappa,
                                 tau = 1/reciprocal_tau, deriv = 1)
      S[-d.index, d.index] <- t(S[d.index, -d.index])

      # Covariance updates
      E.ind <- c(1:4)
      Obs.ind <- -E.ind
      Bt <- solve(S[E.ind, E.ind], S[E.ind, Obs.ind, drop = FALSE])
      Sigma_i <- S[Obs.ind,Obs.ind] - S[Obs.ind, E.ind] %*% Bt
      diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2

      # Cache Cholesky decomposition
      R_i <- base::chol(Sigma_i)

      # Sigma_iB <- backsolve(R_i, forwardsolve(t(R_i), t(Bt)))
      Sigma_iB <- solve(R_i, forwardsolve(t(R_i), t(Bt)))

      BtSinvB <- Bt %*% Sigma_iB

      E <- precomputed_data$graph$E[e, ]
      if (E[1] == E[2]) {
        warning("Circle not implemented")
      }

      BtSinvB <- BtSinvB[c(3,1,4,2), c(3,1,4,2)]
      Qpmu[4 * (e - 1) + 1:4] <- Qpmu[4 * (e - 1) + 1:4] +
        (t(Sigma_iB) %*% y_i)[c(3, 1, 4, 2)]

      # Efficiently add precision matrix entries - use direct indexing instead of loop
      # This is more efficient than using expand.grid in a loop
      idx <- seq(all_count + 1, all_count + 16)

      # Lower edge u diagonal
      all_i[idx[1]] <- 4 * (e - 1) + 1
      all_j[idx[1]] <- 4 * (e - 1) + 1
      all_x[idx[1]] <- BtSinvB[1, 1]

      # Lower edge u' diagonal
      all_i[idx[2]] <- 4 * (e - 1) + 2
      all_j[idx[2]] <- 4 * (e - 1) + 2
      all_x[idx[2]] <- BtSinvB[2, 2]

      # Upper edge u diagonal
      all_i[idx[3]] <- 4 * (e - 1) + 3
      all_j[idx[3]] <- 4 * (e - 1) + 3
      all_x[idx[3]] <- BtSinvB[3, 3]

      # Upper edge u' diagonal
      all_i[idx[4]] <- 4 * (e - 1) + 4
      all_j[idx[4]] <- 4 * (e - 1) + 4
      all_x[idx[4]] <- BtSinvB[4, 4]

      # Lower edge (u, u')
      all_i[idx[5]] <- 4 * (e - 1) + 1
      all_j[idx[5]] <- 4 * (e - 1) + 2
      all_x[idx[5]] <- BtSinvB[1, 2]

      all_i[idx[6]] <- 4 * (e - 1) + 2
      all_j[idx[6]] <- 4 * (e - 1) + 1
      all_x[idx[6]] <- BtSinvB[1, 2]

      # Upper edge (u, u')
      all_i[idx[7]] <- 4 * (e - 1) + 3
      all_j[idx[7]] <- 4 * (e - 1) + 4
      all_x[idx[7]] <- BtSinvB[3, 4]

      all_i[idx[8]] <- 4 * (e - 1) + 4
      all_j[idx[8]] <- 4 * (e - 1) + 3
      all_x[idx[8]] <- BtSinvB[3, 4]

      # Lower u, upper u
      all_i[idx[9]] <- 4 * (e - 1) + 1
      all_j[idx[9]] <- 4 * (e - 1) + 3
      all_x[idx[9]] <- BtSinvB[1, 3]

      all_i[idx[10]] <- 4 * (e - 1) + 3
      all_j[idx[10]] <- 4 * (e - 1) + 1
      all_x[idx[10]] <- BtSinvB[1, 3]

      # Lower u, upper u'
      all_i[idx[11]] <- 4 * (e - 1) + 1
      all_j[idx[11]] <- 4 * (e - 1) + 4
      all_x[idx[11]] <- BtSinvB[1, 4]

      all_i[idx[12]] <- 4 * (e - 1) + 4
      all_j[idx[12]] <- 4 * (e - 1) + 1
      all_x[idx[12]] <- BtSinvB[1, 4]

      # Lower u', upper u
      all_i[idx[13]] <- 4 * (e - 1) + 2
      all_j[idx[13]] <- 4 * (e - 1) + 3
      all_x[idx[13]] <- BtSinvB[2, 3]

      all_i[idx[14]] <- 4 * (e - 1) + 3
      all_j[idx[14]] <- 4 * (e - 1) + 2
      all_x[idx[14]] <- BtSinvB[2, 3]

      # Lower u', upper u'
      all_i[idx[15]] <- 4 * (e - 1) + 2
      all_j[idx[15]] <- 4 * (e - 1) + 4
      all_x[idx[15]] <- BtSinvB[2, 4]

      all_i[idx[16]] <- 4 * (e - 1) + 4
      all_j[idx[16]] <- 4 * (e - 1) + 2
      all_x[idx[16]] <- BtSinvB[2, 4]

      all_count <- all_count + 16

      # Compute quadratic form directly with Cholesky
      v_i <- backsolve(R_i, forwardsolve(t(R_i), y_i))
      quad_form <- sum(y_i * v_i)

      # Update log likelihood
      loglik <- loglik - 0.5 * quad_form - sum(log(diag(R_i)))
    }
  }

  # Build sparse matrix just once after collecting all entries
  if(all_count > 0) {
    BtSB <- Matrix::sparseMatrix(i = all_i[1:all_count],
                               j = all_j[1:all_count],
                               x = all_x[1:all_count],
                               dims = dim(Q))
    Qp <- Q + BtSB
    Qp <- precomputed_data$Tc %*% Qp %*% t(precomputed_data$Tc)
    R_count <- Matrix::Cholesky(forceSymmetric(Qp), LDL = FALSE, perm = TRUE)
    det_R_count <- Matrix::determinant(R_count, sqrt=TRUE)$modulus[1]

    for(i in seq_along(precomputed_data$u_repl)) {
      curr_repl <- precomputed_data$u_repl[i]
      repl_name <- paste0("repl_", curr_repl)

      loglik <- loglik - det_R_count

      v <- c(as.matrix(Matrix::solve(R_count, Matrix::solve(R_count, precomputed_data$Tc%*%Qpmu, system = 'P'),
                                     system='L')))

      loglik <- loglik + 0.5 * t(v) %*% v
    }
  }

  # Add constant term
  loglik <- loglik - 0.5 * n_obs_total * log(2 * pi)

  return(as.numeric(loglik))
}


