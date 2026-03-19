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


matern.p.precision <- function(loc,kappa,p,equally_spaced = FALSE, alpha = 1) {
  
  n <- length(loc)
  if(alpha==1) {
    Q <- exp_precision(loc,kappa)/(2*kappa)
    A <- Diagonal(n)
    
    return(list(Q=Q,A=A))
  } else {
    if(alpha%%1 == 0) {
      fa <- alpha
    } else {
      fa <- floor(alpha) + 1    
    }
    
    if(fa == 1) {
      N <- n  + n - 1 
    } else {
      N <- n*fa^2 + (n-1)*fa^2 - n*fa*(fa -1)/2    
    }
    
    ii <- numeric(N)
    jj <- numeric(N)
    val <- numeric(N)
    
    if(equally_spaced){
      Q.1 <- solve(rbind(cbind(matern.p.joint(loc[1],loc[1],kappa,p,alpha), 
                               matern.p.joint(loc[1],loc[1+1],kappa,p,alpha)),
                         cbind(matern.p.joint(loc[1+1],loc[1],kappa,p,alpha), 
                               matern.p.joint(loc[1+1],loc[1+1],kappa,p,alpha))))
      
      Q.i <- solve(rbind(cbind(matern.p.joint(loc[1],loc[1],kappa,p,alpha), 
                               matern.p.joint(loc[1],loc[2],kappa,p,alpha),
                               matern.p.joint(loc[1],loc[3],kappa,p,alpha)),
                         cbind(matern.p.joint(loc[2],loc[1],kappa,p,alpha),
                               matern.p.joint(loc[2],loc[2],kappa,p,alpha),
                               matern.p.joint(loc[2],loc[3],kappa,p,alpha)),
                         cbind(matern.p.joint(loc[3],loc[1],kappa,p,alpha),
                               matern.p.joint(loc[3],loc[2],kappa,p,alpha),
                               matern.p.joint(loc[3],loc[3],kappa,p,alpha))))[-c(1:fa),-c(1:fa)]
    }
    
    
    for(i in 1:max((n-1),1)){
      if(i==1){
        if(!equally_spaced){
          Q.1 <- solve(rbind(cbind(matern.p.joint(loc[i],loc[i],kappa,p,alpha),
                                   matern.p.joint(loc[i],loc[i+1],kappa,p,alpha)),
                             cbind(matern.p.joint(loc[i+1],loc[i],kappa,p,alpha),
                                   matern.p.joint(loc[i+1],loc[i+1],kappa,p,alpha))))
          
        } 
        counter <- 1
        for(ki in 1:fa) {
          for(kj in ki:(2*fa)) {
            ii[counter] <- ki
            jj[counter] <- kj
            val[counter] <- Q.1[ki,kj]
            counter <- counter + 1
          }
        }
      } else {
        if(!equally_spaced){
          Q.i <- solve(rbind(cbind(matern.p.joint(loc[i-1],loc[i-1],kappa,p,alpha),
                                   matern.p.joint(loc[i-1],loc[i],kappa,p,alpha),
                                   matern.p.joint(loc[i-1],loc[i+1],kappa,p,alpha)),
                             cbind(matern.p.joint(loc[i],loc[i-1],kappa,p,alpha),
                                   matern.p.joint(loc[i],loc[i],kappa,p,alpha),
                                   matern.p.joint(loc[i],loc[i+1],kappa,p,alpha)),
                             cbind(matern.p.joint(loc[i+1],loc[i-1],kappa,p,alpha),
                                   matern.p.joint(loc[i+1],loc[i],kappa,p,alpha),
                                   matern.p.joint(loc[i+1],loc[i+1],kappa,p,alpha))))[-c(1:fa),-c(1:fa)]
        } 
        # Q[(2*n-1):(2*n), (2*n-3):(2*n)] = Q.i[3:4,]
        for(ki in 1:fa){
          for(kj in ki:(2*fa)){
            ii[counter] <- fa*(i-1) + ki
            jj[counter] <- fa*(i-1) + kj
            val[counter] <-Q.i[ki,kj]
            counter <- counter + 1
          }
        }     
      }
    }
    if(n<=2){
      Q.i <- Q.1
    }
    for(ki in 1:fa){
      for(kj in ki:fa){
        ii[counter] <- fa*(n-1) + ki
        jj[counter] <- fa*(n-1) + kj
        val[counter] <-Q.i[ki+fa,kj+fa]
        counter <- counter + 1
      }
    }    
    Q <- Matrix::sparseMatrix(i   = ii,
                              j    = jj,
                              x    = val,
                              dims = c(fa*n, fa*n), symmetric=TRUE)
    
    A <-  Matrix::sparseMatrix(i   = 1:n,
                               j    = seq(from=1,to=n*fa,by=fa),
                               x    = rep(1,n),
                               dims = c(n, fa*n))
    return(list(Q=Q,A=A))    
  }
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

matern.covariance <- function(h, kappa, nu, sigma) {
  if (nu == 1 / 2) {
    C <- sigma^2 * exp(-kappa * abs(h))
  } else {
    C <- (sigma^2 / (2^(nu - 1) * gamma(nu))) *
      ((kappa * abs(h))^nu) * besselK(kappa * abs(h), nu)
    C[h == 0] <- sigma^2
  }
  return(C)
}


matern.derivative <- function(h, kappa, nu, sigma, deriv = 1) 
{
  if(deriv == 0) {
    C = matern.covariance(h, kappa = kappa, nu = nu, sigma = sigma)
    return(C)
  } else if (deriv == 1) {
    C = h * matern.covariance(h, kappa = kappa, nu = nu - 1, sigma = sigma)
    C[h == 0] = 0
    return(-(kappa^2/(2 * (nu - 1))) * C)
  }
  else if (deriv == 2) {
    C = matern.covariance(h, kappa = kappa, nu = nu - 1, sigma = sigma) + 
      h * matern.derivative(h, kappa = kappa, nu = nu - 1, 
                            sigma = sigma, deriv = 1)
    return(-(kappa^2/(2 * (nu - 1))) * C)
  }
  else {
    C = (deriv - 1) * matern.derivative(h, kappa = kappa, 
                                        nu = nu - 1, sigma = sigma, 
                                        deriv = deriv - 2) + 
      h * matern.derivative(h, kappa = kappa, nu = nu - 
                              1, sigma = sigma, deriv = deriv - 1)
    return(-(kappa^2/(2 * (nu - 1))) * C)
  }
  
}


matern.p <- function(s,t,kappa,p,alpha){
  h <- s-t
  if(p==0){
    return(matern.covariance(h, kappa = kappa, nu = alpha - 1/2, sigma = 1))
  } else {
    ca <- gamma(alpha)/gamma(alpha-0.5)
    fa <- floor(alpha)
    kp <- kappa*sqrt(1-p)
    out <- matern.covariance(h, kappa = kp, nu = 1/2, 
                             sigma = sqrt(ca*sqrt(pi)/sqrt(1-p)))
    if(alpha < 1) {
      return(out)
    } else {
      
      for(j in 1:fa) {
        out <- out - matern.covariance(h, kappa = kappa, nu = j-1/2, 
                                       sigma = sqrt(ca*gamma(j-0.5)/gamma(j)))/p^(1 - j)
      }
      out <- out/p^fa
      return(out)    
    }
  }
}

matern.p.deriv <- function(s,t,kappa,p,alpha,deriv = 0){
  h <- s-t
  if(deriv ==0){
    return(matern.p(s,t,kappa,p,alpha))
  } else {
    if(p==0){
      return(matern.derivative(h, kappa = kappa, nu = alpha - 1/2, 
                               sigma = 1, deriv = deriv))
    } else {
      ca <- gamma(alpha)/gamma(alpha-0.5)
      fa <- floor(alpha)
      kp <- kappa*sqrt(1-p)
      out <- matern.derivative(h, kappa = kp, nu =  1/2, 
                               sigma = sqrt(ca*sqrt(pi)/sqrt(1-p)), deriv = deriv)
      if(alpha < 1) {
        return(out)
      } else {
        out <- out/p^fa
        for(j in 1:fa) {
          out <- out - matern.derivative(h, kappa = kappa, nu = j-1/2, 
                                         sigma = sqrt(ca*gamma(j-0.5)/gamma(j)), 
                                         deriv = deriv)/p^(fa + 1 - j)
        }
      }
      return(out)
    }    
  }
  
}

matern.p.joint <- function(s,t,kappa,p, alpha = 1){
  
  if(alpha%%1 == 0) {
    fa <- alpha
  } else {
    fa <- floor(alpha) + 1    
  }
  mat <- matrix(0, nrow = fa, ncol = fa)
  for(i in 1:fa) {
    for(j in i:fa) {
      if(i==j) {
        mat[i,i] <- ((-1)^(i-1))*matern.p.deriv(s,t,kappa,p,alpha, deriv = 2*(i-1))
      } else {
        tmp <- matern.p.deriv(s,t,kappa,p,alpha, deriv = i-1 + j - 1)
        mat[i,j] <- (-1)^(j-1)*tmp
        mat[j,i] <- (-1)^(i-1)*tmp
      }
    }
  }
  
  return(mat)
}


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
  
  # compute parameters
  fa <- floor(alpha)
  ca <- ceiling(alpha)
  
  nu <- alpha - 1/2
  sigma <- sqrt(gamma(nu) / (tau^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))
  c_alpha <- gamma(alpha)/gamma(alpha - 0.5)
  c_1 <- gamma(fa)/gamma(fa - 0.5)
  
  # get edge lengths
  L_e <- graph$edge_lengths
  
  # initialize Qtilde_i, a list containing block diagonal matrices with blocks Qtilde_{i,e} for each i
  Qtilde_i <- list() 
  for(i in 1:m){
    
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
    for(e in 1:length(L_e)){
      
      # compute Q_{i,e}
      Q_e <- matern.p.precision(
        loc = c(0, L_e[e]),
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
  
  # --------------------------------------------------
  # CASE i = 0
  # --------------------------------------------------
  
  factor_0 <- c_1/(2 * k * c_alpha * kappa * sigma^2 * tau^2)
  
  Qtilde_0_star_UU <- MetricGraph:::Qalpha1(
    theta = c(tau, kappa), 
    graph = graph, 
    BC = 3000, 
    build = TRUE) * factor_0
  
  A_0 <- graph$.__enclos_env__$private$A()
  
  # --------------------------------------------------
  # CASE i = 1,...,m
  # --------------------------------------------------
  
  # build conditioning matrix
  graph$buildC(alpha = 2, edge_constraint = TRUE) # should always be TRUE
  COND_i <- graph$CoB
  Tc <- COND_i$T[-c(1:length(COND_i$S)), ]
  
  factor_i <- (2 * kappa^(2 * alpha - 1) * c_alpha * sqrt(pi) * tau^2)/r
  Qtilde_i_star_UU <- purrr::map2(
    Qtilde_i, 
    factor_i, 
    function(Q, x) Tc %*% Q %*% t(Tc) * x)
  
  index.obs_i <- gives.indices(graph = graph, factor = 4, constant = 3)
  A_i <- t(Tc)[index.obs_i, ] 
  
  # Build matrix A and Q_UU
  A <- cbind(A_0, do.call(cbind, rep(list(A_i), m)))
  Q_UU <- bdiag(Qtilde_0_star_UU, bdiag(Qtilde_i_star_UU))
  # Return Sigma
  Sigma <- A %*% solve(Q_UU, t(A)) 
  return(Sigma)
}


