matern.p <- function(s,t,kappa,p,tau,deriv=0){
  h <- s-t
  if(p==0){
   if(deriv==0){
      return(matern.covariance(h,kappa=kappa,nu=3/2,sigma=1/tau))
    } else {
      return(matern.derivative(h, kappa=kappa, nu=3/2, sigma=1/tau,deriv=deriv))
    }
  } else {
    a <- -1/(2*kappa*p)
    b <- sqrt(kappa^2-p)
    if(deriv==0){
      return(tau^(-2)*a*(exp(-kappa*abs(h))-(kappa/b)*exp(-b*abs(h))))
    } else if(deriv==1){
      return(-tau^(-2)*a*kappa*sign(h)*(-exp(-kappa*abs(h)) + exp(-b*abs(h))))
    } else if(deriv==2){
      return(tau^(-2)*a*kappa*(kappa*exp(-kappa*abs(h))-b*exp(-b*abs(h))))
     } 
    else{
      stop("only deriv=0,1,2,3,4 allowed")
    }
  }
}

match_order <- function(coord, V) {
  matched_indices <- integer(nrow(V))
  
  for (i in 1:nrow(V)) {
    match_index <- which(apply(coord, 1, function(x) all.equal(x, V[i, ])) == TRUE)
    if (length(match_index) == 0) {
      matched_indices[i] <- NA
    } else {
      matched_indices[i] <- match_index[1]
    }
  }
  
  return(matched_indices)
}

graph_Qalpha1 <- function(kappa,tau,graph,BC = 0){
i_ <- j_ <- x_ <- rep(0, graph$nE * 4)
  count <- 0

for (i in 1:graph$nE) {
  l_e <- graph$edge_lengths[i]
  Q_edge<-((kappa*tau^2)/(exp(2*kappa*l_e)-1))*
          matrix(c(exp(2*kappa*l_e)+1, -2*exp(kappa*l_e),
          -2*exp(kappa*l_e), exp(2*kappa*l_e)+1),2,2)

if (graph$E[i, 1] == graph$E[i, 2]) {
      warning("Circular edges are not implemented")
    }

      #lower edge precision u
      i_[count + 1] <- 2 * (i - 1) + 1
      j_[count + 1] <- 2 * (i - 1) + 1
      x_[count + 1] <- Q_edge[1, 1]

      #lower edge  u'
      i_[count + 2] <- 2 * (i - 1) + 2
      j_[count + 2] <- 2* (i - 1) + 2
      x_[count + 2] <- Q_edge[2, 2]

      i_[count + 3] <- 2 * (i - 1) + 1
      j_[count + 3] <- 2 * (i - 1) + 2
      x_[count + 3] <- Q_edge[1, 2]

      i_[count + 4] <- 2 * (i - 1) + 2
      j_[count + 4] <- 2 * (i - 1) + 1
      x_[count + 4] <- Q_edge[1, 2]

      count <- count + 4
}
 Q <- Matrix::sparseMatrix(i   = i_[1:count],
                            j    = j_[1:count],
                            x    = x_[1:count],
                            dims = c(2 * graph$nE, 2 * graph$nE))
return(Q)

}

conditioning<-function(graph,alpha=1){
      i_  =  rep(0, 2 * graph$nE)
      j_  =  rep(0, 2 * graph$nE)
      x_  =  rep(0, 2 * graph$nE)

      count_constraint <- 0
      count <- 0
      for (v in 1:graph$nV) {
        lower_edges  <- which(graph$E[, 1] %in% v)
        upper_edges  <- which(graph$E[, 2] %in% v)
        n_e <- length(lower_edges) + length(upper_edges)
        if (n_e > 1) {
          if (length(upper_edges) == 0) {
            edges <- cbind(lower_edges, 1)
          } else if(length(lower_edges) == 0){
            edges <- cbind(upper_edges, 2)
          }else{
            edges <- rbind(cbind(lower_edges, 1),
                           cbind(upper_edges, 2))
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
        # return(K)
                             
        CB <- MetricGraph:::c_basis2(K)
      
}

conditioningalpha2 = function(alpha = 2, self, edge_constraint = TRUE) {

    if(alpha==2){
      i_  =  rep(0, 2 * self$nE)
      j_  =  rep(0, 2 * self$nE)
      x_  =  rep(0, 2 * self$nE)

      count_constraint <- 0
      count <- 0
      for (v in 1:self$nV) {
        lower_edges  <- which(self$E[, 1] %in% v)
        upper_edges  <- which(self$E[, 2] %in% v)
        n_e <- length(lower_edges) + length(upper_edges)

        #derivative constraint
        if ((edge_constraint & n_e ==1) | n_e > 1) {
          i_[count + 1:n_e] <- count_constraint + 1
          j_[count + 1:n_e] <- c(4 * (lower_edges-1) + 2, 4 * (upper_edges-1) + 4)
          x_[count + 1:n_e] <- c(rep(1,length(lower_edges)),
                                 rep(-1,length(upper_edges)))
          count <- count + n_e
          count_constraint <- count_constraint + 1
        }
        if (n_e > 1) {
          if (length(upper_edges) == 0) {
            edges <- cbind(lower_edges, 1)
          } else if(length(lower_edges) == 0){
            edges <- cbind(upper_edges, 3)
          }else{
            edges <- rbind(cbind(lower_edges, 1),
                           cbind(upper_edges, 3))
          }
          for (i in 2:n_e) {
            i_[count + 1:2] <- count_constraint + 1
            j_[count + 1:2] <- c(4 * (edges[i-1,1] - 1) + edges[i-1, 2],
                                 4 * (edges[i,1]   - 1) + edges[i,   2])
            x_[count + 1:2] <- c(1,-1)
            count <- count + 2
            count_constraint <- count_constraint + 1
          }
        }
      }
      C <- Matrix::sparseMatrix(i = i_[1:count],
                                j = j_[1:count],
                                x = x_[1:count],
                                dims = c(count_constraint, 4*self$nE))
      # C<-rbind(C,c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0))

      CB <- MetricGraph:::c_basis2(C)
    }else{
      error("only alpha=2 implemented")
    }
  }

precision_alpha1 <- function(kappa,tau,nu,graph,m,BC = 1, type_rational_approx = "brasil"){

alpha=nu+1/2
coeff <- rSPDE:::interp_rational_coefficients(order = m, 
                                              type_rational_approx = type_rational_approx, 
                                              type_interp = "spline", 
                                              alpha = alpha)
        r <- coeff$r
        p <- coeff$p
        k <- coeff$k
  
I <- Diagonal(graph$nV)
Q <- I/k
A <- Matrix::Diagonal(graph$nV)[graph$PtV, ]

for (i in 1:length(p))
{
Q1<-graph_Qalpha1(sqrt(kappa^2*(1-p[i])),tau,graph,BC = 0) # Unconstrained precision for alpha=1  
Q1<-Q1/(r[i]*kappa^2)
cbmat<-conditioning(graph,alpha=1)
nc1<-1:length(cbmat$S)
Qtilde1<- t(cbmat$T) %*% Q1 %*% (cbmat$T)
W<-Diagonal(dim(Q1)[1])
W <- W[,-nc1]
Qtilde1_uu<-t(W)%*% Qtilde1%*%W

A1 <- (cbmat$T[,-nc1])

index.obs1 <- sapply(graph$PtV, function(i){idx_temp <- i == graph$E[,1]
                                                                      idx_temp <- which(idx_temp)
                                                                      return(idx_temp[1])})
    index.obs1 <- (index.obs1-1)*2+1
    index.obs2 <- NULL
    na_obs1 <- is.na(index.obs1)
    if(any(na_obs1)){
          idx_na <- which(na_obs1)
          PtV_NA <- graph$PtV[idx_na]
          index.obs2 <- sapply(PtV_NA, function(i){idx_temp <- i == graph$E[,2]
                                                                      idx_temp <- which(idx_temp)
                                                                      return(idx_temp[1])})
          index.obs1[na_obs1] <- (index.obs2-1)*2 + 2                                                                      
          }

A1 <- A1[index.obs1,]

Q<-bdiag(Q,Qtilde1_uu)
A<-cbind(A,A1)
}
Q<- Q*kappa^(2*alpha)
return(list(Q=Q, A = A))
}


precision_alpha1_direct <- function(kappa,tau,nu,graph,m,BC = 1){

alpha=nu+1/2
coeff <- rSPDE:::interp_rational_coefficients(order = m, 
                                              type_rational_approx = "chebfun", 
                                              type_interp = "spline", 
                                              alpha = alpha)
        r <- coeff$r
        p <- coeff$p
        k <- coeff$k
  
I <- Diagonal(graph$nV)
Q <- I/k
A <- Matrix::Diagonal(graph$nV)[graph$PtV, ]

for (i in 1:length(p)){
Q<-bdiag(Q, spde_precision(sqrt(kappa^2*(1-p[i])),tau,alpha=1,graph, BC = BC)/(r[i]*kappa^2))
A = cbind(A, Matrix::Diagonal(graph$nV)[graph$PtV, ])
}

Q <- Q*kappa^(2*alpha)
return(list(Q= Q, A = A))
}



precision_alpha2 <- function(kappa,tau,nu,graph,m,BC = 0){

alpha=nu+1/2

coeff <- rSPDE:::interp_rational_coefficients(order = m, 
                                              type_rational_approx = "chebfun", 
                                              type_interp = "spline", 
                                              alpha = alpha)
        r <- coeff$r
        p <- coeff$p
        k <- coeff$k
  
Q1<-graph_Qalpha1(kappa,tau,graph,BC = 0) # Unconstrained precision for alpha=1
Q1<-Q1/(k*kappa^2)
cbmat<-conditioning(graph,alpha=1)
nc1<-1:length(cbmat$S)
Qtilde1<- t(cbmat$T) %*% Q1 %*% (cbmat$T)
W<-Diagonal(dim(Q1)[1])
W <- W[,-nc1]
Qtilde1_uu<-t(W)%*% Qtilde1%*%W
# Qtilde1_uu2<-Qtilde1[-nc1, -nc1]

A <- (cbmat$T[,-nc1])

#  index.obs1 <- (2 * graph$get_PtE()[,1] - 1) * (graph$get_PtE()[,2] < 1e-14) +
#                    2 * (graph$get_PtE()[,1]) * (graph$get_PtE()[,2] > 1e-14)


index.obs1 <- sapply(graph$PtV, function(i){idx_temp <- i == graph$E[,1]
                                                                      idx_temp <- which(idx_temp)
                                                                      return(idx_temp[1])})
    index.obs1 <- (index.obs1-1)*2+1
    index.obs4 <- NULL
    na_obs1 <- is.na(index.obs1)
    if(any(na_obs1)){
          idx_na <- which(na_obs1)
          PtV_NA <- graph$PtV[idx_na]
          index.obs4 <- sapply(PtV_NA, function(i){idx_temp <- i == graph$E[,2]
                                                                      idx_temp <- which(idx_temp)
                                                                      return(idx_temp[1])})
          index.obs1[na_obs1] <- (index.obs4-1)*2 + 2                                                                      
          }

 A <- A[index.obs1,] #A matrix for alpha=1

for (ii in 1:length(p)){

i_ <- j_ <- x_ <- rep(0, graph$nE * 16)
  count <- 0

  R_00 <- matrix(c( matern.p(0,0, kappa = kappa, kappa^2*p[ii], tau = tau, deriv = 0),
                   -matern.p(0,0, kappa = kappa, kappa^2*p[ii], tau = tau, deriv = 1),
                    matern.p(0,0, kappa = kappa, kappa^2*p[ii], tau = tau, deriv = 1),
                   -matern.p(0,0, kappa = kappa, kappa^2*p[ii], tau = tau, deriv = 2)), 2, 2)
  R_node <- rbind(cbind(R_00, matrix(0, 2, 2)),
                  cbind(matrix(0, 2, 2), R_00))
  R00i <- solve(R_00)
  Ajd <- -1 * rbind(cbind(0.5 * R00i, matrix(0, 2, 2)),
                      cbind(matrix(0, 2, 2), (1-0.5)*R00i))
  for (i in 1:graph$nE) {

    l_e <- graph$edge_lengths[i]
    r_0l <-   matern.p(l_e,0, kappa = kappa, kappa^2*p[ii], tau = tau, deriv = 0)
    r_11 <- - matern.p(l_e,0, kappa = kappa, kappa^2*p[ii], tau = tau, deriv = 2)
    # order by node not derivative
     R_01 <- matrix(c(r_0l, -matern.p(l_e, 0, kappa = kappa, kappa^2*p[ii], tau = tau, deriv = 1),
                     matern.p(l_e, 0, kappa = kappa, kappa^2*p[ii], tau = tau, deriv = 1), r_11), 2, 2)

    R_node[1:2, 3:4] <- R_01
    R_node[3:4, 1:2] <- t(R_01)
    Q_adj <- solve(R_node) + Ajd
    if (graph$E[i, 1] == graph$E[i, 2]) {
      warning("Circular edges are not implemented")
    }

      #lower edge precision u
      i_[count + 1] <- 4 * (i - 1) + 1
      j_[count + 1] <- 4 * (i - 1) + 1
      x_[count + 1] <- Q_adj[1, 1]

      #lower edge  u'
      i_[count + 2] <- 4 * (i - 1) + 2
      j_[count + 2] <- 4 * (i - 1) + 2
      x_[count + 2] <- Q_adj[2, 2]

      #upper edge  u
      i_[count + 3] <- 4 * (i - 1) + 3
      j_[count + 3] <- 4 * (i - 1) + 3
      x_[count + 3] <- Q_adj[3, 3]

      #upper edge  u'
      i_[count + 4] <- 4 * (i - 1) + 4
      j_[count + 4] <- 4 * (i - 1) + 4
      x_[count + 4] <- Q_adj[4, 4]

      #lower edge  u, u'
      i_[count + 5] <- 4 * (i - 1) + 1
      j_[count + 5] <- 4 * (i - 1) + 2
      x_[count + 5] <- Q_adj[1, 2]
      i_[count + 6] <- 4 * (i - 1) + 2
      j_[count + 6] <- 4 * (i - 1) + 1
      x_[count + 6] <- Q_adj[1, 2]

      #upper edge  u, u'
      i_[count + 7] <- 4 * (i - 1) + 3
      j_[count + 7] <- 4 * (i - 1) + 4
      x_[count + 7] <- Q_adj[3, 4]
      i_[count + 8] <- 4 * (i - 1) + 4
      j_[count + 8] <- 4 * (i - 1) + 3
      x_[count + 8] <- Q_adj[3, 4]

      #lower edge  u, upper edge  u,
      i_[count + 9]  <- 4 * (i - 1) + 1
      j_[count + 9]  <- 4 * (i - 1) + 3
      x_[count + 9]  <- Q_adj[1, 3]
      i_[count + 10] <- 4 * (i - 1) + 3
      j_[count + 10] <- 4 * (i - 1) + 1
      x_[count + 10] <- Q_adj[1, 3]

      #lower edge  u, upper edge  u',
      i_[count + 11] <- 4 * (i - 1) + 1
      j_[count + 11] <- 4 * (i - 1) + 4
      x_[count + 11] <- Q_adj[1, 4]
      i_[count + 12] <- 4 * (i - 1) + 4
      j_[count + 12] <- 4 * (i - 1) + 1
      x_[count + 12] <- Q_adj[1, 4]

      #lower edge  u', upper edge  u,
      i_[count + 13] <- 4 * (i - 1) + 2
      j_[count + 13] <- 4 * (i - 1) + 3
      x_[count + 13] <- Q_adj[2, 3]
      i_[count + 14] <- 4 * (i - 1) + 3
      j_[count + 14] <- 4 * (i - 1) + 2
      x_[count + 14] <- Q_adj[2, 3]

      #lower edge  u', upper edge  u',
      i_[count + 15] <- 4 * (i - 1) + 2
      j_[count + 15] <- 4 * (i - 1) + 4
      x_[count + 15] <- Q_adj[2, 4]
      i_[count + 16] <- 4 * (i - 1) + 4
      j_[count + 16] <- 4 * (i - 1) + 2
      x_[count + 16] <- Q_adj[2, 4]

      count <- count + 16

  }
  if(BC == 1){
    #Vertices with of degree 1
    i.table <- table(c(graph$E))
    index <- as.integer(names(which(i.table == 1)))
    #for this vertices locate position
    lower.edges <- which(graph$E[, 1] %in% index)
    upper.edges <- which(graph$E[, 2] %in% index)
    for (le in lower.edges) {
      ind <- c(4 * (le - 1) + 1, 4 * (le - 1) + 2)

      i_ <- c(i_, ind)
      j_ <- c(j_, ind)
      x_ <- c(x_, 0.5*c(1 / R_00[1, 1], 1 / R_00[2, 2]))
      count <- count + 2
    }
    for (ue in upper.edges) {
      ind <- c(4 * (ue - 1) + 3, 4 * (ue - 1) + 4)
      i_ <- c(i_, ind)
      j_ <- c(j_, ind)
      x_ <- c(x_, 0.5 * c(1 / R_00[1, 1], 1 / R_00[2, 2]))
      count <- count + 2
    }
  }
  #if (build) {
    Q0 <- Matrix::sparseMatrix(i   = i_[1:count],
                              j    = j_[1:count],
                              x    = x_[1:count],
                              dims = c(4 * graph$nE, 4 * graph$nE))


   Q0<-Q0/(r[ii]*kappa^4)


############################################# Use this part for my A and  K

     index.obs2 <- sapply(graph$PtV, function(i){idx_temp <- i == graph$E[,1]
                                                                      idx_temp <- which(idx_temp)
                                                                      return(idx_temp[1])})
    index.obs2 <- (index.obs2-1)*4+1
    index.obs3 <- NULL
    na_obs1 <- is.na(index.obs2)
    if(any(na_obs1)){
          idx_na <- which(na_obs1)
          PtV_NA <- graph$PtV[idx_na]
          index.obs3 <- sapply(PtV_NA, function(i){idx_temp <- i == graph$E[,2]
                                                                      idx_temp <- which(idx_temp)
                                                                      return(idx_temp[1])})
          index.obs2[na_obs1] <- (index.obs3-1)*4 + 3                                                                      
    }

    cbmat2<-conditioningalpha2(alpha=2, graph)   
    nc2<-1:length(cbmat2$S)
    Qtilde<- t(cbmat2$T) %*% Q0 %*% (cbmat2$T)
    W2<-Diagonal(dim(Q0)[1])
    W2 <- W2[,-nc2]
    Qtilde_uu<-t(W2)%*% Qtilde%*%W2
    A1<-(cbmat2$T[,-nc2])
    A1 <- A1[index.obs2, ] # A matrix for alpha=2
    Qtilde1_uu<-bdiag(Qtilde1_uu,Qtilde_uu)
    A <- cbind(A,A1)

}

   Qtilde1_uu <- Qtilde1_uu*kappa^(2*alpha)
   return(list(Q = Qtilde1_uu, A = A))
}



shiftedrational_Qalpha1 <- function(kappa,tau,nu,graph,m,BC = 1){

alpha=nu+1/2

  coeff <- rSPDE:::interp_rational_coefficients(order = m, 
                                              type_rational_approx = "chebfun", 
                                              type_interp = "spline", 
                                              alpha = alpha)
        r <- coeff$r
        p <- coeff$p
        k <- coeff$k
  

I <- Diagonal(graph$nV)     # in case of alpha=1 gives Q with conditioning.
Q<-I/k
A<-I

for (i in 1:length(p))
{
Q1<-graph_Qalpha1(sqrt(kappa^2*(1-p[i])),tau,graph,BC = 0) # Unconstrained precision for alpha=1  
Q1<-Q1/(r[i]*kappa^2)
cbmat<-conditioning(graph,alpha=1)
nc1<-1:length(cbmat$S)
Qtilde1<- t(cbmat$T) %*% Q1 %*% (cbmat$T)
W<-Diagonal(dim(Q1)[1])
W <- W[,-nc1]
Qtilde1_uu<-t(W)%*% Qtilde1%*%W
# Qtilde1_uu2<-Qtilde1[-nc1, -nc1]

A1 <- (cbmat$T[,-nc1])
#  index.obs1 <- (2 * graph$get_PtE()[,1] - 1) * (graph$get_PtE()[,2] < 1e-14) +
#                    2 * (graph$get_PtE()[,1]) * (graph$get_PtE()[,2] > 1e-14)

  index.obs1 <- sapply(graph$PtV, function(i){idx_temp <- i == graph$E[,1]
                                                                      idx_temp <- which(idx_temp)
                                                                      return(idx_temp[1])})
    index.obs1 <- (index.obs1-1)*2+1
    index.obs2 <- NULL
    na_obs1 <- is.na(index.obs1)
    if(any(na_obs1)){
          idx_na <- which(na_obs1)
          PtV_NA <- graph$PtV[idx_na]
          index.obs2 <- sapply(PtV_NA, function(i){idx_temp <- i == graph$E[,2]
                                                                      idx_temp <- which(idx_temp)
                                                                      return(idx_temp[1])})
          index.obs1[na_obs1] <- (index.obs2-1)*2 + 2                                                                      
          }
A1 <- A1[index.obs1,]

Q<-bdiag(Q,Qtilde1_uu)
A<-cbind(A,A1)
}
Q<- Q*kappa^(2*alpha)
sigma=A%*%solve(Q,t(A)) 
return(sigma)
}

shiftedrational_Qalpha2 <- function(kappa,tau,nu,graph,m,BC = 0){

alpha=nu+1/2

coeff <- rSPDE:::interp_rational_coefficients(order = m, 
                                              type_rational_approx = "chebfun", 
                                              type_interp = "spline", 
                                              alpha = alpha)
        r <- coeff$r
        p <- coeff$p
        k <- coeff$k
   

Q1<-graph_Qalpha1(kappa,tau,graph,BC = 0) # Unconstrained precision for alpha=1
Q1<-Q1/(k*kappa^2)
cbmat<-conditioning(graph,alpha=1)
nc1<-1:length(cbmat$S)
Qtilde1<- t(cbmat$T) %*% Q1 %*% (cbmat$T)
W<-Diagonal(dim(Q1)[1])
W <- W[,-nc1]
Qtilde1_uu<-t(W)%*% Qtilde1%*%W


A <- (cbmat$T[,-nc1])

#  index.obs1 <- (2 * graph$get_PtE()[,1] - 1) * (graph$get_PtE()[,2] < 1e-14) +
#                    2 * (graph$get_PtE()[,1]) * (graph$get_PtE()[,2] > 1e-14)

index.obs1 <- sapply(graph$PtV, function(i){idx_temp <- i == graph$E[,1]
                                                                      idx_temp <- which(idx_temp)
                                                                      return(idx_temp[1])})
    index.obs1 <- (index.obs1-1)*2+1
    index.obs4 <- NULL
    na_obs1 <- is.na(index.obs1)
    if(any(na_obs1)){
          idx_na <- which(na_obs1)
          PtV_NA <- graph$PtV[idx_na]
          index.obs4 <- sapply(PtV_NA, function(i){idx_temp <- i == graph$E[,2]
                                                                      idx_temp <- which(idx_temp)
                                                                      return(idx_temp[1])})
          index.obs1[na_obs1] <- (index.obs4-1)*2 + 2                                                                      
          }

 A <- A[index.obs1,] #A matrix for alpha=1

for (ii in 1:length(p)){

i_ <- j_ <- x_ <- rep(0, graph$nE * 16)
  count <- 0

  R_00 <- matrix(c( matern.p(0,0, kappa = kappa, kappa^2*p[ii], tau = tau, deriv = 0),
                   -matern.p(0,0, kappa = kappa, kappa^2*p[ii], tau = tau, deriv = 1),
                    matern.p(0,0, kappa = kappa, kappa^2*p[ii], tau = tau, deriv = 1),
                   -matern.p(0,0, kappa = kappa, kappa^2*p[ii], tau = tau, deriv = 2)), 2, 2)

  R_node <- rbind(cbind(R_00, matrix(0, 2, 2)),
                  cbind(matrix(0, 2, 2), R_00))
  R00i <- solve(R_00)
  Ajd <- -1 * rbind(cbind(0.5 * R00i, matrix(0, 2, 2)),
                      cbind(matrix(0, 2, 2), (1-0.5)*R00i))
  for (i in 1:graph$nE) {

    l_e <- graph$edge_lengths[i]
    r_0l <-   matern.p(l_e,0, kappa = kappa, kappa^2*p[ii], tau = tau, deriv = 0)
    r_11 <- - matern.p(l_e,0, kappa = kappa, kappa^2*p[ii], tau = tau, deriv = 2)
  
     R_01 <- matrix(c(r_0l, -matern.p(l_e, 0, kappa = kappa, kappa^2*p[ii], tau = tau, deriv = 1),
                     matern.p(l_e, 0, kappa = kappa, kappa^2*p[ii], tau = tau, deriv = 1), r_11), 2, 2)
    

    R_node[1:2, 3:4] <- R_01
    R_node[3:4, 1:2] <- t(R_01)
    Q_adj <- solve(R_node) + Ajd
    if (graph$E[i, 1] == graph$E[i, 2]) {
      warning("Circular edges are not implemented")
    }

      #lower edge precision u
      i_[count + 1] <- 4 * (i - 1) + 1
      j_[count + 1] <- 4 * (i - 1) + 1
      x_[count + 1] <- Q_adj[1, 1]

      #lower edge  u'
      i_[count + 2] <- 4 * (i - 1) + 2
      j_[count + 2] <- 4 * (i - 1) + 2
      x_[count + 2] <- Q_adj[2, 2]

      #upper edge  u
      i_[count + 3] <- 4 * (i - 1) + 3
      j_[count + 3] <- 4 * (i - 1) + 3
      x_[count + 3] <- Q_adj[3, 3]

      #upper edge  u'
      i_[count + 4] <- 4 * (i - 1) + 4
      j_[count + 4] <- 4 * (i - 1) + 4
      x_[count + 4] <- Q_adj[4, 4]

      #lower edge  u, u'
      i_[count + 5] <- 4 * (i - 1) + 1
      j_[count + 5] <- 4 * (i - 1) + 2
      x_[count + 5] <- Q_adj[1, 2]
      i_[count + 6] <- 4 * (i - 1) + 2
      j_[count + 6] <- 4 * (i - 1) + 1
      x_[count + 6] <- Q_adj[1, 2]

      #upper edge  u, u'
      i_[count + 7] <- 4 * (i - 1) + 3
      j_[count + 7] <- 4 * (i - 1) + 4
      x_[count + 7] <- Q_adj[3, 4]
      i_[count + 8] <- 4 * (i - 1) + 4
      j_[count + 8] <- 4 * (i - 1) + 3
      x_[count + 8] <- Q_adj[3, 4]

      #lower edge  u, upper edge  u,
      i_[count + 9]  <- 4 * (i - 1) + 1
      j_[count + 9]  <- 4 * (i - 1) + 3
      x_[count + 9]  <- Q_adj[1, 3]
      i_[count + 10] <- 4 * (i - 1) + 3
      j_[count + 10] <- 4 * (i - 1) + 1
      x_[count + 10] <- Q_adj[1, 3]

      #lower edge  u, upper edge  u',
      i_[count + 11] <- 4 * (i - 1) + 1
      j_[count + 11] <- 4 * (i - 1) + 4
      x_[count + 11] <- Q_adj[1, 4]
      i_[count + 12] <- 4 * (i - 1) + 4
      j_[count + 12] <- 4 * (i - 1) + 1
      x_[count + 12] <- Q_adj[1, 4]

      #lower edge  u', upper edge  u,
      i_[count + 13] <- 4 * (i - 1) + 2
      j_[count + 13] <- 4 * (i - 1) + 3
      x_[count + 13] <- Q_adj[2, 3]
      i_[count + 14] <- 4 * (i - 1) + 3
      j_[count + 14] <- 4 * (i - 1) + 2
      x_[count + 14] <- Q_adj[2, 3]


      #lower edge  u', upper edge  u',
      i_[count + 15] <- 4 * (i - 1) + 2
      j_[count + 15] <- 4 * (i - 1) + 4
      x_[count + 15] <- Q_adj[2, 4]
      i_[count + 16] <- 4 * (i - 1) + 4
      j_[count + 16] <- 4 * (i - 1) + 2
      x_[count + 16] <- Q_adj[2, 4]

      count <- count + 16

  }
  if(BC == 1){
    #Vertices with of degree 1
    i.table <- table(c(graph$E))
    index <- as.integer(names(which(i.table == 1)))
    #for this vertices locate position
    lower.edges <- which(graph$E[, 1] %in% index)
    upper.edges <- which(graph$E[, 2] %in% index)
    for (le in lower.edges) {
      ind <- c(4 * (le - 1) + 1, 4 * (le - 1) + 2)

      i_ <- c(i_, ind)
      j_ <- c(j_, ind)
      x_ <- c(x_, 0.5*c(1 / R_00[1, 1], 1 / R_00[2, 2]))
      count <- count + 2
    }
    for (ue in upper.edges) {
      ind <- c(4 * (ue - 1) + 3, 4 * (ue - 1) + 4)
      i_ <- c(i_, ind)
      j_ <- c(j_, ind)
      x_ <- c(x_, 0.5 * c(1 / R_00[1, 1], 1 / R_00[2, 2]))
      count <- count + 2
    }
  }
  #if (build) {
    Q0 <- Matrix::sparseMatrix(i   = i_[1:count],
                              j    = j_[1:count],
                              x    = x_[1:count],
                              dims = c(4 * graph$nE, 4 * graph$nE))


   Q0<-Q0/(r[ii]*kappa^4)

  index.obs2 <- sapply(graph$PtV, function(i){idx_temp <- i == graph$E[,1]
                                                                      idx_temp <- which(idx_temp)
                                                                      return(idx_temp[1])})
    index.obs2 <- (index.obs2-1)*4+1
    index.obs3 <- NULL
    na_obs1 <- is.na(index.obs2)
    if(any(na_obs1)){
          idx_na <- which(na_obs1)
          PtV_NA <- graph$PtV[idx_na]
          index.obs3 <- sapply(PtV_NA, function(i){idx_temp <- i == graph$E[,2]
                                                                      idx_temp <- which(idx_temp)
                                                                      return(idx_temp[1])})
          index.obs2[na_obs1] <- (index.obs3-1)*4 + 3                                                                      
    }

    cbmat2<-conditioningalpha2(alpha=2, graph)   
    nc2<-1:length(cbmat2$S)
    Qtilde<- t(cbmat2$T) %*% Q0 %*% (cbmat2$T)
    W2<-Diagonal(dim(Q0)[1])
    W2 <- W2[,-nc2]
    Qtilde_uu<-t(W2)%*% Qtilde%*%W2
    A1<-(cbmat2$T[,-nc2])
    A1 <- A1[index.obs2, ] # A matrix for alpha=2
    Qtilde1_uu<-bdiag(Qtilde1_uu,Qtilde_uu)
    A <- cbind(A,A1)

  
}

    Qtilde1_uu <- Qtilde1_uu*kappa^(2*alpha)
    sigma=A%*%solve(Qtilde1_uu,t(A)) 
    return(sigma)
}

