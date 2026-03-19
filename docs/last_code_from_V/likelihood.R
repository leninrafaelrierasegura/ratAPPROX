library(rSPDE)
library(rdist)
library(MetricGraph)
library(Matrix)
library(sp)
library(latex2exp)
library(tidyr)

source("auxiliary_functions.R")

likelihood_rational <- function(theta, y, X_cov, graph, m, repl, BC, parameterization) {


#  repl = NULL
#  parametrization = "matern"

  repl_vec <- graph$.__enclos_env__$private$data[[".group"]]
  #  y <- graph$.__enclos_env__$private$data[["y"]]


  if(is.null(repl)){
    repl <- unique(repl_vec)
  }

  sigma_e <- exp(theta[1])
  sigma <- exp(theta[2])
  # nu <- exp(theta[4])
  nu <- 2.5 * exp(theta[4])/(1+exp(theta[4]))

  if(parameterization == "matern"){
    kappa = sqrt(8 * nu) / exp(theta[3])
  } else{
    kappa = exp(theta[3])
  }

  tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) * (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))

  
  #build Q
  # graph$observation_to_vertex(mesh_warning=FALSE)

  if (nu>0 && nu<0.5){

  Q_tmp <- precision_alpha1(kappa = kappa, tau = tau,
                      nu = nu, graph = graph, m)
  }else if (nu>0.5 && nu<1.5){

 Q_tmp <- precision_alpha2(kappa = kappa, tau = tau,
                      nu = nu, graph = graph, m)
  }else if (nu>1.5 && nu<2.5){
  Q_tmp <- shiftedrational_Qalpha3(kappa = kappa, tau = tau,
                      nu = nu, graph = graph, m)
  }else{
   print("not implemented")
  }

  Q = Q_tmp$Q

  # Q <- Q_tmp$A %*% Q_tmp$Q %*% t(Q_tmp$A)

  if(is.null(graph$PtV)){
    stop("No observation at the vertices! Run observation_to_vertex().")
  }

  # R <- chol(Q)
  # R <- Matrix::chol(Q)

  R <- Matrix::Cholesky(Q)

  l <- 0

  for(i in repl){
      # A <- Matrix::Diagonal(graph$nV)[graph$PtV, ]
      A = Q_tmp$A
      ind_tmp <- (repl_vec %in% i)
      y_tmp <- y[ind_tmp]

      if(ncol(X_cov) == 0){
        X_cov_tmp <- 0
      } else {
        X_cov_tmp <- X_cov[ind_tmp,,drop=FALSE]
      }

      na_obs <- is.na(y_tmp)

      y_ <- y_tmp[!na_obs]
      n.o <- length(y_)
      Q.p <- Q  + t(A[!na_obs,]) %*% A[!na_obs,]/sigma_e^2
      # R.p <- Matrix::chol(Q.p)
      R.p <- Matrix::Cholesky(Q.p)

    #    l <- l + sum(log(diag(R))) - sum(log(diag(R.p))) - n.o*log(sigma_e)

       l <- l + determinant(R, logarithm = TRUE, sqrt = TRUE)$modulus - determinant(R.p, logarithm = TRUE, sqrt = TRUE)$modulus - n.o * log(sigma_e)

      v <- y_

      if(ncol(X_cov) != 0){
        X_cov_tmp <- X_cov_tmp[!na_obs, ]
        v <- v - X_cov_tmp %*% theta[5:(4+ncol(X_cov))]
      }

      # if(!is.null(X_cov)){
      #     n_cov <- ncol(X_cov)
      #     if(n_cov == 0){
      #       X_cov_repl <- 0
      #     } else{
      #       X_cov_repl <- X_cov[graph$.__enclos_env__$private$data[[".group"]] == repl[i], , drop=FALSE]
      #       X_cov_repl <- X_cov_repl[!na_obs, , drop = FALSE]
      #       v <- v - X_cov_repl %*% theta[5:(4+n_cov)]
      #     }
      # }

        # mu.p <- solve(Q.p,as.vector(t(A[!na_obs,]) %*% v / sigma_e^2))

        mu.p <- solve(R.p, as.vector(t(A[!na_obs,]) %*% v / sigma_e^2), system = "A")

      v <- v - A[!na_obs,]%*%mu.p

      l <- l - 0.5*(t(mu.p) %*% Q %*% mu.p + t(v) %*% v / sigma_e^2) -
        0.5 * n.o * log(2*pi)

  }

  return(-as.double(l))
}

fit_rational = function (formula, theta0, graph, m, method){

graph$observation_to_vertex(mesh_warning=FALSE)

  data <- graph$.__enclos_env__$private$data

  y_term <- stats::terms(formula)[[2]]

  y_graph <- eval(y_term, envir = data, enclos = parent.frame())
  y_graph <- as.numeric(y_graph)

  cov_term <- stats::delete.response(terms(formula))

  X_cov <- stats::model.matrix(cov_term, data)

  if(ncol(X_cov)>0){
      model_matrix <- cbind(y_graph, X_cov)
    } else{
      model_matrix <- y_graph
    }

  # print("Model Matrix:")
  # print(model_matrix)
 likelihood_new <- function(theta, graph, y, X_cov, m, repl, BC, parameterization) {
      l_tmp <- tryCatch(likelihood_rational(theta, graph=graph,m=m, repl=repl, y = y_graph, X_cov = X_cov, BC=BC, parameterization = parameterization),
        error = function(e) {
          return(NULL)
        }
      )
      if (is.null(l_tmp)) {
        return(10^100)
      }
      return(l_tmp)
    }
par <- optim(theta0, likelihood_new, method = method, graph = graph, X_cov = X_cov, 
              y = y_graph, m = m, repl= NULL, BC=0, parameterization = "matern")

# par <- optim(theta0, likelihood_rational, method = method, graph = graph, X_cov = X_cov, 
#               y = y_graph, m = m, repl= NULL, BC=0, parameterization = "matern")

  return (list(par = par, model_matrix = model_matrix))
}


