

graph_lme_minimal <- function(formula, graph, BC = 0,
                              model = list(type = "WhittleMatern", alpha = 1, version = 1)) {

  call_graph_lme <- match.call()

  # ── fixed defaults ──────────────────────────────────────────────────────────
  model_type    <- "whittlematern"
  optim_method  <- "L-BFGS-B"
  optim_controls <- list()
  improve_hessian <- FALSE
  which_repl    <- unique(graph$.__enclos_env__$private$data[[".group"]])
  model[["fem"]]         <- FALSE
  model[["directional"]] <- 0
  model[["version"]]     <- 1
  par_names <- c("tau", "kappa")
  fix_vec   <- NULL
  fix_v_val <- NULL

  # ── data ────────────────────────────────────────────────────────────────────
  graph_bkp <- graph$clone()
  data      <- graph$.__enclos_env__$private$data

  all_vars      <- all.vars(formula)
  positional_vars <- c(".edge_number", ".distance_on_edge", ".group", ".coord_x", ".coord_y")
  names_temp    <- c(all_vars, positional_vars)

  filtered_data <- list()
  for (name in names_temp) {
    if (!is.null(data[[name]])) filtered_data[[name]] <- data[[name]]
  }
  data <- filtered_data

  filtered_data2 <- list()
  for (name in names_temp) {
    if (!is.null(graph_bkp$.__enclos_env__$private$data[[name]]))
      filtered_data2[[name]] <- graph_bkp$.__enclos_env__$private$data[[name]]
  }
  graph_bkp$.__enclos_env__$private$data <- filtered_data2

  # ── response & design matrix ───────────────────────────────────────────────
  y_term    <- stats::terms(formula)[[2]]
  cov_term  <- stats::delete.response(terms(formula))
  X_cov     <- stats::model.matrix(cov_term, data)
  y_graph   <- as.numeric(eval(y_term, envir = data, enclos = parent.frame()))

  if (all(dim(X_cov) == c(0, 1))) {
    X_cov <- matrix(1, nrow = length(y_graph))
    colnames(X_cov) <- "(Intercept)"
  }
  n_fixed <- ncol(X_cov)

  # ── starting values ──────────────────────────────────────────────────────────
  start_fixed_values <- graph_starting_values(
    graph        = graph_bkp,
    model        = "alpha1",
    manual_data  = unlist(y_graph),
    log_scale    = TRUE,
    range_par    = FALSE,
    model_options = list(),
    rec_tau      = TRUE
  )
  start_values  <- start_fixed_values$start_values
  fixed_values  <- start_fixed_values$fixed_values

  fix_v <- create_fix_vec_val(fixed_values)
  fix_vec   <- fix_v$fix_vec
  fix_v_val <- fix_v$fix_v_val

  # ── precomputation ───────────────────────────────────────────────────────────
  precomp_data <- precompute_alpha1(
    graph      = graph_bkp,
    manual_y   = y_graph,
    data_name  = NULL,
    X_cov      = X_cov,
    repl       = which_repl
  )

  # ── likelihood ───────────────────────────────────────────────────────────────
  likelihood <- function(theta) {
    n_cov <- ncol(X_cov)
    fix_v_val_full <- c(fix_v_val, rep(NA,  n_cov))
    fix_vec_full   <- c(fix_vec,   rep(FALSE, n_cov))
    new_theta <- fix_v_val_full
    new_theta[!fix_vec_full] <- theta
    return(-likelihood_alpha1_precompute(
      theta           = new_theta,
      graph           = graph_bkp,
      precomputeddata = precomp_data,
      data_name       = NULL,
      manual_y        = y_graph,
      X_cov           = X_cov,
      repl            = which_repl,
      BC              = BC,
      parameterization = "spde"
    ))
  }

  likelihood_new <- function(theta) {
    l_tmp <- tryCatch(likelihood(theta), error = function(e) NULL)
    if (is.null(l_tmp) || is.nan(l_tmp)) return(10^100)
    return(l_tmp)
  }

  # ── optimisation ─────────────────────────────────────────────────────────────
  res <- withCallingHandlers(
    tryCatch(
      optim(start_values, likelihood_new, method = optim_method,
            control = optim_controls, hessian = TRUE),
      error = function(e) NA
    ),
    warning = function(w) invokeRestart("muffleWarning")
  )

  if (is.na(res[1])) stop("Optimisation failed.")

  observed_fisher <- res$hessian

  # ── back-transform coefficients ──────────────────────────────────────────────
  n_random <- length(fix_vec) - 1   # = 2  (tau, kappa)

  if (sum(!fix_vec) > 0) {
    coeff <- c(exp(res$par[1:sum(!fix_vec)]), res$par[-c(1:sum(!fix_vec))])
  } else {
    coeff <- res$par
  }

  # Jacobian correction on Fisher matrix
  tmp_vec    <- c(exp(-res$par[1:sum(!fix_vec)]), rep(1, n_fixed))
  par_change <- if (length(tmp_vec) > 1) diag(tmp_vec) else tmp_vec
  observed_fisher <- par_change %*% observed_fisher %*% par_change

  # ── recover full tmp_coeff (sigma_e, tau, kappa) ─────────────────────────────
  tmp_coeff <- rep(NA, 3)
  tmp_coeff[!fix_vec] <- coeff[1:sum(!fix_vec)]
  tmp_coeff[fix_vec]  <- exp(fix_v_val[!is.na(fix_v_val)])
  tmp_coeff[2]        <- 1 / tmp_coeff[2]   # invert tau

  # Jacobian for tau inversion
  grad_tmp        <- diag(c(c(1, -1 / (tmp_coeff[2]^2), 1)[!fix_vec], rep(1, n_fixed)))
  observed_fisher <- grad_tmp %*% observed_fisher %*% grad_tmp

  # ── Matern re-parameterisation (tau,kappa) → (sigma, range) ──────────────────
  fix_vec_full  <- c(fix_vec, rep(FALSE, n_fixed))
  observed_tmp  <- matrix(NA, nrow = length(fix_vec_full), ncol = length(fix_vec_full))
  observed_tmp[!fix_vec_full, !fix_vec_full] <- observed_fisher

  coeff_tmp        <- tmp_coeff[2:3]
  new_observed_fisher <- observed_tmp[2:3, 2:3]

  change_par <- change_parameterization_graphlme(
    model[["alpha"]] - 0.5,
    coeff_tmp,
    hessian = new_observed_fisher,
    fix_vec[2:3]
  )
  matern_coeff <- list(random_effects = setNames(change_par$coeff, c("sigma", "range")),
                       std_random     = change_par$std_random)

  # ── standard errors ──────────────────────────────────────────────────────────
  inv_fisher  <- tryCatch(solve(observed_fisher),
                          error = function(e) matrix(NA, nrow(observed_fisher), ncol(observed_fisher)))
  std_err     <- sqrt(diag(inv_fisher))
  std_err_tmp <- rep(NA, length(fix_vec_full))
  std_err_tmp[!fix_vec_full] <- std_err
  std_err <- std_err_tmp

  coeff_random <- setNames(tmp_coeff[2:(1 + n_random)], par_names)
  coeff_meas   <- setNames(tmp_coeff[1], "std. dev")
  std_random   <- std_err[2:(1 + n_random)]
  std_meas     <- std_err[1]

  coeff_fixed <- std_fixed <- NULL
  if (n_fixed > 0) {
    coeff_fixed <- res$par[-c(1:sum(!fix_vec))]
    std_fixed   <- std_err[(2 + n_random):length(fix_vec_full)]
  }

  # ── assemble output ───────────────────────────────────────────────────────────
  object <- list(
    coeff            = list(measurement_error = coeff_meas,
                            fixed_effects     = coeff_fixed,
                            random_effects    = coeff_random),
    std_errors       = list(std_meas   = std_meas,
                            std_fixed  = std_fixed,
                            std_random = std_random),
    call             = call_graph_lme,
    formula          = formula,
    loglik           = -res$value,
    BC               = BC,
    latent_model     = model,
    matern_coeff     = matern_coeff,
    optim_method     = optim_method,
    which_repl       = which_repl,
    nobs             = sum(graph$.__enclos_env__$private$data[[".group"]] %in% which_repl),
    niter            = res$counts,
    response_var     = y_term,
    fix_vec          = fix_vec,
    fix_v_val        = fix_v_val,
    start_values     = start_values,
    fixed_values     = fixed_values,
    graph            = graph$clone(),
    df.residual      = sum(graph$.__enclos_env__$private$data[[".group"]] %in% which_repl) -
      (1 + length(coeff_fixed) + length(coeff_random)),
    lik_fun          = likelihood_new,
    mle_par_orig     = res$par,
    precomp_data = precomp_data
  )

  class(object) <- "graph_lme"
  return(object)
}



eval_likelihood_alpha1 <- function(sigma_e, tau, kappa, Y, graph,
                                   precomp_data, BC = 0) {

  # pack into theta as the original function expects
  theta <- c(log(sigma_e), log(1/tau), log(kappa))

  # ── precision matrix Q ───────────────────────────────────────────────────────
  Q.list <- spde_precision(kappa = kappa, tau = tau, alpha = 1,
                           graph = graph, build = FALSE, BC = BC)

  Qp <- Matrix::sparseMatrix(i    = Q.list$i,
                             j    = Q.list$j,
                             x    = Q.list$x,
                             dims = Q.list$dims)

  R       <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE)
  det_R   <- Matrix::determinant(R, sqrt = TRUE)$modulus[1]

  # ── loop setup ───────────────────────────────────────────────────────────────
  obs.edges <- precomp_data$obs.edges
  nV        <- nrow(graph$V)

  repl_vec  <- graph$.__enclos_env__$private$data[[".group"]]
  u_repl    <- unique(repl_vec)

  i_ <- j_ <- x_ <- rep(0, 4 * length(obs.edges))

  loglik        <- 0
  det_R_count   <- NULL
  n.o           <- 0

  # ── replicate loop ───────────────────────────────────────────────────────────
  for (j in seq_along(u_repl)) {

    curr_repl  <- u_repl[j]
    repl_name  <- paste0("repl_", curr_repl)

    loglik <- loglik + det_R
    count  <- 0
    Qpmu   <- numeric(nV)

    # ── edge loop ──────────────────────────────────────────────────────────────
    for (i in seq_along(obs.edges)) {

      e         <- obs.edges[i]
      edge_name <- paste0("edge_", e)

      y_i <- precomp_data$y[[repl_name]][[edge_name]]
      if (is.null(y_i) || length(y_i) == 0) next

      n.o      <- n.o + length(y_i)
      D_matrix <- precomp_data$D_matrix[[repl_name]][[edge_name]]

      # ── edge-level covariance ────────────────────────────────────────────────
      S <- r_1(D_matrix, kappa = kappa, tau = tau)

      E.ind   <- c(1:2)
      Obs.ind <- -E.ind

      Bt      <- solve(S[E.ind, E.ind, drop = FALSE],
                       S[E.ind, Obs.ind, drop = FALSE])
      Sigma_i <- S[Obs.ind, Obs.ind, drop = FALSE] -
        S[Obs.ind, E.ind, drop = FALSE] %*% Bt
      diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2

      R_edge      <- base::chol(Sigma_i)
      Sigma_iB    <- backsolve(R_edge, forwardsolve(t(R_edge), t(Bt)))
      BtSinvB     <- Bt %*% Sigma_iB

      E <- graph$E[e, ]

      if (E[1] == E[2]) {
        Qpmu[E[1]] <- Qpmu[E[1]] + sum(as.vector(t(Sigma_iB) %*% y_i))
        i_[count + 1] <- E[1];  j_[count + 1] <- E[1]
        x_[count + 1] <- sum(BtSinvB)
        count <- count + 1
      } else {
        y_prod    <- as.vector(t(Sigma_iB) %*% y_i)
        Qpmu[E]   <- Qpmu[E] + y_prod
        idx       <- count + 1:4
        i_[idx]   <- c(E[1], E[1], E[2], E[2])
        j_[idx]   <- c(E[1], E[2], E[1], E[2])
        x_[idx]   <- c(BtSinvB[1,1], BtSinvB[1,2],
                       BtSinvB[1,2], BtSinvB[2,2])
        count <- count + 4
      }

      # ── edge log-likelihood contribution ─────────────────────────────────────
      log_det  <- sum(log(diag(R_edge)))
      v_i      <- backsolve(R_edge, forwardsolve(t(R_edge), y_i))
      quad_form <- sum(y_i * v_i)

      loglik <- loglik - 0.5 * quad_form - log_det
    }

    # ── posterior precision Q + B'Σ⁻¹B ───────────────────────────────────────
    if (is.null(det_R_count)) {
      i_all <- c(Q.list$i, i_[1:count])
      j_all <- c(Q.list$j, j_[1:count])
      x_all <- c(Q.list$x, x_[1:count])

      Qp_count <- Matrix::sparseMatrix(i    = i_all,
                                       j    = j_all,
                                       x    = x_all,
                                       dims = Q.list$dims)

      R_count     <- Matrix::Cholesky(Qp_count, LDL = FALSE, perm = TRUE)
      det_R_count <- Matrix::determinant(R_count, sqrt = TRUE)$modulus[1]
    }

    loglik <- loglik - det_R_count

    # ── quadratic form t(mu) Q_post mu ───────────────────────────────────────
    v      <- c(as.matrix(Matrix::solve(R_count,
                                        Matrix::solve(R_count, Qpmu, system = "P"),
                                        system = "L")))
    loglik <- loglik + 0.5 * sum(v^2) - 0.5 * n.o * log(2 * pi)
  }

  return(loglik)
}













graph_lme_minimal_v2 <- function(formula, graph, BC = 0,
                                 model = list(type = "WhittleMatern", alpha = 1, version = 2)) {

  call_graph_lme <- match.call()

  # ── fixed defaults ──────────────────────────────────────────────────────────
  model_type    <- "whittlematern"
  optim_method  <- "L-BFGS-B"
  optim_controls <- list()
  which_repl    <- unique(graph$.__enclos_env__$private$data[[".group"]])
  model[["fem"]]         <- FALSE
  model[["directional"]] <- 0
  model[["version"]]     <- 2
  par_names <- c("tau", "kappa")
  fix_vec   <- NULL
  fix_v_val <- NULL

  # ── key difference vs v1: snap observations to vertices ─────────────────────
  graph_bkp <- graph$clone()
  graph_bkp$observation_to_vertex(mesh_warning = FALSE)
  data <- graph_bkp$.__enclos_env__$private$data

  # ── filter data columns ──────────────────────────────────────────────────────
  all_vars        <- all.vars(formula)
  positional_vars <- c(".edge_number", ".distance_on_edge", ".group", ".coord_x", ".coord_y")
  names_temp      <- c(all_vars, positional_vars)

  filtered_data <- list()
  for (name in names_temp) {
    if (!is.null(data[[name]])) filtered_data[[name]] <- data[[name]]
  }
  data <- filtered_data

  filtered_data2 <- list()
  for (name in names_temp) {
    if (!is.null(graph_bkp$.__enclos_env__$private$data[[name]]))
      filtered_data2[[name]] <- graph_bkp$.__enclos_env__$private$data[[name]]
  }
  graph_bkp$.__enclos_env__$private$data <- filtered_data2

  # ── response & design matrix ─────────────────────────────────────────────────
  y_term   <- stats::terms(formula)[[2]]
  cov_term <- stats::delete.response(terms(formula))
  X_cov    <- stats::model.matrix(cov_term, data)
  y_graph  <- as.numeric(eval(y_term, envir = data, enclos = parent.frame()))

  if (all(dim(X_cov) == c(0, 1))) {
    X_cov <- matrix(1, nrow = length(y_graph))
    colnames(X_cov) <- "(Intercept)"
  }
  n_fixed <- ncol(X_cov)

  # ── starting values ──────────────────────────────────────────────────────────
  start_fixed_values <- graph_starting_values(
    graph         = graph_bkp,
    model         = "alpha1",
    manual_data   = unlist(y_graph),
    log_scale     = TRUE,
    range_par     = FALSE,
    model_options = list(),
    rec_tau       = TRUE
  )
  start_values <- start_fixed_values$start_values
  fixed_values <- start_fixed_values$fixed_values

  fix_v <- create_fix_vec_val(fixed_values)
  fix_vec   <- fix_v$fix_vec
  fix_v_val <- fix_v$fix_v_val

  # ── key difference vs v1: no precomputation, likelihood_alpha1_v2 ────────────
  likelihood <- function(theta) {
    n_cov <- ncol(X_cov)
    fix_v_val_full <- c(fix_v_val, rep(NA,    n_cov))
    fix_vec_full   <- c(fix_vec,   rep(FALSE,  n_cov))
    new_theta <- fix_v_val_full
    new_theta[!fix_vec_full] <- theta
    return(-likelihood_alpha1_v2(
      theta            = new_theta,
      graph            = graph_bkp,   # graph already has obs at vertices
      X_cov            = X_cov,
      y                = y_graph,
      repl             = which_repl,
      BC               = BC,
      parameterization = "spde"
    ))
  }

  likelihood_new <- function(theta) {
    l_tmp <- tryCatch(likelihood(theta), error = function(e) NULL)
    if (is.null(l_tmp) || is.nan(l_tmp)) return(10^100)
    return(l_tmp)
  }

  # ── optimisation ─────────────────────────────────────────────────────────────
  res <- withCallingHandlers(
    tryCatch(
      optim(start_values, likelihood_new, method = optim_method,
            control = optim_controls, hessian = TRUE),
      error = function(e) NA
    ),
    warning = function(w) invokeRestart("muffleWarning")
  )

  if (is.na(res[1])) stop("Optimisation failed.")

  observed_fisher <- res$hessian

  # ── back-transform coefficients ──────────────────────────────────────────────
  n_random <- length(fix_vec) - 1   # = 2  (tau, kappa)

  if (sum(!fix_vec) > 0) {
    coeff <- c(exp(res$par[1:sum(!fix_vec)]), res$par[-c(1:sum(!fix_vec))])
  } else {
    coeff <- res$par
  }

  # Jacobian correction on Fisher matrix
  tmp_vec    <- c(exp(-res$par[1:sum(!fix_vec)]), rep(1, n_fixed))
  par_change <- if (length(tmp_vec) > 1) diag(tmp_vec) else tmp_vec
  observed_fisher <- par_change %*% observed_fisher %*% par_change

  # ── recover full tmp_coeff (sigma_e, tau, kappa) ─────────────────────────────
  tmp_coeff <- rep(NA, 3)
  tmp_coeff[!fix_vec] <- coeff[1:sum(!fix_vec)]
  tmp_coeff[fix_vec]  <- exp(fix_v_val[!is.na(fix_v_val)])
  tmp_coeff[2]        <- 1 / tmp_coeff[2]   # invert tau

  # Jacobian for tau inversion
  grad_tmp        <- diag(c(c(1, -1/(tmp_coeff[2]^2), 1)[!fix_vec], rep(1, n_fixed)))
  observed_fisher <- grad_tmp %*% observed_fisher %*% grad_tmp

  # ── Matern re-parameterisation (tau,kappa) → (sigma, range) ──────────────────
  fix_vec_full  <- c(fix_vec, rep(FALSE, n_fixed))
  observed_tmp  <- matrix(NA, nrow = length(fix_vec_full), ncol = length(fix_vec_full))
  observed_tmp[!fix_vec_full, !fix_vec_full] <- observed_fisher

  coeff_tmp           <- tmp_coeff[2:3]
  new_observed_fisher <- observed_tmp[2:3, 2:3]

  change_par <- change_parameterization_graphlme(
    model[["alpha"]] - 0.5,
    coeff_tmp,
    hessian  = new_observed_fisher,
    fix_vec[2:3]
  )
  matern_coeff <- list(
    random_effects = setNames(change_par$coeff, c("sigma", "range")),
    std_random     = change_par$std_random
  )

  # ── standard errors ──────────────────────────────────────────────────────────
  inv_fisher  <- tryCatch(solve(observed_fisher),
                          error = function(e) matrix(NA, nrow(observed_fisher),
                                                     ncol(observed_fisher)))
  std_err     <- sqrt(diag(inv_fisher))
  std_err_tmp <- rep(NA, length(fix_vec_full))
  std_err_tmp[!fix_vec_full] <- std_err
  std_err <- std_err_tmp

  coeff_random <- setNames(tmp_coeff[2:(1 + n_random)], par_names)
  coeff_meas   <- setNames(tmp_coeff[1], "std. dev")
  std_random   <- std_err[2:(1 + n_random)]
  std_meas     <- std_err[1]

  coeff_fixed <- std_fixed <- NULL
  if (n_fixed > 0) {
    coeff_fixed <- res$par[-c(1:sum(!fix_vec))]
    std_fixed   <- std_err[(2 + n_random):length(fix_vec_full)]
  }

  # ── assemble output ───────────────────────────────────────────────────────────
  object <- list(
    coeff        = list(measurement_error = coeff_meas,
                        fixed_effects     = coeff_fixed,
                        random_effects    = coeff_random),
    std_errors   = list(std_meas   = std_meas,
                        std_fixed  = std_fixed,
                        std_random = std_random),
    call         = call_graph_lme,
    formula      = formula,
    loglik       = -res$value,
    BC           = BC,
    latent_model = model,
    matern_coeff = matern_coeff,
    optim_method = optim_method,
    which_repl   = which_repl,
    nobs         = sum(graph$.__enclos_env__$private$data[[".group"]] %in% which_repl),
    niter        = res$counts,
    response_var = y_term,
    fix_vec      = fix_vec,
    fix_v_val    = fix_v_val,
    start_values = start_values,
    fixed_values = fixed_values,
    graph        = graph$clone(),
    df.residual  = sum(graph$.__enclos_env__$private$data[[".group"]] %in% which_repl) -
      (1 + length(coeff_fixed) + length(coeff_random)),
    lik_fun      = likelihood_new,
    mle_par_orig = res$par
  )

  class(object) <- "graph_lme"
  return(object)
}




eval_likelihood_alpha1_v2 <- function(theta, graph, y, BC, parameterization) {

  kappa = exp(theta[3])
  sigma_e <- exp(theta[1])
  reciprocal_tau <- exp(theta[2])

  #build Q
  Q <- spde_precision(kappa = kappa, tau = 1/reciprocal_tau,
                      alpha = 1, graph = graph, BC=BC)


  R <- Matrix::Cholesky(Q)

  l <- 0

    A <- Matrix::Diagonal(graph$nV)[graph$PtV, ]


    n.o <- length(y)
    Q.p <- Q  + t(A) %*% A/sigma_e^2

    R.p <- Matrix::Cholesky(Q.p)


    l <- l + determinant(R, logarithm = TRUE, sqrt = TRUE)$modulus -
      determinant(R.p, logarithm = TRUE, sqrt = TRUE)$modulus -
      n.o * log(sigma_e)

    v <- y


    mu.p <- solve(R.p, as.vector(t(A) %*% v / sigma_e^2), system = "A")

    v <- v - A%*%mu.p

    l <- l - 0.5*(t(mu.p) %*% Q %*% mu.p + t(v) %*% v / sigma_e^2) -
      0.5 * n.o * log(2*pi)


  return(as.double(l))
}





eva_likelihood_alpha1_v2_no_simplified <- function(theta, graph, X_cov, y, repl, BC, parameterization) {

  repl_vec <- graph$.__enclos_env__$private$data[[".group"]]

  if(is.null(repl)){
    repl <- unique(repl_vec)
  }

  if(parameterization == "matern"){
    kappa = sqrt(8 * 0.5) / exp(theta[3])
  } else{
    kappa = exp(theta[3])
  }

  sigma_e <- exp(theta[1])
  reciprocal_tau <- exp(theta[2])
  #build Q
  Q <- spde_precision(kappa = kappa, tau = 1/reciprocal_tau,
                      alpha = 1, graph = graph, BC=BC)
  if(is.null(graph$PtV)){
    stop("No observation at the vertices! Run observation_to_vertex().")
  }

  # R <- chol(Q)
  # R <- Matrix::chol(Q)

  R <- Matrix::Cholesky(Q)

  l <- 0

  for(i in repl){
    A <- Matrix::Diagonal(graph$nV)[graph$PtV, ]
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

    # l <- l + sum(log(diag(R))) - sum(log(diag(R.p))) - n.o*log(sigma_e)

    l <- l + determinant(R, logarithm = TRUE, sqrt = TRUE)$modulus - determinant(R.p, logarithm = TRUE, sqrt = TRUE)$modulus - n.o * log(sigma_e)

    v <- y_

    if(ncol(X_cov) != 0){
      X_cov_tmp <- X_cov_tmp[!na_obs, , drop=FALSE]
      v <- v - X_cov_tmp %*% theta[4:(3+ncol(X_cov))]
    }

    # mu.p <- solve(Q.p,as.vector(t(A[!na_obs,]) %*% v / sigma_e^2))

    mu.p <- solve(R.p, as.vector(t(A[!na_obs,]) %*% v / sigma_e^2), system = "A")

    v <- v - A[!na_obs,]%*%mu.p

    l <- l - 0.5*(t(mu.p) %*% Q %*% mu.p + t(v) %*% v / sigma_e^2) -
      0.5 * n.o * log(2*pi)

  }

  return(as.double(l))
}




naked_precompute_alpha1 <- function(graph,
                                    data_name = NULL,
                                    manual_y = NULL,
                                    X_cov = NULL,
                                    repl = NULL){

  PtE <- graph$get_PtE() # gets the data location as the matrix [edge_number, distance on edge]
  obs.edges <- unique(PtE[, 1]) # edges where data is

  repl_vec <- graph$.__enclos_env__$private$data[[".group"]] # vector with groups indices

  u_repl <- 1

  y_resp <- graph$.__enclos_env__$private$data[[data_name]]

  # Cache some values used in the loop
  nV <- nrow(graph$V)

  precomputeddata <- list(y = list(),obs.edges=obs.edges,
                          D_matrix = list(),
                          x = list(),
                          u_repl = u_repl)


    curr_repl <- u_repl
    # Use character names for replicate indices
    repl_name <- paste0("repl_", curr_repl)

    # Pre-compute replicate membership only once
    ind_repl_curr <- (repl_vec == curr_repl)
    y_reply <- y_resp[ind_repl_curr]

    precomputeddata$y[[repl_name]] <- list()
    precomputeddata$x[[repl_name]] <- list()
    precomputeddata$D_matrix[[repl_name]] <- list()
    for (i in seq_along(obs.edges)) {
      e <- obs.edges[i]
      # Use character names for edge indices
      edge_name <- paste0("edge_", e)

      # Use pre-computed replicate indices
      obs.id <- PtE[,1] == e
      y_i <- y_reply[obs.id]

      idx_na <- is.na(y_i)

      y_i <- y_i[!idx_na]
      precomputeddata$y[[repl_name]][[edge_name]] <- y_i


      l <- graph$edge_lengths[e]

      PtE_temp <- PtE[obs.id, 2]
      PtE_temp <- PtE_temp[!idx_na]

      # Compute and store time points and distance matrix
      t <- c(0, l, l*PtE_temp)
      precomputeddata$D_matrix[[repl_name]][[edge_name]] <- outer(t, t, `-`)
    }
  return(precomputeddata)
}





naked_eval_likelihood_alpha1 <- function(sigma_e,
                                         tau,
                                         kappa,
                                         Y,
                                         graph,
                                         precomp_data,
                                         BC = 0) {


  # ── precision matrix Q ───────────────────────────────────────────────────────
  Q.list <- spde_precision(kappa = kappa, tau = tau, alpha = 1,
                           graph = graph, build = FALSE, BC = BC)

  # This just builds the precision matrix Q
  Qp <- Matrix::sparseMatrix(i    = Q.list$i,
                             j    = Q.list$j,
                             x    = Q.list$x,
                             dims = Q.list$dims)

  R       <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE)
  det_R   <- Matrix::determinant(R, sqrt = TRUE)$modulus[1]

  # det_R = 0.5*log|Q|

  # ── loop setup ───────────────────────────────────────────────────────────────
  obs.edges <- precomp_data$obs.edges # where data is located
  nV        <- nrow(graph$V) # of the modified graph

  repl_vec  <- graph$.__enclos_env__$private$data[[".group"]]
  u_repl    <- unique(repl_vec)

  i_ <- j_ <- x_ <- rep(0, 4 * length(obs.edges))

  loglik        <- 0
  det_R_count   <- NULL
  n.o           <- 0 # n observations accumulator

  # ── replicate loop ───────────────────────────────────────────────────────────
  for (j in seq_along(u_repl)) {

    curr_repl  <- u_repl[j]
    repl_name  <- paste0("repl_", curr_repl)

    loglik <- loglik + det_R # this adding 0.5*log|Q| n_repl times
    count  <- 0
    Qpmu   <- numeric(nV)

    # ── edge loop ──────────────────────────────────────────────────────────────
    for (i in seq_along(obs.edges)) {

      e         <- obs.edges[i]
      edge_name <- paste0("edge_", e)

      y_i <- precomp_data$y[[repl_name]][[edge_name]]
      if (is.null(y_i) || length(y_i) == 0) next # this almost never could happen

      n.o      <- n.o + length(y_i) # keep count of the data
      D_matrix <- precomp_data$D_matrix[[repl_name]][[edge_name]]
      # recall how D_matrix is computed
      # t <- c(0, l, l*PtE_temp)
      # precomputeddata$D_matrix[[repl_name]][[edge_name]] <- outer(t, t, `-`)

      # ── edge-level covariance ────────────────────────────────────────────────
      S <- r_1(D_matrix, kappa = kappa, tau = tau) # this is just \varrho_M(h) when nu = 1/2, so the exponential covariance
      # > r_1
      # function(D, kappa, tau) {
      #   return((1 / (2 * kappa * tau^2)) * exp(-kappa * abs(D)))
      # }
      E.ind   <- c(1:2) # endpoints of t, that is , 0 and ell
      Obs.ind <- -E.ind

      # S[E.ind, E.ind, drop = FALSE] = Cov between endpoints and endpoints
      # S[E.ind, Obs.ind, drop = FALSE] = Cov between endpoints and observations
      Bt      <- solve(S[E.ind, E.ind, drop = FALSE], S[E.ind, Obs.ind, drop = FALSE]) # this is B^\top

      # S[Obs.ind, Obs.ind, drop = FALSE] = Cov between observations
      # S[Obs.ind, E.ind, drop = FALSE] = Cov between observations and endpoints
      Sigma_i <- S[Obs.ind, Obs.ind, drop = FALSE] - S[Obs.ind, E.ind, drop = FALSE] %*% Bt

      diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2
      # So far I think Sigma_i is Sigma_e

      R_edge      <- base::chol(Sigma_i)
      Sigma_iB    <- backsolve(R_edge, forwardsolve(t(R_edge), t(Bt))) # this computes \Sigma_e^{-1}B
      BtSinvB     <- Bt %*% Sigma_iB # This is B^\top \Sigma_e^{-1}B

      E <- graph$E[e, ] # this gives [underline{e} overline{e}]

      if (E[1] == E[2]) { # for a loop
        Qpmu[E[1]] <- Qpmu[E[1]] + sum(as.vector(t(Sigma_iB) %*% y_i))
        i_[count + 1] <- E[1];  j_[count + 1] <- E[1]
        x_[count + 1] <- sum(BtSinvB)
        count <- count + 1
      } else {
        y_prod    <- as.vector(t(Sigma_iB) %*% y_i) # B^\top\Sigma_e^{-1} y
        Qpmu[E]   <- Qpmu[E] + y_prod
        idx       <- count + 1:4
        i_[idx]   <- c(E[1], E[1], E[2], E[2])
        j_[idx]   <- c(E[1], E[2], E[1], E[2])
        x_[idx]   <- c(BtSinvB[1,1], BtSinvB[1,2],
                       BtSinvB[1,2], BtSinvB[2,2])
        count <- count + 4
      }

      # ── edge log-likelihood contribution ─────────────────────────────────────
      log_det  <- sum(log(diag(R_edge))) # 0.5*log|Sigma_e|
      v_i      <- backsolve(R_edge, forwardsolve(t(R_edge), y_i)) # this computes \Sigma_e^{-1}y
      quad_form <- sum(y_i * v_i) # this computes y^\top\Sigma_e^{-1}y

      loglik <- loglik - 0.5 * quad_form - log_det #here we add -0.5*y^\top\Sigma_e^{-1}y - 0.5*log|Sigma_e|
      # note that we are not adding the sum, just one term in \sum_{e\in\mathcal{E}}
    }
    # finish edge loop

    # ── posterior precision Q + B'Σ⁻¹B ───────────────────────────────────────
    if (is.null(det_R_count)) {
      i_all <- c(Q.list$i, i_[1:count])
      j_all <- c(Q.list$j, j_[1:count])
      x_all <- c(Q.list$x, x_[1:count])

      Qp_count <- Matrix::sparseMatrix(i    = i_all,
                                       j    = j_all,
                                       x    = x_all,
                                       dims = Q.list$dims)

      R_count     <- Matrix::Cholesky(Qp_count, LDL = FALSE, perm = TRUE)
      det_R_count <- Matrix::determinant(R_count, sqrt = TRUE)$modulus[1]
    }

    loglik <- loglik - det_R_count
    # here we add -0.5*log|hat Q_p|, note that it is only computed once but added n_repl times as it does not deppends on data y

    # ── quadratic form t(mu) Q_post mu ───────────────────────────────────────
    v      <- c(as.matrix(Matrix::solve(R_count,
                                        Matrix::solve(R_count, Qpmu, system = "P"),
                                        system = "L")))
    loglik <- loglik + 0.5 * sum(v^2) - 0.5 * n.o * log(2 * pi)
    # here we are adding 0.5*v^\top v = 0.5*\mu^\top Q_p mu
  }

  return(loglik)
}


r_bridge <- function(t1, t2, kappa, alpha, sigma, ell){
  nu <- alpha - 1/2
  r1_t1_0 <- sigma^2*rSPDE:::matern.p.joint(s = t1, t = 0, kappa = kappa, p = 0, alpha = alpha)[1,]
  r1_t1_ell <- sigma^2*rSPDE:::matern.p.joint(s = t1, t = ell, kappa = kappa, p = 0, alpha = alpha)[1,]

  r_0_0 <- sigma^2*rSPDE:::matern.p.joint(s = 0, t = 0, kappa = kappa, p = 0, alpha = alpha)
  r_0_ell <- sigma^2*rSPDE:::matern.p.joint(s = 0, t = ell, kappa = kappa, p = 0, alpha = alpha)
  r_ell_0 <- sigma^2*rSPDE:::matern.p.joint(s = ell, t = 0, kappa = kappa, p = 0, alpha = alpha)
  r_ell_ell <- sigma^2*rSPDE:::matern.p.joint(s = ell, t = ell, kappa = kappa, p = 0, alpha = alpha)

  r1_0_t2 <- sigma^2*rSPDE:::matern.p.joint(s = 0, t = t2, kappa = kappa, p = 0, alpha = alpha)[1,]
  r1_ell_t2 <- sigma^2*rSPDE:::matern.p.joint(s = ell, t = t2, kappa = kappa, p = 0, alpha = alpha)[1,]

  r1_h <- matrix(c(r1_t1_0, r1_t1_ell), nrow = 1, ncol = 2*alpha)
  r1_v <- matrix(c(r1_0_t2, r1_ell_t2), nrow = 2*alpha, ncol = 1)

  cov <- matern.covariance(h = t1 - t2, kappa = kappa, nu = nu, sigma = sigma) - r1_h %*% solve(rbind(cbind(r_0_0, r_0_ell), cbind(r_ell_0, r_ell_ell)), r1_v)
  return(as.numeric(cov))
}



S_mat <-  function(t1, kappa, alpha, sigma, ell){

  nu <- alpha - 1/2

  r1_t1_0 <- sigma^2*rSPDE:::matern.p.joint(s = t1, t = 0, kappa = kappa, p = 0, alpha = alpha)[1,]
  r1_t1_ell <- sigma^2*rSPDE:::matern.p.joint(s = t1, t = ell, kappa = kappa, p = 0, alpha = alpha)[1,]

  r_0_0 <- sigma^2*rSPDE:::matern.p.joint(s = 0, t = 0, kappa = kappa, p = 0, alpha = alpha)
  r_0_ell <- sigma^2*rSPDE:::matern.p.joint(s = 0, t = ell, kappa = kappa, p = 0, alpha = alpha)
  r_ell_0 <- sigma^2*rSPDE:::matern.p.joint(s = ell, t = 0, kappa = kappa, p = 0, alpha = alpha)
  r_ell_ell <- sigma^2*rSPDE:::matern.p.joint(s = ell, t = ell, kappa = kappa, p = 0, alpha = alpha)

  r1_h <- matrix(c(r1_t1_0, r1_t1_ell), nrow = 1, ncol = 2*alpha)
  r1_v <- Matrix::Diagonal(n = 2*alpha)

  cov <-  r1_h %*% solve(rbind(cbind(r_0_0, r_0_ell), cbind(r_ell_0, r_ell_ell)), r1_v)
  return(cov)
}



Simga_e <- function(t_vector, kappa, alpha, sigma, sigma_e, ell){

  n_obs <- length(t_vector)
  Sigma <- matrix(NA, nrow = n_obs, ncol = n_obs)
  for (i in 1:n_obs) {
    for (j in 1:n_obs) {
      Sigma[i,j] = r_bridge(t1 = t_vector[i],
                            t2 = t_vector[j],
                            kappa = kappa,
                            alpha = alpha,
                            sigma = sigma,
                            ell = ell)
    }
  }
  diag(Sigma) <- diag(Sigma) + sigma_e^2
  return(Sigma)
}













graph_lme_alpha2 <- function(formula, graph, BC = 0) {

  # ── 1. Basic checks ──────────────────────────────────────────────────────────
  if (!inherits(graph, "metric_graph")) stop("graph must be a metric_graph")

  model      <- list(type = "WhittleMatern", alpha = 2, fem = FALSE,
                     directional = 0, version = 1)
  which_repl <- unique(graph$.__enclos_env__$private$data[[".group"]])

  # ── 2. Data ──────────────────────────────────────────────────────────────────
  # alpha=2, version=1 → plain data, no observation_to_vertex
  data <- graph$.__enclos_env__$private$data

  all_vars      <- all.vars(formula)
  positional    <- c(".edge_number", ".distance_on_edge", ".group",
                     ".coord_x", ".coord_y")
  names_keep    <- c(all_vars, positional)

  data <- Filter(Negate(is.null),
                 setNames(lapply(names_keep, function(n) data[[n]]), names_keep))

  # ── 3. Design matrix & response ──────────────────────────────────────────────
  y_term   <- stats::terms(formula)[[2]]
  cov_term <- stats::delete.response(stats::terms(formula))
  X_cov    <- stats::model.matrix(cov_term, data)          # y ~ -1 → intercept-free
  y_graph  <- as.numeric(eval(y_term, envir = data, enclos = parent.frame()))

  # ── 4. Graph backup (alpha2 precompute needs a clone) ────────────────────────
  graph_bkp <- graph$clone()

  filtered <- Filter(Negate(is.null),
                     setNames(lapply(names_keep,
                                     function(n) graph_bkp$.__enclos_env__$private$data[[n]]),
                              names_keep))
  graph_bkp$.__enclos_env__$private$data <- filtered

  # ── 5. Starting values ───────────────────────────────────────────────────────
  start_fixed <- graph_starting_values(
    graph        = graph_bkp,
    model        = "alpha2",          # ← alpha = 2
    manual_data  = unlist(y_graph),
    log_scale    = TRUE,
    range_par    = FALSE,
    model_options = list(),
    rec_tau      = TRUE
  )

  start_values  <- start_fixed$start_values
  fixed_values  <- start_fixed$fixed_values

  fix_v     <- create_fix_vec_val(fixed_values)
  fix_vec   <- fix_v$fix_vec
  fix_v_val <- fix_v$fix_v_val

  # ── 6. Precompute (alpha = 2 branch) ─────────────────────────────────────────
  precomp_data <- precompute_alpha2(
    graph    = graph_bkp,
    manual_y = y_graph,
    data_name = NULL,
    X_cov    = X_cov,
    repl     = which_repl
  )

  # ── 7. Likelihood ─────────────────────────────────────────────────────────────
  likelihood <- function(theta) {
    n_cov          <- ncol(X_cov)
    fix_v_val_full <- c(fix_v_val, rep(NA,  n_cov))
    fix_vec_full   <- c(fix_vec,   rep(FALSE, n_cov))
    new_theta                   <- fix_v_val_full
    new_theta[!fix_vec_full]    <- theta
    -likelihood_alpha2_precompute(           # ← alpha2 likelihood
      theta            = new_theta,
      precomputed_data = precomp_data,
      BC               = BC,
      parameterization = "spde"
    )
  }

  likelihood_new <- function(theta) {
    val <- tryCatch(likelihood(theta), error = function(e) NULL)
    if (is.null(val) || is.nan(val)) return(1e100)
    val
  }

  # ── 8. Optimisation (L-BFGS-B, no parallel) ──────────────────────────────────
  res <- withCallingHandlers(
    tryCatch(
      optim(start_values, likelihood_new,
            method  = "L-BFGS-B",
            hessian = TRUE),
      error = function(e) stop("optim failed: ", conditionMessage(e))
    ),
    warning = function(w) invokeRestart("muffleWarning")
  )

  # ── 9. Coefficients & standard errors ────────────────────────────────────────
  n_fixed  <- ncol(X_cov)
  n_random <- length(fix_vec) - 1          # tau, kappa  (2 params)

  # Back-transform from log scale
  coeff <- if (sum(!fix_vec) > 0) {
    c(exp(res$par[1:sum(!fix_vec)]), res$par[-c(1:sum(!fix_vec))])
  } else {
    res$par
  }

  # Delta-method Jacobian for log → original scale
  tmp_vec     <- c(exp(-res$par[1:sum(!fix_vec)]), rep(1, n_fixed))
  par_change  <- if (length(tmp_vec) > 1) diag(tmp_vec) else tmp_vec
  obs_fisher  <- par_change %*% res$hessian %*% par_change

  # Rebuild full coefficient vector (including any fixed params)
  tmp_coeff          <- rep(NA, 3)
  tmp_coeff[!fix_vec] <- coeff[1:sum(!fix_vec)]
  tmp_coeff[fix_vec]  <- exp(fix_v_val[!is.na(fix_v_val)])
  tmp_coeff[2]        <- 1 / tmp_coeff[2]   # kappa stored as 1/kappa internally

  # ── 10. Matern parameterisation (sigma, range) ────────────────────────────────
  fix_vec_full <- c(fix_vec, rep(FALSE, n_fixed))
  obs_tmp      <- matrix(NA, length(fix_vec_full), length(fix_vec_full))
  obs_tmp[!fix_vec_full, !fix_vec_full] <- obs_fisher

  # Jacobian correction for kappa → 1/kappa
  grad_tmp   <- diag(c(c(1, -1/(tmp_coeff[2]^2), 1)[!fix_vec], rep(1, n_fixed)))
  obs_fisher <- grad_tmp %*% obs_fisher %*% grad_tmp

  # matern_coeff <- change_parameterization_graphlme(
  #   nu      = model$alpha - 0.5,           # nu = 1.5 for alpha = 2
  #   coeff   = tmp_coeff[2:3],
  #   hessian = obs_tmp[2:3, 2:3],
  #   fix_vec = fix_vec[2:3]
  # )

  # Standard errors
  inv_fisher            <- tryCatch(solve(obs_fisher),
                                    error = function(e) matrix(NA, nrow(obs_fisher), ncol(obs_fisher)))
  std_err               <- rep(NA, length(fix_vec_full))
  std_err[!fix_vec_full] <- sqrt(diag(inv_fisher))

  coeff_random <- tmp_coeff[2:(1 + n_random)]
  names(coeff_random) <- c("tau", "kappa")
  coeff_meas   <- tmp_coeff[1]
  names(coeff_meas) <- "std. dev"

  # ── 11. Assemble output ───────────────────────────────────────────────────────
  object <- list(
    coeff = list(
      measurement_error = coeff_meas,
      fixed_effects     = NULL,
      random_effects    = coeff_random
    ),
    std_errors = list(
      std_meas   = std_err[1],
      std_fixed  = NULL,
      std_random = std_err[2:(1 + n_random)]
    ),
    # matern_coeff  = list(
    #   random_effects = setNames(matern_coeff$coeff,  c("sigma", "range")),
    #   std_random     = matern_coeff$std_random
    # ),
    loglik        = -res$value,
    formula       = formula,
    latent_model  = model,
    BC            = BC,
    which_repl    = which_repl,
    nobs          = length(y_graph),
    graph         = graph$clone(),
    optim_method  = "L-BFGS-B",
    niter         = res$counts,
    fix_vec       = fix_vec,
    fix_v_val     = fix_v_val,
    start_values  = start_values,
    fixed_values  = fixed_values,
    lik_fun       = likelihood_new,
    mle_par_orig  = res$par,
    precomp_data  = precomp_data
  )

  class(object) <- "graph_lme"
  return(object)
}


eval_likelihood_alpha2 <- function(sigma_e, tau, kappa, precomp, BC = 0) {

  # ── 1. Build theta vector (log scale, spde parameterization) ─────────────────
  theta <- c(
    log(sigma_e),   # theta[1]: log(sigma_e)
    log(1 / tau),   # theta[2]: log(reciprocal_tau)
    log(kappa)      # theta[3]: log(kappa)
  )

  # ── 2. Evaluate ───────────────────────────────────────────────────────────────
  likelihood_alpha2_precompute(
    theta            = theta,
    precomputed_data = precomp,
    BC               = BC,
    parameterization = "spde"
  )
}
