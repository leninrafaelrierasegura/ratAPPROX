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























