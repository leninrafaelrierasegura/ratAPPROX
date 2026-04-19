
# 1 Object `error_with_h0.001.RData`

This object was computed using this code


```
alphaVector <- c(seq(0.51,0.59,by=0.01),
           seq(0.6,0.9,by=0.05), 
           seq(0.91,1.09,by=0.01),
           seq(1.1,1.9,by=0.05), 
           seq(1.91,2.09,by=0.01),
           seq(2.1,2.9,by=0.05), 
           seq(2.91,3,by=0.01)) 
mVector <- 1:7
rhoVector <- c(0.05, 0.1)

sigma <- 1
h <- 0.001
nOverkill <- 1000
type_rational_approxVector <- c("chebfun", "brasil", "chebfunLB")
type_interp <- "spline"


dat <- expand_grid(
  alpha = alphaVector,
  m = mVector,
  rho = rhoVector,
  type_rational_approx = type_rational_approxVector
) |> 
  mutate(
    L2ErrorRational = NA_real_, 
    LinfErrorRational = NA_real_, 
    timeElapsedRational = NA_real_, 
    L2ErrorFEM = NA_real_,
    LinfErrorFEM = NA_real_,
    timeElapsedFEM = NA_real_,
    timeElapsedTrue = NA_real_,
    error_msg = NA_character_,
    status = "not_computed"   
  )


graphFEM <- gets.graph.tadpole(flip_edge = FALSE)
graph <- graphFEM$clone()

graphFEM$build_mesh(h = h)

mesh_locations <- graphFEM$get_mesh_locations() %>% 
  as.data.frame() %>% 
  mutate(dummy = 1) %>%
  rename(edge_number = V1, distance_on_edge = V2)

# add mesh locations as observations
graph$add_observations(
  data = mesh_locations,
  normalized = TRUE)

graph$observation_to_vertex()

unique_params <- unique(dat[, c("alpha", "rho")])
trueSigma_list <- list()
timeTrue_list <- list()
for (j in 1:nrow(unique_params)) {
  alpha_j <- unique_params$alpha[j]
  rho_j <- unique_params$rho[j]
  
  nu_j <- alpha_j - 1/2
  kappa_j <- sqrt(8 * nu_j) / rho_j
  tau_j <- sqrt(gamma(nu_j) / (sigma^2 * kappa_j^(2*nu_j) * (4*pi)^(1/2) * gamma(nu_j + 1/2)))
  
  key <- sprintf("%.6f_%.6f", alpha_j, rho_j)
  
  tryCatch({
  timeTrue_list[[key]] <- system.time({
    trueSigma_list[[key]] <- gets_true_cov_mat(
      graph = graphFEM,
      kappa = kappa_j,
      tau = tau_j,
      alpha = alpha_j,
      n.overkill = nOverkill
    )
  })["elapsed"]
  }, error = function(err){
    warning(paste0("Error occurred while computing true covariance for alpha = ", alpha_j, 
                   " and rho = ", rho_j, 
                   ": ", conditionMessage(err)))
    slackr_msg(text = paste0("Error occurred while computing true covariance for alpha = ", alpha_j, 
                   " and rho = ", rho_j, 
                   ": ", conditionMessage(err)),
               channel = "#research")
    trueSigma_list[[key]] <- NA_real_
    timeTrue_list[[key]] <- NA_real_
  })
  slackr_msg(text = paste0("Finished computing true covariance for alpha = ", alpha_j, 
                          " and rho = ", rho_j, 
                          " , iteration ", j, 
                          " out of ", nrow(unique_params)), 
             channel = "#research")
}

save(trueSigma_list, timeTrue_list, file = "~/Desktop/trueSigma_list.RData")


load("~/Desktop/trueSigma_list.RData")
for (i in 1:nrow(dat)) {
  dat_i <- dat[i, ]
  alpha <- dat_i$alpha
  m <- dat_i$m
  rho <- dat_i$rho
  type_rational_approx <- dat_i$type_rational_approx
  
  nu <- alpha - 1/2
  kappa <- sqrt(8*nu)/rho
  tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))  
  
  res <- tryCatch({
            key <- sprintf("%.6f_%.6f", alpha, rho)

            trueSigma <- trueSigma_list[[key]]
            timeElapsedTrue <- timeTrue_list[[key]]
            
            timeElapsedFEM <- system.time({
            FEMSigma <- covariance_mesh(
              matern.operators(
                alpha = alpha, 
                kappa = kappa, 
                tau = tau,
                m = m, 
                graph = graphFEM,
                type = "covariance",
                type_rational_approximation = type_rational_approx))
            })["elapsed"]
            
            timeElapsedRational <- system.time({
            ratSigma <- rat_covariance(
              graph = graph, 
              kappa = kappa, 
              tau = tau, 
              alpha = alpha, 
              m = m,
              type_rational_approx = type_rational_approx,
              type_interp = type_interp,
              build_cov = TRUE)
            })["elapsed"]
            
            # computing the errors
            LinfErrorRational <- max(abs(trueSigma - ratSigma))
            L2ErrorRational   <- sqrt(as.double(t(graphFEM$mesh$weights)%*%(trueSigma - ratSigma)^2%*%graphFEM$mesh$weights))
            LinfErrorFEM      <- max(abs(trueSigma - FEMSigma))
            L2ErrorFEM        <- sqrt(as.double(t(graphFEM$mesh$weights)%*%(trueSigma - FEMSigma)^2%*%graphFEM$mesh$weights))
            
            list(
                LinfErrorRational =  LinfErrorRational,
                L2ErrorRational   =  L2ErrorRational,
                LinfErrorFEM      =  LinfErrorFEM,
                L2ErrorFEM        =  L2ErrorFEM,
                timeElapsedRational = timeElapsedRational,
                timeElapsedFEM      = timeElapsedFEM,
                timeElapsedTrue     = timeElapsedTrue,
                status = "ok"
                )
            }, error = function(err){
              warning(paste0("Error occurred at iteration rho = ", rho, 
                            ", kappa = ", kappa,
                            ", m = ", m, 
                            ", alpha = ", alpha,
                            ", nu = ", nu, 
                            ", Error:", conditionMessage(err)))
              slackr_msg(text = paste0("Error occurred at iteration rho = ", rho, 
                          ", kappa = ", kappa,
                          ", m = ", m, 
                          ", alpha = ", alpha,
                          ", nu = ", nu, 
                          ", type = ", type_rational_approx,
                          ", Error:", conditionMessage(err)),
                         channel = "#research")
              list(
                LinfErrorRational = NA_real_,
                L2ErrorRational   = NA_real_,
                LinfErrorFEM      = NA_real_,
                L2ErrorFEM        = NA_real_,
                timeElapsedRational = NA_real_,
                timeElapsedFEM      = NA_real_,
                timeElapsedTrue     = NA_real_,
                error_msg = conditionMessage(err),
                status = "error"
                )
            })
  dat[i, names(res)] <- res
  slackr_msg(text = paste0("Current parameters are rho = ", rho, 
                          ", kappa = ", kappa,
                          ", m = ", m, 
                          ", alpha = ", alpha,
                          ", nu = ", nu,
                          ", type = ", type_rational_approx,
                          " , iteration ", i, 
                          " out of ", nrow(dat)), 
             channel = "#research")
}
save(dat, file = here::here("data_files/error_comparison_results.RData"))
```



# 2 Object `error_with_h0.01.RData`

This object was computed using this code


```
alphaVector <- c(seq(0.51,0.59,by=0.01),
           seq(0.6,0.9,by=0.05), 
           seq(0.91,1.09,by=0.01),
           seq(1.1,1.9,by=0.05), 
           seq(1.91,2.09,by=0.01),
           seq(2.1,2.9,by=0.05), 
           seq(2.91,3,by=0.01)) 
mVector <- 1:7
rhoVector <- c(0.001, 0.05, 0.1)

sigma <- 1
h <- 0.01
nOverkill <- 1000
type_rational_approxVector <- c("chebfun", "brasil", "chebfunLB")
type_interp <- "spline"


dat <- expand_grid(
  alpha = alphaVector,
  m = mVector,
  rho = rhoVector,
  type_rational_approx = type_rational_approxVector
) |> 
  mutate(
    L2ErrorRational = NA_real_, 
    LinfErrorRational = NA_real_, 
    timeElapsedRational = NA_real_, 
    L2ErrorFEM = NA_real_,
    LinfErrorFEM = NA_real_,
    timeElapsedFEM = NA_real_,
    timeElapsedTrue = NA_real_,
    error_msg = NA_character_,
    status = "not_computed"   
  )


graphFEM <- gets.graph.tadpole(flip_edge = FALSE)
graph <- graphFEM$clone()

graphFEM$build_mesh(h = h)

mesh_locations <- graphFEM$get_mesh_locations() %>% 
  as.data.frame() %>% 
  mutate(dummy = 1) %>%
  rename(edge_number = V1, distance_on_edge = V2)

# add mesh locations as observations
graph$add_observations(
  data = mesh_locations,
  normalized = TRUE)

graph$observation_to_vertex()

unique_params <- unique(dat[, c("alpha", "rho")])
trueSigma_list <- list()
timeTrue_list <- list()
for (j in 1:nrow(unique_params)) {
  alpha_j <- unique_params$alpha[j]
  rho_j <- unique_params$rho[j]
  
  nu_j <- alpha_j - 1/2
  kappa_j <- sqrt(8 * nu_j) / rho_j
  tau_j <- sqrt(gamma(nu_j) / (sigma^2 * kappa_j^(2*nu_j) * (4*pi)^(1/2) * gamma(nu_j + 1/2)))
  
  key <- sprintf("%.6f_%.6f", alpha_j, rho_j)
  
  tryCatch({
  timeTrue_list[[key]] <- system.time({
    trueSigma_list[[key]] <- gets_true_cov_mat(
      graph = graphFEM,
      kappa = kappa_j,
      tau = tau_j,
      alpha = alpha_j,
      n.overkill = nOverkill
    )
  })["elapsed"]
  }, error = function(err){
    warning(paste0("Error occurred while computing true covariance for alpha = ", alpha_j, 
                   " and rho = ", rho_j, 
                   ": ", conditionMessage(err)))
    slackr_msg(text = paste0("Error occurred while computing true covariance for alpha = ", alpha_j, 
                   " and rho = ", rho_j, 
                   ": ", conditionMessage(err)),
               channel = "#research")
    trueSigma_list[[key]] <- NA_real_
    timeTrue_list[[key]] <- NA_real_
  })
  slackr_msg(text = paste0("Finished computing true covariance for alpha = ", alpha_j, 
                          " and rho = ", rho_j, 
                          " , iteration ", j, 
                          " out of ", nrow(unique_params)), 
             channel = "#research")
  print(j)
}

save(trueSigma_list, timeTrue_list, file = "~/Desktop/trueSigma_list.RData")

load("~/Desktop/trueSigma_list.RData")
for (i in 1:nrow(dat)) {
  dat_i <- dat[i, ]
  alpha <- dat_i$alpha
  m <- dat_i$m
  rho <- dat_i$rho
  type_rational_approx <- dat_i$type_rational_approx
  
  nu <- alpha - 1/2
  kappa <- sqrt(8*nu)/rho
  tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))  
  
  res <- tryCatch({
            key <- sprintf("%.6f_%.6f", alpha, rho)

            trueSigma <- trueSigma_list[[key]]
            timeElapsedTrue <- timeTrue_list[[key]]
            
            timeElapsedFEM <- system.time({
            FEMSigma <- covariance_mesh(
              matern.operators(
                alpha = alpha, 
                kappa = kappa, 
                tau = tau,
                m = m, 
                graph = graphFEM,
                type = "covariance",
                type_rational_approximation = type_rational_approx))
            })["elapsed"]
            
            timeElapsedRational <- system.time({
            ratSigma <- rat_covariance(
              graph = graph, 
              kappa = kappa, 
              tau = tau, 
              alpha = alpha, 
              m = m,
              type_rational_approx = type_rational_approx,
              type_interp = type_interp,
              build_cov = TRUE)
            })["elapsed"]
            
            # computing the errors
            LinfErrorRational <- max(abs(trueSigma - ratSigma))
            L2ErrorRational   <- sqrt(as.double(t(graphFEM$mesh$weights)%*%(trueSigma - ratSigma)^2%*%graphFEM$mesh$weights))
            LinfErrorFEM      <- max(abs(trueSigma - FEMSigma))
            L2ErrorFEM        <- sqrt(as.double(t(graphFEM$mesh$weights)%*%(trueSigma - FEMSigma)^2%*%graphFEM$mesh$weights))
            
            list(
                LinfErrorRational =  LinfErrorRational,
                L2ErrorRational   =  L2ErrorRational,
                LinfErrorFEM      =  LinfErrorFEM,
                L2ErrorFEM        =  L2ErrorFEM,
                timeElapsedRational = timeElapsedRational,
                timeElapsedFEM      = timeElapsedFEM,
                timeElapsedTrue     = timeElapsedTrue,
                status = "ok"
                )
            }, error = function(err){
              warning(paste0("Error occurred at iteration rho = ", rho, 
                            ", kappa = ", kappa,
                            ", m = ", m, 
                            ", alpha = ", alpha,
                            ", nu = ", nu, 
                            ", Error:", conditionMessage(err)))
              slackr_msg(text = paste0("Error occurred at iteration rho = ", rho, 
                          ", kappa = ", kappa,
                          ", m = ", m, 
                          ", alpha = ", alpha,
                          ", nu = ", nu, 
                          ", type = ", type_rational_approx,
                          ", Error:", conditionMessage(err)),
                         channel = "#research")
              list(
                LinfErrorRational = NA_real_,
                L2ErrorRational   = NA_real_,
                LinfErrorFEM      = NA_real_,
                L2ErrorFEM        = NA_real_,
                timeElapsedRational = NA_real_,
                timeElapsedFEM      = NA_real_,
                timeElapsedTrue     = NA_real_,
                error_msg = conditionMessage(err),
                status = "error"
                )
            })
  dat[i, names(res)] <- res
  slackr_msg(text = paste0("Current parameters are rho = ", rho, 
                          ", kappa = ", kappa,
                          ", m = ", m, 
                          ", alpha = ", alpha,
                          ", nu = ", nu,
                          ", type = ", type_rational_approx,
                          " , iteration ", i, 
                          " out of ", nrow(dat)), 
             channel = "#research")
  print(i)
}
save(dat, file = here::here("data_files/error_comparison_results.RData"))
```


# 3 Object `error_with_h0.001_and_rho0.005.RData`

This object was computed using this code


```{r}
alphaVector <- c(seq(0.51,0.59,by=0.01),
           seq(0.6,0.9,by=0.05), 
           seq(0.91,1.09,by=0.01),
           seq(1.1,1.9,by=0.05), 
           seq(1.91,2.09,by=0.01),
           seq(2.1,2.9,by=0.05), 
           seq(2.91,3,by=0.01)) 
mVector <- 1:7
rhoVector <- c(0.005)

sigma <- 1
h <- 0.001
nOverkill <- 1000
type_rational_approxVector <- c("brasil", "chebfunLB")
type_interp <- "spline"


dat <- expand_grid(
  alpha = alphaVector,
  m = mVector,
  rho = rhoVector,
  type_rational_approx = type_rational_approxVector
) |> 
  mutate(
    L2ErrorRational = NA_real_, 
    LinfErrorRational = NA_real_, 
    timeElapsedRational = NA_real_, 
    L2ErrorFEM = NA_real_,
    LinfErrorFEM = NA_real_,
    timeElapsedFEM = NA_real_,
    timeElapsedTrue = NA_real_,
    error_msg = NA_character_,
    status = "not_computed"   
  )


graphFEM <- gets.graph.tadpole(flip_edge = FALSE)
graph <- graphFEM$clone()

graphFEM$build_mesh(h = h)

mesh_locations <- graphFEM$get_mesh_locations() %>% 
  as.data.frame() %>% 
  mutate(dummy = 1) %>%
  rename(edge_number = V1, distance_on_edge = V2)

# add mesh locations as observations
graph$add_observations(
  data = mesh_locations,
  normalized = TRUE)

graph$observation_to_vertex()


unique_params <- unique(dat[, c("alpha", "rho")])
trueSigma_list <- list()
timeTrue_list <- list()
for (j in 1:nrow(unique_params)) {
  alpha_j <- unique_params$alpha[j]
  rho_j <- unique_params$rho[j]
  
  nu_j <- alpha_j - 1/2
  kappa_j <- sqrt(8 * nu_j) / rho_j
  tau_j <- sqrt(gamma(nu_j) / (sigma^2 * kappa_j^(2*nu_j) * (4*pi)^(1/2) * gamma(nu_j + 1/2)))
  
  key <- sprintf("%.6f_%.6f", alpha_j, rho_j)
  
  tryCatch({
  timeTrue_list[[key]] <- system.time({
    trueSigma_list[[key]] <- gets_true_cov_mat(
      graph = graphFEM,
      kappa = kappa_j,
      tau = tau_j,
      alpha = alpha_j,
      n.overkill = nOverkill
    )
  })["elapsed"]
  }, error = function(err){
    warning(paste0("Error occurred while computing true covariance for alpha = ", alpha_j, 
                   " and rho = ", rho_j, 
                   ": ", conditionMessage(err)))
    slackr_msg(text = paste0("Error occurred while computing true covariance for alpha = ", alpha_j, 
                   " and rho = ", rho_j, 
                   ": ", conditionMessage(err)),
               channel = "#research")
    trueSigma_list[[key]] <- NA_real_
    timeTrue_list[[key]] <- NA_real_
  })
  slackr_msg(text = paste0("Finished computing true covariance for alpha = ", alpha_j, 
                          " and rho = ", rho_j, 
                          " , iteration ", j, 
                          " out of ", nrow(unique_params)), 
             channel = "#research")
  print(j)
}

save(trueSigma_list, timeTrue_list, file = "~/Desktop/trueSigma_list.RData")

load("~/Desktop/trueSigma_list.RData")
for (i in 1:nrow(dat)) {
  dat_i <- dat[i, ]
  alpha <- dat_i$alpha
  m <- dat_i$m
  rho <- dat_i$rho
  type_rational_approx <- dat_i$type_rational_approx
  
  nu <- alpha - 1/2
  kappa <- sqrt(8*nu)/rho
  tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu + 1/2)))  
  
  res <- tryCatch({
            key <- sprintf("%.6f_%.6f", alpha, rho)

            trueSigma <- trueSigma_list[[key]]
            timeElapsedTrue <- timeTrue_list[[key]]
            
            timeElapsedFEM <- system.time({
            FEMSigma <- covariance_mesh(
              matern.operators(
                alpha = alpha, 
                kappa = kappa, 
                tau = tau,
                m = m, 
                graph = graphFEM,
                type = "covariance",
                type_rational_approximation = type_rational_approx))
            })["elapsed"]
            
            timeElapsedRational <- system.time({
            ratSigma <- rat_covariance(
              graph = graph, 
              kappa = kappa, 
              tau = tau, 
              alpha = alpha, 
              m = m,
              type_rational_approx = type_rational_approx,
              type_interp = type_interp,
              build_cov = TRUE)
            })["elapsed"]
            
            # computing the errors
            LinfErrorRational <- max(abs(trueSigma - ratSigma))
            L2ErrorRational   <- sqrt(as.double(t(graphFEM$mesh$weights)%*%(trueSigma - ratSigma)^2%*%graphFEM$mesh$weights))
            LinfErrorFEM      <- max(abs(trueSigma - FEMSigma))
            L2ErrorFEM        <- sqrt(as.double(t(graphFEM$mesh$weights)%*%(trueSigma - FEMSigma)^2%*%graphFEM$mesh$weights))
            
            list(
                LinfErrorRational =  LinfErrorRational,
                L2ErrorRational   =  L2ErrorRational,
                LinfErrorFEM      =  LinfErrorFEM,
                L2ErrorFEM        =  L2ErrorFEM,
                timeElapsedRational = timeElapsedRational,
                timeElapsedFEM      = timeElapsedFEM,
                timeElapsedTrue     = timeElapsedTrue,
                status = "ok"
                )
            }, error = function(err){
              warning(paste0("Error occurred at iteration rho = ", rho, 
                            ", kappa = ", kappa,
                            ", m = ", m, 
                            ", alpha = ", alpha,
                            ", nu = ", nu, 
                            ", Error:", conditionMessage(err)))
              slackr_msg(text = paste0("Error occurred at iteration rho = ", rho, 
                          ", kappa = ", kappa,
                          ", m = ", m, 
                          ", alpha = ", alpha,
                          ", nu = ", nu, 
                          ", type = ", type_rational_approx,
                          ", Error:", conditionMessage(err)),
                         channel = "#research")
              list(
                LinfErrorRational = NA_real_,
                L2ErrorRational   = NA_real_,
                LinfErrorFEM      = NA_real_,
                L2ErrorFEM        = NA_real_,
                timeElapsedRational = NA_real_,
                timeElapsedFEM      = NA_real_,
                timeElapsedTrue     = NA_real_,
                error_msg = conditionMessage(err),
                status = "error"
                )
            })
  dat[i, names(res)] <- res
  slackr_msg(text = paste0("Current parameters are rho = ", rho, 
                          ", kappa = ", kappa,
                          ", m = ", m, 
                          ", alpha = ", alpha,
                          ", nu = ", nu,
                          ", type = ", type_rational_approx,
                          " , iteration ", i, 
                          " out of ", nrow(dat)), 
             channel = "#research")
  print(i)
}
save(dat, file = here::here("data_files/error_comparison_results.RData"))
```










