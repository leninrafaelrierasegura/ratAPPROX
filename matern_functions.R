## ---------------------------------------------------------------------------------------------------------------------------------------
# remotes::install_github("davidbolin/rspde", ref = "devel")
# remotes::install_github("davidbolin/metricgraph", ref = "devel")
library(rSPDE)
library(MetricGraph)
library(grateful)
library(Matrix)

library(ggplot2)
library(reshape2)
library(plotly)


## ---------------------------------------------------------------------------------------------------------------------------------------
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


## ---------------------------------------------------------------------------------------------------------------------------------------
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


## ---------------------------------------------------------------------------------------------------------------------------------------
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


## ---------------------------------------------------------------------------------------------------------------------------------------
#Joint covariance of process and derivative for shifted Matern
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


## ---------------------------------------------------------------------------------------------------------------------------------------
mat1 <- matern.p.joint(s = 0.6, t = 0.5, kappa = 10, p = 0.3, alpha = 2)
mat2 <- matern.p.joint(s = 0.6, t = 0.5, kappa = 10, p = 0.3, alpha = 3)
mat1
mat2

## ---------------------------------------------------------------------------------------------------------------------------------------
exp_precision <- function(loc, kappa, boundary = "free") {
    n <- length(loc)
    l <- diff(loc)
    
    # Precompute reusable terms
    a <- exp(-2 * kappa * l)
    b <- 1 - a
    
    # Initialize matrix Q efficiently
    Q <- Matrix(0, nrow = n, ncol = n)
    
    # Set diagonal elements
    diag(Q)[1:(n-1)] <- 1/2 + a/b
    diag(Q)[2:n] <- diag(Q)[2:n] + 1/2 + a/b
    
    # Set off-diagonal elements
    #  indices <- c(which(row(Q) == col(Q) + 1), which(row(Q) == col(Q) - 1))
    # off_diag_values <- -exp(-kappa * l) /b
    # Q[indices] <- off_diag_values
    
    
    row_indices <- seq_len(n - 1)
    col_indices <- row_indices + 1  # Offsets for off-diagonal elements
    
    # Compute the off-diagonal values
    off_diag_values <- -exp(-kappa * l)/b
    
    # Set the off-diagonal elements in matrix Q
    Q[cbind(row_indices, col_indices)] <- off_diag_values
    Q[cbind(col_indices, row_indices)] <- off_diag_values
    
    # Adjust boundary conditions if required
    if (boundary == "free") {
        Q[1, 1] <- Q[1, 1] + 0.5
        Q[n, n] <- Q[n, n] + 0.5
    }
    
    return(2 * kappa * Q[])
}


## ---------------------------------------------------------------------------------------------------------------------------------------
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

