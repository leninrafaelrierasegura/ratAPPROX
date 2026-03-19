library(rSPDE)
library(rdist)
library(MetricGraph)
library(Matrix)
library(sp)
library(latex2exp)

source("auxiliary_functions.R")

# line1 <- Line(rbind(c(1,0),c(0,0)))
# line2 <- Line(rbind(c(0,1/(1+pi/4)),c(0,0)))
# line3 <- Line(rbind(c(-1/(1+pi/4),1/(1+pi/4)),c(0,1/(1+pi/4))))
# theta <- seq(from=pi,to=3*pi/2,length.out = 20000)
# line4 <- Line(cbind(sin(theta)/(1+pi/4),(1+ cos(theta))/(1+pi/4)))
# Lines = SpatialLines(list(Lines(list(line1),ID="1"),
#                           Lines(list(line2),ID="2"),
#                           Lines(list(line3),ID="3"),
#                           Lines(list(line4),ID="4")))

edge1 <- rbind(c(0,0),c(1,0))
theta <- seq(from=-pi,to=pi,length.out = 10000)
edge2 <- cbind(1+1/pi+cos(theta)/pi,sin(theta)/pi)
Lines = list(edge1, edge2)

graph <- metric_graph$new(lines = Lines)
graph$prune_vertices()
graph$plot()

graph$build_mesh(h =.01)

graph_tmp <- graph$clone()
graph_tmp$plot(mesh=TRUE)
PtE_tmp <- graph_tmp$mesh$VtE
df_temp <- data.frame(y = 0, edge_number = PtE_tmp[,1],
                        distance_on_edge = PtE_tmp[,2])
graph_tmp$add_observations(data = df_temp, normalized = TRUE) 
graph_tmp$observation_to_vertex()
graph$plot(mesh=TRUE)

coord = graph_tmp$coordinates(PtE = graph_tmp$get_PtE())

ord <- match_order(graph$mesh$V,coord)

d<-1
range = 0.5
# nu.v<-seq(from=0.1,to=1.45,length.out=10)
nu.v = 1.3

# nu.v<-seq(from=0.1,to=1.45,by = .01)
#  idx <- (nu.v + 0.5)%%1 > 1e-10
#  nu.v <- nu.v[idx]

alpha.v = nu.v+0.5
m.v <- 4
err.L2<- err.L2_fem <-matrix(0,length(nu.v),length(m.v))
err.sup<- err.sup_fem <-matrix(0,length(nu.v),length(m.v))

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

for(i in 1:length(nu.v)){

  cat(i/length(nu.v),"\n")

  kappa = sqrt(8*nu.v[i])/range
  tau<-sqrt(gamma(nu.v[i]) / (kappa^(2 * nu.v[i]) *
    (4 * pi)^(d / 2) * gamma(nu.v[i] + d / 2)))

############ use this for true solution

Sigma.kl <- matrix(0,nrow = dim(graph$mesh$V)[1],ncol = dim(graph$mesh$V)[1])

for(ii in 0:25000){
  phi <- tadpole.eig(ii,graph)$phi
  Sigma.kl <- Sigma.kl + (1/(kappa^2 + (ii*pi/2)^2))^(alpha.v[i])*phi%*%t(phi)
  if(ii>0 && (ii %% 2)==0){
    psi <- tadpole.eig(ii,graph)$psi
    Sigma.kl <- Sigma.kl + (1/(kappa^2 + (ii*pi/2)^2))^(alpha.v[i])*psi%*%t(psi)
  }
}
  Sigma.kl <- Sigma.kl/tau^2
  graph$compute_fem()

 for(j in 1:length(m.v)){


rspde.order <- m.v[j]

############## use this for reference FEM solution

op <- matern.operators(alpha = alpha.v[i], kappa = kappa, tau = tau,
                         parameterization = "spde",
                         m = rspde.order, graph = graph, type_rational_approximation = "chebfun")   
sigma_fem = op$covariance_mesh()
sigma_fem = sigma_fem[ord,ord]


#########################################

  
  if (alpha.v[i]>0 && alpha.v[i]<1){
           sigma_rat<- shiftedrational_Qalpha1(kappa,tau,nu.v[i],graph_tmp,m=rspde.order,BC = 0) # For alpha b/w 0 and 1
       }else if(alpha.v[i]>1 && alpha.v[i]<2){
           sigma_rat<- shiftedrational_Qalpha2(kappa,tau,nu.v[i],graph_tmp,m=rspde.order,BC = 0) # For alpha b/w 1 and 2
       } else{
         Q2<- shiftedrational_Qalpha3(kappa,tau,nu.v[i],graph_tmp,m=rspde.order,BC = 0) # For alpha b/w 2 and 3
       }

err.L2[i,j]<- sqrt(as.double(t(graph$mesh$weights) %*% (Sigma.kl - sigma_rat)^2 %*% graph$mesh$weights))
err.L2_fem[i,j]<- sqrt(as.double(t(graph$mesh$weights) %*% (sigma_rat - sigma_fem)^2 %*% graph$mesh$weights))

err.sup[i,j] <- max(abs(Sigma.kl-sigma_rat))
err.sup_fem[i,j] <- max(abs(sigma_fem-sigma_rat))

}
}
# err.L2[nu.v==0.5] = 1e-16

png("error_tad.png", 840,480)
par (mfrow = c(1,2))
plot(nu.v,err.L2[,1],type="l",ylim=c(min(err.L2),max(err.L2)),log="y", col=1,xlab=expression(nu), ylab=expression(L[2]~error),
       lwd=2,cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.1)
 lines(nu.v,err.L2[,2],col=2)
 lines(nu.v,err.L2[,3],col=3)
 lines(nu.v,err.L2[,4],col=4)
 lines(nu.v,err.L2[,5],col=5)
  lines(nu.v,err.L2[,6],col=6)
 legend("bottomleft", legend=c("m=1", "m=2", "m=3", "m=4","m=5","m=6"), 
 fill = c("1","2", "3", "4","5","6"), cex=0.8)
  plot(nu.v,err.sup[,1],type="l",ylim=c(min(err.sup),max(err.sup)),log="y", col=1, xlab=expression(nu), ylab=expression(L[infinity]~error), lwd=2,
       cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.1)
 lines(nu.v,err.sup[,2],col=2, lwd=2)
 lines(nu.v,err.sup[,3],col=3,lwd=2)
 lines(nu.v,err.sup[,4],col=4,lwd=2)
  lines(nu.v,err.sup[,5],col=5,lwd=2)
   lines(nu.v,err.sup[,6],col=6,lwd=2)
legend("bottomleft", legend=c("m=1", "m=2", "m=3", "m=4", "m=5", "m=6"),  fill = c("1","2","3","4","5","6"), cex=0.8)
dev.off()

p <- graph$plot_function(Sigma.kl[92,], vertex_size = 2, plotly = TRUE, edge_width = 3)
p<- graph$plot_function(sigma_rat [92,], vertex_size = 2, plotly = TRUE, p = p, line_color = "red", edge_width = 3)
p

error_tadpole_range0.5 = list(err.L2 = err.L2, err.sup = err.sup)

save(error_tadpole_range0.5, file = "error_tadpole_range0.5_N100_alpha0to2.RData")
