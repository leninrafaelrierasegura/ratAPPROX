library(rSPDE)
library(rdist)
library(MetricGraph)
library(Matrix)
library(sp)
library(latex2exp)

source("auxiliary_functions.R")

d<-1
range = 0.2
sigma<-1

nu.v<-seq(from=0.1,to=1.45,length.out = 10)

# nu.v<-seq(from=0.1,to=1.45,by = .05)
# idx <- (nu.v + 0.5)%%1 > 1e-10
# nu.v <- nu.v[idx]

alpha.v <- nu.v+0.5
m.v <- 1:6

err.L2<- matrix(0,length(nu.v),length(m.v))
err.sup<- matrix(0,length(nu.v),length(m.v))

line1 <- Line(rbind(c(0,0),c(0.5,0)))
line2 <- Line(rbind(c(1,0),c(0.5,0)))
Lines = SpatialLines(list(Lines(list(line1),ID="1"),
                          Lines(list(line2),ID="2")))
                         

graph <- metric_graph$new(lines = Lines)

graph$build_mesh(h =.01)

graph_tmp <- graph$clone()
graph_tmp$plot(mesh=TRUE)
PtE_tmp <- graph_tmp$mesh$VtE 
df_temp <- data.frame(y = 0, edge_number = PtE_tmp[,1],
                        distance_on_edge = PtE_tmp[,2])
graph_tmp$add_observations(data = df_temp, normalized = TRUE) 
graph_tmp$observation_to_vertex() 

graph_tmp$plot(direction=TRUE)

coord = graph_tmp$coordinates(PtE = graph_tmp$get_PtE())
ord = match_order(graph$mesh$V, coord)


for(i in 1:length(nu.v)){
   cat(i/length(nu.v),"\n")

   nu = nu.v[i]
   alpha = alpha.v[i]
  kappa = sqrt(8*nu)/range
  tau <- sqrt(gamma(nu) / (kappa^(2 * nu) *
 (4 * pi)^(d / 2) * gamma(nu + d / 2)))

  
c.true_mat<-folded_matern_cov(graph$mesh$V[,1],kappa,nu,sigma,N=10,L=1)
c.true_mat = c.true_mat[ord,ord]

  for(j in 1:length(m.v)){

  rspde.order = m.v[j]

if (alpha>0 && alpha<1){
        sigma_rat<- shiftedrational_Qalpha1(kappa,tau,nu,graph_tmp,m=rspde.order,BC = 0) # For alpha b/w 0 and 1
       }else if(alpha>1 && alpha<2){
        sigma_rat<- shiftedrational_Qalpha2(kappa,tau,nu,graph_tmp,m=rspde.order,BC = 0) # For alpha b/w 1 and 2
       } else{
         Q2<- shiftedrational_Qalpha3(kappa,tau,nu,graph_tmp,m=rspde.order,BC = 0) # For alpha b/w 2 and 3
       }


graph$compute_fem() # needed to run graph$mesh$weights
err.L2[i,j]<- sqrt(as.double(t(graph$mesh$weights) %*% (c.true_mat - sigma_rat)^2 %*% graph$mesh$weights))
err.sup[i,j] <- max(abs(sigma_rat-c.true_mat))

}
}


par (mfrow = c(1,2))
plot(nu.v,err.L2[,1],type="l",ylim=c(min(err.L2),max(err.L2)),log="y", col=1, xlab=expression(nu), ylab=expression(L[2]~error), lwd=2,
       cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.1)
 lines(nu.v,err.L2[,2],col=2, lwd=2)
 lines(nu.v,err.L2[,3],col=3,lwd=2)
 lines(nu.v,err.L2[,4],col=4,lwd=2)
  lines(nu.v,err.L2[,5],col=5,lwd=2)
  lines(nu.v,err.L2[,6],col=6,lwd=2)
legend("bottomleft", legend=c("m=1", "m=2", "m=3", "m=4", "m=5", "m=6"),  fill = c("1","2","3","4","5","6"), cex=1.2)

plot(nu.v,err.sup[,1],type="l",ylim=c(min(err.sup),max(err.sup)),log="y", col=1, xlab=expression(nu), ylab=expression(L[infinity]~error), lwd=2,
       cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.1)
 lines(nu.v,err.sup[,2],col=2, lwd=2)
 lines(nu.v,err.sup[,3],col=3,lwd=2)
 lines(nu.v,err.sup[,4],col=4,lwd=2)
  lines(nu.v,err.sup[,5],col=5,lwd=2)
  lines(nu.v,err.sup[,6],col=6,lwd=2)
legend("bottomleft", legend=c("m=1", "m=2", "m=3", "m=4", "m=5", "m=6"),  fill = c("1","2","3","4","5","6"), cex=1.2)

