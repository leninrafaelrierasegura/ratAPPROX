buildKirchooffConditioningMatrixCaseAlphaEqual1 <- function(graph) {
  E <- graph$E
  degrees <- graph$get_degrees()
  numberOfConstraints <- sum(degrees[degrees > 1] - 1)
  K <- matrix(0, nrow = numberOfConstraints, ncol = 2*nrow(E))
  onesCounter <- 0
  
  edgeMatrixFlatten <- c(t(E))
  for (vertex in seq_along(degrees)) {
    degreeOfVertex <- degrees[vertex]
    if (degreeOfVertex < 2) next
    positionOfVertex <-  which(edgeMatrixFlatten == vertex)
    whereToPutOne <- positionOfVertex[1]
    whereToPutMinusOne <- positionOfVertex[-1]
    howManyOnes <- length(whereToPutMinusOne)
    K[cbind(c(1:howManyOnes) + onesCounter, rep(whereToPutOne,howManyOnes))] <- 1
    K[cbind(c(1:howManyOnes) + onesCounter, whereToPutMinusOne)] <- -1
    onesCounter <- onesCounter + howManyOnes
  }
  return(K)
}

buildKirchooffConditioningMatrixCaseAlphaEqualOneNonSparse <- function(graph) {
  E <- graph$E
  degrees <- graph$get_degrees()
  numberOfConstraints <- sum(degrees[degrees > 1] - 1)
  K <- matrix(0, nrow = numberOfConstraints, ncol = 2*nrow(E))
  onesCounter <- 0
  
  edgeMatrixFlatten <- c(t(E))
  for (vertex in seq_along(degrees)) {
    degreeOfVertex <- degrees[vertex]
    if (degreeOfVertex < 2) next
    indicesOfVertex <-  which(edgeMatrixFlatten == vertex)
    colIndicesForOne <- indicesOfVertex[-degreeOfVertex]
    colIndicesForMinusOne <- indicesOfVertex[-1]
    howManyContinuityConditions <- degreeOfVertex - 1
    rowIndicesForBothOneAndMinusOne <- c(1:howManyContinuityConditions) + onesCounter
    K[cbind(rowIndicesForBothOneAndMinusOne, colIndicesForOne)] <- 1
    K[cbind(rowIndicesForBothOneAndMinusOne, colIndicesForMinusOne)] <- -1
    onesCounter <- onesCounter + howManyContinuityConditions
  }
  return(K)
}


gets_graph_circle <- function(n){
  r = 1/(pi)
  theta <- seq(from=-pi,to=pi,length.out = 10000)
  edge <- cbind(1+r+r*cos(theta),r*sin(theta))
  edges = list(edge)
  graph <- metric_graph$new(edges = edges)
  graph$set_manual_edge_lengths(edge_lengths = 2)
  graph$build_mesh(n = n)
  return(graph)
}

gets_graph_double_circle <- function(n){
  r = 1/(pi)
  theta1 <- seq(from=-pi,to=pi,length.out = 1000)
  theta2 <- seq(from=0,to=2*pi,length.out = 1000)
  edge1 <- cbind(1+r+r*cos(theta1),r*sin(theta1))
  edge2 <- cbind(1-r+r*cos(theta2),r*sin(theta2))  
  edges = list(edge1,edge2)
  graph <- metric_graph$new(edges = edges)
  graph$set_manual_edge_lengths(edge_lengths = c(2,2))
  graph$build_mesh(n = n)
  return(graph)
}




graph <- gets_graph_double_circle(10)
graph$plot(direction = TRUE)

graph$E

buildMatrixAWhichMapsUToUv(graph, alpha)
buildKirchooffConditioningMatrixCaseAlphaEqualOne(graph)
  
  
K <- buildKirchooffConditioningMatrixCaseAlphaEqualOne(graph)

MetricGraph:::c_basis2(K)
graph$buildC()
graph$C

























