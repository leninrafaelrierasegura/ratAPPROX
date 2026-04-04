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
