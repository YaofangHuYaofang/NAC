Ac = function(X, n){
  z = c(2 * X %*% matrix(1, n, 1), sum(diag(X)))
  return(z)
}

Acs = function(z, n){
  #mu = z[, 1][1:n]
  mu = z[1:n]
  U = matrix(rep(mu, n), byrow = TRUE, nrow = n)
  V = matrix(rep(mu, n), byrow = FALSE, nrow = n)
  #Z = U + V + z[, 1][n + 1]*diag(n)
  Z = U + V + z[n + 1]*diag(n)
  return(Z)
}

Pinv = function(z, n){
  library(reshape2)
  #mu = z[, 1][1:n]
  mu = z[1:n]
  #nu = melt(z)$value[(n+1): length(melt(z)$value)]
  nu = z[(n+1): length(z)]
  x1 = 1/2/n * (diag(n) - (n-2)/n/(2*n-2)*matrix(1, n, n)) %*% mu 
  x2 = 1/n/(2-2*n)*matrix(1, n, 1) %*% nu 
  x3 = - 1/n/(2*n-2)*matrix(1, 1, n) %*% mu 
  x4 = nu/(n - 1)
  X = c(x1 + x2, x3 + x4)
  #X = rbind(matrix(rep(x1, dim(x2)[2]), nrow = 2) + x2, rep(x3, length(x4))+ x4)
  return(X)
}

projAXb = function(X0, k, n){
  # k is the trace of X
  b = c(2*rep(1, n), k)
  X = X0 - Acs(Pinv(Ac(X0, n) - b, n), n)
  return(X)
}
