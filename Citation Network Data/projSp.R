symmetrize = function(X){
  Z = 0.5*(X + t(X))
  return(Z)
}

projSp = function(X0){
  n = dim(X0)[1]
  tryCatch({temp = eigen(symmetrize(X0)); U = temp$vectors; D = diag(temp$values)}, error = function(x)
      {temp  = svd(symmetrize(X0)); U = temp$u; V = temp$v;S = diag(temp$d); D = diag(diag(t(U) %*% X0 %*% U))})
  idx = as.integer(which(diag(D) >= 0))
  X = as.matrix(U[, idx]) %*% as.matrix(D[idx,idx]) %*% as.matrix(t(U[, idx])); 
  return(X)
}
