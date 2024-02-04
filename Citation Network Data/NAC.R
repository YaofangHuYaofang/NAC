# library(pracma)
# library(Matrix)
# 

  #Regularity check
  netrow = rowSums(A);
  covrow = rowSums(abs(X));
  ind_reg = which(netrow != 0 | covrow != 0)
  length(ind_reg) #3232 all nodes
  A = A[ind_reg, ind_reg];
  X = X[ind_reg, ];
  
  ##Algorithm
  lambda = log(n)/(degree + log(n));
  alpha = mean(degree)/2
  
  Newmat = (A + alpha*diag(lambda))%*%X;
  zz1 = Newmat%*%t(Newmat);
  nac_eig = eigen(zz1)$vectors

  covmean = colMeans(X, na.rm = TRUE);
  beta = sum(covmean^2);
  print(beta)
  AA = A%*%A;
  zz2 = zz1 +  beta*n*AA;
  nagc_eig = eigen(zz2)$vectors
  

NACeig <- function(vec, K, itermax = 100, startn = 10){

    vecn = vec[,1:K]/apply(vec[,1:K], 1, Norm);
    result = kmeans(vecn, K, iter.max = itermax, nstart = startn);
    if (result$ifault==4) { result = kmeans(X, K,  iter.max = itermax, nstart = startn, algorithm="Lloyd"); }
    est = as.factor(result$cluster);

    return(est)
  }




