library(Matrix)
library(igraph)
library(pracma)
source("Various Functions.R")
source("Mtable.R")

############ Parameter Set  ############

# Simulate the Network
#n = 1000; # number of nodes
#K1 = 4; # number of communities revealed by adjacency matrix


# Set the labels
set.seed(2019)
l = sample(1:K1, n, replace=TRUE); # node labels

Pi = matrix(0, n, K1) # label matrix
for (k in 1:K1){
  Pi[l == k, k] = 1
}

# use theta from "ExpSets.R" and the generated labels for the dense community and sparse communities
theta[l < 3] = theta[l < 3]*0.3 + 0.3; #dense community
theta[l >= 3] = theta[l >= 3]*0.03 + 0.03; #sparse communities
Theta = diag(theta); # node degree heterogeneity

# Get the expected adjacency matrix
Omega = Theta %*% Pi %*% P %*% t(Pi) %*% Theta;

############ End Parameter Set  ############
errormat = matrix(0, nrow = repetition, ncol = 12) #store the error rate in each repitition
colnames(errormat) = c("New", "New2", "CASC", "ADMM", "Net_based", "Cov_based", "New,small", "New2,small",
                       "CASC,small", "ADMM,small", "Cov_based, small", "giant/all")

#Calculation of the error 
bigp = 600; smallp = 20;

for (ii in 1: repetition){
  set.seed(ii+9999)
  
  #Set up E[X] matrix 
    # Suppose 20% of the covariates are signals
    # to make sure $|mu|/sigma_x$ is the same for both p = 20 and p = 600, #signals should be the same
  
#    mu = 0.2;
  Qsmall = 0.1*matrix(sign(runif(smallp*K2) - 0.5), nrow = smallp); 
  for(smalli in 1:K2){
    Qsmall[(smalli-1)*(smallp/K2)+(1:(smallp/K2)), smalli] = mu + 0.1; #remark. has a change here
  }
#    signal_big = rbinom(bigp, 1, 2*sqrt(bigp)/bigp);
  signal_big = rbinom(bigp, 1, 0.1);
  Qbig = mu*matrix(rep(signal_big,K2),nrow=bigp,byrow=FALSE)*matrix(sign(runif(bigp*K2) - 0.5), nrow = bigp);
    W = matrix(0, nrow = n, ncol = K2); 
    
    # use prob1 from "ExpSets.R"
    for(jj in 1:n) {
      pp = rep(1/(K2-1), K2); pp[l[jj]] = 0;
      if(runif(1) <= prob1) {W[jj, 1:K1] = Pi[jj, ];}
      else 
        W[jj, sample(K2, 1, prob = pp)] = 1;
    }
    W = t(W)
    
    Dsmall0 = Qsmall %*% W # expectation of covariate matrix 
    Dbig0 = Qbig %*% W # expectation of covariate matrix 
    
  
  
  Dsmall = matrix(0, n, smallp)
  Dbig = matrix(0, n, bigp)
  #the parameter in multivariate distribution 
  for (i in 1:n){
      Dsmall[i,] = rnorm(smallp, mean = Dsmall0[,i], sd = 1);
      Dbig[i,] = rnorm(bigp, mean = Dbig0[,i], sd = 1);
  }

  A = matrix(runif(n*n, 0, 1), nrow = n);
  A = Omega - A;
  A = 1*(A >= 0)
  diag(A) = 0
  A <- as.matrix(forceSymmetric(A))
  is.igraph(A) # [1] FALSE
  ix = components(graph.adjacency(A))
  componentLabel = ix$membership
  giantLabel = which(componentLabel == which.max(ix$csize))
  Giant = A[giantLabel, giantLabel]
  errormat[ii, "giant/all"] = length(giantLabel)/n;


  ######### This is the section for large p ###################
  # New approach: NAC
  degree = colSums(A); alpha = mean(degree)/2;
  est_nac = NAC(A, Dbig, K1, alpha = alpha, beta = 0)
  comp = table(est_nac, l);
  err_nac = cluster(comp)$error;
  errormat[ii, "New"] = err_nac;

  # New approach: NAC-generalized
  degree = colSums(A); alpha = mean(degree)/2;
  est_nagc = NAC(A, Dbig, K1, alpha = alpha, beta = 1)
  comp = table(est_nagc, l);
  err_nagc = cluster(comp)$error;
  errormat[ii, "New2"] = err_nagc;
  
  # CA-clustering
  est_casc = CAclustering(A, Dbig, K1)
  comp = table(est_casc, l);
  err_casc = cluster(comp)$error;
  errormat[ii, "CASC"] = err_casc;

  # admm
  est_admm = admm(A, C = Dbig, K = K1, lambda = 0.5, alpha = 2, rho = 0.5, TT = 100, tol = 100)
  comp = table(est_admm, l);
  errormat[ii, "ADMM"] = cluster(comp)$error;
  
  # Network-based method RSC
  est_net <- Net_based(Adj = A, K = K1)
  comp = table(est_net, l);
  errormat[ii, "Net_based"] = cluster(comp)$error;
  
  # Covariant-based method: SpectralGem
  est_cov <- Cov_based(Covariate = Dbig, K = K1)
  comp = table(est_cov, l);
  errormat[ii, "Cov_based"] = cluster(comp)$error;

  ######### This is the section for small p ###################
  # New approach: NAC
  degree = colSums(A); alpha = mean(degree)/2;
  est_nac = NAC(A, Dsmall, K1, alpha = alpha, beta = 0)
  comp = table(est_nac, l);
  err_nac = cluster(comp)$error;
  errormat[ii, "New,small"] = err_nac;
  
  # New approach: NAC-generalized
  degree = colSums(A); alpha = mean(degree)/2;
  est_nagc = NAC(A, Dsmall, K1, alpha = alpha, beta = 1)
  comp = table(est_nagc, l);
  err_nagc = cluster(comp)$error;
  errormat[ii, "New2,small"] = err_nagc;
  
  # CA-clustering
  est_casc = CAclustering(A, Dsmall, K1)
  comp = table(est_casc, l);
  err_casc = cluster(comp)$error;
  errormat[ii, "CASC,small"] = err_casc;

  # admm
  est_admm = admm(A, C = Dsmall, K = K1, lambda = 0.5, alpha = 2, rho = 0.5, TT = 100, tol = 100)
  comp = table(est_admm, l);
  errormat[ii, "ADMM,small"] = cluster(comp)$error;
  
  # Network-based method
  #est_net <- Net_based(Adj = A, K = K1)
  #errormat[ii, "Net_based, small"] = err_net;
  
  # Covariant-based method
  est_cov <- Cov_based(Covariate = Dsmall, K = K1)
  comp = table(est_cov, l);
  errormat[ii, "Cov_based, small"] = cluster(comp)$error;
  
  print(ii)
#  print(errormat[ii,])
}
print(colMeans(errormat))
