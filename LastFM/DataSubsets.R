### run algorithm on the countries and compare ###
allcomp <- function(A, X, label, class){
  A1 = A; X1 = X;
  
  result = rep(0,6)
  names(result) = c("NAC", "NAGC",
                    "CASC",
                    "SDP",
                    "Net-Based", "Cov-Based")
  timeresult = result;
  NMIresult = result;
  error = result; 
  
  kk = length(unique(label))
  
  # NAC
  t1 = Sys.time()
  d = rowSums(A1);
  est_nac = NAC(A1, X1, kk, alpha = mean(d)/2, beta = 0, startn = 50)
  est_nac[is.na(est_nac)] = which.max(summary(as.factor(est_nac)))
  timeresult["NAC"] = Sys.time() - t1
  result["NAC"] = NMI(est_nac, label)
  if(kk <= 7) error["NAC"] = cluster(table(est_nac, label))$error; 
  
  # NAGC
  t1 = Sys.time()
  d = rowSums(A1);
  est_nagc = NAC(A1, X1, kk, alpha = mean(d)/2, beta = 1, startn = 50)
  est_nagc[is.na(est_nagc)] = which.max(summary(as.factor(est_nac)))
  timeresult["NAGC"] = Sys.time() - t1
  result["NAGC"] = NMI(est_nagc, label)
  if(kk <= 7) error["NAGC"] = cluster(table(est_nagc, label))$error; 
  
    # CAclustering
  t1 = Sys.time()
  est_casc = CAclustering(A1, X1, kk)
  timeresult["CASC"] = Sys.time() - t1
  est_casc[is.na(est_casc)] = which.max(summary(as.factor(est_casc)))
  result["CASC"] = NMI(est_casc, label)
  if(kk <= 7) error["CASC"] = cluster(table(est_casc, label))$error;

  # SDP
  t1 = Sys.time()
  est_sdp = admm(A1, X1, 0.5, kk, 2, 0.5, 100, 5)
  timeresult["SDP"] = Sys.time() - t1
  est_sdp[is.na(est_sdp)] = which.max(summary(as.factor(est_sdp)))
  result["SDP"] = NMI(est_sdp, label)
  if(kk <= 7) error["SDP"] = cluster(table(est_sdp, label))$error;


#  est_casc = 0; est_sdp = 0;
  
  # Net-Based (RSC)
  t1 = Sys.time()
  est_rsc = Net_based(A1, kk, tau = log(n), startn = 50)
  timeresult["Net-Based"] = Sys.time() - t1
  result["Net-Based"] = NMI(est_rsc, label)
  if(kk <= 7) error["Net-Based"] = cluster(table(est_rsc, label))$error;
  
  # Cov-Based (SpectralGem)
  t1 = Sys.time()
  est_spec = Cov_based(X1, kk, startn = 50)
  timeresult["Cov-Based"] = Sys.time() - t1
  result["Cov-Based"] = NMI(est_spec, label)
  if(kk <= 7) error["Cov-Based"] = cluster(table(est_spec, label))$error;
  
  return(
    result = list(
      class = class,
      error = error, 
      NMI = result,
      time = timeresult,
      Sizes = dim(A1)[1],
      est = list(nac = est_nac, casc = est_casc, 
                 sdp =est_sdp, net = est_rsc, cov = est_spec),  
      label = label
    )
  )
}

###Select the corresponding nodes 
ind.select = which(label %in% as.vector(class_select))
Aselect = A[ind.select, ind.select]; Xselect = X[ind.select,]; 
labelselect = label[ind.select];

n = dim(Aselect)[1]; p = dim(Xselect)[2]
kk = length(unique(labelselect))

# Select the regional popular artists##
dartist = colSums(Xselect);
prop = dartist/dartist_all;
prob = 1 - round(min(n/2, 600)/p, 4);
Xselect = Xselect[, prop > quantile(prop, probs = prob, na.rm = TRUE)]

## Remove the users who like no artists
dnode = rowSums(Xselect)
like = which(dnode > 0)
Xselect = Xselect[like, ]
Aselect = Aselect[like, like]
labelselect = labelselect[like]
dselect = rowSums(Aselect)
print(mean(d))

kk = length(unique(labelselect))
dcomm1 = rep(0, kk); dcomm2 = rep(0, kk)
size.com = rep(0, kk);
for(ii in 1:kk){
  commi = which(labelselect == class_select[ii]);
  dcomm1[ii] = mean(d[commi])
  dcomm2[ii] = sum(Aselect[commi, commi])/length(commi)
  size.com[ii] = length(commi)
}
print("Average degree:")
print(dcomm1)
print("Community size and In-community degree:")
print(paste("size:", size.com, "mean:", dcomm2, "\n"))
n = dim(Aselect)[1]; p = dim(Xselect)[2]

#run results
result <- allcomp(Aselect, Xselect, labelselect, class_select)
