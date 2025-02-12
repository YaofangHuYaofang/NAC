library(Matrix)
library(igraph)
library(pracma)
source("Various Functions.R")
source("Mtable.R")


##################### Simulation Studies  ######################
# For all experiments, we consider the Gaussian covariates. 
# We consider p = 20 versus p = 600, which gives small and larger number of covariates. 
# In supplementary materials, we include the result about Bernoulli random variables


# We always consider 3 communities, one is dense and the other two are sparse. 
# Experiment 1. explore what happens when the diagonal and off-diagonal of P changes


# Experiment 2. explore what happens when sigma_x in X changes


# Experiment 3. explore what happens when the mis-specification between A and X changes.
# the result should be based on the sparse nodes only. 




##################### Off-diagnals in P ######################
# Experiment 1. see how the connection pattern will affect the result
# Change a from 0.1 to 0.9, where a is the off-diagonal in P


caseno = 1;
n = 1200;
repetition = 50;


K1 = 4; K2 = 5; prob1 = 0.8; #there are 10% nodes mis-matched
aseq = seq(0.1, 0.9, by = 0.1);
#aseq = c(0.1, 0.4, 0.7, 0.9)
mu = 0.5;
errorall = matrix(0, nrow = length(aseq), ncol = 12) #store the error rate in each repetition
colnames(errorall) = c("New", "New2", "CASC", "ADMM", "Net_based", "Cov_based", "New,small", "New2,small",
                       "CASC,small", "ADMM,small", "Cov_based, small", "giant/all")
errsd = errorall;


#set up the difference between the different communities
for(ajj in 1:length(aseq)){
  a = aseq[ajj]
  P = diag(rep(a, K1)) + (1-a)*ones(K1);
  theta = runif(n); 
  source('Simulation.R')
  filename = paste('Simulation Results/Exp1_a_', a*10, '.Rdata', sep = "");
  save(Dsmall0, Dbig0, errormat, P, Pi, Qsmall, Qbig, Theta, W, caseno, aseq, n, l, K1, K2, mu, a, repetition, theta, file = filename)
  errorall[ajj,] = colMeans(errormat); errsd[ajj,] = apply(errormat, 2, sd)
  print(paste("a = ", a));
}
print(errorall)


save(errorall, errsd, K1, K2, n, repetition, Theta, P, Pi, aseq, prob1, mu, 
     file = "Simulation Results/Exp1.Rdata")


##################### Signal Strength in X ######################
# Experiment 2. see how the signal strength will affect the result
# Let mu changes from 0.1 to 0.5; a stronger mu indicates stronger signal strength


caseno = 2;
n = 1200;
repetition = 50;


K1 = 4; K2 = 5; prob1 = 0.8; #there are 10% nodes mis-matched
a = 0.4; 
museq = seq(0.3, 0.8, by = 0.05);
#museq = c(0.3, 0.5, 0.7);
errorall = matrix(0, nrow = length(museq), ncol = 12) #store the error rate in each repetition
colnames(errorall) = c("New", "New2", "CASC", "ADMM", "Net_based", "Cov_based", "New,small", "New2,small",
                       "CASC,small", "ADMM,small", "Cov_based, small", "giant/all")
errsd = errorall;


#set up the difference between the dense community and sparse communities
for(mujj in 1:length(museq)){
  mu = museq[mujj];
  P = diag(rep(a, K1)) + (1-a)*ones(K1);
  theta = runif(n); 
  source('Simulation.R')
  filename = paste('Simulation Results/Exp2_mu_', mu*10, '.Rdata', sep = "");
  save(Dsmall0, Dbig0, errormat, P, Pi, Qsmall, Qbig, Theta, W, caseno, n, l, K1, K2, mu, a, repetition, theta, file = filename)
  errorall[mujj,] = colMeans(errormat); errsd[mujj,] = apply(errormat, 2, sd)
  print(paste("mu = ", mu));
}


print(errorall)
save(errorall, errsd, K1, K2, n, repetition, Theta, P, Pi, prob1, museq, a, 
     file = "Simulation Results/Exp2.Rdata")


##################### Mis-specification ######################
# Experiment 3. See how the mis-specification rate will impact the result
# change prob1 from 0.3 to 1, where larger prob1 means less mis-specification


caseno = 3;
n = 1200;
repetition = 50;


K1 = 4; K2 = 5; 
a = 0.4; mu = 0.5; #this change comes from the change in the simulation file
#prob1seq = c(0.3, 0.5, 0.7, 1);
prob1seq = seq(0.5, 1, by = 0.05);
errorall = matrix(0, nrow = length(prob1seq), ncol = 12) #store the error rate in each repetition
colnames(errorall) = c("New", "New2", "CASC", "ADMM", "Net_based", "Cov_based", "New,small", "New2,small",
                       "CASC,small", "ADMM,small", "Cov_based, small", "giant/all")
errsd = errorall;


#set up the difference between the dense community and sparse communities
for(probjj in 1:length(prob1seq)){
  prob1 = prob1seq[probjj];
  P = diag(rep(a, K1)) + (1-a)*ones(K1);
  theta = runif(n); 
  source('Simulation.R')
  filename = paste('Simulation Results/Exp3_prob1_', prob1*10, '.Rdata', sep = "");
  save(Dsmall0, Dbig0, errormat, P, Pi, Qsmall, Qbig, Theta, W, caseno, n, l, K1, K2, mu, a, prob1, repetition, theta, file = filename)
  errorall[probjj,] = colMeans(errormat); errsd[probjj,] = apply(errormat, 2, sd)
  print(paste("prob1 = ", prob1));
}
print(errorall)


save(errorall, errsd, K1, K2, n, repetition, Theta, P, Pi, prob1seq, mu, a, 
     file = "Simulation Results/Exp3.Rdata")






####### Plot Section ############
#no color, no points
setEPS()                                             # Set postscript arguments
postscript("Simulation.eps", height = 6.5, width =5)  
#pdf(file = "Simulation.pdf", height = 9, width =6)
par(mfrow = c(3,2))
load("Simulation Results/Exp1.Rdata")
ltyvec = 1:10; pchvec = 1:10; colvec = 1:10; 
plot(aseq, errorall[,1], type = "l", lty = ltyvec[1], col = colvec[1],
     xlab = expression(paste("1 - Between-community Intensity ", alpha)), ylab = "Error Rate", 
     main = "(a) High-dim covariates, varying graph", ylim = c(0, 0.7),
     cex.lab = 1, cex.main = 1, lwd = 1.5)
lines(aseq, errorall[, 2], lty = "99")
for(i in 2:4){
  lines(aseq, errorall[, i+1], lty = ltyvec[i], col = colvec[i])
}
lines(aseq, errorall[, 6], lty = "62121212", col = colvec[5])
#legend("topright", legend = c("CASCORE", "CASC", "SDP", "Net-based", "Cov-based"), lty = 1:5, cex = 0.5)


plot(aseq, errorall[,7], type = "l", lty = ltyvec[1],col = colvec[1],
     xlab = expression(paste("1 - Between-community Intensity ", alpha)), ylab = "Error Rate", main = "(b) Low-dim covariates, varying graph", ylim = c(0, 0.7),
     cex.lab = 1, cex.main = 1, lwd = 1.5)
lines(aseq, errorall[, 8], lty = "99")
for(i in 2:3){
  lines(aseq, errorall[, i+7], lty = ltyvec[i], col = colvec[i])
}
lines(aseq, errorall[, 5], lty = ltyvec[4], col = colvec[4])
lines(aseq, errorall[, 11], lty = "62121212", col = colvec[5])


#legend("topright", legend = c("CASCORE", "CASC", "SDP", "Net-based", "Cov-based"), lty = 1:5, cex = 0.5)




load("Simulation Results/Exp2.Rdata")
colvec = 1:10; pchvec = 1:10; ltyvec = 1:10;
plot(museq, errorall[,1], type = "l", lty = ltyvec[1], col =colvec[1],
     xlab = expression(paste("Covariates Signal Strength ", mu[1])), ylab = "Error Rate", main = "(c) High-dim covariates, varying covariates", ylim = c(0, 0.7),
     cex.lab = 1, cex.main = 1, lwd = 1.5)
lines(museq, errorall[, 2], lty = "99")
for(i in 2:4){
  lines(museq, errorall[, i+1], lty = ltyvec[i], col = colvec[i])
}
lines(museq, errorall[, 6], lty = "62121212", col = colvec[5])




museq = museq + 0.1;
plot(museq, errorall[,7], type = "l", lty = ltyvec[1],col =colvec[1],
     xlab = expression(paste("Covariates Signal Strength ", mu[2])), ylab = "Error Rate", main = "(d) Low-dim covariates, varying covariates", ylim = c(0, 0.7),
     cex.lab = 1, cex.main = 1, lwd = 1.5)
lines(museq, errorall[, 8], lty = "99")
for(i in 2:3){
  lines(museq, errorall[, i+7], lty = ltyvec[i], col =colvec[i])
}
lines(museq, errorall[, 5], lty = ltyvec[4], col =colvec[4])
lines(museq, errorall[, 11], lty = "62121212", col =colvec[5])
#legend("topright", legend = c("CASCORE", "CASC", "SDP", "Net-based", "Cov-based"), lty = 1:5, cex = 0.5)




load("Simulation Results/Exp3.Rdata")
colvec = 1:10; pchvec = 1:10;
plot(prob1seq, errorall[,1], type = "l", lty = ltyvec[1], col =colvec[1],
     xlab = expression(paste("1 - Mis-specification Rate ", gamma)), ylab = "Error Rate", main = "(e) High-dim covariates, mis-specification", ylim = c(0, 0.7),
     cex.lab = 1, cex.main = 1, lwd = 1.5)
lines(prob1seq, errorall[, 2], lty = "99")
for(i in 2:4){
  lines(prob1seq, errorall[, i+1], lty = ltyvec[i], col =colvec[i])
}
lines(prob1seq, errorall[, 6], lty = "62121212", col =colvec[5])


plot(prob1seq, errorall[,7], type = "l", lty = ltyvec[1], col =colvec[1],
     xlab = expression(paste("1 - Mis-specification Rate ", gamma)), ylab = "Error Rate", main = "(f) Low-dim covariates, mis-specification", ylim = c(0, 0.7),
     cex.lab = 1, cex.main = 1, lwd = 1.5)
lines(prob1seq, errorall[, 8], lty = "99")
for(i in 2:3){
  lines(prob1seq, errorall[, i+7], lty = ltyvec[i], col =colvec[i])
}
lines(prob1seq, errorall[, 5], lty = ltyvec[4], col =colvec[4])
lines(prob1seq, errorall[, 11], lty = "62121212", col =colvec[5])


dev.off()
