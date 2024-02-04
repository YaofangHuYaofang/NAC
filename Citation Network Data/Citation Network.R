#setwd("~/Project/Proj-Social Network/New")
library(pracma)
library(Matrix)
library(ggplot2)
library(cowplot)
library(gridGraphics)
source("Various Functions.R")

A = as.matrix(read.table("connection.txt", header=F))
# X = as.matrix(read.table("Nfreq.txt", header=T))
load("Nfreq.Rdata")
colnames(A) = NULL
n = dim(A)[1]; p = dim(X)[2]

degree = colSums(A)

newind = which(degree >= 50)
length(newind) #326

require(igraph)
graph1 = graph.adjacency(A[newind, newind], mode = "max")
summary(graph1)
set.seed(2015)
L = layout_nicely(graph1, dim = 2)


# choose the number of communities
#Generate a figure about the eigenvalues of A
Aeig = eig(A)
plot(Aeig[1:6])

######## Step 1. Apply Net-based Method and see the results #############
#In this section, we consider two net-based method: Regularized community detection and SCORE
#for each method, we try K from 3 to 6

#SCORE: we calculate the giant component and the eigenvectors first, then try options of K
#giant component
ix = components(graph.adjacency(A))
componentLabel = ix$membership
GiantLabel = which(componentLabel == which.max(ix$csize))
length(GiantLabel) #2179
Giant = A[GiantLabel, GiantLabel]
gianteig = eigen(Giant)$vectors;

SCOREeig = function(vec, K, itermax = 100, startn = 10){
  n = dim(vec)[1]
  R = vec[, -1]
  R = R[, 1: (K-1)]
  R = R / vec[, 1]
  R[R > sqrt(log(n))] = sqrt(log(n));
  R[R < -1*sqrt(log(n))] = -1*sqrt(log(n));
  
  # apply Kmeans to assign nodes into communities
  result = kmeans(R, K, iter.max = itermax, nstart = startn) #apply kmeans on ratio matrix

  est = as.factor(result$cluster)
  return(est)
}

print("SCORE sizes:")
for(kk in 3:6){
  est_score = SCOREeig(gianteig, kk, startn = 50);
  print(summary(as.factor(est_score)))
}

#Net-based: Regularized method
n = dim(A)[1]
tau = log(n)  
A_tau = A + tau * matrix(1, n, n)/n
s = rowSums(A_tau)
s =  s^(-1/2)
S = diag(s)
Z = S %*% A_tau %*% S
Zeig = eigen(Z)$vectors

Net_based_eig <- function(vec, K, itermax = 100, startn = 10){
  R = vec
  R = R[, 1: K]
  R <- t(apply(R, 1, function(x) x/sqrt(sum(x^2))))
  
  # apply Kmeans to assign nodes into communities
  result = kmeans(R, K, iter.max = itermax, nstart = startn) #apply kmeans on ratio matrix
  
  est = as.factor(result$cluster)
  return(est)
}

print("Net-based method sizes:")
for(kk in 3:6){
  est_net = Net_based_eig(Zeig, kk, startn = 50);
  print(summary(as.factor(est_net)))
}


####### Step 2. Apply Our method and Choose K = 5 #############
source("NAC.R")
est_nac = matrix(0, n, nrow = 4, ncol = n); 
est_nagc = matrix(0, n, nrow = 4, ncol = n); 

for(kk in 3:6){  
  est_nac[kk-2,] = NACeig(nac_eig, kk, startn = 50)
  est_nagc[kk-2,] = NACeig(nagc_eig, kk, startn = 50)
}

#NAC: Match the labels for different K
est_nac_large = est_nac[,newind]; 
for(kk in 4:6){
  tmp = table(est_nac_large[kk-3,], est_nac_large[kk-2, ])
  tmplabel = est_nac[kk - 2,]; 
  recordj = NULL;
  for(ii in 1:nrow(tmp)){
    j = which.max(tmp[ii,]); 
    est_nac[kk - 2, tmplabel == j] = ii; 
    recordj = c(recordj, j); 
  }
  j = setdiff(1:kk, recordj); 
  est_nac[kk-2, tmplabel == j] = kk;
  est_nac_large = est_nac[,newind]; 
}

# Plot 4 figures with K = 3:6 (Figure in Supplementary Material)
#pdf(file = "CitationK.pdf", height = 8, width =8)
setEPS()  
postscript("CitationK_NAC.eps", height = 7, width = 5.5)  
par(mfrow = c(2,2), mai = c(0.1, 0.1, 0.2, 0.1))
for(kk in 3:6){
plot(graph1, layout = L, vertex.color = est_nac[kk - 2, newind], vertex.label = NA, 
     vertex.size = 5, main = paste("NAC, K = ", kk))
}
dev.off()


#match nac and nagc for K = 3 first
est_nac_large = est_nac[,newind]; 
est_nagc_large = est_nagc[,newind]; 

tmp = table(est_nac_large[1,], est_nagc_large[1, ])
tmplabel = est_nac[1,]; 
recordj = NULL;
for(ii in 1:nrow(tmp)){
  j = which.max(tmp[ii,]); 
  est_nac[1, tmplabel == j] = ii; 
  recordj = c(recordj, j); 
}


#NAGC: Match the labels for different K
est_nagc_large = est_nagc[,newind]; 

for(kk in 4:6){
  tmp = table(est_nagc_large[kk-3,], est_nagc_large[kk-2, ])
  tmplabel = est_nagc[kk - 2,]; 
  recordj = NULL;
  for(ii in 1:nrow(tmp)){
    j = which.max(tmp[ii,]);  
    est_nagc[kk - 2, tmplabel == j] = ii; 
    recordj = c(recordj, j); 
  }
  j = setdiff(1:kk, recordj); 
  est_nagc[kk-2, tmplabel == j] = kk;
  est_nagc_large = est_nagc[,newind]; 
}

# Plot 4 figures with K = 3:6 (Figure in Supplementary Material)
#pdf(file = "CitationK_NAGC.pdf", height = 8, width =8)
setEPS()  
postscript("CitationK_NAGC.eps", height = 7, width = 5.5)  
par(mfrow = c(2,2), mai = c(0.1, 0.1, 0.2, 0.1))
for(kk in 3:6){
  plot(graph1, layout = L, vertex.color = est_nagc[kk - 2, newind], vertex.label = NA, 
       vertex.size = 5, main = paste("NAGC, K = ", kk))
}
dev.off()


####### Step 3. Given K = 5, interpret the communities and compare with other methods #############

# Given K = 5, we want to compare with other methods
# consider NAC, NAGC, CASC and SDP
K = 5;
est_casc = CAclustering(A, X, 5, alphan = 5, startn = 50)
est_sdp = admm (A, C = X, K = 5, lambda = 0.5, alpha = 2, rho = 0.5, TT = 100, tol = 100)  
est_net = Net_based_eig(Zeig, K = 5, startn = 50);

#Find the NMI matrix first
nmi_matrix = matrix(0, nrow = 5, ncol = 5)
method_matrix = cbind(est_nac[3,], est_sdp, est_nagc[3,], est_casc, est_net);
for(i in 1:5){
  for(j in 1:i)
  nmi_matrix[i,j] = NMI(method_matrix[,i], method_matrix[,j])
}
colnames(nmi_matrix) <- c("New", "SDP", "New(g)", "CAL", "Net-based")
rownames(nmi_matrix) <- c("New", "SDP", "New(g)", "CAL", "Net-based")

#Draw the heatmap of NMI
nmi_df <- as.data.frame(as.table(nmi_matrix))

heatmap_plot <- ggplot(nmi_df, aes(Var1, Var2, fill = Freq, label = ifelse(Freq != 0, round(Freq, 2), ""))) +
  geom_tile() + 
  geom_text(color = "black", size = 2, vjust = 1) +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() + 
  theme(legend.position = "none", # Change legend position to the right
         legend.key.width = unit(0.5, "cm"), # Adjust the legend key width
         legend.key.height = unit(2, "cm"),
        legend.text = element_text(size = 5), 
         axis.title = element_blank(), 
        axis.text = element_text(size = 5)) + 
  theme()

#Find the NMI matrix first
nmi_matrix_giant = matrix(0, nrow = 5, ncol = 5)
method_matrix = cbind(est_nac[3,], est_sdp, est_nagc[3,], est_casc, est_net);
for(i in 1:5){
  for(j in 1:i)
    nmi_matrix_giant[i,j] = NMI(method_matrix[GiantLabel,i], method_matrix[GiantLabel,j])
}
colnames(nmi_matrix_giant) <- c("New", "SDP", "New(g)", "CAL", "Net-based")
rownames(nmi_matrix_giant) <- c("New", "SDP", "New(g)", "CAL", "Net-based")

#Draw the heatmap of NMI
nmi_df_giant <- as.data.frame(as.table(nmi_matrix_giant))

heatmap_plot_giant <- ggplot(nmi_df_giant, aes(Var1, Var2, fill = Freq, label = ifelse(Freq != 0, round(Freq, 2), ""))) +
  geom_tile() + 
  geom_text(color = "black", size = 2, vjust = 1) +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() + 
  theme(legend.position = "none", # Change legend position to the right
        legend.key.width = unit(0.5, "cm"), # Adjust the legend key width
        legend.key.height = unit(2, "cm"),
        legend.text = element_text(size = 5), 
        axis.title = element_blank(), 
        axis.text = element_text(size = 5)) + 
  theme()

plot_grid(heatmap_plot, heatmap_plot_giant, nrow = 1)
ggsave(file = "NMI.pdf", height = 2, width = 4)  

###The main difference comes from the low degree part 
kk = 5
for(ii in 1:5){
  allsize = summary(as.factor(method_matrix[,ii]))
  print(colnames(nmi_matrix)[ii])
  print(allsize)
  giantsize = summary(as.factor(method_matrix[GiantLabel,ii]))
#  print(giantsize)
  print(allsize - giantsize)
}
#CAL, Net-based and NAGC turns to cluster all the low-degree nodes into the same cluster (the largest cluster)
#NAC and SDP are more uniform
#SDP is very uniform, all clusters have very similar sizes. 


# Given K = 5, we want to interpret the communities
#Find the hot corpus
K = 5;
for(kk in 1:5){
  subclass = which(est_nac[K - 2, ] == kk);
  print(names(sort(apply(X[subclass,], 2, sum), decreasing = TRUE)[1:10]))
  subclass = which(est_nagc[K - 2, ] == kk);
  print(names(sort(apply(X[subclass,], 2, sum), decreasing = TRUE)[1:10]))
}

library(gridBase)
library(grid)

setEPS()  
postscript("Citation.eps", height = 4, width = 6)
par(mfrow=c(1, 2),  mai = c(0.1, 0.1, 0, 0.1))
## the last one is the current plot
plot.new()              ## suggested by @Josh
vps <- baseViewports()
pushViewport(vps$figure) ##   I am in the space of the autocorrelation plot
vp1 <-plotViewport(c(4,1,5,1)) ## create new vp with margins, you play with this values 
print(heatmap_plot_giant,vp = vp1)        ## suggested by @bpatiste

plot(graph1, layout = L, xlim = c(-1, 1), vertex.color = newlabel[est_nac[kk - 2, newind]], vertex.label = NA, 
     vertex.size = 5)
legendnames = c("Variable Selection \n (Regression)\n", "Large-Scale \n Multiple Testing", "Biostatistics", 
                "Bayesian\n", "Variable Selection \n (Semiparametric)\n");
legend(-1, -0.2, legend = legendnames, col = categorical_pal(kk), pch = 16, 
       cex = 0.4, y.intersp = 0.5)

dev.off()
### End of Plot ####

#find the most popular papers in each community
paperList = read.table("paperList.txt", sep=",", stringsAsFactors=F, header=T, 
                       na.strings = "")
for(ii in 1:kk){
  dcomm = degree*0;
  dcomm[est_cascore[kk-2, ] == ii] = degree[est_cascore[kk-2, ] == ii]; 
  indcomm = order(dcomm, decreasing = TRUE);
  print(paperList[indcomm[1:10], 3])
}

#number of nodes in each community
allsize = summary(as.factor(newlabel[est_cascore[kk-2,]]))
print(allsize)

giantsize = summary(as.factor(newlabel[est_cascore[kk-2, GiantLabel]]))
print(giantsize)

print(allsize - giantsize)
## Consider an example isolated paper
nodes = c(2893, 2481);
degree[nodes]

paperList[2893, 3]
#Bayesian pseudo-empirical-likelihood intervals for complex surveys
legendnames[newlabel[est_cascore[kk-2, 2893]]]
#Bayesian

paperList[2481, 3]
#Testing dependence among serially correlated multicategory variables
legendnames[newlabel[est_cascore[kk-2, 2481]]]
#Large-Scale Multiple Testing

write.table(newlabel[est_cascore[kk-2, ]], file = "Result.txt", row.names = paperList[,3])
save.image(file = "Citation.Rdata");
