##### Description
# This code is created to explore more possibilities of LastFM dataset, which is believed to contain sub-graph that is suitable for CASCORE.

##### Read libraries / functions
#setwd("~/Project/Proj-Social Network/New")
library(igraph)
library(Matrix)
library(pracma)
library(aricode)
source("Various Functions.R")
source("Mtable.R")

############ ############ ############ ############ ############
##### Read LastFM Asian Dataset and Pre-processing ############
############ ############ ############ ############ ############

### attribute matrix 
load("CovariatesMatrix.Rdata")
X = as.matrix(sparseX)
n = dim(X)[1]; p = dim(X)[2]

### adjacency matrix
edgelist = read.csv("lastfm_asia_edges.csv")
edgelist = edgelist + 1
A = matrix(0, n, n)
for (i in 1: dim(edgelist)[1]){
  A[edgelist[i, 1], edgelist[i, 2]] = 1
}
A = A + t(A)

### country label
label = read.csv('lastfm_asia_target.csv')
label = label$target
label = label + 1
unique(label)
table(label)


## delete country with 17 users only
ind5 = which(label != 5);
A = A[ind5, ind5]; X = X[ind5,]; label = label[ind5];
X = X[, colSums(X) > 0]; 

labelnew = label;
label[labelnew > 5] = label[labelnew > 5] - 1;
rm(labelnew); 

n = dim(A)[1]; p = dim(X)[2]
dartist_all = colSums(X); 

############ ############ ############ ############ 
############## End of Reading Data ##############
############ ############ ############ ############ 


## We only consider the comparable countries
sizes = summary(as.factor(label))

# Curret Degrees
d = rowSums(A);
k = length(unique(label))
dave = rep(0, k)
for(i in 1:length(unique(label))){
  dave[i] = mean(d[label == i])
}


print("Small Countries")
class_select = names(sizes)[sizes < 100 & sizes > 50]
source("DataSubsets.R")
print(result$NMI)
print(result$error)
print(dim(Xselect))
result_small = result; 

print("Medium Countries")
class_select = names(sizes)[sizes < 300 & sizes > 100]
source("DataSubsets.R")
print(result$NMI)
print(result$error)
print(dim(Xselect))
result_medium = result; 

print("Large Countries")
class_select = names(sizes)[sizes > 300 & sizes < 1000]
source("DataSubsets.R")
print(result$NMI)
print(result$error)
print(dim(Xselect))
result_large = result;

print("Giant Countries")
class_select = names(sizes)[sizes > 1000]
source("DataSubsets.R")
print(result$NMI)
print(result$error)
print(dim(Xselect))
result_giant = result;

rm(Aselect)
rm(Xselect)
rm(edgelist)
save.image(file = "LastFM_Results.Rdata")
