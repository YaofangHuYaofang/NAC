cl2mat = function(labels){
  n = length(labels);
  k = length(unique(labels));
  mat = matrix(0, n, k)
  for (j in 1: k){
    mat[, j] =  1*(labels == j);
  }
  return(mat)
}