cluster <- function(T){
n = sum(T); k = dim(T)[1]; truek = dim(T)[2];
if(k > truek) stop("Error: Number of estimated clusters is larger than truth")
match = 0;

if(sum(T) == 0){
    return(list(error = 0, flagrecord = 1:k))
}   

if(k == 1) {match = max(T);  flagrecord = which.max(T);}  
if(k == 2 & truek == 2) {
	m = c(sum(diag(T)),  n - sum(diag(T)) );
	match = max(m); flagrecord = c(which.max(m), 3 - which.max(m));
}
	
if(k == 2 & truek == 3){
    for(i1 in 1:3){
        m = 0; flag = numeric(3);
        flag(i1) = 1; 
        m = m + T(1, i1);
        m = m + max(T[2, flag == 0]);
        if(match < m)
            {match = m; flagrecord = c(i1, which(T[2, ] == max(T[2, flag == 0])) );}
    }
}
if(k == 3 & truek == 3){
    for(i1 in 1:3){
        m = 0; flag = numeric(3);
        flag[i1] = 1; 
        m = m + T[1, i1];
        for(i2 in 1:3){
            if(flag[i2]==0) 
                {flag[i2] = 1; m = m + T[2, i2];}
            else
                next;
            m = m + T[3, flag == 0];
            if(m > match)
                match = m; flagrecord = c(i1, i2, which(flag==0));
            m = T[1, i1];
            flag[i2] = 0;
        }
    }
}

if(truek == 4){
     for(i1 in 1:4){
        flag = numeric(4);
        flag[i1] = 1;
        T2 = T[2:k, flag==0];
        result = cluster(T2); 
        m = (1-result$error)*sum(T2) + T[1, i1]; flag = result$flagrecord;
        flag[flag >= i1] = flag[flag >= i1] + 1;
        if(match < m) {match = m; flagrecord = c(i1, flag);}
     }
}

if(truek == 5){
     for(i1 in 1:5){
        flag = numeric(5);
        flag[i1] = 1;
        T2 = T[2:k, flag==0];
        result = cluster(T2); 
        m = (1-result$error)*sum(T2) + T[1, i1]; flag = result$flagrecord;
        flag[flag >= i1] = flag[flag >= i1] + 1;
        if(match < m) {match = m; flagrecord = c(i1, flag);}
     }
}

if(truek == 6){
    for(i1 in 1:6){
        flag = numeric(6);
        flag[i1] = 1;
        T2 = T[2:k, flag==0];
        result = cluster(T2); 
        m = (1-result$error)*sum(T2) + T[1, i1]; flag = result$flagrecord;
        flag[flag >= i1] = flag[flag >= i1] + 1;
        if(match < m) {match = m; flagrecord = c(i1, flag);}
    }
}

if(truek == 7){
    for(i1 in 1:7){
        flag = numeric(7);
        flag[i1] = 1;
        T2 = T[2:k, flag==0];
        result = cluster(T2); 
        m = (1-result$error)*sum(T2) + T[1, i1]; flag = result$flagrecord;
        flag[flag >= i1] = flag[flag >= i1] + 1;
        if(match < m) {match = m; flagrecord = c(i1, flag);}
    }
}
    error = 1 - match/n;    
    return(list(error = error, flagrecord = flagrecord))
}
