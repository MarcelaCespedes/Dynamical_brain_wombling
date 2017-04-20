#####################################################3
# Sunday 22 Jan
# this code corresponds with GenSpatTemp_Wombling_T4.r
# Here we generate different starting values for 
# gamma matrices and W

diff.gamma.start<- function(K, I, R, age){
  W_young = matrix(c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0,
                     1, 0, 1, 1, 1, 0, 0, 0, 0, 0,
                     1, 1, 0, 1, 1, 0, 0, 0, 0, 0,
                     1, 1, 1, 0, 1, 0, 0, 0, 0, 0,
                     1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
                     0, 0, 0, 0, 0, 1, 0, 1, 1, 1,
                     0, 0, 0, 0, 0, 1, 1, 0, 1, 1,
                     0, 0, 0, 0, 0, 1, 1, 1, 0, 1,
                     0, 0, 0, 0, 0, 1, 1, 1, 1, 0), K, K)

  W_old<- matrix(c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
                   1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                   0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 
                   0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 
                   0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
                   0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 
                   0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 
                   0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 
                   0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 
                   0, 0, 0, 0, 0, 0, 0, 0, 1, 0), K, K)
  
  W_young<- W_old
  
  gamma_0<- gamma_1<- matrix(0, K,K)
  g0.high<-2
  g0.low<- -2
  
  g1.high<- 2
  g1.low<- -2
  for(r in 1:(K-1)){
    for(c in r:K){
      
      if(W_young[r, c] == 1){
        gamma_0[r, c]<- gamma_0[c, r]<- g0.high #rnorm(1, mean = 3, sd = 0.1)
      }else{
        gamma_0[r, c]<- gamma_0[c, r]<- g0.low #rnorm(1, mean = -3, sd = 0.1)
      }
      
      if(W_old[r, c] == 1){
        gamma_1[r, c]<- gamma_1[c, r]<- g1.high #rnorm(1, mean = 5, sd = 0.1)
      }else{
        gamma_1[r, c]<- gamma_1[c, r]<- g1.low #rnorm(1, mean = -5, sd = 0.1)
      }
    }
  }
  
  diag(gamma_0)<- 0; diag(gamma_1)<- 0
  # check
  W_young
  round(gamma_0, 2)
  
  W_old
  round(gamma_1, 2)
  ###########################################
  return(list(init.g0 = gamma_0, 
              init.g1 = gamma_1, W1 = W_young, w2 = W_old))
}