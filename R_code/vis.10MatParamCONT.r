##################################################
# This is similar to vis.MatParam2_T11a.r
# only here it's for continuous age groups
# as per meeting with JM (10.3.2017) here we do not want to view the W matrices any more
# for continuous age groups -- too complex
# rather we aim to view posterior probability matrices from the logit^-1(.)
# with respect to age

vis.10MatParamCONT<- function(g.acc.rej, #acc.rate.MAT,  
                         age,
                         mcmc,  mat.g0, mat.g1, rand.w,
                         gamma_0, gamma_1, mat.W, K = 10){
  
  library(raster)
  library(ggplot2)
  
  off.diag = K*(K-1)/2
  I = dim(age)[1]
  
  ##
  ## acceptance rates & random walk ********************************************
  ## 
  acc.r<-randWALK<- raster(xmn = 0, xmx =K, ymn = 0, ymx = K, nrows = K, ncols = K)
  acc.r[]<- as.vector(g.acc.rej/mcmc)
  
  diag(rand.w)<- 1
  randWALK[]<- as.vector(rand.w)
  
  #x11()
  #par(mfrow = c(1,2))
  #plot(acc.r, main = "Acceptance rates for gamma")
  #plot(randWALK, main = "Random walk variance mat")
  
  ## 
  ## view trace plot of acceptance rates - to check adaptive MH *************************
  ##   Note- for large runs (eg 50K) this takes up a lot of RAM
#    acc.rate.MAT = acc.rate.MAT[,,2:mcmc] # at t starts from t = 2
#    acc.Rate.trace<- matrix(0, nrow = dim(acc.rate.MAT)[3], ncol = off.diag)
#    cOuNt = 1
  
#    for(row in 1:(K-1) ){  # loop thru upper triangular
#      for(col in (row + 1): K){
#        acc.Rate.trace[,cOuNt]<- acc.rate.MAT[row,col,]
#        cOuNt = cOuNt + 1
#      }
#    }
  
#    acc.R.t<- data.frame(acc.Rate.trace)
#    acc.R.t$iter<- seq(2:mcmc)
#    acc.R.t1<- melt(acc.R.t, id.vars = "iter")
#    acc.R.t2<- acc.R.t1[acc.R.t1$value !=0, ] # now remove those instances where the acceptance rate was 0, 
  # i.e on the t^th iteration is did not accept a proposal
  #x11()
  #  acc.trace<- ggplot(acc.R.t2, aes(x = iter, y = value)) + geom_line() + facet_wrap(~variable) + theme_bw() + 
  #    ggtitle("Acceptance rate trace plots") +
  #    geom_hline(yintercept = c(0.2, 0.4), colour = "red")
  
  ##
  ## Posterior Mean: gamma0 and gamma1 *******************************************
  ## 
  g0.m<- apply(mat.g0, 1:2, mean)
  g1.m<- apply(mat.g1, 1:2, mean)
  
  r.g0<- r.g0.sol<- r.g1.sol<- r.g1<-raster(xmn = 0, xmx =K, ymn = 0, ymx = K, nrows = K, ncols = K)
  r.g0[]<- as.vector(g0.m); r.g0.sol[]<- as.vector(gamma_0); r.g1.sol[]<- as.vector(gamma_1)
  r.g1[]<- as.vector(g1.m)
  
  #x11()
  #par(mfrow = c(2,2))
  #plot(r.g0, main = "G0 posterior mean")
  #plot(r.g0.sol, main = "G0 solution")
  #plot(r.g1, main= "G1 posterior mean")
  #plot(r.g1.sol, main = "G1 solution")
  
  ##
  ## trace plots for G0 and G1 ****************************************
  ##
  G0<- G1<- matrix(0, nrow = dim(mat.g0)[3], ncol = off.diag) # take the off-diagonals of the matrices pkeep and Wkeep
  couNt= 1
  for(row in 1:(K-1) ){  # loop thru upper triangular
    for(col in (row + 1): K){
      G0[,couNt]<- mat.g0[row,col,]
      G1[, couNt]<- mat.g1[row, col,]
      couNt = couNt + 1
    }
  }
  
  G0<- data.frame(G0); G0$iter<- seq(1:mcmc); G0.1<- melt(G0, id.vars = "iter")
  G1<- data.frame(G1); G1$iter<- seq(1:mcmc); G1.1<- melt(G1, id.vars = "iter")
  
#  x11()
  tr.G0<- ggplot(G0.1, aes(x = iter, y = value)) + geom_line() + facet_wrap(~variable) + 
    theme_bw() + ggtitle("Gamma 0 trace plot") + geom_hline(yintercept = c(g0.upp, g0.low), colour = "red")
  
  #x11()
  tr.G1<- ggplot(G1.1, aes(x = iter, y = value)) + geom_line() + facet_wrap(~variable) + theme_bw() + 
    ggtitle("Gamma 1 trace plot") + geom_hline(yintercept =c(g1.upp, g1.low), colour = "red")

  ##  
  ## initial values 
  ##
  G0.s<- G1.s<- matrix(0, nrow = dim(mat.g0)[3], ncol = 45)
  
  init.g0<- init.g1<-raster(xmn = 0, xmx =K, ymn = 0, ymx = K, nrows = K, ncols = K)
  init.g0[]<- as.vector(mat.g0[,,1]) #as.vector(d.gamma$init.g0)
  init.g1[]<- as.vector(mat.g1[,,1]) #as.vector(d.gamma$init.g1)
  
  #x11()
  #par(mfrow = c(1,2))
  #plot(init.g0, main = "initial G0")
  #plot(init.g1, main = "initial G1")
  
  ##
  ## binary matrices to denote if gamma_0 and gamma_1 solutions are in 95% CI intervals
  ##
  
  low.g0<- apply(mat.g0, 1:2, quantile, prob=0.025)
  low.g1<- apply(mat.g1, 1:2, quantile, prob=0.025)
  
  up.g0<- apply(mat.g0, 1:2, quantile, prob=0.975)
  up.g1<- apply(mat.g1, 1:2, quantile, prob=0.975)
  
  g1.inSol<- ifelse(low.g0< gamma_0 & gamma_0 < up.g0, 1, 0)
  g0.inSol<- ifelse(low.g1< gamma_1 & gamma_1 < up.g1, 1, 0)
  
  g1.InSol<- g0.InSol<-raster(xmn = 0, xmx =K, ymn = 0, ymx = K, nrows = K, ncols = K)
  g1.InSol[]<-as.vector(g1.inSol)
  g0.InSol[]<- as.vector(g0.inSol)
  
#  x11()
#  par(mfrow = c(1,2))
#  plot(g0.InSol,legend=F, main = "Solution for G0 in 95% CI")
#  plot(g1.InSol, legend=F, main = "Solution for G1 in 95% CI")
  
  ##
  ## here to plot the posterior probability matrices (from the output of mat.g0 and mat.g1) ***********************
  ## 
  
  # to divi up continuous age into 4 distinct age groups
  p.age1<- p.age2<- p.age3<- p.age4<- p.age5<- array(0, dim = c(K, K, mcmc))
  
  age.r2<- quantile(seq(from=min(age), to=max(age), by=0.1), probs = c(0, 0.33, 0.66, 1))
  for(t in 1:mcmc){
  
    for(r in 1:(K-1)){
      for(c in (r+1):K){
        p.age1[r,c,t]<- p.age1[c,r,t]<- exp(mat.g0[r,c,t] + mat.g1[r,c,t]*age.r2[1])/
                        (1 + exp(mat.g0[r,c,t] + mat.g1[r,c,t]*age.r2[1]))
        
        p.age2[r,c,t]<- p.age2[c,r,t]<- exp(mat.g0[r,c,t] + mat.g1[r,c,t]*age.r2[2])/
                        (1 + exp(mat.g0[r,c,t] + mat.g1[r,c,t]*age.r2[2]))
        
        p.age3[r,c,t]<- p.age3[c,r,t]<- exp(mat.g0[r,c,t] + mat.g1[r,c,t]*age.r2[3])/
                        (1 + exp(mat.g0[r,c,t] + mat.g1[r,c,t]*age.r2[3]))
        
        p.age4[r,c,t]<- p.age4[c,r,t]<- exp(mat.g0[r,c,t] + mat.g1[r,c,t]*age.r2[4])/
          (1 + exp(mat.g0[r,c,t] + mat.g1[r,c,t]*age.r2[4]))
        
        #p.age5[r,c,t]<- p.age4[c,r,t]<- exp(mat.g0[r,c,t] + mat.g1[r,c,t]*age.r2[5])/
        #  (1 + exp(mat.g0[r,c,t] + mat.g1[r,c,t]*age.r2[5]))
        }
    }
    
  }
  p.m1<- apply(p.age1, 1:2, mean)  # <=== need to find the mean and a way to plot this!
  p.m2<- apply(p.age2, 1:2, mean)
  p.m3<- apply(p.age3, 1:2, mean)
  p.m4<- apply(p.age4, 1:2, mean)
#  p.m5<- apply(p.age5, 1:2, mean)
  
  p1<- p2<- p3<- p4<-raster(xmn = 0, xmx =K, ymn = 0, ymx = K, nrows = K, ncols = K)
  p1[]<- as.vector(p.m1)
  p2[]<- as.vector(p.m2)
  p3[]<- as.vector(p.m3)
  p4[]<- as.vector(p.m4)
  
#  x11()
#  par(mfrow=c(2,2))
#  plot(p1, main = paste("Posterior prob for age ",round(age.r2[1],2), sep = ""))
#  plot(p2, main = paste("Posterior prob for age ",round(age.r2[2],2), sep = ""))
#  plot(p3, main = paste("Posterior prob for age ",round(age.r2[3],2), sep = ""))
#  plot(p4, main = paste("Posterior prob for age ",round(age.r2[4],2), sep = ""))
  
  ##########################################################
  return(list(acc.r = acc.r, randWALK = randWALK, #acc.trace=acc.trace, 
              r.g0=r.g0, r.g1=r.g1, r.g1.sol = r.g1.sol, r.g0.sol = r.g0.sol,tr.G0=tr.G0, tr.g1=tr.G1,
              init.g0= init.g0, init.g1=init.g1, p1=p1, p2=p2, p3=p3, p4=p4, age.r2=age.r2,
              g1.InSol=g1.InSol, g0.InSol=g0.InSol))
  
}