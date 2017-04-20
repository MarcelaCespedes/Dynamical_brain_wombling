##################################################################
# Monday 3rd April 2017

# This code is for UNLABANCED Spatio-Temporal (dynamical) wombling.
# profiled to run faster by running sections written in C++
# ... soo many of my full conditionals are the same as the original wombling algorithm
# ... so this contains a lot of the samplers (b0, s2 and s2s) as the original code

# useful links for starting on Rcpp
# http://adv-r.had.co.nz/Rcpp.html
# commands on how to use C++
# http://arma.sourceforge.net/docs.html

library(RcppArmadillo)
library(Rcpp)
library(installr) 
install.Rtools()  # <== test for Rcpp
evalCpp("1+1") 

# note: mvnfast has trouble working in the HPC (Linux systems - works well on Windows)
library(mvnfast) # <-- faster alternative to mvtnorm, note: MASS also has MVN random sim and density
#library(MASS)  # <-- to draw MVN random samples
library(ggplot2)
library(mefa)
library(reshape2)
library(raster)
rm(list = ls())
#########################################################
# Now generate Age at discrete values 
# there are now 4 age groups

I = 200  # <==== as we will now have people in between the whole (-2, 1) spectrum, I need to have
# enough groups of people with ages similar to each other

no.female = 109

av.1<- c(4,7,3,2,5,3,4,6,7,3)
R.vec<- rep(av.1, times = I/10)
length(R.vec)  # <== this should  == I

K=10
no.Obs<- K*sum(R.vec)
no.Obs

# Binary matrices to derive gamma_0 and gamma_1 matrices <=== NOT young or old matrices
W_g1<- matrix(c(1,1,1,1,1,1,1,1,0,0,
                1,1,1,1,1,1,1,1,0,0,
                1,1,1,1,1,1,1,1,0,0,
                1,1,1,1,1,1,1,1,0,0,
                1,1,1,1,1,1,1,1,0,0,
                1,1,1,1,1,1,1,1,0,0,
                1,1,1,1,1,1,1,1,0,0,
                1,1,1,1,1,1,1,1,0,0,
                0,0,0,0,0,0,0,0,1,1,
                0,0,0,0,0,0,0,0,1,1), K, K)
diag(W_g1)<- 0
W_g1

W_g2<- matrix(c(1,1,1,1,1,1,0,0,0,0,
                1,1,1,1,1,1,0,0,0,0,
                1,1,1,1,1,1,0,0,0,0,
                1,1,1,1,1,1,0,0,0,0,
                1,1,1,1,1,1,0,0,0,0,
                1,1,1,1,1,1,0,0,0,0,
                0,0,0,0,0,0,1,1,1,1,
                0,0,0,0,0,0,1,1,1,1,
                0,0,0,0,0,0,1,1,1,1,
                0,0,0,0,0,0,1,1,1,1), K, K)
diag(W_g2)<- 0
W_g2

##
## Age matrix (for baseline only)

age<- c()
for(i in 1:I){
  age[i]<- runif(n=1, min = -2, max=2)
}
head(age, 50)

# Take a peek as the age vector divi into 3 groups
age.g1<- ifelse(age < -0.5, 1, 0)
sum(age.g1)

age.g2<- ifelse(-0.5 < age & age < 0.5, 1, 0)
sum(age.g2)

age.g3<- ifelse(0.5< age, 1,0)
sum(age.g3)

gamma_0<- gamma_1<- matrix(0, K,K)
g0.upp<- 2   #5
g0.low<- -2 #-5

g1.upp<- 3 #10
g1.low<- -3 #-10

for(r in 1:(K-1)){
  for(c in (r+1):K){
    
    if(W_g1[r, c] == 1){
      gamma_0[r, c]<- gamma_0[c, r]<- g0.upp #rnorm(1, mean = g0.upp, sd = 0.1)
    }else{
      gamma_0[r, c]<- gamma_0[c, r]<- g0.low #rnorm(1, mean = g0.low, sd = 0.1)
      
    }
    
    if(W_g2[r, c] == 1){
      gamma_1[r, c]<- gamma_1[c, r]<- g1.upp  #rnorm(1, mean = g1.upp, sd = 0.1)
    }else{
      gamma_1[r, c]<- gamma_1[c, r]<- g1.low  #rnorm(1, mean = g1.low, sd = 0.1)
    }
  }
}

diag(gamma_0)<- 0; diag(gamma_1)<- 0
# check
W_g1
round(gamma_0, 2)

W_g2
round(gamma_1, 2)

################## now to generate data according to the model
#P.kj<- W.kj<- array(0, dim = c(K, K, I*R))
couNt<- 1
rho= 0.9
s2s = 1; s2 = 0.5
b0 = 3
b1 = 0.5
X.mat<- matrix(1, nrow = I, ncol = 2)
X.mat[, 2]<-c(rep(1, no.female), rep(0, I-no.female))
e.1<- rep(1, K)

sol.mat.W<- array(0, dim=c(K, K, I))
W_g1_2<- W_g1-W_g2
W.t<- ifelse(W_g2 == 1, 0,1)
diag(W.t)<- 0


##
## Set mat.W for each person at baseline age
##

for(i in 1:I){   # <== initialise W at solution
  
  for(r in 1:(K-1)){
    for(c in (r+1):K){
      p.kj.i<- exp(gamma_0[r,c] + gamma_1[r,c]*age[i])/(1 + exp(gamma_0[r,c] + gamma_1[r,c]*age[i]))
      sol.mat.W[r, c, i]<- sol.mat.W[c, r, i]<- rbinom(1,1, prob = p.kj.i)  
      
    }
  }
}

spatDat<- matrix(0, nrow = I, ncol = K)
for(i in 1:I){ # <=== this only works when people don't change age groups
  
  Q = rho*(diag(K)*rowSums(sol.mat.W[,,i]) - sol.mat.W[,,i]) + (1-rho)*diag(K)
  spatDat[i, ]<- rmvn(1, mu = rep(0, K), sigma = s2s*solve(Q)) # <== create unique set of spatial random effects for all I
}

# spatial random effects in IR format
spatDat1<- as.matrix(rep(data.frame(spatDat),R.vec)) # <== repeat for each replicate

# try a different way to generate data
vec.data<- as.vector(t(spatDat1))
x1<- rep(X.mat[, 2], K*R.vec)

#x1 = c(rep(1, K*sum(R.vec[1:no.female])), rep(0,K*sum(R.vec[(1 + no.female):I]) )   )
y = b0 + b1*x1 + vec.data + rnorm(no.Obs, mean = 0, sd = sqrt(s2))
y.vec = y
tR<- matrix(y, nrow = K, ncol = sum(R.vec))
y.mat<- t(tR)

########################################################################################
########################################################################################
# Initialise MCMC chains for dynamical wombling

# first prep MCMC
mcmc = 10000

sigma2<- sigma2.s<- beta0<- beta1<- rep(0, mcmc)
sigma2[1] = runif(1, min = 0.0001, max = 1)
sigma2.s[1] = runif(1, min = 0.0001, max = 3)
beta0[1] = rnorm(1, mean = b0, sd= 0.5)
beta1[1] = rnorm(1, mean = b1, sd = 0.5)

# priors and covariates for these parameters
X.mat.long = as.matrix(rep(data.frame(X.mat), (K*R.vec)))

#library(invgamma)  # <=====these priors have been set for optimal values -DO NOT CHANGE
#x<- seq(from = 0.001, to=5, by = 0.01)

sh.s2s = 1; rat.s2s = 2 # prior for s2s
#y.s2s = dinvgamma(x, shape = sh.s2s, rate = rat.s2s)

# prior for s2
sh.s2 = 1; rat.s2 = 1
#y.s2 <- dinvgamma(x, shape = sh.s2, rate = rat.s2)

#x11() # this is to check priors for variance terms are well supported by the prior
#ggplot(data.frame(x =c(x,x), y = c(y.s2s,y.s2), s= factor(c(rep("s2s", length(x)), rep("s2", length(x))))), 
#       aes(x=x, y=y, colour = s)) + geom_line() + ggtitle("Inverse gamma prior")+ geom_vline(xintercept = s2s, colour = "blue") +
#       geom_vline(xintercept = s2, colour = "red")

Sigma.o.inv<- diag(2)*0.1  # priors for beta.vec
Sigma.o.inv
mu.o<- c(b0, b1) 


###################################################
###################################################
# need to set priors for gamma_0 and gamma_1

mat.W<- array(0, dim = c(K, K, I, mcmc))

mu_g<- c(0, 0)
Sigma.g <- diag(2)*4 # <=== Although it seems that we have informative priors on the 
                     # covariance structure of Gamma priors, it comes out to be non-informative
                     # in the logit scale.

mat.g0<- mat.g1<- array(0, dim = c(K, K, mcmc))

# initiate different starting values for gamma's ::: same data, just different starting values
source("diff.gamma.start.r")
d.gamma<- diff.gamma.start(I=I, K=K, R=R, age=age)

mat.g0[,,1]<- d.gamma$init.g0
mat.g1[,,1]<- d.gamma$init.g1

#mat.g0[,,1]<- gamma_0  # <== again for now start at solution
#mat.g1[,,1]<- gamma_1

for(i in 1:I){  # initialise mat.W
  for(r in 1:(K-1)){
    for(c in (r+1):K){
      p.kj.i<- exp(mat.g0[r,c,1] + mat.g1[r,c,1]*age[i])/(1 +exp(mat.g0[r,c,1] + mat.g1[r,c,1]*age[i]))  
      mat.W[r,c,i,1]<- mat.W[c,r,i,1]<- rbinom(1,1,prob=p.kj.i)
    }
  }
}

g.acc.rej<- matrix(0, K, K)
acc.rate.MAT<- array(0, dim = c(K,K,mcmc)); acc.rate.counter<- 1
rand.w<- matrix(5, K, K)  # random walk parameter

##
## set random effects array

b.mat<- matrix(0, nrow = I, ncol = K)
b.mat.m<- array(0, dim = c(I, K, mcmc))
b.mat.m[,,1]<- spatDat  # <== initialise at solution
e.1<- rep(1, K)

head(age)

# note: snipper below only works for balanced design
#n<- R  # sum columns-wise every 3 elements (the spatial random effects) from y.mat
#sum.i.rep<-apply(y.mat, 2, function(x) tapply(x, ceiling(seq_along(x)/n), sum))

# As we have unbalanced simulations: use code below
# sum over all the replicates for person i - note this doesn't change... can do this external to the code
sum.i.rep<- matrix(0, nrow = I, ncol = K)
for(i in 1:I){
  
  if(R.vec[i] != 1){
    if(i == 1){
      sum.i.rep[i,] <- colSums(y.mat[ (i + (i-1)*R.vec[i]) :R.vec[i], ]) # checked :)
    }else{
      sum.i.rep[i, ]<- colSums(y.mat[ (sum(R.vec[1:i-1])+1)  :sum(R.vec[1:i]), ])
    }
  }else{
    sum.i.rep[i,]<- y.mat[i,]
  }
  
}

########################################################################################
# start adaptive Metropolis Hastings within Gibbs for dynamical wombling
#
########################################################################################

# load up C++ snippets
Rcpp::sourceCpp('SampleG0G1.cpp')
Rcpp::sourceCpp('SampleWs10.cpp')

#----------------------------------------------------------------------
tiMe<- proc.time() 
pb <- txtProgressBar(style = 3)

for(t in 2:mcmc){
  
  ##################################################################################
  # b_ir  --- as W changes with respect to age at the rth replicate
  # b.mat<- spatDat
  
  beta.vec <- c(beta0[t-1], beta1[t-1])
  for(i in 1:I){
    Q.inv = rho*(diag(K)*rowSums(mat.W[,,i,t-1]) -mat.W[,,i,t-1]) + (1-rho)*diag(K)
    Omega = solve(R.vec[i]/sigma2[t-1]*diag(K) + 1/sigma2.s[t-1]*Q.inv)
    mu = Omega%*%(1/sigma2[t-1]*sum.i.rep[i,]- R.vec[i]/sigma2[t-1]*X.mat[i,]%*%beta.vec*e.1)
    b.mat[i,]<- rmvn(1, mu=mu, sigma=Omega)
  }  
  b.mat.m[,,t]<- b.mat
  
  #round(b.mat[1:5, ], 2)
  #round(spatDat[1:5,], 2)
  
  ###############################################
  # sigma2.s  ---- Metropolis Hastings scheme ----------------
  # sigma2.s[t] <- sigma2.s[t-1]
  
  shape.s2s = (I*K + 2*sh.s2s)/2    # from Gibbs attempt with earlier model (no replicates)
  #shape.s2s                         # the one which keeps track of the number of islands s2s^(K-G)/2
  
  sum.re<- c()
  for(i in 1:I){
    inv.Q = rho*(diag(K)*rowSums(mat.W[,,i,t-1]) -mat.W[,,i,t-1]) + (1-rho)*diag(K)
    sum.re[i]<- b.mat[i,]%*%inv.Q%*%b.mat[i,]
  }
  
  rate.s2s = sum(sum.re)/2 + rat.s2s
  #rate.s2s
  sigma2.s[t]<- 1/rgamma(1, shape = shape.s2s, rate = rate.s2s)
  #head(sigma2.s)
  
  ################################################
  # sigma2  
  #  sigma2[t] <- sigma2[t-1]
  
  shape.s2 = no.Obs/2 + sh.s2  # note the shape for this distribution won't change
  #shape.s2
  
  linear.pred<- c(); count = 1
  b.rep<- as.matrix(rep(data.frame(b.mat), R.vec)) # repeat each b.mat column R times
  b.vec<- as.vector(t(b.rep))
  
  for(obs in 1:no.Obs){
    linear.pred[obs]<- (y.vec[obs] - X.mat.long[obs,]%*%c(beta0[t-1], beta1[t-1]) - b.vec[obs])^2
  }
  
  rate.s2 = sum(linear.pred)/2 + rat.s2
  #rate.s2
  sigma2[t]<- 1/rgamma(1, shape = shape.s2, rate = rate.s2)
  head(sigma2)
  
  ################################################
  # [beta0, beta1] --- same as old wombling
  #beta0[t] <- beta0[t-1]; beta1[t] <- beta1[t-1]
  
  Omega.b<- solve(t(X.mat.long)%*%X.mat.long/sigma2[t] + Sigma.o.inv)
  b.rep<- as.matrix(rep(data.frame(b.mat), R.vec)) # repeat each b.mat column R.vec times
  b.vec<- as.vector(t(b.rep))
  
  mu = t(t(X.mat.long)%*%(y.vec - b.vec)/sigma2[t] +  t(mu.o%*%Sigma.o.inv))%*%Omega.b
  beta.vec<- rmvn(1, mu=mu, sigma=Omega.b)
  #beta.vec
  beta0[t]<- beta.vec[1]
  beta1[t]<- beta.vec[2]
  
  ##########################################################################################
  # Update gamma_0, gamma_1 
  #mat.g0[,,t]<-mat.g0[,,t-1]
  #mat.g1[,,t]<- mat.g1[,,t-1]
  
  op2<- SampleG0G1(matW = mat.W[,,,t-1], age=age, G0 = mat.g0[,,t-1], G1 = mat.g1[,,t-1], 
                   mu_g=mu_g, Sigma_g=Sigma.g,
                   g_acc_rej=g.acc.rej, rand_w= rand.w, I=I)
  #op2
  
  mat.g0[,,t] <- op2[,,1]
  mat.g1[,,t] <- op2[,,2]
  g.acc.rej <- op2[,,3]
  
  # old fashion way: run with R
  #for(r in 1:(K-1)){
  #  for(c in (r+1):K){
      
  #    gamm.vec.curr<- c(mat.g0[r, c, t-1], mat.g1[r, c, t-1])
  #    gamm.vec.curr
  #    gamm.vec.star<- gamm.vec.curr + rmvn(1, mu = c(0, 0), sigma =diag(2)*rand.w[r, c])  
  #    gamm.vec.star
      
  #    sum.i.kj.curr<- sum.i.kj.star<- c()
  #    for(i in 1:I){
  #      sum.i.kj.curr[i]<- mat.W[r, c, i,t-1]*(gamm.vec.curr[1] + gamm.vec.curr[2]*age[i]) -
  #        log(1 + exp(gamm.vec.curr[1] + gamm.vec.curr[2]*age[i]))
        
  #      sum.i.kj.star[i]<- mat.W[r, c, i,t-1]*(gamm.vec.star[1] + gamm.vec.star[2]*age[i]) -
  #        log(1 + exp(gamm.vec.star[1] + gamm.vec.star[2]*age[i]))
        
  #    }
  #    log.g.curr<- sum(sum.i.kj.curr) + dmvn(gamm.vec.curr, mu = mu_g, sigma = Sigma.g, log=TRUE)
  #    log.g.star<- sum(sum.i.kj.star) + dmvn(gamm.vec.star, mu = mu_g, sigma = Sigma.g, log=TRUE)
      
  #    alpha = exp(log.g.star - log.g.curr)
  #    alpha
      
  #    mat.g0[r, c, t]<- mat.g0[c, r, t]<- gamm.vec.curr[1]
  #    mat.g1[r, c, t]<- mat.g1[c, r, t]<-gamm.vec.curr[2]
      
  #    if(runif(1)< alpha){
  #      mat.g0[r, c, t]<- mat.g0[c, r, t]<- gamm.vec.star[1]
  #      mat.g1[r, c, t]<- mat.g1[c, r, t]<-gamm.vec.star[2]
  #      g.acc.rej[r, c]<- g.acc.rej[c, r]<- g.acc.rej[r, c] + 1
        # below is just for plotting purposes - for alrge MCMC this also take up a lot of memory
        #acc.rate.MAT[row, col, t]<- acc.rate.MAT[col, row, t]<- g.acc.rej[row, col]/t
  #    }
  #  }
  #} 
  
  #mat.g0[,,t-1]
  #round(mat.g0[,,t],3)
  
  #round(mat.g1[,,t], 2)
  #mat.g1[,,1]
  
  #****************************************************************************************
  # adaptive random walk MH step: to be implemented after the first 50 MCMC iterations
  for(row in 1:(K-1)){
    for(col in (row + 1):K){
      
      if(g.acc.rej[row, col]/t < 0.2 & g.acc.rej[row, col]/t > 0.01 & t > 100 & rand.w[row, col] > 0.11){
        rand.w[row, col]<- rand.w[row,col] - 0.1 # <== make random walk smaller if acceptance rate too low
      }
      
      if(g.acc.rej[row, col]/t > 0.4 & t > 100 & rand.w[row, col] < 100){
        rand.w[row, col]<- rand.w[row, col] + 0.1 # <== make random walk bigger if acceptance rate high
      }
    }
  }
  
  ########################################################################################
  # sample W_ir  -- this is taken from GenSpatTemp_Wombling_T4.r
  
  # mat.W[,,1,t]<- mat.W[,,1,t-1]
  # mat.W[,,2,t]<- mat.W[,,2,t-1]
  
  for(i in 1:I){
    
    #mat.W[,,i,t]<- mat.W[,,i, t-1]  #<== leave it at solution
    op1<- SampleWs10(w = mat.W[,,i,t-1], age=age[i], rho=rho, s2s=sigma2.s[t], b_mat=b.mat[i,],
                     G0 = mat.g0[,,t], G1 = mat.g1[,,t])
    #op1
    mat.W[,,i,t]<- op1
  }
  
  # old fashion way: run with R
  #for(i in 1:I){
  #  for(row in 1:(K-1)){
  #    for(col in (row + 1):K){
        
  #      m1<- m0<- mat.W[,,i,t-1]
  #      m1[row, col]<- m1[col, row]<- 1
  #      m0[row, col]<- m0[col, row]<-0
        
  #      Q.inv.1 = rho*(diag(K)*rowSums(m1) - m1) + (1-rho)*diag(K)
  #      Q.inv.0 = rho*(diag(K)*rowSums(m0) - m0) + (1-rho)*diag(K)
        
  #      l.p0<- 1/2*log(det(Q.inv.0)) - 1/(2*sigma2.s[t])*b.mat[i,]%*%Q.inv.0%*%b.mat[i,] + 
  #        m0[row, col]*(mat.g0[row, col,t-1] + mat.g1[row, col,t-1]*age[i]) - 
  #        log(1 + exp(mat.g0[row, col,t-1] + mat.g1[row, col,t-1]*age[i])) 
        
  #      l.p1<- 1/2*log(det(Q.inv.1)) - 1/(2*sigma2.s[t])*b.mat[i,]%*%Q.inv.1%*%b.mat[i,] + 
  #        m1[row, col]*(mat.g0[row, col,t-1] + mat.g1[row, col,t-1]*age[i]) - 
  #        log(1 + exp(mat.g0[row, col,t-1] + mat.g1[row, col,t-1]*age[i])) 
        
  #      prob.1<- 1/(1 + exp(l.p0 - l.p1))
  #      mat.W[row, col, i, t]<- mat.W[col, row, i, t]<- rbinom(1,1,prob = prob.1)
  #    }   
  #  }
  #}
  
  #mat.W[,,1,1]
  #mat.W[,,1,2]
  ####################################################################
  prog<- t/mcmc 
  setTxtProgressBar(pb, prog) 
}

(proc.time() - tiMe)/60 # time in minutes


#startOFFSol<- list(mat.W=mat.W, beta=beta0, beta1=beta1, sigma2=sigma2, sigma2.s=sigma2.s, mat.g0=mat.g0,
#                    mat.g1=mat.g1, b.mat.m=b.mat.m, y.vec=y.vec, X.mat=X.mat, X.mat.long=X.mat.long, mcmc=mcmc,
#                    gamma_0=gamma_0, gamma_1=gamma_1, age=age, g.acc.rej=g.acc.rej, rand.w=rand.w, 
#                    note=  "This took xx hours", spatDat=spatDat)

#save(startOFFSol, file = "startOFFSol.Rdata")

source("multiplot.r")
#############################################################################
#############################################################################
# process fixed effect and variance parameters

source("vis.FEParametets.r")
op.p<- vis.FEParametets(mcmc=mcmc, sigma2=sigma2, sigma2.s=sigma2.s, beta0=beta0, beta1=beta1,
                        s2=s2, s2s=s2s, b0=b0, b1=b1)

x11() # <== trace plots
multiplot(op.p$trace.s2s, op.p$trace.s2, op.p$trace.b0, op.p$trace.b1, cols = 2)

x11() # <== density plots
multiplot(op.p$d.s2s, op.p$d.s2, op.p$d.b0, op.p$d.b1, cols = 2)

x11() # <=== autocorrelation
multiplot(op.p$s2s.ac, op.p$s2.ac, op.p$b0.ac, op.p$b1.ac, cols = 2)




############################################################################################
source("vis.10MatParamCONT.r") # <== note: these functions do not burnIn or thin chainc - only visualise!

op.m<- vis.10MatParamCONT(g.acc.rej = g.acc.rej, #acc.rate.MAT=acc.rate.MAT, 
                          mat.g0=mat.g0, mat.W = mat.W, mcmc=mcmc,rand.w=rand.w,age=age,
                          mat.g1=mat.g1, gamma_0=gamma_0, gamma_1=gamma_1)

x11()# <=== diagnositcs for G0 and G1
par(mfrow = c(2,2))
plot(op.m$r.g0, main = "G0 posterior mean")
plot(op.m$r.g0.sol, main = "G0 solution")
plot(op.m$r.g1, main= "G1 posterior mean")
plot(op.m$r.g1.sol, main = "G1 solution")

x11()
par(mfrow = c(1,2))
plot(op.m$g0.InSol,legend=F, main = "Solution for G0 in 95% CI")
plot(op.m$g1.InSol, legend=F, main = "Solution for G1 in 95% CI")


x11()
par(mfrow = c(1,2))
plot(op.m$acc.r, main = "Acceptance rates for gamma")
plot(op.m$randWALK, main = "Random walk variance mat")

#x11() # <== acceptance rate trace plots. !! when MCMC is large this 
#plot(op.m$acc.trace) # will require a lot of RAM

##### trace plots for G0 and G1 ****************************************
x11()
plot(op.m$tr.G0)

X11()
plot(op.m$tr.g1)

x11() # <=== initial values
par(mfrow = c(1,2))
plot(op.m$init.g0, main = "initial G0")
plot(op.m$init.g1, main = "initial G1")

x11() # <==== posterior mean for probabilities at quantiles (0, 0.33, 0.66, 1) 
par(mfrow =c(2,2))
plot(op.m$p1, main = paste("Posterior prob for age ",round(op.m$age.r2[1],2), sep = ""))
plot(op.m$p2, main = paste("Posterior prob for age ",round(op.m$age.r2[2],2), sep = ""))
plot(op.m$p3, main = paste("Posterior prob for age ",round(op.m$age.r2[3],2), sep = ""))
plot(op.m$p4, main = paste("Posterior prob for age ",round(op.m$age.r2[4],2), sep = ""))

#############################################################
##############################################################
pdf.name<- "DynamicalWombling_10K.pdf"

pdf(pdf.name, onefile = T)
plot(1:20, type = 'n', xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
text(10, 16, paste("Note: processed chains, burnIn=2K, thin=20. All started away from sol"))
text(10, 14, paste("I ppl: ", I, sep = ""))
text(10, 13, paste("R unbalanced ", sep = ""))
text(10, 12, paste("K regions: ", K, sep = ""))
text(10, 11, paste("No mcmc: ", mcmc, sep = ""))

multiplot(op.p$trace.s2s, op.p$trace.s2, op.p$trace.b0, op.p$trace.b1, cols = 2) # <== trace plots
multiplot(op.p$d.s2s, op.p$d.s2, op.p$d.b0, op.p$d.b1, cols = 2)  # <== density plots
multiplot(op.p$s2s.ac, op.p$s2.ac, op.p$b0.ac, op.p$b1.ac, cols = 2) # <=== autocorrelation

par(mfrow = c(1,2))
plot(op.m$acc.r, main = "Acceptance rates for gamma")
plot(op.m$randWALK, main = "Random walk variance mat")
#plot(op.m$acc.trace)

par(mfrow = c(2,2))
plot(op.m$r.g0, main = "G0 posterior mean")
plot(op.m$r.g0.sol, main = "G0 solution")
plot(op.m$r.g1, main= "G1 posterior mean")
plot(op.m$r.g1.sol, main = "G1 solution")

par(mfrow = c(1,2)) # <== binary matrices which determine if solution in 95% CI for each value mat.G0 and mat.G1
plot(op.m$g0.InSol,legend=F, main = "Solution for G0 in 95% CI")
plot(op.m$g1.InSol, legend=F, main = "Solution for G1 in 95% CI") 

# trace plot G0 and G1
plot(op.m$tr.G0)
plot(op.m$tr.g1)

#par(mfrow = c(1,2)) # <=== initial values
#plot(op.m$init.g0, main = "initial G0")
#plot(op.m$init.g1, main = "initial G1")

par(mfrow =c(2,2)) #<== assess posterior probability matrices for 4 age quantiles 
plot(op.m$p1, main = paste("Posterior prob for age ",round(op.m$age.r2[1],2), sep = ""))
plot(op.m$p2, main = paste("Posterior prob for age ",round(op.m$age.r2[2],2), sep = ""))
plot(op.m$p3, main = paste("Posterior prob for age ",round(op.m$age.r2[3],2), sep = ""))
plot(op.m$p4, main = paste("Posterior prob for age ",round(op.m$age.r2[4],2), sep = ""))
dev.off()

