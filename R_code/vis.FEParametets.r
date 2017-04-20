################################################################
# Thursday 16.2.2017: similar to vis.MatParameters.r
# here we visualise (NOT burn in or thinning) the fixed effect parameters

# mebe adapt this to process chains

vis.FEParametets<- function(mcmc, sigma2.s, sigma2, beta0, beta1, b0, b1, s2, s2s){
  
  # trace -----------------------------------------------------
  x = seq(1:mcmc)
  
  trace.s2s<- ggplot(data.frame(x=x, y = sigma2.s), aes(x=x, y=y)) + geom_line() + geom_hline(yintercept = s2s, colour = "red") + 
    theme_bw() + theme(legend.position = "none") + ggtitle("s2s")
  
  trace.s2<- ggplot(data.frame(x=x, y = sigma2), aes(x=x, y=y)) + geom_line() + geom_hline(yintercept = s2, colour = "red") + 
    theme_bw() + theme(legend.position = "none") + ggtitle("s2")
  
  trace.b0<- ggplot(data.frame(x = x, y = beta0), aes(x=x, y=y)) + geom_line() + geom_hline(yintercept = b0, colour = "red") + 
    theme_bw() + theme(legend.position = "none") + ggtitle("b0")
  
  trace.b1<- ggplot(data.frame(x = x, y=beta1), aes(x=x,y=y)) + geom_line() + geom_hline(yintercept = b1, colour = "red") + 
    theme_bw() + theme(legend.position = "none") + ggtitle("b1")
  
  # density ---------------------------------------------------
  d.s2s<- ggplot(data.frame(x = sigma2.s), aes(x=x)) + geom_density() + geom_vline(xintercept = s2s, colour = "red") + 
    theme_bw() + theme(legend.position = "none") + ggtitle("s2s")
  
  d.s2<- ggplot(data.frame(x = sigma2), aes(x=x)) + geom_density() + geom_vline(xintercept = s2, colour = "red") + 
    theme_bw() + theme(legend.position = "none") + ggtitle("s2")
  
  d.b0<- ggplot(data.frame(x = beta0), aes(x=x)) + geom_density() + geom_vline(xintercept = b0, colour = "red") + 
    theme_bw() + theme(legend.position = "none") + ggtitle("b0")
  
  d.b1<- ggplot(data.frame(x = beta1), aes(x=x)) + geom_density() + geom_vline(xintercept = b1, colour = "red") + 
    theme_bw() + theme(legend.position = "none") + ggtitle("b1")
  
  # autocorr ---------------------------------------------------
  s2s.auto.corr<- with(acf(sigma2.s, plot=FALSE), data.frame(lag, acf))
 
  s2s.ac<- ggplot(data = s2s.auto.corr, aes(x=lag, y=acf)) + geom_hline(aes(yintercept = 0)) + theme_bw() +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + ggtitle("sigma2.s")
  
  s2.auto.corr<- with(acf(sigma2, plot=FALSE), data.frame(lag, acf))
  s2.ac<- ggplot(data = s2.auto.corr, aes(x=lag, y=acf)) + geom_hline(aes(yintercept = 0)) + theme_bw() +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + ggtitle("sigma2")
  
  b1.auto.corr<- with(acf(beta1, plot=FALSE), data.frame(lag, acf))
  b1.ac<- ggplot(data = s2.auto.corr, aes(x=lag, y=acf)) + geom_hline(aes(yintercept = 0)) + theme_bw() +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + ggtitle("beta1")
  
  b0.auto.corr<- with(acf(beta0, plot=FALSE), data.frame(lag, acf))
  b0.ac<- ggplot(data = s2.auto.corr, aes(x=lag, y=acf)) + geom_hline(aes(yintercept = 0)) + theme_bw() +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + ggtitle("beta0")
  
  ##################################################
   return(list(trace.s2s=trace.s2s, trace.s2=trace.s2, trace.b1=trace.b1, trace.b0=trace.b0,
               d.s2s=d.s2s, d.s2=d.s2, d.b1=d.b1, d.b0=d.b0,
               s2s.ac= s2s.ac, s2.ac=s2.ac, b1.ac= b1.ac, b0.ac=b0.ac))
}