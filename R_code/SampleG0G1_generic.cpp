// Code in C++ to replace sampling  mat.g0 and mat.G1

//#include <Rcpp.h>
  #include <RcppArmadillo.h>
  // [[Rcpp::depends(RcppArmadillo)]]

//using namespace Rcpp;
using namespace std;
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::cube SampleG0G1_generic(cube matW, NumericVector age, NumericMatrix G0, NumericMatrix G1, NumericVector mu_g, 
                      NumericMatrix Sigma_g, NumericMatrix g_acc_rej, NumericMatrix rand_w, int I, int K){
  
  colvec g_curr(2);
  colvec g_star(2);
  
  colvec g_curr_mean(2);
  colvec g_star_mean(2);
  
  double log_g_curr;
  double log_g_star;
  double alpha;
  double Urand;
  double det_Sigma_g;
  //  mat acc_rej1;
  //  acc_rej1.zeros(10, 10);
  
  mat Sigma_g1 =Rcpp::as<arma::mat>(Sigma_g);
  det_Sigma_g = det(Sigma_g1);
  //  acc_rej1 = Rcpp::as<arma::mat>(g_acc_rej);
  mat rand_w1 = Rcpp::as<arma::mat>(rand_w);
  colvec mu_g1 = Rcpp::as<arma::colvec>(mu_g);
  mat g_Acc_rej = Rcpp::as<arma::mat>(g_acc_rej);
  
  int nrow = K;
  colvec sum_i_kj_curr(I);
  colvec sum_i_kj_star(I);
  
  colvec temp1(2);
  temp1.fill(1);
  mat var1 = diagmat(temp1);
  mat tempMat(K, K);
  mat invSigma = Sigma_g1.i();
  
  mat g0 = Rcpp::as<arma::mat>(G0);
  mat g1 = Rcpp::as<arma::mat>(G1);
  
  for(int r = 0; r < (nrow-1); r ++){
    for(int c = (r+1); c < nrow; c++){
      //int r=0;
      //int c=1;
      g_curr[0] = g0(r,c);
      g_curr[1] = g1(r,c);
      
      var1(0,0) = rand_w1(r,c);
      var1(1,1) = rand_w1(r,c);
      
      g_star = g_curr + trans(randn(1,2)*chol(var1)); // generate a multivariate normal
      g_curr_mean = g_curr - mu_g1;
      g_star_mean = g_star - mu_g1;
      
      for(int i = 0; i < I; i++){
        tempMat = matW.slice(i);
        sum_i_kj_curr[i] = tempMat(r,c)*(g_curr[0] + g_curr[1]*age[i]) - log(1 + exp(g_curr[0] + g_curr[1]*age[i]));
        sum_i_kj_star[i] = tempMat(r,c)*(g_star[0] + g_star[1]*age[i]) -log(1 + exp(g_star[0] + g_star[1]*age[i]));     
      } 
      
      log_g_curr = sum(sum_i_kj_curr) - K/2*log(2*datum::pi) - 0.5*log(det_Sigma_g) - 0.5*as_scalar(trans(g_curr_mean)*invSigma*g_curr_mean);
      log_g_star = sum(sum_i_kj_star) - K/2*log(2*datum::pi) - 0.5*log(det_Sigma_g) - 0.5*as_scalar(trans(g_star_mean)*invSigma*g_star_mean);                                                              
      
      alpha = exp(log_g_star - log_g_curr);
      
      //      g0(r,c) = g_curr[0]; 
      //      g0(c,r) = g_curr[0];
      //      g1(r,c) = g_curr[1];
      //      g1(c,r) = g_curr[1];
      
      Urand = randu();
      if(Urand < alpha){
        g0(r,c) = g_star[0]; 
        g0(c,r) = g_star[0];
        g1(r,c) = g_star[1];
        g1(c,r) = g_star[1];
        g_Acc_rej(r,c) = g_Acc_rej(r,c) + 1;
        g_Acc_rej(c,r) = g_Acc_rej(c,r) + 1;
      }
      
    }
  }
  
  cube reSults(K, K, 3, fill::zeros);
  reSults.slice(0) = g0;
  reSults.slice(1) = g1;
  reSults.slice(2) = g_Acc_rej;
  
  return(reSults );
  //return(Urand);
  //return(wrap(g0));
}  // end function