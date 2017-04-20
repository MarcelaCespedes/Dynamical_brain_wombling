// Code in C++ to replace samplign W's

//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//using namespace Rcpp;
using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericMatrix SampleWs10(NumericMatrix w, double age, double rho, double s2s, 
                         NumericVector b_mat, NumericMatrix G0, NumericMatrix G1){
  
  int nrow = w.nrow();
  colvec rhoVec1(10);
  rhoVec1.fill(1-rho);

  double DeT_0;
  double DeT_1;
  double l_p0;
  double l_p1;
  double InvLogit;
  double proB;
  
  int Wdraws;
//  double test1;
//  double test2;
//  rowvec probs(45);
//  int couNt=1;
  //NumericVector Wdraws;
  
  mat Q1;
  mat Q0;
  mat Rho_a;
  mat D1;
  mat D0;
  mat temp1;
  mat temp0;
  NumericMatrix matWfinal(10, 10);
  
  colvec t1;  // with help from Ace on getting this working 15.3.2017
  colvec t0;
  colvec Ones_p(10);
  Ones_p.fill(1);
  
  colvec Bmat;
  Bmat = b_mat;
  
  mat mat_g0 = Rcpp::as<arma::mat>(G0);
  mat mat_g1 = Rcpp::as<arma::mat>(G1);
  
  Rho_a = diagmat(rhoVec1);
  for(int r = 0; r < (nrow-1); r ++){
    for(int c = (r+1); c < nrow; c++){
//int r=0;
//int c=1;
    
      mat w1 = Rcpp::as<arma::mat>(w);
      mat q0 = Rcpp::as<arma::mat>(w);
      q0(r, c) = 0; 
      q0(c, r) = 0;
      w1(r, c) = 1;
      w1(c, r) = 1;
      
      t1 = w1*Ones_p;  // row sum of w1
      D1 = diagmat(t1);  // convert it into a diagonal matrix
      temp1 = D1 - w1;
      Q1 = rho*temp1 + Rho_a;
      
      t0 = q0*Ones_p;
      D0 = diagmat(t0);
      temp0 = D0 - q0;
      Q0 = rho*temp0 + Rho_a;  
      
      DeT_0 = det(Q0);
      DeT_1 = det(Q1);
      
      InvLogit = mat_g0(r,c) + mat_g1(r,c)*age;
      //InvLogit = G0(r,c) + G1(r,c)*age;
//      test1 = as_scalar(trans(Bmat)*Q0*Bmat );
//      test2 = as_scalar(trans(Bmat)*Q1*Bmat);


//        l.p0<- 1/2*log(det(Q.inv.0)) - 1/(2*sigma2.s[t])*b.mat[i,]%*%Q.inv.0%*%b.mat[i,] + 
//          m0[row, col]*(mat.g0[row, col,t-1] + mat.g1[row, col,t-1]*age[i]) - 
//          log(1 + exp(mat.g0[row, col,t-1] + mat.g1[row, col,t-1]*age[i])) 

      l_p0 = 0.5*log(DeT_0) - 1/(2*s2s)*as_scalar(trans(Bmat)*Q0*Bmat ) + q0(r,c)*InvLogit - log(1 + exp(InvLogit));
      l_p1 = 0.5*log(DeT_1) - 1/(2*s2s)*as_scalar(trans(Bmat)*Q1*Bmat ) + w1(r,c)*InvLogit - log(1 + exp(InvLogit));
     
      proB = 1/(1 + exp(l_p0 - l_p1));
//      probs[couNt] = proB;
//      couNt = couNt + 1;
      //Wdraws = as<double>(rbinom(1, 1, proB));  // this results in a lot of ones at t-> increases
      Wdraws = R::rbinom(1,proB);
      matWfinal(r,c) = Wdraws;
      matWfinal(c,r) = Wdraws;
    }
  }
  //return(wrap(matWfinal));
  return(wrap(matWfinal));
  //return();
}