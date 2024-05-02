#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;


// [[Rcpp::export]]
double dmvtnorm_diag_arma(arma::vec y, arma::vec mu, arma::vec sigma_sq){
  int d = y.n_rows;
  return( pow(2*arma::datum::pi, -0.5 * d) * pow(arma::prod(sigma_sq),-0.5) * exp(-0.5 * arma::dot(y - mu, arma::inv_sympd(diagmat(sigma_sq)) * (y-mu)) ));
}

// [[Rcpp::export]]
double log_dmvtnorm_diag_arma(arma::vec y, arma::vec mu, arma::vec sigma_sq){
  int d = y.n_rows;
  return( -0.5* d * log(2 * arma::datum::pi) - 0.5 * arma::sum(log(sigma_sq)) -0.5 * arma::dot(y - mu, arma::inv_sympd(diagmat(sigma_sq)) * (y-mu)) );
}

// [[Rcpp::export]]
int sample_arma(arma::vec probs) {
  int K = probs.n_rows;
  IntegerVector clusts = Range(1,K);
  IntegerVector samp = RcppArmadillo::sample(clusts, 1, TRUE, probs);
  int s = samp(0);
  return(s);
}

// [[Rcpp::export]]
arma::vec update_c_marginal(arma::mat y, arma::vec c, 
                            arma::mat mu, arma::mat sigma_sq, 
                            double alpha_lambda) {
  int n = y.n_rows;
  int M = mu.n_rows;
  int n_minus_i = 0;
  arma::vec prbs = arma::zeros(M,1);
  arma::vec eprbs = arma::zeros(M,1);
  for (int i=0;i<n;i++) {
    for (int m=0; m<M; m++) {
      // compute cell counts without ith entry
      if (c(i)==(m+1)) {
        n_minus_i = sum(c==(m+1))-1;
      } else {
        n_minus_i = sum(c==(m+1));
      }
      prbs(m) = log(n_minus_i + alpha_lambda/(M*1.0)) + log_dmvtnorm_diag_arma(y.row(i).t(), mu.row(m).t(), sigma_sq.row(m).t());
    }
    eprbs = arma::exp(prbs - arma::max(prbs));
    c(i) = sample_arma(eprbs/sum(eprbs));
  }
  return(c);
}
