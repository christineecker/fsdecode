################################################################################
#                                                                              #
#                          Compile fast cpp functions                          #
#                                                                              #
################################################################################

############################ cpp version of cor    #############################

# Note. can be used to correlate two matrices, as well as a matrix with a vector
# returns correlation matrix with dimensions n.cols.x x n.cols.y

cor_eigen <- '
Eigen::MatrixXd cor_eigen(Eigen::Map<Eigen::MatrixXd> & X,
                          Eigen::Map<Eigen::MatrixXd> & Y) {

  // Handle degenerate cases
  if (X.cols() == 0 || Y.cols() == 0) {
    return Eigen::MatrixXd::Constant(0, 0, 0);
  } else if (X.rows() == 0) { // && X.cols() > 0 && Y.cols() > 0 implicit
    return Eigen::MatrixXd::Constant(X.cols(), Y.cols(),
                                     Rcpp::NumericVector::get_na());
  }

  // Computing degrees of freedom
  // n - 1 is the unbiased estimate whereas n is the MLE
  const int df = X.rows() - 1; // Subtract 1 by default

  // Centering matrices
  Y.rowwise() -= Y.colwise().mean();
  X.rowwise() -= X.colwise().mean();

  // The covariance matrix
  Eigen::MatrixXd cor = X.transpose() * Y / df;

  // Compute 1 over the standard deviations of X and Y
  Eigen::VectorXd inv_sds_X = (X.colwise().norm()/sqrt(df)).array().inverse();
  Eigen::VectorXd inv_sds_Y = (Y.colwise().norm()/sqrt(df)).array().inverse();

  // Scale the covariance matrix
  cor = cor.cwiseProduct(inv_sds_X * inv_sds_Y.transpose());
  return cor;
}
'

Rcpp::cppFunction(code = cor_eigen,
                  depends = "RcppEigen",
                  includes = "#include <RcppEigen.h>",
                  echo = FALSE)


############################ cpp version of cor with arma    #############################

# xcorArma <- '
# arma::mat xcorArma(const arma::mat &X, const arma::mat &Y){
#   return arma::cor(X,Y);
# }
# '
#
# Rcpp::cppFunction(code = xcorArma,
#                   depends = "RcppArmadillo",
#                   includes = "#include <RcppArmadillo.h>")
