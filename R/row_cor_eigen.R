################################################################################
#                                                                              #
#                          Compile fast cpp functions                          #
#                                                                              #
################################################################################

############################ cpp version of cor    #############################

# Note. can be used to correlate two matrices, as well as a matrix with a vector
# returns correlation matrix with dimensions n.rows.x x n.rows.y
# this is synonymous to cor_eigen(t(X), t(Y))

row_cor_eigen <- '
Eigen::MatrixXd row_cor_eigen(Eigen::Map<Eigen::MatrixXd> & X1,
                              Eigen::Map<Eigen::MatrixXd> & Y1) {

  Eigen::MatrixXd X = X1;
  Eigen::MatrixXd Y = Y1;

  // Handle degenerate cases
  if (X.rows() == 0 || Y.rows() == 0) {
    return Eigen::MatrixXd::Constant(0, 0, 0);
  } else if (X.cols() == 0) {
    return Eigen::MatrixXd::Constant(X.rows(), Y.rows(),
                                     Rcpp::NumericVector::get_na());
  }

  // Computing degrees of freedom
  // n - 1 is the unbiased estimate whereas n is the MLE
  const int df = X.cols() - 1; // Subtract 1 by default

  // Centering matrices
  Y.colwise() -= Y.rowwise().mean();
  X.colwise() -= X.rowwise().mean();

  // The covariance matrix
  Eigen::MatrixXd cor = X * Y.transpose() / df;

  // Compute 1 over the standard deviations of X and Y
  Eigen::VectorXd inv_sds_X = (X.rowwise().norm()/sqrt(df)).array().inverse();
  Eigen::VectorXd inv_sds_Y = (Y.rowwise().norm()/sqrt(df)).array().inverse();

  // Scale the covariance matrix
  cor = cor.cwiseProduct(inv_sds_X * inv_sds_Y.transpose());
  return cor;
}
'

Rcpp::cppFunction(code = row_cor_eigen,
                  depends = "RcppEigen",
                  includes = "#include <RcppEigen.h>",
                  echo = FALSE)
