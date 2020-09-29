betaest <- function(betaest_object, ...){
  siginv_xo <- betaest_object$siginv_o[, seq.int(1, ncol(xo))]
  siginv_yo <- betaest_object$siginv_o[, ncol(xo) + 1]
  cov_beta_hat <- chol2inv(chol(t(betaest_object$xo) %*% siginv_xo))
  beta_hat <- cov_beta_hat %*% t(xo) %*% siginv_yo
  return(list(beta_hat = beta_hat, cov_beta_hat = cov_beta_hat))
}
