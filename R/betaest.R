betaest <- function(xo, sigmainv_xyo, diag_tol, return_estlist = FALSE){
  ncol_xo <- ncol(xo)
  xo_dims <- seq.int(1, ncol_xo)
  sigmainv_xo <- sigmainv_xyo[, xo_dims, drop = FALSE]
  sigmainv_yo <- sigmainv_xyo[, ncol_xo + 1, drop = FALSE]
  invcov_betahat <- t(xo) %*% sigmainv_xo
  diag(invcov_betahat) <- diag(invcov_betahat) + diag_tol
  chol_invcov_betahat <- chol(invcov_betahat)
  cov_betahat <- chol2inv(chol_invcov_betahat)
  betahat <- cov_betahat %*% t(xo) %*% sigmainv_yo
  betaest_output <- list(betahat = betahat, cov_betahat = cov_betahat)
  if (return_estlist) {
    betaest_output$estlist <- list(ldet_cicb = 2 * sum(log(diag(chol_invcov_betahat))),
                                          sigmainv_xo = sigmainv_xo, sigmainv_yo = sigmainv_yo,
                                   p = ncol_xo)
  }
  return(betaest_output)
}
