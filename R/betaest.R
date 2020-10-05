betaest <- function(xo, sigmainv_xyo, diag_tol){
  xo_dims <- seq.int(1, ncol(xo))
  yo_dims <- seq.int(ncol(xo) + 1, ncol(sigmainv_xyo))
  sigmainv_xo <- sigmainv_xyo[, xo_dims, drop = FALSE]
  sigmainv_yo <- sigmainv_xyo[, yo_dims, drop = FALSE]
  invcov_betahat <- t(xo) %*% sigmainv_xo
  diag(invcov_betahat) <- diag(invcov_betahat) + diag_tol
  cov_betahat <- chol2inv(chol(invcov_betahat))
  betahat <- cov_betahat %*% t(xo) %*% sigmainv_yo
  return(list(betahat = betahat, cov_betahat = cov_betahat))
}
