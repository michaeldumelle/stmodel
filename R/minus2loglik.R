minus2loglik <- function(invert_object, invert_output){
  UseMethod("minus2loglik", object = invert_object)
}

minus2loglik.reml <- function(invert_object, invert_output){
  n <- nrow(invert_output$sigmainv_o)
  betaest_output <- betaest(xo = invert_object$xo, sigmainv_xyo = invert_output$sigmainv_o,
                            diag_tol = invert_object$diag_tol, return_estlist = TRUE)
  r_siginv_r <- t(invert_object$yo - invert_object$xo %*% betaest_output$betahat) %*%
    (betaest_output$estlist$sigmainv_yo - betaest_output$estlist$sigmainv_xo %*% betaest_output$betahat)
  m2ll <- invert_output$logdet + (n - betaest_output$estlist$p) * log(r_siginv_r) +
    betaest_output$estlist$ldet_cicb + (n - betaest_output$estlist$p) * (1 + log(2 * pi / (n - betaest_output$estlist$p)))
  return(m2ll)
}








minus2loglik <- function(dim_xo, dim_yo, siginv_xyc_o, logdet_o, diag_tol, estmethod){
  if (!(estmethod %in% c("reml", "ml"))) {
    stop("invalid estimation method")
  }
  n <- nrow(siginv_xyc_o)
  betaest_output
  siginv_xo <- siginv_xyc_o[, dim_xo, drop = FALSE]
  siginv_yo <- siginv_xyc_o[, dim_yo, drop = FALSE]
  xo_siginv_xo <- t(xo) %*% siginv_xo
  diag(xo_siginv_xo) <- diag(xo_siginv_xo) + diag_tol
  chol_xo_siginv_xo <- chol(xo_siginv_xo)
  beta <- chol2inv(chol_xo_siginv_xo) %*% t(xo) %*% siginv_yo
  r_siginv_r <- t(yo - xo %*% beta) %*% (siginv_yo - siginv_xo %*% beta)

  m2ll <- logdet_o + n * log(r_siginv_r)
  if (estmethod == "reml") {
    m2ll <- m2ll - dim_xo * log(r_siginv_r) + 2 * sum(log(diag(chol_xo_siginv_xo)))
  }
  return(m2ll)
}
