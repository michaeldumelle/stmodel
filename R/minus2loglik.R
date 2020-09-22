minus2loglik <- function(dim_xo, dim_yo, siginv_xyc_o, logdet_o, diag_tol, estmethod){
  if (!(estmethod %in% c("reml", "ml"))) {
    stop("invalid estimation method")
  }
  n <- nrow(siginv_xyc_o)
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
