make_beta_object <- function(xo, sigmainv_xyo, diag_tol, stcov){
  xo_dims <- seq.int(1, ncol(xo))
  yo_dims <- seq.int(ncol(xo) + 1, ncol(sigmainv_xyo))
  sigmainv_xo <- sigmainv_xyo[, xo_dims, drop = FALSE]
  sigmainv_yo <- sigmainv_xyo[, yo_dims, drop = FALSE]
  beta_object <- structure(list(xo = xo, diag_tol = diag_tol, sigmainv_xo = sigmainv_xo,
                                sigmainv_yo = sigmainv_yo),
                           class = stcov)
  return(beta_object)
}
