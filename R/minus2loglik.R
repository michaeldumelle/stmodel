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
