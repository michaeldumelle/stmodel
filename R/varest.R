varest <- function(invert_object, invert_output){
  UseMethod("varest", object = invert_object)
}

varest.reml <- function(invert_object, invert_output){
  n <- nrow(invert_output$sigmainv_o)
  betaest_output <- betaest(xo = invert_object$xo, sigmainv_xyo = invert_output$sigmainv_o,
                            diag_tol = invert_object$diag_tol, return_estlist = TRUE)
  r_siginv_r <- t(invert_object$yo - invert_object$xo %*% betaest_output$betahat) %*%
    (betaest_output$estlist$sigmainv_yo - betaest_output$estlist$sigmainv_xo %*% betaest_output$betahat)
  ws_l2 <- as.vector(r_siginv_r / (n - betaest_output$estlist$p))
  return(ws_l2) # thanks russ and john!
}
