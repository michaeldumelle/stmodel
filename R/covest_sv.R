covest_sv <- function(covest_object, ...){
  UseMethod("covest_sv", object = covest_object)
}

covest_sv.productsum <- function(covest_object, ...) {

  covest_output <- optim(par = covest_object$lop_initial, fn = sv_fn.productsum,
        covest_object = covest_object, ...)
  covest_output$r_par <- exp(covest_output$par) / (1 + exp(covest_output$par))
  covest_output$r_par <- p2r_sv.productsum(par = covest_output$r_par, covest_object = covest_object)
  return(covest_output)
}
