covest_sv <- function(est_object, ...){
  UseMethod("covest_sv", object = est_object)
}

covest_sv.productsum <- function(est_object, ...) {

  opt_output <- optim(par = est_object$lo_initial, fn = sv_fn.productsum,
        est_object = est_object, method = "Nelder-Mead", ...)
  opt_output$par <- exp(opt_output$par) / (1 + exp(opt_output$par))
  opt_output$par <- p2r_sv.productsum(par = opt_output$par, est_object = est_object)
  return(opt_output)
}
