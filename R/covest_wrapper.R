covest_wrapper <- function(data_object, optim_defaults, ...){
  UseMethod("covest_wrapper", object = data_object)
}

covest_wrapper.svwls <- function(data_object, optim_defaults, ...){

  covest_output <- optim(par = data_object$initial_plo, fn = covest.svwls,
                         data_object = data_object,
                         method = optim_defaults$method,
                         control = optim_defaults$control,
                         ...)
  covest_output$par_r <- plo2r.svwls(par = covest_output$par, data_object = data_object)

  # need to write a wrapper to give optim defaults
  if (covest_output$convergence != 0) {
    warning("covariance parameter convergence may not have been achieved - consider
            setting new initial values, lowering the relative tolerance, or increasing
            the maximum iterations")
  }

  return(covest_output)


}
