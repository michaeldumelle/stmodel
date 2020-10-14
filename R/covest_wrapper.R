covest_wrapper <- function(covest_object, data_object){
  UseMethod("covest_wrapper", object = covest_object)
}

covest_wrapper.svwls <- function(covest_object, data_object){



  covest_output <- optim(par = covest_object$initial_plo, fn = covest.svwls,
                         covest_object = covest_object,
                         data_object = data_object,
                         method = optim_options$method,
                         control = optim_options$control)
  covest_output$par_r <- plo2r.svwls(par = covest_output$par, covest_object = covest_object)

  # need to write a wrapper to give optim defaults
  if (covest_output$convergence != 0) {
    warning("covariance parameter convergence may not have been achieved - consider
            setting new initial values, lowering the relative tolerance, or increasing
            the maximum iterations")
  }

  return(covest_output)


}
