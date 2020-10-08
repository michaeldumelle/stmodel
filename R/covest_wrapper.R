covest_wrapper <- function(data_object, optim_options, ...){
  UseMethod("covest_wrapper", object = data_object)
}

covest_wrapper.svwls <- function(data_object, optim_options, ...){

  if (is.null(optim_options)){
    optim_options <- list(method = "Nelder-Mead", control = list(reltol = 1e-8, maxit = 10000))
  }

  covest_output <- optim(par = data_object$initial_plo, fn = covest.svwls,
                         data_object = data_object,
                         method = optim_options$method,
                         control = optim_options$control,
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
