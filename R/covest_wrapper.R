covest_wrapper <- function(covest_object, data_object){
  UseMethod("covest_wrapper", object = covest_object)
}

covest_wrapper.svwls <- function(covest_object, data_object){



  covest_output <- optim(par = covest_object$initial_plo, fn = covest.svwls,
                         covest_object = covest_object,
                         data_object = data_object,
                         method = covest_object$optim_options$method,
                         control = covest_object$optim_options$control)
  covest_output$par_r <- plo2r.svwls(par = covest_output$par, covest_object = covest_object)

  # need to write a wrapper to give optim defaults
  if (covest_output$convergence != 0) {
    warning("covariance parameter convergence may not have been achieved - consider
            setting new initial values, lowering the relative tolerance, or increasing
            the maximum iterations")
  }

  return(covest_output)


}

covest_wrapper.reml <- function(covest_object, data_object){

  initial_plo_noclass <- covest_object$initial_plo
  class(covest_object$initial_plo) <- class(covest_object)
  invert_object <- make_invert_object(covparam_object = covest_object$initial_plo,
                                      chol = covest_object$chol, co = NULL,
                                      diag_tol = covest_object$diag_tol,
                                      h_s_large = data_object$h_s_large,
                                      h_t_large = data_object$h_t_large,
                                      h_s_small = data_object$h_s_small,
                                      h_t_small = data_object$h_t_small,
                                      logdet = covest_object$logdet,
                                      m_index = data_object$m_index,
                                      o_index = data_object$o_index,
                                      sp_cor = covest_object$sp_cor,
                                      t_cor = covest_object$t_cor,
                                      xo = data_object$ordered_xo,
                                      yo = data_object$ordered_yo)

  covest_output <- optim(par = initial_plo_noclass, fn = covest.reml,
                         covest_object = covest_object,
                         invert_object = invert_object,
                         method = covest_object$optim_options$method,
                         control = covest_object$optim_options$control)

  invert_object$covparams <- plo2r.reml(covest_output$par, covest_object = covest_object, ov_var = 1)
  invert_output <- invert(invert_object)
  ov_var <- varest.reml(invert_object = invert_object, invert_output = invert_output)
  covest_output$par_r <- plo2r.reml(covest_output$par, covest_object = covest_object, ov_var = ov_var)
  # get overall variance
  #getovvar
  #covest_output$par_r <- plo2r.reml(par = covest_output$par, covest_object = covest_object, #ov var)

  # need to write a wrapper to give optim defaults
  if (covest_output$convergence != 0) {
    warning("covariance parameter convergence may not have been achieved - consider
            setting new initial values, lowering the relative tolerance, or increasing
            the maximum iterations")
  }

  return(covest_output)


}


