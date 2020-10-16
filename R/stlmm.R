stlmm <- function(data, ...){
  UseMethod("stlmm", object = data)
}

stlmm.data.frame <- function(formula, data, xcoord, ycoord = NULL, tcoord, stcov,
                             estmethod, sp_cor, t_cor, chol = FALSE, diag_tol = 1e-4,
                             logdet = FALSE, weights = NULL, initial = NULL,
                             optim_options = NULL, h_options = NULL,
                             max_options = NULL, stempsv_options = NULL, ...){

  # make the data object
  data_object <- make_data_object(formula = formula, xcoord = xcoord, ycoord = ycoord,
                                  tcoord = tcoord, data = data, h_options = h_options)

  # make the covest object
  covest_object <- make_covest_object(initial = initial, estmethod = estmethod,
                                        stcov = stcov, data_object = data_object, diag_tol = diag_tol, chol = chol,
                                        sp_cor = sp_cor, t_cor = t_cor, weights = weights, max_options = max_options,
                                        optim_options = optim_options, stempsv_options = stempsv_options)




  # estimate the profiled covariance parameters
  covest_output <- covest_wrapper(covest_object = covest_object, data_object = data_object)
  class(covest_output$par_r) <- class(covest_object)

  # invert object
  invert_object <- make_invert_object(covparam_object = covest_output$par_r,
                                      chol = chol,
                                      diag_tol = diag_tol,
                                      h_s_large = data_object$h_s_large, h_t_large = data_object$h_t_large,
                                      h_s_small = data_object$h_s_small, h_t_small = data_object$h_t_small,
                                      logdet = logdet,
                                      m_index = data_object$m_index,
                                      o_index = data_object$o_index,
                                      sp_cor = sp_cor,
                                      t_cor = t_cor,
                                      xo = data_object$ordered_xo,
                                      yo = data_object$ordered_yo)

  invert_output <- invert(invert_object)


  # estimate the fixed effects
  betaest_output <- betaest(xo = data_object$ordered_xo, sigmainv_xyo = invert_output$sigmainv_o,
                            diag_tol = diag_tol, return_estlist = FALSE)

  # return the relevant output
  Coefficients <-  betaest_output$betahat
  names(Coefficients) <- colnames(data_object$original_xo)
  stlmm_object <- structure(list(CovarianceParameters = covest_output$par_r, Coefficients = betaest_output$betahat,
                                 NamesCoefficients = names(Coefficients), CovCoefficients = betaest_output$cov_betahat,
                                 ObjectiveFn = covest_output$value,
                                 CovarianceForms = c(stcov = stcov, sp_cor = sp_cor, t_cor = t_cor),
                                 formula = formula,
                                 model = list(FixedDesignMatrix = data_object$original_xo, Response = data_object$original_yo)),
                            class = "stlmm")

  # computing the residuals
  stlmm_object$Residuals <- residuals(stlmm_object = stlmm_object)

  # finally returning the overall output
  return(stlmm_object)
}

