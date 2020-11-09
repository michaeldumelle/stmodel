stlmm <- function(data, formula, ...){
  UseMethod("stlmm", object = data)
}

stlmm.data.frame <- function(data, formula, xcoord, ycoord = NULL, tcoord, stcov,
                             estmethod = "reml", s_cor = "exponential", t_cor = "exponential", chol = FALSE, condition = 1e-4,
                             logdet = FALSE, weights = "cressie", initial = NULL,
                             optim_options = NULL, h_options = NULL,
                             max_options = NULL, stempsv_options = NULL, ...){

  # make the data object
  data_object <- make_data_object(formula = formula, xcoord = xcoord, ycoord = ycoord,
                                  tcoord = tcoord, data = data, h_options = h_options)

  # make the covest object
  covest_object <- make_covest_object(initial = initial, estmethod = estmethod,
                                        stcov = stcov, data_object = data_object, condition = condition, chol = chol,
                                        s_cor = s_cor, t_cor = t_cor, weights = weights, max_options = max_options,
                                        optim_options = optim_options, stempsv_options = stempsv_options)




  # estimate the profiled covariance parameters
  covest_output <- covest_wrapper(covest_object = covest_object, data_object = data_object)
  class(covest_output$par_r) <- class(covest_object)

  # invert object
  invert_object <- make_invert_object(covparam_object = covest_output$par_r,
                                      chol = chol,
                                      condition = condition,
                                      h_s_large = data_object$h_s_large, h_t_large = data_object$h_t_large,
                                      h_s_small = data_object$h_s_small, h_t_small = data_object$h_t_small,
                                      logdet = logdet,
                                      m_index = data_object$m_index,
                                      o_index = data_object$o_index,
                                      s_cor = s_cor,
                                      t_cor = t_cor,
                                      xo = data_object$ordered_xo,
                                      yo = data_object$ordered_yo)

  invert_output <- invert(invert_object)


  # estimate the fixed effects
  betaest_output <- betaest(xo = data_object$ordered_xo, sigmainv_xyo = invert_output$sigmainv_o,
                            condition = condition, return_estlist = FALSE)

  # return the relevant output
  Coefficients <-  betaest_output$betahat
  names(Coefficients) <- colnames(data_object$original_xo)
  stlmm_object <- structure(list(CovarianceParameters = covest_output$par_r, Coefficients = betaest_output$betahat,
                                 NamesCoefficients = names(Coefficients), CovCoefficients = betaest_output$cov_betahat,
                                 ObjectiveFn = covest_output$value,
                                 CovarianceForms = c(stcov = stcov, s_cor = s_cor, t_cor = t_cor),
                                 formula = formula,
                                 model = list(FixedDesignMatrix = data_object$original_xo, Response = data_object$original_yo),
                                 data_object = data_object, invert_output = invert_output,
                                 coordnames = list(xcoord = xcoord, ycoord = ycoord, tcoord = tcoord),
                                 coords = list(xcoord = data_object$ordered_data_o[[xcoord]],
                                               ycoord = data_object$ordered_data_o[[ycoord]],
                                               tcoord = data_object$ordered_data_o[[tcoord]]),
                                 h_options = data_object$h_options,
                                 stempsv_options = covest_object$stempsv_options,
                                 optim_options = covest_object$optim_options,
                                 max_options = covest_object$max_options,
                                 chol = chol,
                                 condition = condition),
                            class = "stlmm")

  # computing the residuals
  stlmm_object$Residuals <- residuals(stlmm_object = stlmm_object)

  # finally returning the overall output
  return(stlmm_object)
}

