#' Fit a spatio-temporal linear mixed model
#'
#' @param formula
#' @param xcoord
#' @param ycoord
#' @param tcoord
#' @param stcov
#' @param data
#' @param estmethod
#' @param sp_cor
#' @param t_cor
#' @param weights
#' @param initial
#' @param chol
#' @param diag_tol
#' @param max_v
#' @param max_s_range
#' @param max_t_range
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
stlmm <- function(data, ...){
  UseMethod("stlmm", object = data)
}

stlmm.data.frame <- function(formula, xcoord, ycoord = NULL, tcoord, stcov, data,
                             estmethod, sp_cor, t_cor, chol, diag_tol = 1e-4,
                             logdet = FALSE, optim_options = NULL,
                             h_options = NULL,
                             weights = NULL, initial = NULL,
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
  # invert_object <- make_invert_object(stcov = stcov,
  #                                     chol = chol, co = NULL,
  #                                     covparams = covest_output$par_r,
  #                                     diag_tol = diag_tol,
  #                                     h_s_large = data_object$h_s_large, h_t_large = data_object$h_t_large,
  #                                     h_s_small = data_object$h_s_small, h_t_small = data_object$h_t_small,
  #                                     n_s = data_object$n_s, n_t = data_object$n_t,
  #                                     logdet = logdet,
  #                                     m_index = data_object$m_index,
  #                                     o_index = data_object$o_index,
  #                                     sp_cor = sp_cor,
  #                                     t_cor = t_cor,
  #                                     xo = data_object$ordered_xo,
  #                                     yo = data_object$ordered_yo)
  # compute the inverse
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
                                 CovarianceForms = c(st_cov = stcov, sp_cor = sp_cor, t_cor = t_cor),
                                 formula = formula,
                                 model = list(FixedDesignMatrix = data_object$original_xo, Response = data_object$original_yo)),
                            class = "stlmm")

  # computing the residuals
  stlmm_object$Residuals <- residuals(stlmm_object = stlmm_object)

  # finally returning the overall output
  return(stlmm_object)
}








# stlmm.data.frame <- function(formula, xcoord, ycoord = NULL, tcoord, stcov, data,
#                   estmethod, sp_cor, t_cor, chol, diag_tol = 1e-4,
#                   logdet = FALSE, optim_options = NULL,
#                   h_options = NULL,
#                   weights = NULL, initial = NULL,
#                   max_v = NULL, max_s_range = NULL,
#                   max_t_range = NULL, ...){
#
#   # NOTES: MAKE THIS A GENERIC AT SOME POINT - stlmm.data.frame, stlmm.tibble, stlmm.sf, stlmm.sp, etc.
#
#   #can do the switch call here for the data object
#   data_object <- make_data_object(formula = formula, xcoord = xcoord, ycoord = ycoord,
#                                   tcoord = tcoord, stcov = stcov, data = data,
#                                   estmethod = estmethod, sp_cor = sp_cor, t_cor = t_cor,
#                                   weights = weights, initial = initial, chol = chol,
#                                   diag_tol = diag_tol, max_v = max_v, max_s_range = max_s_range,
#                                   max_t_range = max_t_range, h_options = h_options, logdet = logdet, ...)
#
#
#
#  # estimate the profiled covariance parameters
#   covest_output <- covest_wrapper(data_object = data_object, optim_options = optim_options, ...)
#
#
#   # invert object
#   invert_object <- make_invert_object(stcov = data_object$stcov,
#                                       chol = data_object$chol, co = NULL,
#                                       covparams = covest_output$par_r,
#                                       diag_tol = data_object$diag_tol,
#                                       h_s_large = data_object$h_s_large, h_t_large = data_object$h_t_large,
#                                       h_s_small = data_object$h_s_small, h_t_small = data_object$h_t_small,
#                                       n_s = data_object$n_s, n_t = data_object$n_t,
#                                       logdet = data_object$logdet,
#                                       m_index = data_object$m_index,
#                                       o_index = data_object$o_index,
#                                       sp_cor = data_object$sp_cor,
#                                       t_cor = data_object$t_cor,
#                                       xo = data_object$ordered_xo,
#                                       yo = data_object$ordered_yo)
#   # compute the inverse
#   inverse <- invert(invert_object)
#
#
#   # estimate the fixed effects
#   betaest_output <- betaest(xo = data_object$ordered_xo, sigmainv_xyo = inverse$siginv_o,
#                             diag_tol = data_object$diag_tol)
#
#   # return the relevant output
#   Coefficients <-  betaest_output$betahat
#   names(Coefficients) <- colnames(data_object$original_xo)
#   stlmm_object <- structure(list(CovarianceParameters = covest_output$par_r, Coefficients = betaest_output$betahat,
#               NamesCoefficients = names(Coefficients), CovCoefficients = betaest_output$cov_betahat,
#               ObjectiveFn = covest_output$value,
#               CovarianceForms = c(st_cov = stcov, sp_cor = sp_cor, t_cor = t_cor),
#               formula = formula,
#               model = list(FixedDesignMatrix = data_object$original_xo, Response = data_object$original_yo)),
#               class = "stlmm")
#
#   # computing the residuals
#   stlmm_object$Residuals <- residuals(stlmm_object = stlmm_object)
#
#   # finally returning the overall output
#   return(stlmm_object)
# }
#
#
#
