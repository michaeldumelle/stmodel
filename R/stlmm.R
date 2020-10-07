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
#' @param max_srange
#' @param max_trange
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
# stlmm <- function(formula, xcoord, ycoord = NULL, tcoord, stcov, data,
#                  estmethod = c("reml", "svwls"), sp_cor, t_cor, weights, initial = NULL, chol,
#                  diag_tol = 1e-4, max_v = NULL, max_srange = NULL, max_trange = NULL, ...){
#
#   # to display the possible estimation methods
#   estmethod <- match.arg(estmethod)
#
#   # running the appropriate stlm function
#   output <- switch(estmethod,
#          "svwls" = stlm_svwls(formula = formula, xcoord = xcoord, ycoord = ycoord,
#                               tcoord = tcoord, stcov = stcov, data = data,
#                               sp_cor = sp_cor, t_cor = t_cor, weights = weights,
#                               initial = initial, chol = chol, diag_tol = diag_tol,
#                               max_v = max_v, max_srange = max_srange, max_trange = max_trange,
#                               ...),
#          stop("choose valid estimation method")) # error if they don't choose the proper one
#   # storing the class of the output for use in generics
#   class(output) <- "stlmm"
#
#   # computing the residuals
#   output$Residuals <- residuals(output, ...)
#
#   # finally returning the overall output
#   return(output)
# }

stlmm <- function(formula, xcoord, ycoord = NULL, tcoord, stcov, data,
                  estmethod, sp_cor, t_cor, weights = "cressie", initial = NULL, chol,
                  diag_tol = 1e-4, max_v = NULL, max_srange = NULL,
                  max_trange = NULL, logdet = FALSE, optim_defaults, ...){

  #can do the switch call here for the data object
  data_object <- make_data_object(formula = formula, xcoord = xcoord, ycoord = ycoord,
                                  tcoord = tcoord, stcov = stcov, data = data,
                                  estmethod = estmethod, sp_cor = sp_cor, t_cor = t_cor,
                                  weights = weights, initial = initial, chol = chol,
                                  diag_tol = diag_tol, max_v = max_v, max_srange = max_srange,
                                  max_trange = max_trange)


  # making the covariance parameter estimation object with the appropriate class
  # covest_object <- make_covest_object(initial = initial, max_srange = max_srange,
  #                              max_trange = max_trange, max_v = max_v,
  #                              sp_cor = sp_cor, sv = sv, t_cor = t_cor, weights = weights)
  # returning the parameter vector used for profiling (which does not help optimization here
  # but does give us the abilit to easily set maxes on the overall variance and ranges)

  # estimate the profiled covariance parameters
  covest_output <- covest_wrapper(data_object = data_object, optim_defaults = optim_defaults, ...)
  # covest_output <- optim(par = data_object$initial_plo, fn = covest,
  #                        data_object = data_object,
  #                        method = optim_defaults$method,
  #                        control = optim_defaults$control,
  #                        ...)
  #
  # # need to write a wrapper to give optim defaults
  # if (covest_output$convergence != 0) {
  #   warning("covariance parameter convergence may not have been achieved - consider
  #           setting new initial values, lowering the relative tolerance, or increasing
  #           the maximum iterations")
  # }
  # covest_output$par_r <- plo2r(par = covest_output$par, data_object = data_object)



  # invert object
  invert_object <- make_invert_object(stcov = data_object$stcov,
                                      chol = data_object$chol, co = NULL,
                                      covparams = covest_output$par_r,
                                      diag_tol = data_object$diag_tol,
                                      f_s = data_object$f_s, f_t = data_object$f_t,
                                      h_s = data_object$h_s, h_t = data_object$h_t,
                                      logdet = data_object$logdet,
                                      m_index = data_object$m_index,
                                      o_index = data_object$o_index,
                                      sp_cor = data_object$sp_cor,
                                      t_cor = data_object$t_cor,
                                      xo = data_object$ordered_xo,
                                      yo = data_object$ordered_yo)
  # compute the inverse
  inverse <- invert(invert_object)


  # estimate the fixed effects
  betaest_output <- betaest(xo = data_object$ordered_xo, sigmainv_xyo = inverse$siginv_o,
                            diag_tol = data_object$diag_tol)

  # return the relevant output
  Coefficients <-  betaest_output$betahat
  names(Coefficients) <- colnames(xo)
  stlmm_object <- structure(list(CovarianceParameters = covest_output$par_r, Coefficients = betaest_output$betahat,
              NamesCoefficients = colnames(raw_xo), CovCoefficients = betaest_output$cov_betahat,
              ObjectiveFn = covest_output$value,
              CovarianceForms = c(st_cov = stcov, sp_cor = sp_cor, t_cor = t_cor),
              formula = formula, model = list(FixedDesignMatrix = raw_xo, Response = raw_yo)),
              class = "stlmm")

  # computing the residuals
  stlmm_object$Residuals <- residuals(stlmm_object = stlmm_object)

  # finally returning the overall output
  return(stlmm_object)
}


# covest_sv_optim <- function(par, fn, covest_object, method = "Nelder-Mead",
#                          control = list(reltol = 1e-10, maxit = 5000), ...) {
#   return(optim(par = par, fn = fn,
#                covest_object = covest_object, method = method,
#                control = control, ...))
# }
