stlm <- function(){

}

stlm_svwls <- function(formula, xcoord, ycoord = NULL, tcoord, data, stcov,
                       sp_cor, t_cor, weights, initial, chol,
                       diag_tol = 1e-4, max_v = NULL, max_srange = NULL, max_trange = NULL, ...){

  # order the data by space within time
  spint <- order_spint(data = data, xcoord = xcoord,
                                  ycoord = ycoord, tcoord = tcoord, chol = FALSE)



  # create the model frame using the provided formula
  stmodel_frame <- model.frame(formula, spint$data_o,
                               na.action = stats::na.omit)

  # creating the fixed design matrix
  xo <- model.matrix(formula, stmodel_frame)

  # creating the response column
  yo <- model.response(stmodel_frame)

  # find the linear model residuals
  lmod_r <- as.vector((yo - xo %*% (chol2inv(chol(t(xo) %*% xo)) %*% (t(xo) %*% yo))))

  # find the sample variance
  lmod_s2 <- sum(lmod_r^2) / (nrow(xo) - ncol(xo))

  # make the semivariogram
  sv <- st_empsv(response = lmod_r, xcoord = data_o[[xcoord]], ycoord = data_o[[ycoord]],
                   tcoord = data_o[[tcoord]], ...)

  # provide default value for the maximum possible variance
  if (is.null(max_v)){
    max_v <- 4 * lmod_s2
  }
  # provide default value for the maximum possible spatial range
  if (is.null(max_srange)){
    max_srange <- 4 * max(spint$h_s)
  }
  # provide default value for the maximum possible temporal range
  if (is.null(max_trange)){
    max_trange <- 4 * max(spint$h_t)
  }

  # making the covariance parameter estimation object with the appropriate class
  covest_object <- make_covest_object(object_type = "covest", stcov = stcov,
                               initial = initial, max_srange = max_srange,
                               max_trange = max_trange, max_v = max_v,
                               sp_cor = sp_cor, sv = sv, t_cor = t_cor, weights = weights)
  # returning the parameter vector used for profiling (which does not help optimization here
  # but does give us the abilit to easily set maxes on the overall variance and ranges)

  # estimate the profiled covariance parameters
  covest_output <- optim(par = covest_object$initial_plo, fn = sv_fn,
                         covest_object = covest_object, ...)
  covest_output$par_r <- plo2r_sv(par = covest_output$par, covest_object = covest_object)



  # invert object
  invert_object <- make_invert_object(stcov = stcov,
                                      chol = chol, co = NULL,
                                      covparams = covest_output$par_r,
                                      diag_tol = diag_tol,
                                      f_s = spint$f_s, f_t = spint$f_t,
                                      h_s = spint$h_s, h_t = spint$h_t,
                                      logdet = logdet,
                                      m_index = spint$m_index,
                                      o_index = spint$o_index,
                                      sp_cor = sp_cor,
                                      t_cor = t_cor,
                                      xo = xo,
                                      yo = yo)
  inverse <- invert(invert_object = est_object)
  est_object$siginv_o <- inverse$siginv_o
  beta <- betaest(est_object)

  return(list(covparams = est_object$covparams, beta = beta$beta_hat, cov_beta = beta$cov_beta_hat))
}


