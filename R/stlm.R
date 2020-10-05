stlm <- function(formula, xcoord, ycoord = NULL, tcoord, stcov, data,
                 estmethod = c("reml", "svwls"), sp_cor, t_cor, weights, initial = NULL, chol,
                 diag_tol = 1e-4, max_v = NULL, max_srange = NULL, max_trange = NULL, ...){
  estmethod <- match.arg(estmethod)
  switch(estmethod,
         "svwls" = stlm_svwls(formula = formula, xcoord = xcoord, ycoord = ycoord,
                              tcoord = tcoord, stcov = stcov, data = data,
                              sp_cor = sp_cor, t_cor = t_cor, weights = weights,
                              initial = initial, chol = chol, diag_tol = diag_tol,
                              max_v = max_v, max_srange = max_srange, max_trange = max_trange,
                              ...),
         stop("choose valid estimation method"))
}

stlm_svwls <- function(formula, xcoord, ycoord, tcoord, stcov, data,
                       sp_cor, t_cor, weights, initial, chol,
                       diag_tol, max_v, max_srange, max_trange, ...){

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
  sv <- st_empsv(response = lmod_r, xcoord = spint$data_o[[xcoord]], ycoord = spint$data_o[[ycoord]],
                   tcoord = spint$data_o[[tcoord]], ...)

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


  # setting initial values if there are none specified
  if (is.null(initial)){
    initial <- make_covparam_object(s_de = 1, s_ie = 1, t_de = 1,
                                    t_ie = 1, st_de = 1, st_ie = 1,
                                    s_range = max_srange / 8,
                                    t_range = max_trange / 8,
                                    stcov = stcov)
    vparm_names <- c("s_de", "s_ie", "t_de", "t_ie",
    "st_de", "st_ie")
    numparams <- sum(vparm_names %in% names(initial))
    initial[names(initial) %in% vparm_names] <- lmod_s2 / numparams
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
                         covest_object = covest_object, method = "Nelder-Mead",
                         control = list(reltol = 1e-10, maxit = 5000), ...)
  if (covest_output$convergence != 0) {
    warning("covariance parameter convergence may not have been achieved")
  }
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
  # compute the inverse
  inverse <- invert(invert_object)


  # estimate the fixed effects
  betaest_output <- betaest(xo = xo, sigmainv_xyo = inverse$siginv_o,
                            diag_tol = diag_tol)

  # return the relevant output
  return(list(covparams = covest_output$par_r, betahat = betaest_output$betahat,
              cov_betahat = betaest_output$cov_betahat))
}


