stlm <- function(){

}

stlm_svwls <- function(formula, xcoord, ycoord = NULL, tcoord, data, stcov,
                       sp_cor, t_cor, weights, initial, chol, diagtol = 1e-4, max_v, max_srange, max_trange, ...){

  # order the data
  spint <- order_spint(data = data, xcoord = xcoord,
                                  ycoord = ycoord, tcoord = tcoord)
  o_index <- spint$data$index[spint$data$observed]
  m_index <- spint$data$index[!spint$data$observed]
  data_o <- spint$data[o_index, , drop = FALSE]

  # run the linear model and get residuals
  # creating the model frame
  stmodel_frame <- model.frame(formula, data_o,
                               na.action = stats::na.omit)

  # creating the fixed design matrix
  xo <- model.matrix(formula, stmodel_frame)

  # creating the response
  yo <- model.response(stmodel_frame)

  lmod_r <- (yo - xo %*% (chol2inv(chol(t(xo) %*% xo)) %*% (t(xo) %*% yo)))
  lmod_s2 <- sum(lmod_r^2) / (nrow(xo) - ncol(xo))

  # make the semivariogram
  sv <- st_empsv(response = lmod_r, xcoord = data_o[[xcoord]], ycoord = data_o[[ycoord]],
                   tcoord = data_o[[tcoord]], ...)

  # transform the initial values
  if (missing(max_v)){
    max_v <- 4 * lmod_s2
  }
  if (missing(max_srange)){
    max_srange <- 4 * max(spint$h_s)
  }
  if (missing(max_trange)){
    max_trange <- 4 * max(spint$h_t)
  }


  if (chol) {
    f_s <- h_make(data_o[[xcoord]], data_o[[ycoord]], ...)
    f_t <- h_make(data_o[[tcoord]], ...)
  } else {
    f_s <- NULL
    f_t <- NULL
  }

  est_object <- structure(list(initial = initial, xo = xo, yo = yo,
                               co = NULL, sv = sv, sp_cor = sp_cor,
                               t_cor = t_cor, max_srange = max_srange,
                               max_trange = max_trange, max_v = max_v, weights = weights,
                               h_s = spint$h_s, h_t = spint$h_t, f_s = f_s, f_t = f_t,
                               chol = chol, diag_tol = diag_tol, logdet = FALSE,
                               o_index = o_index, m_index = m_index), class = stcov)
  est_object$xyc_o <- cbind(est_object$xo, est_object$yo, est_object$co)
  est_object$lo_initial <- r2p_sv(est_object)


  output <- covest_sv(est_object) # estimate the parameters
  est_object$covparams <- output$par


  if (chol){
    est_object$rf_s <-  r_make(h = est_object$f_s, range = est_object$covparams[["s_range"]],
                    structure = est_object$sp_cor)
    est_object$rf_t <- r_make(h = est_object$f_t, range = est_object$covparams[["t_range"]],
                   structure = est_object$t_cor)
  } else {
    est_object$r_s <- r_make(h = est_object$h_s, range = est_object$covparams[["s_range"]],
                  structure = est_object$sp_cor)
    est_object$r_t <- r_make(h = est_object$h_t, range = est_object$covparams[["t_range"]],
           structure = est_object$t_cor)
  }
  inverse <- invert(invert_object = est_object)
  est_object$siginv_o <- inverse$siginv_o
  beta <- betaest(est_object)

  return(list(covparams = est_object$covparams, beta = beta$beta_hat, cov_beta = beta$cov_beta_hat))
}


