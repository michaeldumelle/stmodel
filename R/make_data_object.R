make_data_object <- function(formula, xcoord, ycoord = NULL, tcoord, stcov, data,
                             estmethod, sp_cor, t_cor, weights, initial = NULL, chol,
                             diag_tol = 1e-4, max_v = NULL, max_srange = NULL,
                             max_trange = NULL, ...){

  stmodel_frame <- model.frame(formula = formula, data,
                                    na.action = stats::na.omit)

  # creating the fixed design matrix
  xo <- model.matrix(formula, stmodel_frame)

  # creating the response column
  yo <- model.response(stmodel_frame)

  # order the data by space within time
  spint <- order_spint(data = data, xcoord = xcoord,
                       ycoord = ycoord, tcoord = tcoord, chol = chol)



  # create the model frame using the provided formula
  ordered_stmodel_frame <- model.frame(formula, spint$ordered_data_o,
                               na.action = stats::na.omit)

  # creating the fixed design matrix
  ordered_xo <- model.matrix(formula, ordered_stmodel_frame)

  # creating the response column
  ordered_yo <- model.response(ordered_stmodel_frame)

  # find the linear model residuals
  lmod_r <- as.vector((ordered_yo - ordered_xo %*%
                         (chol2inv(chol(t(ordered_xo) %*% ordered_xo)) %*% (t(ordered_xo) %*% ordered_yo))))

  # find the sample variance
  lmod_s2 <- sum(lmod_r^2) / (nrow(ordered_xo) - ncol(ordered_xo))

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
                                    v_s = 0.5, v_t = 0.5,
                                    s_range = max_srange / 8, # 8 chosen so that it is half the max observed distance
                                    t_range = max_trange / 8, # 8 chosen so that it is half the max observed distance
                                    estmethod = estmethod, stcov = stcov)
    vparm_names <- c("s_de", "s_ie", "t_de", "t_ie",
                     "st_de", "st_ie")
    numparams <- sum(vparm_names %in% names(initial))
    initial[names(initial) %in% vparm_names] <- lmod_s2 / numparams
  }

  initial_plo <- r2plo(covparam_object = initial, max_v = max_v, max_srange = max_srange, max_trange = max_trange)

  # make the semivariogram
  if (estmethod == "svwls"){
  sv <- st_empsv(response = lmod_r, xcoord = spint$ordered_data_o[[xcoord]], ycoord = spint$ordered_data_o[[ycoord]],
                 tcoord = spint$ordered_data_o[[tcoord]], ...)
  } else {
    sv <- NULL
  }
  data_object <- structure(list(formula = formula, data = data, xo = xo, yo = yo,
                                xyc_o = cbind(xo, yo), stcov = stcov, estmethod = estmethod,
                                chol = chol, logdet = logdet, diag_tol = diag_tol,
                                ordered_data_dense = spint$ordered_data_dense,
                                ordered_data_o = spint$ordered_data_o,
                                h_s = spint$h_s, h_t = spint$h_t, o_index = spint$o_index,
                                m_index = spint$m_index,
                                sp_cor = sp_cor, t_cor = t_cor,
                                f_s = spint$f_s, f_t = spint$f_t, key_s = spint$key_s, key_t = spint$key_t,
                                ordered_xo = ordered_xo, ordered_yo = ordered_yo,
                                ordered_xyc_o = cbind(ordered_xo, ordered_yo),
                                lmod_s2 = lmod_s2, max_v = max_v, max_srange = max_srange,
                                max_trange = max_trange, initial = initial,
                                initial_plo = initial_plo, sv = sv, weights = weights),
                           class = c(estmethod, stcov))
  return(data_object)
}
