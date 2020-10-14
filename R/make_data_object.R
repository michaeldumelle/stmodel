make_data_object <- function(formula, xcoord, ycoord, tcoord, data, h_options){

  if (is.null(h_options)){
    h_options = list(h_large = TRUE,
                     h_t_distmetric = "euclidean",
                     h_s_distmetric = "euclidean")
  }
  original_data <- data
  # making the original model frame
  original_stmodel_frame <- model.frame(formula, original_data,
                                        na.action = stats::na.omit)

  # creating the fixed design matrix
  original_xo <- model.matrix(formula, original_stmodel_frame)

  # creating the response column
  original_yo <- model.response(original_stmodel_frame)

  # order the data by space within time
  spint <- storder(data = original_data, xcoord = xcoord,
                   ycoord = ycoord, tcoord = tcoord, h_options = h_options)

  # create the model frame using the provided formula
  ordered_stmodel_frame <- model.frame(formula, spint$ordered_data_o,
                                       na.action = stats::na.omit)

  # creating the fixed design matrix
  ordered_xo <- model.matrix(formula, ordered_stmodel_frame)

  # creating the response column
  ordered_yo <- model.response(ordered_stmodel_frame)

  data_object <- list(formula = formula, original_data = original_data,
                                original_xo = original_xo, original_yo = original_yo,
                                ordered_data_dense = spint$ordered_data_dense,
                                ordered_data_o = spint$ordered_data_o,
                                h_s_small = spint$h_s_small, h_t_small = spint$h_t_small,
                                n_s = spint$n_s,n_t = spint$n_t,
                                o_index = spint$o_index, m_index = spint$m_index,
                                h_s_large = spint$h_s_large, h_t_large = spint$h_t_large,
                                key_s = spint$key_s, key_t = spint$key_t,
                                ordered_xo = ordered_xo, ordered_yo = ordered_yo)
  return(data_object)
}


#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# make_data_object <- function(formula, xcoord, ycoord, tcoord, stcov, estmethod, data, ...){
#   data_object <- switch(estmethod,
#          "svwls" = make_data_object_svwls(formula = formula, xcoord = xcoord,
#                                           ycoord = ycoord, tcoord = tcoord,
#                                           stcov = stcov, estmethod = estmethod,
#                                           data = data, ...),
#          stop("Need valid estimation method"))
#   return(data_object)
# }
#
#
#
# make_data_object_svwls <- function(formula, xcoord, ycoord, tcoord, stcov, estmethod, data,
#                              sp_cor, t_cor, chol, initial,
#                              max_v, max_s_range,
#                              max_t_range, weights, h_options, logdet, diag_tol, ...){
#
#
#   if (is.null(h_options)){
#     h_options = list(h_t_distmetric = "euclidean", h_s_distmetric = "euclidean")
#   }
#
#   original_data <- data
#   # making the original model frame
#   original_stmodel_frame <- model.frame(formula, original_data,
#                                     na.action = stats::na.omit)
#
#   # creating the fixed design matrix
#   original_xo <- model.matrix(formula, original_stmodel_frame)
#
#   # creating the response column
#   original_yo <- model.response(original_stmodel_frame)
#
#   # order the data by space within time
#   spint <- storder(data = original_data, xcoord = xcoord,
#                        ycoord = ycoord, tcoord = tcoord, chol = chol, h_options = h_options)
#
#
#
#   # create the model frame using the provided formula
#   ordered_stmodel_frame <- model.frame(formula, spint$ordered_data_o,
#                                na.action = stats::na.omit)
#
#   # creating the fixed design matrix
#   ordered_xo <- model.matrix(formula, ordered_stmodel_frame)
#
#   # creating the response column
#   ordered_yo <- model.response(ordered_stmodel_frame)
#
#   # find the linear model residuals
#   lmod_r <- as.vector((ordered_yo - ordered_xo %*%
#                          (chol2inv(chol(t(ordered_xo) %*% ordered_xo)) %*% (t(ordered_xo) %*% ordered_yo))))
#
#   # find the sample variance
#   lmod_s2 <- sum(lmod_r^2) / (nrow(ordered_xo) - ncol(ordered_xo))
#
#   # provide default value for the maximum possible variance
#   if (is.null(max_v)){
#     max_v <- 4 * lmod_s2
#   }
#   # provide default value for the maximum possible spatial range
#   if (is.null(max_s_range)){
#     max_s_range <- 4 * max(spint$h_s_small)
#   }
#   # provide default value for the maximum possible temporal range
#   if (is.null(max_t_range)){
#     max_t_range <- 4 * max(spint$h_t_small)
#   }
#
#   # setting initial values if there are none specified
#   if (is.null(initial)){
#     initial <- make_covparam_object(s_de = 1, s_ie = 1, t_de = 1,
#                                     t_ie = 1, st_de = 1, st_ie = 1,
#                                     v_s = 0.5, v_t = 0.5,
#                                     s_range = max_s_range / 8, # 8 chosen so that it is half the max observed distance
#                                     t_range = max_t_range / 8, # 8 chosen so that it is half the max observed distance
#                                     estmethod = estmethod, stcov = stcov)
#     vparm_names <- c("s_de", "s_ie", "t_de", "t_ie",
#                      "st_de", "st_ie")
#     numparams <- sum(vparm_names %in% names(initial))
#     initial[names(initial) %in% vparm_names] <- lmod_s2 / numparams
#   }
#
#   initial_plo <- r2plo(covparam_object = initial, max_v = max_v,
#                        max_s_range = max_s_range, max_t_range = max_t_range)
#
#
#   sv <- stempsv(response = lmod_r, xcoord = spint$ordered_data_o[[xcoord]], ycoord = spint$ordered_data_o[[ycoord]],
#                  tcoord = spint$ordered_data_o[[tcoord]], h_options = h_options)
#
#   data_object <- structure(list(formula = formula, original_data = original_data,
#                                 original_xo = original_xo, original_yo = original_yo,
#                                 original_xyc_o = cbind(original_xo, original_yo), stcov = stcov, estmethod = estmethod,
#                                 chol = chol, logdet = logdet, diag_tol = diag_tol,
#                                 ordered_data_dense = spint$ordered_data_dense,
#                                 ordered_data_o = spint$ordered_data_o,
#                                 h_s_small = spint$h_s_small, h_t_small = spint$h_t_small, n_s = spint$n_s,
#                                 n_t = spint$n_t, o_index = spint$o_index,
#                                 m_index = spint$m_index,
#                                 sp_cor = sp_cor, t_cor = t_cor,
#                                 h_s_large = spint$h_s_large, h_t_large = spint$h_t_large, key_s = spint$key_s, key_t = spint$key_t,
#                                 ordered_xo = ordered_xo, ordered_yo = ordered_yo,
#                                 ordered_xyc_o = cbind(ordered_xo, ordered_yo),
#                                 lmod_s2 = lmod_s2, max_v = max_v, max_s_range = max_s_range,
#                                 max_t_range = max_t_range, initial = initial,
#                                 initial_plo = initial_plo, sv = sv, weights = weights),
#                            class = c(estmethod, stcov))
#   return(data_object)
# }
#
#
#
#
#
# make_data_object_reml <- function(formula, xcoord, ycoord, tcoord, stcov, estmethod, data,
#                                    sp_cor, t_cor, chol, initial,
#                                    max_v, max_s_range,
#                                    max_t_range, h_options, ...){
#
#   if (is.null(h_options)){
#     h_options = list(h_t_distmetric = "euclidean", h_s_distmetric = "euclidean")
#   }
#
#   original_data <- data
#   # making the original model frame
#   original_stmodel_frame <- model.frame(formula = formula, data = original_data,
#                                         na.action = stats::na.omit)
#
#   # creating the fixed design matrix
#   original_xo <- model.matrix(formula, original_stmodel_frame)
#
#   # creating the response column
#   original_yo <- model.response(original_stmodel_frame)
#
#   # order the data by space within time
#   spint <- storder(data = original_data, xcoord = xcoord,
#                        ycoord = ycoord, tcoord = tcoord, chol = chol, h_options = h_options)
#
#
#
#   # create the model frame using the provided formula
#   ordered_stmodel_frame <- model.frame(formula, spint$ordered_data_o,
#                                        na.action = stats::na.omit)
#
#   # creating the fixed design matrix
#   ordered_xo <- model.matrix(formula, ordered_stmodel_frame)
#
#   # creating the response column
#   ordered_yo <- model.response(ordered_stmodel_frame)
#
#   # find the linear model residuals
#   lmod_r <- as.vector((ordered_yo - ordered_xo %*%
#                          (chol2inv(chol(t(ordered_xo) %*% ordered_xo)) %*% (t(ordered_xo) %*% ordered_yo))))
#
#   # find the sample variance
#   lmod_s2 <- sum(lmod_r^2) / (nrow(ordered_xo) - ncol(ordered_xo))
#
#   # provide default value for the maximum possible variance
#   if (is.null(max_v)){
#     max_v <- 4 * lmod_s2
#   }
#   # provide default value for the maximum possible spatial range
#   if (is.null(max_s_range)){
#     max_s_range <- 4 * max(spint$h_s_small)
#   }
#   # provide default value for the maximum possible temporal range
#   if (is.null(max_t_range)){
#     max_t_range <- 4 * max(spint$h_t_small)
#   }
#
#   # setting initial values if there are none specified
#   if (is.null(initial)){
#     initial <- make_covparam_object(s_de = 1, s_ie = 1, t_de = 1,
#                                     t_ie = 1, st_de = 1, st_ie = 1,
#                                     v_s = 0.5, v_t = 0.5,
#                                     s_range = max_s_range / 8, # 8 chosen so that it is half the max observed distance
#                                     t_range = max_t_range / 8, # 8 chosen so that it is half the max observed distance
#                                     estmethod = estmethod, stcov = stcov)
#     vparm_names <- c("s_de", "s_ie", "t_de", "t_ie",
#                      "st_de", "st_ie")
#     numparams <- sum(vparm_names %in% names(initial))
#     initial[names(initial) %in% vparm_names] <- lmod_s2 / numparams
#   }
#
#   initial_plo <- r2plo(covparam_object = initial, max_v = max_v,
#                        max_s_range = max_s_range, max_t_range = max_t_range)
#
#
#
#
#   data_object <- structure(list(formula = formula, original_data = original_data,
#                                 original_xo = original_xo, original_yo = original_yo,
#                                 original_xyc_o = cbind(original_xo, original_yo), stcov = stcov, estmethod = estmethod,
#                                 chol = chol, logdet = logdet, diag_tol = diag_tol,
#                                 ordered_data_dense = spint$ordered_data_dense,
#                                 ordered_data_o = spint$ordered_data_o,
#                                 h_s_small = spint$h_s_small, h_t_small = spint$h_t_small, n_s = spint$n_s,
#                                 n_t = spint$n_t, o_index = spint$o_index,
#                                 m_index = spint$m_index,
#                                 sp_cor = sp_cor, t_cor = t_cor,
#                                 h_s_large = spint$h_s_large, h_t_large = spint$h_t_large, key_s = spint$key_s, key_t = spint$key_t,
#                                 ordered_xo = ordered_xo, ordered_yo = ordered_yo,
#                                 ordered_xyc_o = cbind(ordered_xo, ordered_yo),
#                                 lmod_s2 = lmod_s2, max_v = max_v, max_s_range = max_s_range,
#                                 max_t_range = max_t_range, initial = initial,
#                                 initial_plo = initial_plo),
#                            class = c(estmethod, stcov))
#   return(data_object)
# }
