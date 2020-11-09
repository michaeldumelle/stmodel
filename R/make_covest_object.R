# make_covest_object <- function(initial = NULL, estmethod, stcov, sp_cor, t_cor, data_object,
#                                weights = NULL, chol = NULL, condition = NULL,
#                                max_options = NULL, optim_options = NULL, stempsv_options = NULL){
#
#   if (is.null(optim_options)){
#     if (estmethod == "reml"){
#       control = list(reltol = 1e-4, maxit = 2000)
#     } else {
#       control = list(reltol = 1e-8, maxit = 10000)
#     }
#     optim_options <- list(method = "Nelder-Mead", control = control)
#   }
#   if (is.null(max_options)){
#     max_options <- list(max_v = NULL, max_s_range = NULL, max_t_range = NULL)
#   }
#
#   if (is.null(stempsv_options)){
#     if (estmethod == "svwls"){
#       stempsv_options <- list(n_s_lag = 16, n_t_lag = 16, h_s_max = NULL, h_t_max = NULL)
#     } else {
#       stempsv_options <- NULL
#     }
#   }
#
#   # find the linear model residuals
#   lmod_r <- as.vector((data_object$ordered_yo - data_object$ordered_xo %*%
#                          (chol2inv(chol(t(data_object$ordered_xo) %*% data_object$ordered_xo)) %*%
#                             (t(data_object$ordered_xo) %*% data_object$ordered_yo))))
#
#   if (stcov == "svwls" | is.null(initial)){
#     if (is.null(data_object$h_s_large) || is.null(data_object$h_t_large)){
#       stop("h_s_large and h_t_large must be non NULL in data_object: if using stlmm, set h_options$h_large = TRUE")
#     }
#     h_response <- make_h(coord1 = lmod_r, distmetric = "euclidean")^2
#     stempsv <- stempsv(h_response = h_response, h_s_large = data_object$h_s_large,
#                        h_t_large = data_object$h_t_large, stempsv_options = stempsv_options)
#
#   }
#
#   if (is.null(initial)){
#     varinitial <- get_varinitial(stempsv = stempsv, stcov = stcov)
#     initial <- initial(s_de = varinitial[["s_de"]],
#                        s_ie = varinitial[["s_ie"]],
#                        t_de = varinitial[["t_de"]],
#                        t_ie = varinitial[["t_ie"]],
#                        st_de = varinitial[["st_de"]],
#                        st_ie = varinitial[["st_ie"]],
#                        v_s = varinitial[["v_s"]],
#                        v_t = varinitial[["v_t"]],
#                        s_range = max(data_object$h_s_small), # chosen so that correlation is approx zero at boundary
#                        t_range = max(data_object$h_t_small), # chosen so that correlation is approx zero at boundary
#                        estmethod = estmethod,
#                        stcov = stcov)
#   }
#
#   # provide default value for the maximum possible variance
#   if (is.null(max_options$max_v)){
#     # find the sample variance
#     s2 <- sum(initial[!(names(initial) %in% c("s_range", "t_range"))])
#     max_options$max_v <- 4 * s2
#   }
#   # provide default value for the maximum possible spatial range
#   if (is.null(max_options$max_s_range)){
#     max_options$max_s_range <- 4 * max(data_object$h_s_small)
#   }
#   # provide default value for the maximum possible temporal range
#   if (is.null(max_options$max_t_range)){
#     max_options$max_t_range <- 4 * max(data_object$h_t_small)
#   }
#
#   if (!identical(class(initial), c(estmethod, stcov))){
#     stop("class of initial value parameter vector must match estmethod and stcov used in stlmm")
#   }
#
#   initial_plo <- r2plo(covparam_object = initial, max_options = max_options)
#
#   if (estmethod == "svwls"){
#     weights <- weights
#     chol <- NULL
#     logdet <- NULL
#     condition <- NULL
#   } else {
#     stempsv <- NULL
#     weights <- NULL
#     chol <- chol
#     logdet <- TRUE
#     condition <- condition
#   }
#
#   covest_object <- structure(list(chol = chol, condition = condition,
#                                   initial = initial, initial_plo = initial_plo,
#                                   logdet = logdet, max_options = max_options,
#                                   optim_options = optim_options,
#                                   sp_cor = sp_cor, stempsv = stempsv,
#                                   stempsv_options = stempsv_options, t_cor = t_cor,
#                                   weights = weights), class = class(initial))
#   return(covest_object)
# }















make_covest_object <- function(initial = NULL, estmethod, stcov, sp_cor, t_cor, data_object,
                                weights = NULL, chol = NULL, condition = NULL,
                                max_options = NULL, optim_options = NULL, stempsv_options = NULL){

  if (is.null(optim_options)){
    if (estmethod == "reml"){
      control = list(reltol = 1e-4, maxit = 2000)
    } else {
      control = list(reltol = 1e-8, maxit = 10000)
    }
    optim_options <- list(method = "Nelder-Mead", control = control)
  }
  if (is.null(max_options)){
    max_options <- list(max_v = NULL, max_s_range = NULL, max_t_range = NULL)
  }

  if (is.null(stempsv_options)){
    if (estmethod == "svwls"){
      stempsv_options <- list(n_s_lag = 16, n_t_lag = 16, h_s_max = NULL, h_t_max = NULL)
    } else {
      stempsv_options <- NULL
    }
  }

  # find the linear model residuals
  lmod_r <- as.vector((data_object$ordered_yo - data_object$ordered_xo %*%
                         (chol2inv(chol(t(data_object$ordered_xo) %*% data_object$ordered_xo)) %*%
                            (t(data_object$ordered_xo) %*% data_object$ordered_yo))))

  # find the sample variance
  lmod_s2 <- sum(lmod_r^2) / (nrow(data_object$ordered_xo) - ncol(data_object$ordered_xo))

  if (is.null(data_object$h_s_large) || is.null(data_object$h_t_large)){
    stop("h_s_large and h_t_large must be non NULL in data_object: if using stlmm, set h_options$h_large = TRUE")
  }

  if (is.null(initial)){
    initial <- initial(s_de = 1, s_ie = 1, t_de = 1,
                          t_ie = 1, st_de = 1, st_ie = 1,
                          v_s = 0.5, v_t = 0.5,
                          s_range = max(data_object$h_s_small) / 2, # 2 chosen so that it is half the max observed distance
                          t_range = max(data_object$h_t_small) / 2, # 2 chosen so that it is half the max observed distance
                          estmethod = estmethod, stcov = stcov)
    vparm_names <- c("s_de", "s_ie", "t_de", "t_ie",
                     "st_de", "st_ie")
    numparams <- sum(vparm_names %in% names(initial))
    initial[names(initial) %in% vparm_names] <- lmod_s2 / numparams
  }

  # provide default value for the maximum possible variance
  if (is.null(max_options$max_v)){
    # find the sample variance
    s2 <- sum(initial[!(names(initial) %in% c("s_range", "t_range"))])
    max_options$max_v <- 4 * s2
  }
  # provide default value for the maximum possible spatial range
  if (is.null(max_options$max_s_range)){
    max_options$max_s_range <- 4 * max(data_object$h_s_small)
  }
  # provide default value for the maximum possible temporal range
  if (is.null(max_options$max_t_range)){
    max_options$max_t_range <- 4 * max(data_object$h_t_small)
  }

  if (!identical(class(initial), c(estmethod, stcov))){
    stop("class of initial value parameter vector must match estmethod and stcov used in stlmm")
  }

  initial_plo <- r2plo(covparam_object = initial, max_options = max_options)

  if (estmethod == "svwls"){
    h_response <- make_h(coord1 = lmod_r, distmetric = "euclidean")^2
    stempsv <- stempsv(h_response = h_response, h_s_large = data_object$h_s_large,
                       h_t_large = data_object$h_t_large, stempsv_options = stempsv_options)
    weights <- weights
    chol <- NULL
    logdet <- NULL
    condition <- NULL
  } else {
    stempsv <- NULL
    weights <- NULL
    chol <- chol
    logdet <- TRUE
    condition <- condition
  }

  covest_object <- structure(list(chol = chol, condition = condition,
                                   initial = initial, initial_plo = initial_plo,
                                   logdet = logdet, max_options = max_options,
                                   optim_options = optim_options,
                                   sp_cor = sp_cor, stempsv = stempsv,
                                   stempsv_options = stempsv_options, t_cor = t_cor,
                                   weights = weights), class = class(initial))
  return(covest_object)
}





# if (is.null(initial)){
#   initial <- initial(s_de = 1, s_ie = 1, t_de = 1,
#                      t_ie = 1, st_de = 1, st_ie = 1,
#                      v_s = 0.5, v_t = 0.5,
#                      s_range = max(data_object$h_s_small) / 2, # 2 chosen so that it is half the max observed distance
#                      t_range = max(data_object$h_t_small) / 2, # 2 chosen so that it is half the max observed distance
#                      estmethod = estmethod, stcov = stcov)
#   vparm_names <- c("s_de", "s_ie", "t_de", "t_ie",
#                    "st_de", "st_ie")
#   numparams <- sum(vparm_names %in% names(initial))
#   initial[names(initial) %in% vparm_names] <- lmod_s2 / numparams
# }
