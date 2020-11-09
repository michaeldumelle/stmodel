predict.stlmm <- function(object,
                          newdata,
                          interval = c("none", "confidence", "prediction"),
                          se.fit = TRUE,
                          predcov = FALSE,
                          ...) {
  interval <- match.arg(interval)
  pred_output <- switch(interval,
                        "none" = predict.stlmm_none(object = object,
                                                    newdata = newdata,
                                                    ...),
                        "confidence" = predict.stlmm_confidence(object = object,
                                                                newdata = newdata,
                                                                se.fit = se.fit,
                                                                predcov = predcov,
                                                                ...),
                        "prediction" = predict.stlmm_prediction(object = object,
                                                                newdata = newdata,
                                                                se.fit = se.fit,
                                                                predcov = predcov,
                                                                ...),
                        stop("Must choose confidence or prediction interval"))
  pred_output$interval <- interval
  return(pred_output)
}


predict.stlmm_none <- function(object,
                               newdata,
                               ...) {
  newdata_stmodel_frame <- model.frame(object$formula[-2], newdata,
                                       na.action = stats::na.omit)
  # creating the fixed design matrix
  # -2 is to remove response from formula
  newdata_xo <- model.matrix(object$formula[-2], newdata_stmodel_frame)

  fit <- newdata_xo %*% object$Coefficients

  pred_output <- list(fit = fit)
  return(pred_output)
}

predict.stlmm_confidence <- function(object,
                                     newdata,
                                     se.fit,
                                     predcov,
                                     ...) {

  newdata_stmodel_frame <- model.frame(object$formula[-2], newdata,
                                       na.action = stats::na.omit)
  # creating the fixed design matrix
  # -2 is to remove response from formula
  newdata_xo <- model.matrix(object$formula[-2], newdata_stmodel_frame)

  fit <- newdata_xo %*% object$Coefficients

  if (predcov) {
    predcov <- newdata_xo %*% object$CovCoefficients %*% t(newdata_xo)
    if (se.fit) {
      se.fit <- sqrt(diag(predcov))
    } else {
      se.fit <- NULL
    }
  } else {
    predcov <- NULL
    if (se.fit) {
      npred <- nrow(newdata_xo)
      var.fit <- vapply(1:npred, function(x) {
        newdata_xo[x, , drop = FALSE] %*% object$CovCoefficients %*% t(newdata_xo[x, , drop = FALSE])
      }, double(1))
      se.fit <- sqrt(var.fit)
    } else {
      se.fit <- NULL
    }
  }

  pred_output <- list(fit = fit, se.fit = se.fit, predcov = predcov)
  pred_output <- pred_output[!vapply(pred_output, function(x) is.null(x), logical(1))]
  return(pred_output)
}

predict.stlmm_prediction <- function(object,
                                     newdata,
                                     se.fit,
                                     predcov,
                                     ...) {

  newdata_stmodel_frame <- model.frame(object$formula[-2], newdata,
                                       na.action = stats::na.omit)
  # creating the fixed design matrix
  # -2 is to remove response from formula
  newdata_xo <- model.matrix(object$formula[-2], newdata_stmodel_frame)

  # make spatial distance matrices
  if (is.null(object$coordnames$ycoord)){
    newdata_data_h_s_large <- make_newdata_h(newdata_coord1 = newdata[[object$coordnames$xcoord]],
                                             data_coord1 = object$coords$xcoord,
                                             distmetric = object$h_options$h_s_distmetric)
    newdata_h_s_large <- make_h(coord1 = newdata[[object$coordnames$xcoord]],
                                distmetric = object$h_options$h_s_distmetric)
  } else {
    newdata_data_h_s_large <- make_newdata_h(newdata_coord1 = newdata[[object$coordnames$xcoord]],
                                             data_coord1 = object$coords$xcoord,
                                             newdata_coord2 = newdata[[object$coordnames$ycoord]],
                                             data_coord2 = object$coords$ycoord,
                                             distmetric = object$h_options$h_s_distmetric)

    newdata_h_s_large <- make_h(coord1 = newdata[[object$coordnames$xcoord]],
                                coord2 = newdata[[object$coordnames$ycoord]],
                                distmetric = object$h_options$h_s_distmetric)

  }

  # make temporal distance matrices
  newdata_data_h_t_large <- make_newdata_h(newdata_coord1 = newdata[[object$coordnames$tcoord]],
                                           data_coord1 = object$coords$tcoord,
                                           distmetric = object$h_options$h_t_distmetric)
  newdata_h_t_large <- make_h(coord1 = newdata[[object$coordnames$tcoord]],
                              distmetric = object$h_options$h_t_distmetric)


  # make spatio-temopral covariance matrix
  newdata_data_stcovariance <- make_stcovariance(covparam_object = object$CovarianceParameters,
                                                 h_s_large = newdata_data_h_s_large,
                                                 h_t_large = newdata_data_h_t_large,
                                                 s_cor = object$CovarianceForms[["s_cor"]],
                                                 t_cor = object$CovarianceForms[["t_cor"]])


  # make the object to invert and multiply by covariance on the right
  invert_object <- make_invert_object(covparam_object = object$CovarianceParameters,
                                      chol = object$chol,
                                      condition = object$condition,
                                      co = t(newdata_data_stcovariance),
                                      h_s_small = object$data_object$h_s_small,
                                      h_t_small = object$data_object$h_t_small,
                                      h_s_large = object$data_object$h_s_large,
                                      h_t_large = object$data_object$h_t_large,
                                      logdet = FALSE,
                                      m_index = object$data_object$m_index,
                                      o_index = object$data_object$o_index,
                                      s_cor = object$CovarianceForms[["s_cor"]],
                                      t_cor = object$CovarianceForms[["t_cor"]],
                                      xo = NULL,
                                      yo = NULL)

  # compute the inverse
  invert_output <- invert(invert_object)

  invxo <- object$invert_output$sigmainv_o[, 1:(ncol(object$invert_output$sigmainv_o) - 1), drop = FALSE]
  newdata_invxo <- newdata_data_stcovariance %*% invxo


  fit <- newdata_xo %*% object$Coefficients +
    newdata_data_stcovariance %*% object$invert_output$sigmainv_o[, ncol(object$invert_output$sigmainv_o), drop = FALSE] -
    newdata_invxo %*% object$Coefficients



  if (predcov) {
    newdata_stcovariance <- make_stcovariance(covparam_object = object$CovarianceParameters,
                                              h_s_large = newdata_h_s_large,
                                              h_t_large = newdata_h_t_large,
                                              s_cor = object$CovarianceForms[["s_cor"]],
                                              t_cor = object$CovarianceForms[["t_cor"]])
    H <- newdata_xo - newdata_invxo

    predcov <- newdata_stcovariance - newdata_data_stcovariance %*% invert_output$sigmainv_o +
      H %*% object$CovCoefficients %*% t(H)

    if (se.fit) {
      se.fit <- sqrt(diag(predcov))
    } else {
      se.fit <- NULL
    }
  } else {
    predcov <- NULL
    if (se.fit) {
      vparm_names <- c("s_de", "s_ie", "t_de", "t_ie",
                       "st_de", "st_ie")
      varsum <- sum(object$CovarianceParameters[names(object$CovarianceParameters) %in% vparm_names])
      se.fit <- vapply(1:nrow(newdata_xo), function(x) {
        H <- newdata_xo[x, , drop = FALSE] - newdata_invxo[x, , drop = FALSE]
        predvar <- varsum - newdata_data_stcovariance[x, , drop = FALSE] %*% invert_output$sigmainv_o[, x, drop = FALSE] +
          H %*% object$CovCoefficients %*% t(H)
        se.fit <- sqrt(predvar)
      }, double(1))
    } else {
      se.fit <- NULL
    }
  }
  pred_output <- list(fit = fit, se.fit = se.fit, predcov = predcov)
  pred_output <- pred_output[!vapply(pred_output, function(x) is.null(x), logical(1))]
  return(pred_output)
}
