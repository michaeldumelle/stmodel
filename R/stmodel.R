stmodel <- function(formula, stcov, estmethod, spcor, tcor,
                    xcoord, ycoord, tcoord, data,
                    initial = NULL, chol = FALSE, logdet = TRUE, diagtol = 1e-4, ...){

  if (!chol) {
  ordered <- order_sp_in_time(data = data, xcoord = xcoord,
                               ycoord = ycoord, tcoord = tcoord, ...)
  data_ord <- ordered$data
  o_index <- data_ord$index[data_ord$observed]
  m_index <- data_ord$index[!data_ord$observed]

  data_o <- data_ord[o_index, , drop = FALSE]
  } else {
    ordered <- list(data = data, h_s = h_make(data[[xcoord]], data[[ycoord]], ...),
                    h_t = h_make(data[[tcoord]], ...))
    data_o <- ordered$data
  }
  # creating the model frame
  stmodel_frame <- model.frame(formula, data_o,
                               na.action = stats::na.omit)

  # creating the fixed design matrix
  xo <- model.matrix(formula, stmodel_frame)

  # creating the response
  yo <- model.response(stmodel_frame)



  # make a prepare function that uses estmethod and a switch call

  construct_stcov <- function(stcov) {
    structure(list(xo = xo, yo = yo,
                   h_s = ordered$h_s, h_t = ordered$h_t,
                   spcor = spcor, tcor = tcor,
                   logdet = logdet, diagtol = diagtol,
                   chol = chol, estmethod = estmethod),
              class = stcov)
  }

  object <- construct_stcov(stcov)
  return(object)

}
