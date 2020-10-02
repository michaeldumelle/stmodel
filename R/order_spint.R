order_spint <- function(data, xcoord, ycoord = NULL, tcoord, chol = FALSE, ...){
  # make the ordered data frame - do this using single brackets

  # making the temporal ordering
  ut <- unique(data[, tcoord, drop = FALSE])
  ut <- ut[order(ut[[tcoord]]), , drop = FALSE]
  h_t <- h_make(ut[[tcoord]], ...)
  n_t <- nrow(ut)
  ut$tindex <- seq.int(1, n_t)

  us <- unique(data[, c(xcoord, ycoord), drop = FALSE])
  if (is.null(ycoord)) {
    us <- us[order(us[[xcoord]]), , drop = FALSE]
    h_s <- h_make(us[[xcoord]], ...)
  } else {
    us <- us[order(us[[ycoord]], us[[xcoord]]), , drop = FALSE]
    h_s <- h_make(us[[xcoord]], us[[ycoord]], ...)
  }
  n_s <- nrow(us)
  us$sindex <- seq.int(1, n_s)
  data <- merge(merge(data, us), ut)
  full_grid <- expand.grid(sindex = us$sindex, tindex = ut$tindex)
  full_grid$index <- seq.int(1, n_t * n_s)
  data <- merge(full_grid, data, all = TRUE)
  data <- data[order(data$index), , drop = FALSE]
  data$observed <- !(is.na(data[[tcoord]]) & is.na(data[[xcoord]]))


  # find missing and observed indices
  o_index <- data$index[data$observed]
  m_index <- data$index[!data$observed]

  # create a subsetted data frame of the observed values
  data_o <- data[o_index, , drop = FALSE]

  # setting the cholesky distances matrices or NULL
  if (chol) {
    f_s <- h_make(data_o[[xcoord]], data_o[[ycoord]], ...)
    f_t <- h_make(data_o[[tcoord]], ...)
  } else {
    f_s <- NULL
    f_t <- NULL
  }
  return(list(data = data, data_o = data_o, h_s = h_s, h_t = h_t, o_index = o_index, m_index = m_index,
              f_s = f_s, f_t = f_t, us = us, ut = ut))
}


# order_spint_svwls <- function(data, xcoord, ycoord = NULL, tcoord){
#
#   # making the temporal ordering
#   ut <- unique(data[, tcoord, drop = FALSE])
#   ut <- ut[order(ut[[tcoord]]), , drop = FALSE]
#   n_t <- nrow(ut)
#   ut$tindex <- seq.int(1, n_t)
#
#   us <- unique(data[, c(xcoord, ycoord), drop = FALSE])
#   if (is.null(ycoord)) {
#     us <- us[order(us[[xcoord]]), , drop = FALSE]
#   } else {
#     us <- us[order(us[[ycoord]], us[[xcoord]]), , drop = FALSE]
#   }
#   n_s <- nrow(us)
#   us$sindex <- seq.int(1, n_s)
#   data <- merge(merge(data, us), ut)
#   full_grid <- expand.grid(sindex = us$sindex, tindex = ut$tindex)
#   full_grid$index <- seq.int(1, n_t * n_s)
#   data <- merge(full_grid, data, all = TRUE)
#   data <- data[order(data$index), , drop = FALSE]
#   data$observed <- !(is.na(data[[tcoord]]) & is.na(data[[xcoord]]))
#   return(list(data = data))
# }
