order_spint <- function(data, xcoord, ycoord = NULL, tcoord, chol = FALSE, h_options){


  # making the temporal ordering
  ## saving the unique temporal coordinates
  key_t <- unique(data[, tcoord, drop = FALSE])
  ## putting them in order
  key_t <- key_t[order(key_t[[tcoord]]), , drop = FALSE]
  ## making the key distance matrix
  h_t_small <- make_h(coord1 = key_t[[tcoord]], distmetric = h_options$h_t_distmetric)
  ## number of unique temporal locations
  n_t <- nrow(key_t)
  ## index for unique temporal locations
  key_t$tindex <- seq.int(1, n_t)

  ## making the
  key_s <- unique(data[, c(xcoord, ycoord), drop = FALSE])
  if (is.null(ycoord)) {
    key_s <- key_s[order(key_s[[xcoord]]), , drop = FALSE]
    h_s_small <- make_h(coord1 = key_s[[xcoord]], dismetric = h_options$h_s_distmetric)
  } else {
    key_s <- key_s[order(key_s[[ycoord]], key_s[[xcoord]]), , drop = FALSE]
    h_s_small <- make_h(coord1 = key_s[[xcoord]], coord2 = key_s[[ycoord]], distmetric = h_options$h_s_distmetric)
  }
  n_s <- nrow(key_s)
  key_s$sindex <- seq.int(1, n_s)
  data <- merge(merge(data, key_s), key_t)
  full_grid <- expand.grid(sindex = key_s$sindex, tindex = key_t$tindex)
  full_grid$index <- seq.int(1, n_t * n_s)
  data <- merge(full_grid, data, all = TRUE)
  data <- data[order(data$index), , drop = FALSE]
  data$observed <- !(is.na(data[[tcoord]]) & is.na(data[[xcoord]]))


  # find missing and observed indices
  o_index <- data$index[data$observed]
  m_index <- data$index[!data$observed]

  # create a subsetted data frame of the observed values
  ordered_data_o <- data[o_index, , drop = FALSE]

  # setting the cholesky distances matrices or NULL
  if (chol) {
    if (is.null(ycoord)){
      h_s_large <- make_h(coord1 = ordered_data_o[[xcoord]], distmetric = h_options$h_s_distmetric)
      h_t_large <- make_h(coord1 = ordered_data_o[[tcoord]], distmetric = h_options$h_t_distmetric)
    } else{
      h_s_large <- make_h(coord1 = ordered_data_o[[xcoord]], coord2 = ordered_data_o[[ycoord]], distmetric = h_options$h_s_distmetric)
      h_t_large <- make_h(coord1 = ordered_data_o[[tcoord]], distmetric = h_options$h_t_distmetric)
    }
  } else {
    h_s_large <- NULL
    h_t_large <- NULL
  }
  #raw_data is saved becakey_se the data merging mixes up the indices
  return(list(ordered_data_dense = data, ordered_data_o = ordered_data_o,
              h_s_small = h_s_small, h_t_small = h_t_small, n_s = n_s, n_t = n_t, o_index = o_index, m_index = m_index,
              h_s_large = h_s_large, h_t_large = h_t_large, key_s = key_s, key_t = key_t))
}
