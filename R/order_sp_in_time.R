order_sp_in_time <- function(data, xcoord, ycoord = NULL, tcoord, ...){
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
  return(list(data = data, h_s = h_s, h_t = h_t, us = us, ut = ut))
}
