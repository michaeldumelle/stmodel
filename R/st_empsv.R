st_empsv <- function(response, xcoord, ycoord = NULL, tcoord,
                     n_sp_lag = 16, n_t_lag = 16, sp_max = NULL, t_max = NULL,
                     sp_dist = "euclidean", t_dist = "euclidean", r_diff = "euclidean", ...){
  h_spatial <- h_make(coord1 = xcoord, coord2 = ycoord, distmetric = sp_dist)
  h_spatial <- h_spatial[upper.tri(h_spatial, diag = F)]
  h_temporal <- h_make(coord1 = tcoord, distmetric = t_dist)
  h_temporal <- h_temporal[upper.tri(h_temporal, diag = F)]
  sqdif_response <- h_make(response, distmetric = r_diff)^2
  sqdif_response <- sqdif_response[upper.tri(sqdif_response, diag = F)]

  if (is.null(sp_max)) {
    sp_max <- max(h_spatial) / 2
  }

  if (is.null(t_max)) {
    t_max <- max(h_temporal) / 2
  }

  sp_lags_upper <- seq(0, sp_max, length.out = n_sp_lag)
  sp_lags_lower <- c(-0.1, sp_lags_upper[-n_sp_lag])
  sp_lags_upper <- rep(sp_lags_upper, times = n_t_lag)
  sp_lags_lower <- rep(sp_lags_lower, times = n_t_lag)
  sp_h_mid <- pmax(0, rowMeans(cbind(sp_lags_lower, sp_lags_upper)))

  t_lags_upper <- rep(seq(0, t_max, length.out = n_t_lag))
  t_lags_lower <- rep(c(-0.1, t_lags_upper[-n_t_lag]))
  t_lags_upper <- rep(t_lags_upper, each = n_sp_lag)
  t_lags_lower <- rep(t_lags_lower, each = n_sp_lag)
  t_h_mid <- pmax(0, rowMeans(cbind(t_lags_lower, t_lags_upper)))

  output <- pmap_dfr(list(sp_lag_lower = sp_lags_lower, sp_lag_upper = sp_lags_upper,
                          t_lag_lower = t_lags_lower, t_lag_upper = t_lags_upper,
                          sp_h = sp_h_mid, t_h = t_h_mid),
                     .f = compute_sv, h_spatial, h_temporal, sqdif_response)
  # output$t_lower <- pmax(0, t_lags_lower)
  # output$t_upper <- t_lags_upper
  # output$sp_lower <- pmax(0, sp_lags_lower)
  # output$sp_upper <- sp_lags_upper
  return_output <- output[output[["n"]] > 0, ]

  if (nrow(return_output) > 0) {
    return(return_output)
  } else {
    stop("No semivariogram bins meet distance requirements: Choose larger values for sp_max and t_max")
  }
}

compute_sv <- function(sp_lag_lower, sp_lag_upper, t_lag_lower, t_lag_upper,
                       sp_h, t_h, h_spatial, h_temporal, sqdif_response){
  hsp <- (h_spatial > sp_lag_lower) & (h_spatial <= sp_lag_upper)
  avg_hsp <- mean(h_spatial[hsp])
  tsp <- (h_temporal > t_lag_lower) & (h_temporal <= t_lag_upper)
  avg_tsp <- mean(h_temporal[tsp])
  sqdifs <- sqdif_response[hsp & tsp]
  # could multiply by 2 here to match the pairs but doesnt matter for optimization
  n_sqdifs <- length(sqdifs)
  mean_sqdifs <- mean(sqdifs)/2
  return(data.frame(n = n_sqdifs, mean_sqdifs = mean_sqdifs, sp_h = sp_h, avg_hsp = avg_hsp, t_h = t_h, avg_tsp = avg_tsp))
}


# library(gstat)
# library(spacetime)
# set.seed(1)
# x = rep(runif(35, max = 5), times = 30)
# y = rep(1, 35 * 30)
# s_index = rep(1:35, times = 30)
# t = rep(1:30, each = 35)
# response = rnorm(35 * 30)
# h_spatial <- h_make(x)
# h_temporal <- h_make(t)
# sqdif_response <- outer(response, response, sqr_dif)
# s_cutoff = max(h_make(x))/2
# t_cutoff = max(h_make(t))/2
#
#
# #creating a data frame of the observed coordinates / indices / response and arranging
#   #by space within time
#   semivario_data <- data.frame(x = x, y = y, s_index = s_index,
#                                t = t, response = response) %>% dplyr::arrange(t, s_index)
#
#   #creating ordered spatial coordinates
#   s_coords <- unique(semivario_data[, c("x", "y", "s_index")] %>% plyr::arrange(s_index) %>% dplyr::select(x, y))
#   #creating ordered temporal coordinates
#   t_coords <- unique(data.frame(t = semivario_data[, c("t")])  %>% dplyr::arrange(t))
#
#     #creating a spatial points data frame with coordinate information
#     sp_df <- SpatialPoints(coords = s_coords)
#     #creating a time series object with dummy dates spaced by one
#     t_df <- xts::xts(seq(1, nrow(t_coords), 1), order.by = as.Date('0001-01-01') + seq(1, nrow(t_coords), 1))
#     #creating a data frame of the response vector to be used later
#     resp <- data.frame(resp = semivario_data$response)
#     #creating a spatio-temporal index, storing as a matrix
#     index = as.matrix(semivario_data[c("s_index", "t")])
#     #creating a space time data frame that is not dense
#     stdf <- STSDF(sp_df, t_df, data = resp, index = index)
#     #creating a space time semivariogram from the space time data frame - omitting bins with zero observations
#     #renaming the default output to be more standard
#     #bug with new gstat, sp, or spactime - works on gstat_2.0-2, sp_1.3-1, spacetime_1.2-2
#     st_semivariogram <- variogramST(resp ~ 1, stdf, progress = F, cutoff = s_cutoff,
#                                     tlags = seq(from = 0, to = t_cutoff, by = 1), na.omit = T) %>% dplyr::rename(n = np)
#
#     test3 = st_empsv(response, xcoord = x, tcoord = t, n_sp_lag = 16, n_t_lag = floor(29/2) + 1, t_max = floor(29/2))
#
#     View(st_semivariogram)
#     View(test3)
# microbenchmark::microbenchmark(variogramST(resp ~ 1, stdf, progress = F, cutoff = s_cutoff,
#                                            tlags = seq(from = 0, to = t_cutoff, by = 1)), times = 20)
#
#
# microbenchmark::microbenchmark(st_empsv(response, xcoord = x, tcoord = t, n_sp_lag = 16, n_t_lag = 15, t_max = 14), times = 20)

