st_empsv <- function(response, xcoord, ycoord = NULL, tcoord,
                     n_sp_lag, n_t_lag, sp_max = NULL, t_max = NULL,
                     sp_dist = "euclidean", t_dist = "euclidean"){
  h_spatial <- h_make(coord1 = xcoord, coord2 = , distmetric = sp_dist)
  h_temporal <- h_make(coord1 = tcoord, distmetric = t_dist)
  sqdif_response <- outer(response, response, sqr_dif)

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
  output$t_lower <- pmax(0, t_lags_lower)
  output$t_upper <- t_lags_upper
  output$sp_lower <- pmax(0, sp_lags_lower)
  output$sp_upper <- sp_lags_upper
  return(output)
  #return(output[output[["n"]] > 0, ])
}

compute_sv <- function(sp_lag_lower, sp_lag_upper, t_lag_lower, t_lag_upper,
                       sp_h, t_h, h_spatial, h_temporal, sqdif_response){
  sqdifs <- sqdif_response[(h_spatial > sp_lag_lower) &
                             (h_spatial <= sp_lag_upper) &
                             (h_temporal > t_lag_lower) &
                             (h_temporal <= t_lag_upper)]
  n_sqdifs <- length(sqdifs)
  mean_sqdifs <- mean(sqdifs)/2
  return(data.frame(n = n_sqdifs, mean_sqdifs = mean_sqdifs, sp_h = sp_h, t_h = t_h))
}



x = rep(1:4, times = 4)
t = rep(1:4, each = 4)
response = rnorm(16)
h_spatial <- h_make(x)
h_temporal <- h_make(t)
sqdif_response <- outer(response, response, sqr_dif)
