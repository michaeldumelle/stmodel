stempsv <- function(response, xcoord, ycoord = NULL, tcoord, h_options,
                    h_response = NULL, h_s_large = NULL, h_t_large = NULL,
                    stempsv_options = NULL){
  if (is.null(stempsv_options)){
    stempsv_options <- list(n_s_lag = 16, n_t_lag = 16, h_s_max = NULL, h_t_max = NULL)
  }

  if (missing(h_options)){
    h_options = list(h_t_distmetric = "euclidean",
                     h_s_distmetric = "euclidean")
  }

  # make this take distance matrices or coordinates
  ## h_response, h_s_large, h_t_large
  if (is.null(h_s_large)) {
    h_s_large <- make_h(coord1 = xcoord, coord2 = ycoord, distmetric = h_options$h_s_distmetric)
  }
  if (is.null(h_t_large)){
    h_t_large <- make_h(coord1 = tcoord, distmetric = h_options$h_t_distmetric)
  }
  if (is.null(h_response)){
    h_response <- make_h(coord1 = response, distmetric = "euclidean")^2
  }
  h_s_large <- h_s_large[upper.tri(h_s_large, diag = F)]
  h_t_large <- h_t_large[upper.tri(h_t_large, diag = F)]
  h_response <- h_response[upper.tri(h_response, diag = F)]

  if (is.null(stempsv_options$h_s_max)) {
    stempsv_options$h_s_max <- max(h_s_large) / 2
  }

  if (is.null(stempsv_options$h_t_max)) {
    stempsv_options$h_t_max <- max(h_t_large) / 2
  }

  s_lags_upper <- seq(0, stempsv_options$h_s_max, length.out = stempsv_options$n_s_lag)
  s_lags_lower <- c(-0.1, s_lags_upper[-stempsv_options$n_s_lag])
  s_lags_upper <- rep(s_lags_upper, times = stempsv_options$n_t_lag)
  s_lags_lower <- rep(s_lags_lower, times = stempsv_options$n_t_lag)
  h_s_mid <- pmax(0, rowMeans(cbind(s_lags_lower, s_lags_upper)))

  t_lags_upper <- rep(seq(0, stempsv_options$h_t_max, length.out = stempsv_options$n_t_lag))
  t_lags_lower <- rep(c(-0.1, t_lags_upper[-stempsv_options$n_t_lag]))
  t_lags_upper <- rep(t_lags_upper, each = stempsv_options$n_s_lag)
  t_lags_lower <- rep(t_lags_lower, each = stempsv_options$n_s_lag)
  h_t_mid <- pmax(0, rowMeans(cbind(t_lags_lower, t_lags_upper)))

  output <- pmap_dfr(list(s_lag_lower = s_lags_lower, s_lag_upper = s_lags_upper,
                          t_lag_lower = t_lags_lower, t_lag_upper = t_lags_upper,
                          h_s_mid = h_s_mid, h_t_mid = h_t_mid),
                     .f = compute_stempsv, h_s_large, h_t_large, h_response)
  return_output <- output[output[["n"]] > 0, , drop = FALSE]

  if (nrow(return_output) > 0) {
    return(return_output)
  } else {
    stop("No semivariogram bins meet distance requirements: Choose larger values for h_s_max and h_t_max")
  }
}

compute_stempsv <- function(s_lag_lower, s_lag_upper, t_lag_lower, t_lag_upper,
                       h_s_mid, h_t_mid, h_s_large, h_t_large, h_response){
  h_s <- (h_s_large > s_lag_lower) & (h_s_large <= s_lag_upper)
  h_s_avg <- mean(h_s_large[h_s])
  h_t <- (h_t_large > t_lag_lower) & (h_t_large <= t_lag_upper)
  h_t_avg <- mean(h_t_large[h_t])
  sqrdifs <- h_response[h_s & h_t]
  # could multiply by 2 here to match the pairs but doesnt matter for optimization
  n_gammahat <- length(sqrdifs)
  gammahat <- mean(sqrdifs)/2
  return(data.frame(n = n_gammahat, gammahat = gammahat, h_s_mid = h_s_mid,
                    h_s_avg = h_s_avg, h_t_mid = h_t_mid, h_t_avg = h_t_avg))
}
