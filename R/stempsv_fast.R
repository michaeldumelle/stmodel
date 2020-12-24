#' Title
#'
#' @param response
#' @param xcoord
#' @param ycoord
#' @param tcoord
#' @param h_options
#' @param h_response
#' @param h_s_large
#' @param h_t_large
#' @param stempsv_options
#' @import stats
#' @importFrom purrr pmap_dfr
#' @return
#' @export
#'
#' @examples
stempsv_fast <- function(response,
                    xcoord,
                    ycoord = NULL,
                    tcoord,
                    h_options = NULL,
                    h_response = NULL,
                    h_s_large = NULL,
                    h_t_large = NULL,
                    stempsv_options = NULL
                    ) {

  # setting default options if none are given
  if (is.null(stempsv_options)) {
    stempsv_options <- list(n_s_lag = 16, n_t_lag = 16, h_s_max = NULL, h_t_max = NULL)
  }

  # setting default h options if none are given
  if (is.null(h_options)) {
    h_options = list(h_t_distmetric = "euclidean", h_s_distmetric = "euclidean")
  }

  # creating a large spatial distance matrix if not provided
  if (is.null(h_s_large)) {
    h_s_large <- make_h(coord1 = xcoord, coord2 = ycoord, distmetric = h_options$h_s_distmetric)
  }

  # creating a large temporal distance matrix if not provided
  if (is.null(h_t_large)){
    h_t_large <- make_h(coord1 = tcoord, distmetric = h_options$h_t_distmetric)
  }

  # creating a large squared response matrix if not provided
  if (is.null(h_response)) {
    h_response <- make_h(coord1 = response, distmetric = "euclidean")^2
  }

  # only storing the upper triangular portions of the spatial, temporal,
  # and squared distance matrices
  h_s_large <- h_s_large[upper.tri(h_s_large, diag = F)]
  h_t_large <- h_t_large[upper.tri(h_t_large, diag = F)]
  h_response <- h_response[upper.tri(h_response, diag = F)]

  # setting a max for the spatial disance if none is provided
  # (as the max distance over 2)
  if (is.null(stempsv_options$h_s_max)) {
    stempsv_options$h_s_max <- max(h_s_large) / 2
  }

  # setting a max for the temporal disance if none is provided
  # (as the max distance over 2)
  if (is.null(stempsv_options$h_t_max)) {
    stempsv_options$h_t_max <- max(h_t_large) / 2
  }

  hmax_index <- (h_s_large <= stempsv_options$h_s_max) & (h_t_large <= stempsv_options$h_t_max)
  h_s_large <- h_s_large[hmax_index]
  h_t_large <- h_t_large[hmax_index]
  h_response <- h_response[hmax_index]

  s_lags <- c(-0.1, seq(0, stempsv_options$h_s_max, length.out = stempsv_options$n_s_lag))
  h_s_index <- cut(h_s_large, s_lags, include.lowest = T)
  t_lags <- c(-0.1, seq(0, stempsv_options$h_t_max, length.out = stempsv_options$n_t_lag))
  h_t_index <- cut(h_t_large, t_lags, include.lowest = T)
  st_index <- interaction(h_s_index, h_t_index)


  h_s_avg <- tapply(h_s_large, st_index, mean)
  h_t_avg <- tapply(h_t_large, st_index, mean)
  gammahat <- tapply(h_response, st_index, FUN = function(x) mean(x) / 2)
  n <- tapply(h_response, st_index, length)

  return_output <- na.omit(data.frame(gammahat, n, h_s_avg, h_t_avg))
  attr(return_output, "na.action") <- NULL
  row.names(return_output) <- NULL
  return(return_output)
}
