# invert_chol <- function(f_s, f_t, s_de, s_ie, t_de, t_ie, st_de, st_ie,
#                         xyc_o, diag_tol = 1e-4,
#                         st_cov = c("productsum", "product", "sum_with_error"),
#                         log_determinant = TRUE) {
#   st_cov <- match.arg(st_cov)
#   sigma <- switch(st_cov,
#          "productsum" = invert_chol_productsum(f_s, f_t, s_de, s_ie, t_de, t_ie, st_de, st_ie,
#                                           xyc_o, diag_tol, log_determinant),
#          "product" = invert_product(),
#          "sum_with_error" = invert_chol_sum_with_error(f_s, f_t, s_de, s_ie, t_de, t_ie, st_ie = st_ie,
#                                                   xyc_o = xyc_o, diag_tol = diag_tol,
#                                                   log_determinant = log_determinant),
#          stop("choose a valid spatio-temporal covariance structure"))
#   chol_sigma <- chol(sigma)
#   siginv <- chol2inv(chol_sigma)
#   siginv_o <- siginv %*% xyc_o
#   if (log_determinant){
#     logdet <- 2 * sum(log(diag(chol_sigma)))
#   }
#   output <- list(siginv_o = siginv_o, logdet = logdet)
#   output_non_null <- output[!unlist(lapply(output, is.null))]
#   return(output_non_null)
# }

invert_chol_productsum <- function(f_s, f_t, s_de, s_ie, t_de, t_ie, st_de, st_ie,
                       xyc_o, diag_tol, log_determinant) {
  cov_s <- s_de * f_s + s_ie * (f_s == 1)
  cov_t <- t_de * f_t + t_ie * (f_t == 1)
  cov_st <- st_de * f_t * f_s
  diag(cov_st) <- diag(cov_st) + st_ie
  sigma <- cov_s + cov_t + cov_st
  chol_sigma <- chol(sigma)
  siginv <- chol2inv(chol_sigma)
  siginv_o <- siginv %*% xyc_o
  if (logdet){
    logdet <- 2 * sum(log(diag(chol_sigma)))
  } else {
    logdet <- NULL
  }
  output <- list(siginv_o = siginv_o, logdet = logdet)
  output_non_null <- output[!unlist(lapply(output, is.null))]
  return(output_non_null)
}

invert_chol_sum_with_error <- function(f_s, f_t, s_de, s_ie, t_de, t_ie, st_ie,
                                       xyc_o, diag_tol,
                                       log_determinant) {
  cov_s <- s_de * f_s + s_ie * (f_s == 1)
  cov_t <- t_de * f_t + t_ie * (f_t == 1)
  cov_st <- diag(st_ie, nrow = nrow(f_s), ncol = ncol(f_t))
  sigma <- cov_s + cov_t + cov_st
  chol_sigma <- chol(sigma)
  siginv <- chol2inv(chol_sigma)
  siginv_o <- siginv %*% xyc_o
  if (logdet){
    logdet <- 2 * sum(log(diag(chol_sigma)))
  } else {
    logdet <- NULL
  }
  output <- list(siginv_o = siginv_o, logdet = logdet)
  output_non_null <- output[!unlist(lapply(output, is.null))]
  return(output_non_null)
}

invert_chol_separable <- function(f_s, f_t, st_de, st_ie,
                                  xyc_o, diag_tol, log_determinant){
  sigma <- st_de * f_t * f_s
  chol_sigma <- chol(sigma)
  siginv <- chol2inv(chol_sigma)
  siginv_o <- siginv %*% xyc_o
  if (logdet){
    logdet <- 2 * sum(log(diag(chol_sigma)))
  } else {
    logdet <- NULL
  }
  output <- list(siginv_o = siginv_o, logdet = logdet)
  output_non_null <- output[!unlist(lapply(output, is.null))]
  return(output_non_null)
}
