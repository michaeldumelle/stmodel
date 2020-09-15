invert <- function(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de, st_ie,
                   o_index, x_o, y_o, cov_o, diag_tol = 1e-4, st_cov) {
  switch(st_cov,
         "productsum" = invert_productsum(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de, st_ie,
                                          o_index, x_o, y_o, cov_o, diag_tol),
         "product" = invert_product(),
         "sum_with_error" = invert_sum_with_error(),
         stop("choose a valid spatio-temporal covariance structure"))
}

invert_productsum <- function(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de, st_ie,
                              o_index, x_o, y_o, cov_o, diag_tol) {
  full_index <- seq(1, ncol(r_s) * ncol(r_t))
  m_index <- full_index[-o_index]
  # stegle step

  ## adding diagonal tolerances for invertibility stability
  diag(r_s) <- diag(r_s) + diag_tol
  diag(r_t) <- diag(r_t) + diag_tol

  ## finding eigendecompositions
  r_s_eigen <- eigen(r_s)
  r_t_eigen <- eigen(r_t)

  ## creating w matrix
  w <- kronecker(r_t_eigen$vectors, r_s_eigen$vectors)
  vinvroot <- 1 / sqrt(st_de * kronecker(r_t_eigen$values, r_s_eigen$values) + st_ie)
  w_vinvroot <- w %*% diag(vinvroot)






}




