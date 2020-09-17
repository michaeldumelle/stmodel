invert <- function(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de, st_ie,
                   o_index, xyc_o, diag_tol = 1e-4,
                   st_cov = c("productsum", "product", "sum_with_error"),
                   log_determinant = TRUE) {
  n_s <- ncol(r_s)
  n_t <- ncol(r_t)
  n_st <- n_s * n_t
  m_index <- full_index[-o_index]
  dense <- length(o_index) == n_st

  ## adding diagonal tolerances for invertibility stability
  diag(r_s) <- diag(r_s) + diag_tol
  diag(r_t) <- diag(r_t) + diag_tol
  s_ie <- s_ie + diag_tol
  t_ie <- t_ie + diag_tol
  st_ie <- st_ie + diag_tol

  st_cov <- match.arg(st_cov)

  switch(st_cov,
         "productsum" = invert_productsum(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de, st_ie,
                                          o_index, xyc_o, diag_tol),
         "product" = invert_product(),
         "sum_with_error" = invert_sum_with_error(),
         stop("choose a valid spatio-temporal covariance structure"))
}

invert_productsum <- function(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de, st_ie,
                              o_index, xyc_o, diag_tol) {

  ## finding eigendecompositions
  r_s_eigen <- eigen(r_s)
  r_t_eigen <- eigen(r_t)

  ## creating w matrix
  w <- kronecker(r_t_eigen$vectors, r_s_eigen$vectors)
  # creating v matrix
  v <- st_de * kronecker(r_t_eigen$values, r_s_eigen$values) + st_ie
  # creating inverse square root of v matrix
  vinvroot <- 1 / sqrt(v)
  # creating inverse of w * sqrt(v)
  w_vinvroot <- w %*% diag(vinvroot)
  #creating the transpose
  t_w_vinvroot <- t(w_vinvroot)


  # storing the output we will need for the iterative smw
  c_t <- chol(sigma_make(t_de, r_t, t_ie))
  c_s <- chol(sigma_make(s_de, r_t, s_ie))

  ist_zt <- w_vinvroot %*% multiply_z(t_w_vinvroot, "temporal", n_s)
            # st x st %*% (st x st * st x t) = st x t
  tr_ist_zt <- t(ist_zt)
  c_mt <- chol(chol2inv(c_t) + t(multiply_z(mx = t(ist_zt), z_type = "temporal", n_s = n_s)))
            # (t x t + tr(tr(st x t) * st x t) = t x t
  ic_mt <- chol2inv(c_mt)
            # t x t

  istpt_zs <- w_vinvroot %*% multiply_z(mx = t_w_vinvroot, z_type = "spatial", n_s = n_s) -
    ist_zt %*% (ic_mt %*% multiply_z(mx = tr_ist_zt, z_type = "spatial", n_s = n_s))
    # st x st * (st x st * st x s) - st x t (t x t * (tr(st x t) * st x s)) =
    # st x s - st x t * t x s = st x s
  tr_istpt_zs <- t(istpt_zs)
  c_ms <- chol(chol2inv(c_s) + t(multiply_z(mx = tr_istpt_zs, z_type = "spatial", n_s = n_s)))
    # s x s + tr(tr(st x s) * st x s) = s x s
  ic_ms <- chol2inv(c_ms)
    # s x s

  # now to implement algorithm

  if (dense) {
    siginv_o <- w_vinvroot %*% (t_w_vinvroot %*% xyc_o) -
    # st x st * (st x st * st x p) = st x p
      ist_zt %*% (ic_mt %*% (tr_ist_zt %*% xyc_o)) -
      # st x t * (t x t * (t x st * st x p)) = st x p
      istpt_zs %*% (ic_ms %*% (tr_istpt_zs %*% xyc_o))
      # st x s * (s x s * (s x st * st x p)) = st x p
  } else {
    d_oo <- w_vinvroot[o_index, o_index, drop = FALSE] %*% (t_w_vinvroot[o_index, o_index, drop = FALSE] %*% xyc_o) +
      w_vinvroot[o_index, m_index, drop = FALSE] %*% (t_w_vinvroot[m_index, o_index, drop = FALSE] %*% xyc_o) -
      ist_zt[o_index, , drop = FALSE] %*%
      (ic_mt %*% (tr_ist_zt[ , o_index, drop = FALSE] %*% xyc_o)) -
      istpt_zs[o_index, , drop = FALSE] %*%
      (ic_ms %*% (tr_istpt_zs[ , o_index, drop = FALSE] %*% xyc_o))

    d_om <- w_vinvroot[o_index, o_index, drop = FALSE] %*% t_w_vinvroot[o_index, m_index, drop = FALSE]  +
      w_vinvroot[o_index, m_index, drop = FALSE] %*% t_w_vinvroot[m_index, m_index, drop = FALSE] -
      ist_zt[o_index, , drop = FALSE] %*%
      (ic_mt %*% tr_ist_zt[ , m_index, drop = FALSE]) -
      istpt_zs[o_index, , drop = FALSE] %*%
      (ic_ms %*% tr_istpt_zs[ , m_index, drop = FALSE])

    d_mm <- w_vinvroot[m_index, o_index, drop = FALSE] %*% t_w_vinvroot[o_index, m_index, drop = FALSE]  +
      w_vinvroot[m_index, m_index, drop = FALSE] %*% t_w_vinvroot[m_index, m_index, drop = FALSE] -
      ist_zt[m_index, , drop = FALSE] %*%
      (ic_mt %*% tr_ist_zt[ , m_index, drop = FALSE]) -
      istpt_zs[m_index, , drop = FALSE] %*%
      (ic_ms %*% tr_istpt_zs[ , m_index, drop = FALSE])

    # return the correct object
    c_mm <- chol(d_mm)
    siginv_o <- d_oo - d_om %*% (chol2inv(c_mm) %*% (t(d_om) %*% xyc_o))

  }

  if (log_determinant) {
    logdet <- sum(log(v)) +
      2 * sum(log(diag(c_t))) + 2 * sum(log(diag(c_mt))) +
      2 * sum(log(diag(c_s))) + 2 * sum(log(diag(c_ms)))

    if (!dense)
      logdet <- logdet + 2 * sum(log(diag(c_mm)))
  }

  # # smw step - storage 1 right mult by temporal design matrix
  # inv_st_zt <- w_vinvroot %*% multiply_z(t_w_vinvroot, "temporal", n_s)
  # chol_inv_t <- chol(simga_make(t_de, r_t, t_ie + diag_tol, diag(n_t)))
  # chol_smw1_p2mid <- chol(chol2inv(chol_inv_t) + t(multiply_z(inv_st_zt, "temporal", n_s)))
  # smw1_p2 <- inv_st_zs %*% chol2inv(chol_smw1_p2mid) %*% t(inv_st_zt)
  #
  # chol_inv_s <- chol(sigma-make(s_de, r_s, s_ie + diag_tol, diag(n_s)))
  #
  #
  # ## smw step 1.5 - storage 1.5 right mult by fixed effects design matrix
  # smw2_p1 <- w_vinvroot %*% (t_w_vinvroot %*% xyc_o) -
  #   # st x st * (st x st * st x t) = st x t
  #   inv_st_zt %*% chol2inv(chol_smw1_p2mid) %*% (t(inv_st_zt) %*% xyc_o)
  #   # st x t * t x t * (t x st * st x t) = st x t * t x t * t x t = st x t
  #
  # inv_stpt_zs <- w_vinvroot %*% (t_w_vinvroot %*% multiply_z(t_w_vinvroot, "spatial", n_s)) -
  #   # st x st * (st x st * st x t) = st x t
  #   inv_st_zt %*% chol2inv(chol_smw1_p2mid) %*% (t(inv_st_zt) %*% multiply_z(t_w_vinvroot, "spatial", n_s))
  # # st x t * t x t * (t x st * st x t) = st x t * t x t * t x t = st x t
  #
  # chol_smw2_p2mid <- chol(chol2inv(chol_inv_s) + t(multiply_z(iv_stpt_zs, "spatial", n_s)))
  # smw2 <- smw2_p1 - inv_stpt_zs %*% chol2inv(chol_smw2_p2mid) %*%
  #   (t(inv_st_zt) %*% multiply_z(t_w_vinvroot, "spatial", n_s)) %*%
  #   xyc_o
  #
  # chol_inv_smw1_p2mid <- chol2inv(chol_smw1_p2mid)
  # chol_inv_smw2_p2mid <- chol2inv(chol_smw2_p2mid)
#
#   invdense_oo <- w_vinvroot[o_index, o_index, drop = FALSE] %*% (t_w_vinvroot[ob, ob, drop = FALSE] %*% xyc_o) +
#     w_vinvroot[o_index, m_index, drop = FALSE] %*% (t_w_vinvroot[mi, ob, drop = FALSE] %*% xyc_o) -
#     inv_st_zt[o_index, , drop = FALSE] %*%
#     (chol_inv_smw1_p2mid %*% (t(inv_st_zt[o_index, , drop = FALSE]) %*% xyc_o)) -
#     iv_stpt_zs[o_index, , drop = FALSE] %*%
#     (chol_inv_smw2_p2mid %*% (t(inv_stpt_zs[o_index, , drop = FALSE]) %*% xyc_o))
#
#   invdense_om <- w_vinvroot[o_index, o_index, drop = FALSE] %*% w_vinvroot[o_index, m_index, drop = FALSE] +
#     w_vinvroot[o_index, m_index, drop = FALSE] %*% w_vinvroot[m_index, m_index, drop = FALSE] -
#     inv_st_zt[o_index, , drop = FALSE] %*%
#     (chol_inv_smw1_p2mid %*% t(inv_st_zt[m_index, , drop = FALSE])) -
#     iv_stpt_zs[o_index, , drop = FALSE] %*%
#     (chol_inv_smw2_p2mid %*% t(inv_stpt_zs[m_index, , drop = FALSE]))
#
#   invdense_mm <- w_vinvroot[m_index, o_index, drop = FALSE] %*% w_vinvroot[o_index, m_index, drop = FALSE] +
#     w_vinvroot[m_index, m_index, drop = FALSE] %*% w_vinvroot[m_index, m_index, drop = FALSE] -
#     inv_st_zt[m_index, , drop = FALSE] %*%
#     (chol_inv_smw1_p2mid %*% t(inv_st_zt[m_index, , drop = FALSE])) -
#     iv_stpt_zs[m_index, , drop = FALSE] %*%
#     (chol_inv_smw2_p2mid %*% t(inv_stpt_zs[m_index, , drop = FALSE]))
#
#
#   # return the correct object
#   chol_mm <- chol(invdense_mm)
#   sigmainv_o <- invdense_oo - invdense_om %*% (chol2inv(chol_mm) %*% (t(invdense_om) %*% xyc_o))
#
#   # log determinant
#
#   if (log_determinant) {
#     logdet <- sum(log(v)) +
#       2 * sum(log(diag(chol_sigma_t))) + 2 * sum(log(diag(chol_smw1_p2mid))) +
#       2 * sum(log(diag(chol_sigma_s))) + 2 * sum(log(diag(chol_smw2_p2mid)))
#
#     if (dense_check)
#       logdet <- logdet + 2 * sum(log(diag(chol_mm)))
#   }
#
  return(list(siginv_o = siginv_o, logdet = logdet))

}


# x = c(1, 1, 2, 2)
# t = c(1, 2, 1, 2)
# r_s = exp(-h_make(unique(x)))
# r_t = exp(-h_make(unique(t)))
# s_de = t_de = s_ie = t_ie = st_de = st_ie = 1
# o_index = 1:4
# xyc_o = matrix(1, nrow = 4)
# diag_tol = 0
# n_s <- ncol(r_s)
# n_t <- ncol(r_t)
# n_st <- n_s * n_t
# m_index <- full_index[-o_index]
# dense <- length(o_index) == n_st
#
# ## adding diagonal tolerances for invertibility stability
# diag(r_s) <- diag(r_s) + diag_tol
# diag(r_t) <- diag(r_t) + diag_tol
# s_ie <- s_ie + diag_tol
# t_ie <- t_ie + diag_tol
# st_ie <- st_ie + diag_tol
# log_determinant = TRUE
# invert_productsum(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de, st_ie, o_index, xyc_o, diag_tol)
#
#
#
# full_rs = exp(-h_make(x))
# full_rt = exp(-h_make(t))
# z_s = model.matrix(~ as.factor(x) - 1)
# z_t = model.matrix(~ as.factor(t) - 1)
# pscov = s_de * full_rs + s_ie * (z_s %*% diag(2) %*% t(z_s)) +
#   t_de * full_rt + t_ie * (z_t %*% diag(2) %*% t(z_t)) +
#   st_de * full_rs * full_rt + st_ie * diag(4)
# determinant(pscov)
# solve(pscov) %*% xyc_o
#
#
#
#
# x = c(1, 1, 2, 2)
# t = c(1, 2, 1, 2)
# r_s = exp(-h_make(unique(x)))
# r_t = exp(-h_make(unique(t)))
# s_de = t_de = s_ie = t_ie = st_de = st_ie = 1
# o_index = c(1)
# xyc_o = matrix(1, nrow = 3)
# diag_tol = 0
# n_s <- ncol(r_s)
# n_t <- ncol(r_t)
# n_st <- n_s * n_t
# full_index = 1:4
# m_index <- full_index[-o_index]
# dense <- length(o_index) == n_st
#
# ## adding diagonal tolerances for invertibility stability
# diag(r_s) <- diag(r_s) + diag_tol
# diag(r_t) <- diag(r_t) + diag_tol
# s_ie <- s_ie + diag_tol
# t_ie <- t_ie + diag_tol
# st_ie <- st_ie + diag_tol
# log_determinant = TRUE
# xyc_o = matrix(1, nrow = 4)
# xyc_o <- xyc_o[o_index, , drop = FALSE]
# invert_productsum(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de, st_ie, o_index, xyc_o, diag_tol)
#
#
# full_rs = exp(-h_make(x))
# full_rt = exp(-h_make(t))
# z_s = model.matrix(~ as.factor(x) - 1)
# z_t = model.matrix(~ as.factor(t) - 1)
# pscov = s_de * full_rs + s_ie * (z_s %*% diag(2) %*% t(z_s)) +
#   t_de * full_rt + t_ie * (z_t %*% diag(2) %*% t(z_t)) +
#   st_de * full_rs * full_rt + st_ie * diag(4)
# pscov = pscov[o_index, o_index]
# determinant(pscov)
# solve(pscov) %*% xyc_o
#
#
