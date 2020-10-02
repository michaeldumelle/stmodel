invert <- function(invert_object, ...) {
  UseMethod("invert", object = invert_object)
}


# make cov param vector a list within invert object
# can just store r_s and r_t elsewhere for now in the invert object


invert.productsum <- function(invert_object, ...) {



  if (invert_object$chol) {

    # layout the arguments
    rf_s <- invert_object$rf_s
    rf_t <- invert_object$rf_t
    s_de <- invert_object$covparams[["s_de"]]
    s_ie <- invert_object$covparams[["s_ie"]]
    t_de <- invert_object$covparams[["t_de"]]
    t_ie <- invert_object$covparams[["t_ie"]]
    st_de <- invert_object$covparams[["st_de"]]
    st_ie <- invert_object$covparams[["st_ie"]]
    xyc_o <- invert_object$xyc_o
    diag_tol <- invert_object$diag_tol
    logdet <- invert_object$logdet

    cov_s <- s_de * rf_s + s_ie * (rf_s == 1)
    cov_t <- t_de * rf_t + t_ie * (rf_t == 1)
    cov_st <- st_de * rf_t * rf_s + st_ie * (rf_s == 1) * (rf_t == 1)
    sigma <- cov_s + cov_t + cov_st
    diag(sigma) <- diag(sigma) + diag_tol
    chol_sigma <- chol(sigma)
    siginv <- chol2inv(chol_sigma)
    siginv_o <- siginv %*% xyc_o
    if (logdet){
      logdet <- 2 * sum(log(diag(chol_sigma)))
    } else {
      logdet <- NULL
    }
  } else {


    # layout the arguments
    r_s <- invert_object$r_s
    r_t <- invert_object$r_t
    s_de <- invert_object$covparams[["s_de"]]
    s_ie <- invert_object$covparams[["s_ie"]]
    t_de <- invert_object$covparams[["t_de"]]
    t_ie <- invert_object$covparams[["t_ie"]]
    st_de <- invert_object$covparams[["st_de"]]
    st_ie <- invert_object$covparams[["st_ie"]]
    xyc_o <- invert_object$xyc_o
    diag_tol <- invert_object$diag_tol
    logdet <- invert_object$logdet
    o_index <- invert_object$o_index
    m_index <- invert_object$m_index

    n_s <- ncol(r_s)
    n_t <- ncol(r_t)
    n_st <- n_s * n_t
    dense <- length(o_index) == n_st




    ## adding diagonal tolerances for invertibility stability
    diag(r_s) <- diag(r_s) + diag_tol
    diag(r_t) <- diag(r_t) + diag_tol
    s_ie <- s_ie + diag_tol
    t_ie <- t_ie + diag_tol
    st_ie <- st_ie + diag_tol

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
    t_w_vinvroot <- as.vector(vinvroot) * t(w)
    #creating the transpose
    w_vinvroot <- t(t_w_vinvroot)


    # storing the output we will need for the iterative smw
    c_t <- chol(sigma_make(t_de, r_t, t_ie))
    c_s <- chol(sigma_make(s_de, r_s, s_ie))

    ist_zt <- w_vinvroot %*% multiply_z(t_w_vinvroot, "temporal", n_s, n_t, "right")
              # st x st %*% (st x st * st x t) = st x t
    tr_ist_zt <- t(ist_zt)
    c_mt <- chol(chol2inv(c_t) + multiply_z(ist_zt, "temporal", n_s, n_t, "p_left"))
                   #t(multiply_z(mx = t(ist_zt), z_type = "temporal", n_s = n_s)))
              # (t x t + tr(tr(st x t) * st x t) = t x t
    ic_mt <- chol2inv(c_mt)
              # t x t

    istpt_zs <- w_vinvroot %*% multiply_z(mx = t_w_vinvroot, z_type = "spatial", n_s = n_s, n_t = n_t, side = "right") -
      ist_zt %*% (ic_mt %*% multiply_z(mx = tr_ist_zt, z_type = "spatial", n_s = n_s, n_t = n_t, side = "right"))
      # st x st * (st x st * st x s) - st x t (t x t * (tr(st x t) * st x s)) =
      # st x s - st x t * t x s = st x s
    tr_istpt_zs <- t(istpt_zs)
    c_ms <- chol(chol2inv(c_s) + multiply_z(mx = istpt_zs, z_type = "spatial", n_s = n_s, n_t = n_t, side = "p_left"))
                   #t(multiply_z(mx = tr_istpt_zs, z_type = "spatial", n_s = n_s)))
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

    if (logdet) {
      logdet <- sum(log(v)) +
        2 * sum(log(diag(c_t))) + 2 * sum(log(diag(c_mt))) +
        2 * sum(log(diag(c_s))) + 2 * sum(log(diag(c_ms)))

      if (!dense)
        logdet <- logdet + 2 * sum(log(diag(c_mm)))
    } else {
      logdet <- NULL
    }
  }

  output <- list(siginv_o = siginv_o, logdet = logdet)
  output_non_null <- output[!unlist(lapply(output, is.null))]
  return(output_non_null)
}

invert_sum_with_error <- function(r_s, r_t, s_de, s_ie, t_de, t_ie, st_ie,
                              o_index, m_index, xyc_o, diag_tol, logdet) {

  n_s <- ncol(r_s)
  n_t <- ncol(r_t)
  n_st <- n_s * n_t
  dense <- length(o_index) == n_st

  ## adding diagonal tolerances for invertibility stability
  diag(r_s) <- diag(r_s) + diag_tol
  diag(r_t) <- diag(r_t) + diag_tol
  s_ie <- s_ie + diag_tol
  t_ie <- t_ie + diag_tol
  st_ie <- st_ie + diag_tol

  # storing the output we will need for the iterative smw
  c_t <- chol(sigma_make(t_de, r_t, t_ie))
  c_s <- chol(sigma_make(s_de, r_s, s_ie))

  c_mt <- chol(chol2inv(c_t) + multiply_z(z_type = "temporal", n_s = n_s, n_t = n_t, side = "pz_z")/st_ie)
  ic_mt <- chol2inv(c_mt)
  istpt <- - multiply_z(multiply_z(ic_mt, z_type = "temporal", n_s = n_s, n_t = n_t, side = "p_right"),
                        z_type = "temporal", n_s = n_s, n_t = n_t, side = "left") / (st_ie^2)
  diag(istpt) <- diag(istpt) + 1/st_ie


  istpt_zs <- multiply_z(mx = istpt, z_type = "spatial", n_s = n_s, n_t = n_t, side = "right")
  c_ms <- chol(chol2inv(c_s) + multiply_z(mx = istpt_zs, z_type = "spatial", n_s = n_s, n_t = n_t, side = "p_left"))
  ic_ms <- chol2inv(c_ms)


  siginv <- istpt - istpt_zs %*% (ic_ms %*% t(istpt_zs))

  if (dense) {
    siginv_o <- siginv
  } else {
    c_mm <- chol(siginv[m_index, m_index])
    siginv_o <- (siginv[o_index, o_index] - siginv[o_index, m_index] %*% (chol2inv(c_mm) %*% siginv[m_index, o_index])) %*% xyc_o
  }


  if (logdet) {
    logdet <- n_st * (log(st_ie)) +
      2 * sum(log(diag(c_t))) + 2 * sum(log(diag(c_mt))) +
      2 * sum(log(diag(c_s))) + 2 * sum(log(diag(c_ms)))

    if (!dense)
      logdet <- logdet + 2 * sum(log(diag(c_mm)))
  } else {
    logdet <- NULL
  }

  output <- list(siginv_o = siginv_o, logdet = logdet)
  output_non_null <- output[!unlist(lapply(output, is.null))]
  return(output_non_null)
}

invert_product <- function(r_s, r_t, vs_ie, vt_ie, st_de,
                              o_index, m_index, xyc_o, diag_tol, logdet) {

  n_s <- ncol(r_s)
  n_t <- ncol(r_t)
  n_st <- n_s * n_t
  dense <- length(o_index) == n_st

  diag(r_s) <- diag(r_s) + diag_tol
  diag(r_t) <- diag(r_t) + diag_tol
  st_de <- st_de + diag_tol

  scale_r_s <- sigma_make(r_mx = r_s, v_ie = vs_ie, e = 1, scale = TRUE)
  c_scale_r_s  <- chol(scale_r_s)
  scale_r_t <- sigma_make(r_mx = r_t, v_ie = vt_ie, e = 1, scale = TRUE)
  c_scale_r_t  <- chol(scale_r_t)

  siginv <- kronecker(chol2inv(c_scale_r_t), chol2inv(c_scale_r_s)) / st_de

  if (dense) {
    siginv_o <- siginv
  } else {
    c_mm <- chol(siginv[m_index, m_index])
    siginv_o <- (siginv[o_index, o_index] - siginv[o_index, m_index] %*% (chol2inv(c_mm) %*% siginv[m_index, o_index])) %*% xyc_o
  }


  if (logdet){
    logdet <- n_st * log(st_de) +
      n_s * 2 * sum(log(diag(c_scale_r_t))) +
      n_t * 2 * sum(log(diag(c_scale_r_s)))
    if (!dense){
      logdet <- logdet + 2 * sum(log(diag(c_mm)))
    }
  } else {
    logdet <- NULL
  }
  output <- list(siginv_o = siginv_o, logdet = logdet)
  output_non_null <- output[!unlist(lapply(output, is.null))]
  return(output_non_null)
}

#
# x = rep(1:4, times = 4)
# t = rep(1:4, each = 4)
# r_s = exp(-h_make(unique(x)))
# r_t = exp(-h_make(unique(t)))
# # s_de = t_de = s_ie = t_ie = st_ie = 0
# # st_de <- 6
# # vs_ie = vt_ie = 0.5
# s_de = t_de = s_ie = t_ie = st_ie = st_de = 1/6
# full_index = 1:length(x)
# index_sample <- sample(full_index, 2)
# o_index = full_index[-index_sample]
# m_index = full_index[index_sample]
# xyc_o = matrix(1, nrow = length(x))
# diag_tol = 1e-4
# n_s <- ncol(r_s)
# n_t <- ncol(r_t)
# n_st <- n_s * n_t
# # full_index = o_index[-sample(full_index, 10)]
# # m_index <- full_index[-c(1, 100, 400, 500, 600, 725, 900, 1250)]
# dense <- length(o_index) == n_st
#
# ## adding diagonal tolerances for invertibility stability
# diag(r_s) <- diag(r_s) + diag_tol
# diag(r_t) <- diag(r_t) + diag_tol
# s_ie <- s_ie + diag_tol
# t_ie <- t_ie + diag_tol
# st_ie <- st_ie + diag_tol
# log_determinant = TRUE
# logdet = TRUE
# xyc_o = matrix(1, nrow = length(x))
# xyc_o <- xyc_o[o_index, , drop = FALSE]
# xyc_o = cbind(xyc_o, rnorm(length(o_index)))
# test = invert_product(r_s = r_s, r_t = r_t, vs_ie = vs_ie, vt_ie = vt_ie, st_de = st_de,
#                       o_index = o_index, m_index = m_index, xyc_o = xyc_o, diag_tol = diag_tol)
#
#
# full_rs = (1 - vs_ie) * exp(-h_make(x)) + vs_ie * (h_make(x) == 0)
# full_rt = (1 - vt_ie) * exp(-h_make(t)) + vt_ie * (h_make(t) == 0)
# z_s = model.matrix(~ as.factor(x) - 1)
# z_t = model.matrix(~ as.factor(t) - 1)
# pscov =   st_de * full_rs * full_rt
# pscov = pscov[o_index, o_index]
# determinant(pscov)
# psinv <- chol2inv(chol(pscov))
# test2 = psinv %*% xyc_o
#
# testf = cbind(test$siginv_o, test2)
# all.equal(testf[, 1:2], testf[, 3:4])
# log_determinant = FALSE
#
# microbenchmark::microbenchmark(invert_product(r_s = r_s, r_t = r_t, vs_ie = vs_ie, vt_ie = vt_ie, st_de = st_de,
#                o_index = o_index, m_index = m_index, xyc_o = xyc_o, diag_tol = diag_tol, logdet = logdet), times = 15)
# microbenchmark::microbenchmark(invert_sum_with_error(r_s, r_t, s_de, s_ie, t_de, t_ie,
#                                                      st_ie, o_index, m_index,  xyc_o, diag_tol = diag_tol, logdet = logdet), times = 15)
# microbenchmark::microbenchmark(invert_productsum(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de,
#                                                  st_ie, o_index, m_index,  xyc_o, diag_tol, logdet = logdet), times = 15)
# microbenchmark::microbenchmark(chol2inv(chol(pscov)) %*% xyc_o, times = 15)
# test3 = invert_product(r_s, r_t, vs_ie, vt_ie, st_de,
#                o_index, m_index, xyc_o, diag_tol)
















# invert <- function(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de, st_ie,
#                    o_index, m_index, xyc_o, diag_tol = 1e-4,
#                    st_cov = c("productsum", "product", "sum_with_error"),
#                    log_determinant = TRUE) {
#
#   # do i combine this function and the cholesky function?
#
#
#   st_cov <- match.arg(st_cov)
#   switch(st_cov,
#          "productsum" = invert_productsum(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de, st_ie,
#                                           o_index, m_index, xyc_o, diag_tol),
#          "product" = invert_product(),
#          "sum_with_error" = invert_sum_with_error(),
#          stop("choose a valid spatio-temporal covariance structure"))
# }
