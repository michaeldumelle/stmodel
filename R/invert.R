invert <- function(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de, st_ie,
                   o_index, m_index, xyc_o, diag_tol = 1e-4,
                   st_cov = c("productsum", "product", "sum_with_error"),
                   log_determinant = TRUE) {
  n_s <- ncol(r_s)
  n_t <- ncol(r_t)
  n_st <- n_s * n_t
  dense <- length(o_index) == n_st

  st_cov <- match.arg(st_cov)

  # mike - here are the things you need to do on 9/18
  ## for the linear with error covariance you must
  ## 1. make a change multiply_z to be multiply_z_r, multiply_zp_l, and multiply_z_l functions
  ## you only need multiply_z_l to work for z_t on the left in the linear_with_sum covariance
  ## 2. make a multiply_zp_z function
  ## you only need this to work for temporal
  ## 3. store sigma inverse (st + t) - you can store this because sigma inverse st is easy (diagonal)
  ## 4. compute sigma inverse (st + t + s)
  ## 5. use helmert wolf on 4. to compute observed inverse and log det

  ## for the separable covariance you must
  # 1. incorporate nugget proportion variances into sigma_make()
  # 2. make the inversion algorithm for kronecker products
  # 3. use helmert wolf on 3. to compute observed inverse and log det

  # closing

  ## 1. it would be nice to never have to store st x st matrices, but this would require
  ## very detailed code taking a kronecker product multiplied on the right by a matrix and find a
  ## way to do that right multiplication mid kronecker product
  ## right now - that is INFEASIBLE and one st x st matrix will unfortunately have to be stored
  ## for each algorithm (note to self unstore second st x st that is the transpose in ps algorithm)
  ## 2. clean up all commented stuff and make a reasonable test function - matt can you do this?

  switch(st_cov,
         "productsum" = invert_productsum(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de, st_ie,
                                          o_index, m_index, xyc_o, diag_tol),
         "product" = invert_product(),
         "sum_with_error" = invert_sum_with_error(),
         stop("choose a valid spatio-temporal covariance structure"))
}

invert_productsum <- function(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de, st_ie,
                              o_index, m_index, xyc_o, diag_tol) {

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

  if (log_determinant) {
    logdet <- sum(log(v)) +
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

x = rep(1:40, times = 40)
t = rep(1:40, each = 40)
r_s = exp(-h_make(unique(x)))
r_t = exp(-h_make(unique(t)))
s_de = t_de = s_ie = t_ie = st_de = st_ie = 1
full_index = 1:length(x)
index_sample <- sample(full_index, 10)
o_index = full_index[-index_sample]
m_index = full_index[index_sample]
xyc_o = matrix(1, nrow = length(x))
diag_tol = 0
n_s <- ncol(r_s)
n_t <- ncol(r_t)
n_st <- n_s * n_t
# full_index = o_index[-sample(full_index, 10)]
# m_index <- full_index[-c(1, 100, 400, 500, 600, 725, 900, 1250)]
dense <- length(o_index) == n_st

## adding diagonal tolerances for invertibility stability
diag(r_s) <- diag(r_s) + diag_tol
diag(r_t) <- diag(r_t) + diag_tol
s_ie <- s_ie + diag_tol
t_ie <- t_ie + diag_tol
st_ie <- st_ie + diag_tol
log_determinant = TRUE
xyc_o = matrix(1, nrow = length(x))
xyc_o <- xyc_o[o_index, , drop = FALSE]
xyc_o = cbind(xyc_o, rnorm(length(o_index)))
test = invert_productsum(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de, st_ie, o_index, m_index, xyc_o, diag_tol)


full_rs = exp(-h_make(x))
full_rt = exp(-h_make(t))
z_s = model.matrix(~ as.factor(x) - 1)
z_t = model.matrix(~ as.factor(t) - 1)
pscov = s_de * full_rs + s_ie * (z_s %*% diag(40) %*% t(z_s)) +
  t_de * full_rt + t_ie * (z_t %*% diag(40) %*% t(z_t)) +
  st_de * full_rs * full_rt + st_ie * diag(length(x))
pscov = pscov[o_index, o_index]
determinant(pscov)
test2 = solve(pscov) %*% xyc_o

testf = cbind(test$siginv_o, test2)
all.equal(testf[, 1:2], testf[, 3:4])
log_determinant = FALSE
microbenchmark::microbenchmark(invert_productsum(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de,
                                                 st_ie, o_index, m_index,  xyc_o, diag_tol), times = 15)
microbenchmark::microbenchmark(chol2inv(chol(pscov)) %*% xyc_o, times = 15)







invert_sum_with_error <- function(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de, st_ie,
                              o_index, m_index, xyc_o, diag_tol) {

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


  if (log_determinant) {
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


x = rep(1:50, times = 50)
t = rep(1:50, each = 50)
r_s = exp(-h_make(unique(x)))
r_t = exp(-h_make(unique(t)))
s_de = t_de = s_ie = t_ie = st_ie = 1
st_de <- 0
#st_ie = 0
full_index = 1:length(x)
index_sample <- sample(full_index, 10)
o_index = full_index[-index_sample]
m_index = full_index[index_sample]
xyc_o = matrix(1, nrow = length(x))
diag_tol = 0
n_s <- ncol(r_s)
n_t <- ncol(r_t)
n_st <- n_s * n_t
# full_index = o_index[-sample(full_index, 10)]
# m_index <- full_index[-c(1, 100, 400, 500, 600, 725, 900, 1250)]
dense <- length(o_index) == n_st

## adding diagonal tolerances for invertibility stability
diag(r_s) <- diag(r_s) + diag_tol
diag(r_t) <- diag(r_t) + diag_tol
s_ie <- s_ie + diag_tol
t_ie <- t_ie + diag_tol
st_ie <- st_ie + diag_tol
log_determinant = TRUE
xyc_o = matrix(1, nrow = length(x))
xyc_o <- xyc_o[o_index, , drop = FALSE]
xyc_o = cbind(xyc_o, rnorm(length(o_index)))
test = invert_sum_with_error(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de, st_ie, o_index, m_index, xyc_o, diag_tol)


full_rs = exp(-h_make(x))
full_rt = exp(-h_make(t))
z_s = model.matrix(~ as.factor(x) - 1)
z_t = model.matrix(~ as.factor(t) - 1)
pscov = s_de * full_rs + s_ie * (z_s %*% diag(50) %*% t(z_s)) +
  t_de * full_rt + t_ie * (z_t %*% diag(50) %*% t(z_t)) +
  st_de * full_rs * full_rt + st_ie * diag(length(x))
pscov = pscov[o_index, o_index]
determinant(pscov)
#chol(pscov)
test2 = solve(pscov) %*% xyc_o

testf = cbind(test$siginv_o, test2)
all.equal(testf[, 1:2], testf[, 3:4])
log_determinant = FALSE
microbenchmark::microbenchmark(invert_sum_with_error(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de,
                                                 st_ie, o_index, m_index,  xyc_o, diag_tol), times = 15)
microbenchmark::microbenchmark(invert_productsum(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de,
                                                     st_ie, o_index, m_index,  xyc_o, diag_tol), times = 15)
microbenchmark::microbenchmark(chol2inv(chol(pscov)) %*% xyc_o, times = 15)


invert_product <- function(r_s, r_t, vs_ie, vt_ie, st_de,
                              o_index, m_index, xyc_o, diag_tol) {
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


  if (log_determinant) {
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


# x = rep(1:50, times = 50)
# t = rep(1:50, each = 50)
# r_s = exp(-h_make(unique(x)))
# r_t = exp(-h_make(unique(t)))
# s_de = t_de = s_ie = t_ie = st_ie = 0
# st_de <- 6
# vs_ie = vt_ie = 0.5
# full_index = 1:length(x)
# index_sample <- sample(full_index, 10)
# o_index = full_index[-index_sample]
# m_index = full_index[index_sample]
# xyc_o = matrix(1, nrow = length(x))
# diag_tol = 0
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
#                o_index = o_index, m_index = m_index, xyc_o = xyc_o, diag_tol = diag_tol), times = 15)
# microbenchmark::microbenchmark(invert_sum_with_error(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de,
#                                                      st_ie, o_index, m_index,  xyc_o, diag_tol), times = 15)
# microbenchmark::microbenchmark(invert_productsum(r_s, r_t, s_de, s_ie, t_de, t_ie, st_de,
#                                                  st_ie, o_index, m_index,  xyc_o, diag_tol), times = 15)
# microbenchmark::microbenchmark(chol2inv(chol(pscov)) %*% xyc_o, times = 15)
# test3 = invert_product(r_s, r_t, vs_ie, vt_ie, st_de,
#                o_index, m_index, xyc_o, diag_tol)
