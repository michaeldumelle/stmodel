strnorm <- function(object, mu, size, condition, ...){
  UseMethod("strnorm", object = object)
}

strnorm.matrix <- function(object, mu, size, condition = 1e-4){
  n_st <- nrow(object)
  if (length(mu) != 1 & length(mu) != n_st) stop("choose mu as a vector having size nst or a single scalar")
  chol_siginv <- t(chol(object))
  return(vapply(1:size, FUN = function(x) mu + chol_siginv %*% rnorm(n_st), double(n_st)))
}


strnorm.default <- function(object, mu, size, xcoord, ycoord = NULL, tcoord, data,
                            sp_cor, t_cor, chol = FALSE, condition = 1e-4, h_options = NULL){
  if (is.null(h_options)){
    h_options = list(h_large = TRUE, h_t_distmetric = "euclidean", h_s_distmetric = "euclidean")
  }

  if (chol){
    if (is.null(ycoord)){
      h_s_large <- make_h(coord1 = data[[xcoord]],
                          distmetric = h_options$h_s_distmetric)
      h_t_large <- make_h(coord1 = data[[tcoord]], distmetric = h_options$h_t_distmetric)
    } else{
      h_s_large <- make_h(coord1 = data[[xcoord]], coord2 = data[[ycoord]],
                          distmetric = h_options$h_s_distmetric)
      h_t_large <- make_h(coord1 = data[[tcoord]], distmetric = h_options$h_t_distmetric)
    }
    st_covariance <- make_stcovariance(covparam_object = object, h_s_large = h_s_large,
                                                  h_t_large = h_t_large, sp_cor = sp_cor,
                                                  t_cor = t_cor)
    strnorm_sim <- strnorm.matrix(object = st_covariance, mu = mu, size = size, condition = condition)
  } else {
    covparam_object <- object
    data$original_key <- seq.int(1, nrow(data))
    spint <- storder(data = data, xcoord = xcoord, ycoord = ycoord, tcoord = tcoord,
                     h_options = h_options)
    r_s_small <- make_r(h = spint$h_s_small, range = covparam_object[["s_range"]], structure = sp_cor)
    diag(r_s_small) <- diag(r_s_small) + condition


    r_t_small <- make_r(h = spint$h_t_small, range = covparam_object[["t_range"]], structure = t_cor)
    diag(r_t_small) <- diag(r_t_small) + condition



    strnorm_sim <- strnorm_small(covparam_object = object, mu = mu, size = size,
                                 r_s_small = r_s_small, r_t_small = r_t_small)

    strnorm_sim <- strnorm_sim[spint$ordered_data_dense$observed, , drop = FALSE ]
    strnorm_sim <- strnorm_sim[order(spint$ordered_data_o$original_key), , drop = FALSE ]
    return(mu + strnorm_sim)
  }
}

strnorm_small <- function(covparam_object, mu, size, r_s_small, r_t_small){
  UseMethod("strnorm_small", object = covparam_object)
}


strnorm_small.productsum <- function(covparam_object, mu, size, r_s_small, r_t_small){
  chol_r_s_small <- t(chol(r_s_small))
  chol_r_t_small <- t(chol(r_t_small))

  n_s <- nrow(chol_r_s_small)
  n_t <- nrow(chol_r_t_small)
  n_st <- n_s * n_t

  zs_chol_r_s_small <- multiply_z(mx = chol_r_s_small, z_type = "spatial", n_s = n_s,
                                  n_t = n_t, side = "left")
  zs_chol_r_t_small <- multiply_z(mx = chol_r_t_small, z_type = "temporal",
                                  n_s = n_s, n_t = n_t, side = "left")

  chol_r_st <- kronecker(chol_r_t_small, chol_r_s_small)
  strnorm_small_sim <- vapply(1:size, function(x) {
    return(zs_chol_r_s_small %*% rnorm(n_s, sd = sqrt(covparam_object[["s_de"]])) +
             multiply_z(mx = rnorm(n_s, sd = sqrt(covparam_object[["s_ie"]])), z_type = "spatial", n_s = n_s, n_t = n_t, side = "left") +
             zs_chol_r_t_small %*% rnorm(n_t, sd = sqrt(covparam_object[["t_de"]])) +
             multiply_z(mx = rnorm(n_t, sd = sqrt(covparam_object[["t_ie"]])), z_type = "temporal", n_s = n_s, n_t = n_t, side = "left") +
             chol_r_st %*% rnorm(n_st, sd = sqrt(covparam_object[["st_de"]])) +
             rnorm(n_st, sd = sqrt(covparam_object[["st_ie"]])))
  }, double(n_st))
  return(strnorm_small_sim)
}

strnorm_small.sum_with_error <- function(covparam_object, mu, size, r_s_small, r_t_small){
  chol_r_s_small <- t(chol(r_s_small))
  chol_r_t_small <- t(chol(r_t_small))

  n_s <- nrow(chol_r_s_small)
  n_t <- nrow(chol_r_t_small)
  n_st <- n_s * n_t

  zs_chol_r_s_small <- multiply_z(mx = chol_r_s_small, z_type = "spatial", n_s = n_s,
                                  n_t = n_t, side = "left")
  zs_chol_r_t_small <- multiply_z(mx = chol_r_t_small, z_type = "temporal",
                                  n_s = n_s, n_t = n_t, side = "left")


  strnorm_small_sim <- vapply(1:size, function(x) {
    return(zs_chol_r_s_small %*% rnorm(n_s, sd = sqrt(covparam_object[["s_de"]])) +
             multiply_z(mx = rnorm(n_s, sd = sqrt(covparam_object[["s_ie"]])), z_type = "spatial", n_s = n_s, n_t = n_t, side = "left") +
             zs_chol_r_t_small %*% rnorm(n_t, sd = sqrt(covparam_object[["t_de"]])) +
             multiply_z(mx = rnorm(n_t, sd = sqrt(covparam_object[["t_ie"]])), z_type = "temporal", n_s = n_s, n_t = n_t, side = "left") +
             rnorm(n_st, sd = sqrt(covparam_object[["st_ie"]])))
  }, double(n_st))
  return(strnorm_small_sim)
}


strnorm_small.product <- function(covparam_object, mu, size, r_s_small, r_t_small){
  chol_r_s_small <- t(chol(r_s_small))
  chol_r_t_small <- t(chol(r_t_small))

  n_s <- nrow(chol_r_s_small)
  n_t <- nrow(chol_r_t_small)
  n_st <- n_s * n_t

  chol_r_st <- kronecker(chol_r_t_small, chol_r_s_small)

  strnorm_small_sim <- vapply(1:size, function(x) {
    return(chol_r_st %*% rnorm(n_st, sd = sqrt(covparam_object[["st_de"]])))
  }, double(n_st))
  return(strnorm_small_sim)
}
