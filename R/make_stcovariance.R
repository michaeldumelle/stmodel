make_stcovariance <- function(covparam_object, h_s_large, h_t_large,
                               s_cor, t_cor, ...){
  UseMethod("make_stcovariance", object = covparam_object)
}

make_stcovariance.productsum <- function(covparam_object, h_s_large, h_t_large,
                                          s_cor, t_cor){
  r_s_large <- make_r(h = h_s_large, range = covparam_object[["s_range"]],
                structure = s_cor)
  r_t_large <- make_r(h = h_t_large, range = covparam_object[["t_range"]],
                structure = t_cor)
  r_st <- r_s_large * r_t_large

  # make covariance matrices
  sigma_s <- make_sigma(de = covparam_object[["s_de"]], r_mx = r_s_large, ie = covparam_object[["s_ie"]])
  sigma_t <- make_sigma(de = covparam_object[["t_de"]], r_mx = r_t_large, ie = covparam_object[["t_ie"]])
  sigma_st <- make_sigma(de = covparam_object[["st_de"]], r_mx = r_st, ie = covparam_object[["st_ie"]])
  sigma <- sigma_s + sigma_t + sigma_st
  return(sigma)
}

make_stcovariance.sum_with_error <- function(covparam_object, h_s_large, h_t_large,
                                          s_cor, t_cor){
  r_s_large <- make_r(h = h_s_large, range = covparam_object[["s_range"]],
                      structure = s_cor)
  r_t_large <- make_r(h = h_t_large, range = covparam_object[["t_range"]],
                      structure = t_cor)
  r_st <- r_s_large * r_t_large

  # make covariance matrices
  sigma_s <- make_sigma(de = covparam_object[["s_de"]], r_mx = r_s_large, ie = covparam_object[["s_ie"]])
  sigma_t <- make_sigma(de = covparam_object[["t_de"]], r_mx = r_t_large, ie = covparam_object[["t_ie"]])
  sigma_st <- make_sigma(0, r_mx = r_st, ie = covparam_object[["st_ie"]])
  sigma <- sigma_s + sigma_t + sigma_st
  return(sigma)
}



make_stcovariance.product <- function(covparam_object, h_s_large, h_t_large,
                                          s_cor, t_cor){
r_s <- make_r(h = h_s_large, range = covparam_object[["s_range"]], structure = s_cor)
r_t <- make_r(h = h_t_large, range = covparam_object[["t_range"]], structure = t_cor)
r_st <- r_s * r_t

# make covariance matrices
sigma_st <- make_sigma(de = covparam_object[["st_de"]], r_mx = r_st, ie = 0)
sigma <- sigma_st
return(sigma)
}
