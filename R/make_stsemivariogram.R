make_stsemivariogram <- function(covparam_object, h_s_large, h_t_large,
                               sp_cor, t_cor, ...){
  UseMethod("make_stsemivariogram", object = covparam_object)
}

make_stsemivariogram.productsum <- function(covparam_object, h_s_large, h_t_large,
                                          sp_cor, t_cor){
  variances <- c(covparam_object[c("s_de", "s_ie", "t_de", "t_ie", "st_de", "st_ie")])
  gamma <- sum(variances) - make_stcovariance.productsum(covparam_object = covparam_object,
                                                          h_s_large = h_s_large, h_t_large = h_t_large,
                                                          sp_cor = sp_cor, t_cor = t_cor)
  return(gamma)
}

make_stsemivariogram.sum_with_error <- function(covparam_object, h_s_large, h_t_large,
                                             sp_cor, t_cor){
  variances <- c(covparam_object[c("s_de", "s_ie", "t_de", "t_ie", "st_ie")])
  gamma <- sum(variances) - make_stcovariance.sum_with_error(covparam_object = covparam_object,
                                                          h_s_large = h_s_large, h_t_large = h_t_large,
                                                          sp_cor = sp_cor, t_cor = t_cor)
  return(gamma)
}

make_stsemivariogram.product <- function(covparam_object, h_s_large, h_t_large,
                                             sp_cor, t_cor){
  variances <- c(covparam_object[c("st_de")])
  gamma <- sum(variances) - make_stcovariance.product(covparam_object = covparam_object,
                                                          h_s_large = h_s_large, h_t_large = h_t_large,
                                                          sp_cor = sp_cor, t_cor = t_cor)
  return(gamma)
}
