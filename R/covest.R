covest <- function(par, covest_object){
  UseMethod("covest", object = covest_object)
}

covest.svwls <- function(par, covest_object){
  UseMethod("covest.svwls", object = covest_object)
}

covest.svwls.productsum <- function(par, covest_object){
# multiply overall var by exponentiated
# transform profiled to regular
plo2r <- plo2r.svwls.productsum(par, covest_object)
# make correlation matrices
r_s <- r_make(h = covest_object$sv$avg_hsp, range = plo2r[["s_range"]], structure = covest_object$sp_cor)
r_t <- r_make(h = covest_object$sv$avg_tsp, range = plo2r[["t_range"]], structure = covest_object$t_cor)
r_st <- r_s * r_t

# make covariance matrices
sigma_s <- sigma_make(de = plo2r[["s_de"]], r_mx = r_s, ie = plo2r[["s_ie"]])
sigma_t <- sigma_make(de = plo2r[["t_de"]], r_mx = r_t, ie = plo2r[["t_ie"]])
sigma_st <- sigma_make(de = plo2r[["st_de"]], r_mx = r_st, ie = plo2r[["st_ie"]])
sigma <- sigma_s + sigma_t + sigma_st

# make C(0, 0)
vparams <- c(plo2r[["s_de"]], plo2r[["s_ie"]], plo2r[["t_de"]],
             plo2r[["t_ie"]], plo2r[["st_de"]], plo2r[["st_ie"]])

# compute the semivariogram
theo_sv <- sum(vparams) - sigma

# create the weights
wts <- switch(covest_object$weights,
              "cressie" = weight_cressie(sv = covest_object$sv, theo_sv = theo_sv),
              stop("choose valid weights"))

# create the objective function
sumsq <- (covest_object$sv$mean_sqdifs - theo_sv)^2

# return the objective function
return(sum(wts * sumsq))
}

covest.svwls.sum_with_error <- function(par, covest_object, ...){
  # multiply overall var by exponentiated
  # transform profiled to regular
  plo2r <- plo2r.svwls.sum_with_error(par, covest_object)
  # make correlation matrices
  r_s <- r_make(h = covest_object$sv$avg_hsp, range = plo2r[["s_range"]], structure = covest_object$sp_cor)
  r_t <- r_make(h = covest_object$sv$avg_tsp, range = plo2r[["t_range"]], structure = covest_object$t_cor)
  r_st <- r_s * r_t

  # make covariance matrices
  sigma_s <- sigma_make(de = plo2r[["s_de"]], r_mx = r_s, ie = plo2r[["s_ie"]])
  sigma_t <- sigma_make(de = plo2r[["t_de"]], r_mx = r_t, ie = plo2r[["t_ie"]])
  sigma_st <- sigma_make(0, r_mx = r_st, ie = plo2r[["st_ie"]])
  sigma <- sigma_s + sigma_t + sigma_st

  # make C(0, 0)
  vparams <- c(plo2r[["s_de"]], plo2r[["s_ie"]], plo2r[["t_de"]],
               plo2r[["t_ie"]], plo2r[["st_ie"]])

  # compute the semivariogram
  theo_sv <- sum(vparams) - sigma

  # create the weights
  wts <- switch(covest_object$weights,
                "cressie" = weight_cressie(sv = covest_object$sv, theo_sv = theo_sv),
                stop("choose valid weights"))

  # create the objective function
  sumsq <- (covest_object$sv$mean_sqdifs - theo_sv)^2

  # return the objective function
  return(sum(wts * sumsq))
}





covest.svwls.product <- function(par, covest_object, ...){
  # multiply overall var by exponentiated
  # transform profiled to regular
  plo2r <- plo2r.svwls.product(par, covest_object)
  # make correlation matrices
  r_s <- r_make(h = covest_object$sv$avg_hsp, range = plo2r[["s_range"]], structure = covest_object$sp_cor)
  r_t <- r_make(h = covest_object$sv$avg_tsp, range = plo2r[["t_range"]], structure = covest_object$t_cor)
  r_st <- r_s * r_t

  # make covariance matrices
  sigma_st <- sigma_make(de = plo2r[["st_de"]], r_mx = r_st, ie = 0)
  sigma <- sigma_st

  # make C(0, 0)
  vparams <- plo2r[["st_de"]]

  # compute the semivariogram
  theo_sv <- sum(vparams) - sigma

  # create the weights
  wts <- switch(covest_object$weights,
                "cressie" = weight_cressie(sv = covest_object$sv, theo_sv = theo_sv),
                stop("choose valid weights"))

  # create the objective function
  sumsq <- (covest_object$sv$mean_sqdifs - theo_sv)^2

  # return the objective function
  return(sum(wts * sumsq))
}

weight_cressie <- function(sv, theo_sv) {
  wts <- sv$n / theo_sv^2
  return(wts)
}
