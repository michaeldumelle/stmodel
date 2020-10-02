sv_fn <- function(par, covest_object, ...){
  UseMethod("sv_fn", object = covest_object)
}



sv_fn.productsum <- function(par, covest_object, ...){
  # multiply overall var by exponentiated
  invlogit <- exp(par) / (1 + exp(par))
  p2r <- p2r_sv.productsum(par, covest_object)
  r_s <- r_make(h = covest_object$sv$avg_hsp, range = p2r[["s_range"]], structure = covest_object$sp_cor)
  r_t <- r_make(h = covest_object$sv$avg_tsp, range = p2r[["t_range"]], structure = covest_object$t_cor)

  sigma_s <- p2r[["s_de"]] * r_s + p2r[["s_ie"]] * (r_s == 1)
  sigma_t <- p2r[["t_de"]] * r_t + p2r[["t_ie"]] * (r_t == 1)
  sigma_st <- p2r[["st_de"]] * r_s * r_t + p2r[["st_ie"]] * (r_s == 1) * (r_t == 1)
  sigma <- sigma_s + sigma_t + sigma_st
  vparams <- c(p2r[["s_de"]], p2r[["s_ie"]], p2r[["t_de"]],
               p2r[["t_ie"]], p2r[["st_de"]], p2r[["st_ie"]])
  theo_sv <- sum(vparams) - sigma
  wts <- switch(covest_object$weights,
                "cressie" = weight_cressie(sv = sv, theo_sv = theo_sv),
                stop("choose valid weights"))
  sumsq <- (sv$mean_sqdifs - theo_sv)^2
  return(sum(wts * sumsq))
}


weight_cressie <- function(sv, theo_sv) {
  wts <- sv$n / theo_sv^2
  return(wts)
}
