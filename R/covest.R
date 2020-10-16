covest <- function(par, covest_object, ...){
  UseMethod("covest", object = covest_object)
}

covest.svwls <- function(par, covest_object, data_object){
  # multiply overall var by exponentiated
  # transform profiled to regular
  plo2r <- plo2r.svwls(par, covest_object)
  class(plo2r) <- class(covest_object)
  theo_sv <- make_stsemivariogram(covparam_object = plo2r, h_s_large = covest_object$stempsv$h_s_avg,
                        h_t_large = covest_object$stempsv$h_t_avg, sp_cor = covest_object$sp_cor,
                        t_cor = covest_object$t_cor)
  # create the weights
  wts <- switch(covest_object$weights,
                "cressie" = weights_cressie(sv = covest_object$stempsv, theo_sv = theo_sv),
                stop("choose valid weights"))

  # create the objective function
  sumsq <- (covest_object$stempsv$gammahat - theo_sv)^2

  # return the objective function
  return(sum(wts * sumsq))
}

weights_cressie <- function(sv, theo_sv) {
  wts <- sv$n / theo_sv^2
  return(wts)
}

covest.reml <- function(par, covest_object, invert_object){
  plo2r <- plo2r.reml(par, covest_object, ov_var = 1)
  invert_object$covparams <- plo2r
  invert_output <- invert(invert_object)
  m2ll <- minus2loglik.reml(invert_object = invert_object, invert_output = invert_output)
  return(m2ll)
}


