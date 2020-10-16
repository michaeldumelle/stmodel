r2plo <- function(covparam_object, ...){
  UseMethod("r2plo", object = covparam_object)
}

r2plo.svwls <- function(covparam_object, ...) {
  UseMethod("r2plo.svwls", object = covparam_object)
}

r2plo.svwls.productsum <- function(covparam_object, max_options){
  params <- covparam_object
  s_vc <- c(params[["s_de"]], params[["s_ie"]])
  svar <- sum(s_vc)
  t_vc <- c(params[["t_de"]], params[["t_ie"]])
  tvar <- sum(t_vc)
  st_vc <- c(params[["st_de"]], params[["st_ie"]])
  stvar <- sum(st_vc)
  vc <- c(s_vc, t_vc, st_vc)
  lambda <- (svar + tvar) / (svar + tvar + stvar)
  alpha <- svar / (svar + tvar)
  n_s <- params[["s_ie"]] / svar
  n_t <- params[["t_ie"]] / tvar
  n_st <- params[["st_ie"]] / stvar
  pparm <- c(lambda = lambda, alpha = alpha, n_s = n_s, n_t = n_t, n_st = n_st)
  pparm <- c(pparm, var_prop = pmin(1, (svar + tvar + stvar) / max_options$max_v),
             srange_prop = params[["s_range"]] / max_options$max_s_range,
             trange_prop = params[["t_range"]] / max_options$max_t_range)
  pparm <- log(pparm / (1 - pparm))
  return(pparm)
}


r2plo.svwls.sum_with_error <- function(covparam_object, max_options){
  params <- covparam_object
  s_vc <- c(params[["s_de"]], params[["s_ie"]])
  svar <- sum(s_vc)
  t_vc <- c(params[["t_de"]], params[["t_ie"]])
  tvar <- sum(t_vc)
  st_vc <- params[["st_ie"]]
  stvar <- sum(st_vc)
  vc <- c(s_vc, t_vc, st_vc)
  lambda <- (svar + tvar) / (svar + tvar + stvar)
  alpha <- svar / (svar + tvar)
  n_s <- params[["s_ie"]] / svar
  n_t <- params[["t_ie"]] / tvar
  pparm <- c(lambda = lambda, alpha = alpha, n_s = n_s, n_t = n_t)
  pparm <- c(pparm, var_prop = pmin(1, (svar + tvar + stvar) / max_options$max_v),
             srange_prop = params[["s_range"]] / max_options$max_s_range,
             trange_prop = params[["t_range"]] / max_options$max_t_range)
  pparm <- log(pparm / (1 - pparm))
  return(pparm)
}

r2plo.svwls.product <- function(covparam_object, max_options){
  params <- covparam_object
  pparm <- c(v_s = params[["v_s"]], v_t = params[["v_t"]])
  pparm <- c(pparm, var_prop = pmin(1, params[["st_de"]] / max_options$max_v),
             srange_prop = params[["s_range"]] / max_options$max_s_range,
             trange_prop = params[["t_range"]] / max_options$max_t_range)
  pparm <- log(pparm / (1 - pparm))
  return(pparm)
}

r2plo.reml <- function(covparam_object, ...) {
  UseMethod("r2plo.reml", object = covparam_object)
}

r2plo.reml.productsum <- function(covparam_object, max_options){
  params <- covparam_object
  s_vc <- c(params[["s_de"]], params[["s_ie"]])
  svar <- sum(s_vc)
  t_vc <- c(params[["t_de"]], params[["t_ie"]])
  tvar <- sum(t_vc)
  st_vc <- c(params[["st_de"]], params[["st_ie"]])
  stvar <- sum(st_vc)
  vc <- c(s_vc, t_vc, st_vc)
  lambda <- (svar + tvar) / (svar + tvar + stvar)
  alpha <- svar / (svar + tvar)
  n_s <- params[["s_ie"]] / svar
  n_t <- params[["t_ie"]] / tvar
  n_st <- params[["st_ie"]] / stvar
  pparm <- c(lambda = lambda, alpha = alpha, n_s = n_s, n_t = n_t, n_st = n_st)
  pparm <- c(pparm,
             srange_prop = params[["s_range"]] / max_options$max_s_range,
             trange_prop = params[["t_range"]] / max_options$max_t_range)
  pparm <- log(pparm / (1 - pparm))
  return(pparm)
}


r2plo.reml.sum_with_error <- function(covparam_object, max_options){
  params <- covparam_object
  s_vc <- c(params[["s_de"]], params[["s_ie"]])
  svar <- sum(s_vc)
  t_vc <- c(params[["t_de"]], params[["t_ie"]])
  tvar <- sum(t_vc)
  st_vc <- params[["st_ie"]]
  stvar <- sum(st_vc)
  vc <- c(s_vc, t_vc, st_vc)
  lambda <- (svar + tvar) / (svar + tvar + stvar)
  alpha <- svar / (svar + tvar)
  n_s <- params[["s_ie"]] / svar
  n_t <- params[["t_ie"]] / tvar
  pparm <- c(lambda = lambda, alpha = alpha, n_s = n_s, n_t = n_t)
  pparm <- c(pparm,
             srange_prop = params[["s_range"]] / max_options$max_s_range,
             trange_prop = params[["t_range"]] / max_options$max_t_range)
  pparm <- log(pparm / (1 - pparm))
  return(pparm)
}

r2plo.reml.product <- function(covparam_object, max_options){
  params <- covparam_object
  pparm <- c(v_s = params[["v_s"]], v_t = params[["v_t"]])
  pparm <- c(pparm,
             srange_prop = params[["s_range"]] / max_options$max_s_range,
             trange_prop = params[["t_range"]] / max_options$max_t_range)
  pparm <- log(pparm / (1 - pparm))
  return(pparm)
}



