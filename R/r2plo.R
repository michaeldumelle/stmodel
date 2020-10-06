r2plo_sv <- function(covest_object) {
  UseMethod("r2plo_sv", object = covest_object)
}

r2plo_sv.productsum <- function(covest_object){
  params <- covest_object$initial
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
  pparm <- c(pparm, var_prop = pmin(1, (svar + tvar + stvar) / covest_object$max_v),
             srange_prop = params[["s_range"]] / covest_object$max_srange,
             trange_prop = params[["t_range"]] / covest_object$max_trange)
  pparm <- log(pparm / (1 - pparm))
  return(pparm)
}


r2plo_sv.sum_with_error <- function(covest_object){
  params <- covest_object$initial
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
  pparm <- c(pparm, var_prop = pmin(1, (svar + tvar + stvar) / covest_object$max_v),
             srange_prop = params[["s_range"]] / covest_object$max_srange,
             trange_prop = params[["t_range"]] / covest_object$max_trange)
  pparm <- log(pparm / (1 - pparm))
  return(pparm)
}





r2p_productsum <- function(s_de, s_ie, t_de, t_ie, st_de, st_ie, estmethod){
  s_vc <- c(s_de, s_ie)
  svar <- sum(s_vc)
  t_vc <- c(t_de, t_ie)
  tvar <- sum(t_vc)
  st_vc <- c(st_de, st_ie)
  stvar <- sum(st_vc)
  vc <- c(s_vc, t_vc, st_vc)
  lambda <- (svar + tvar) / (svar + tvar + stvar)
  alpha <- svar / (svar + tvar)
  n_s <- s_ie / svar
  n_t <- t_ie / tvar
  n_st <- st_ie / stvar
  pparm <- c(lambda = lambda, alpha = alpha, n_s = n_s, n_t = n_t, n_st = n_st)
  if (estmethod == "sv") {
    pparm <- c(pparm, ov_var = (svar + tvar + stvar))
  }
  return(pparm)
}

r2p_sum_with_error <- function(s_de, s_ie, t_de, t_ie, st_ie, estmethod){
  s_vc <- c(s_de, s_ie)
  svar <- sum(s_vc)
  t_vc <- c(t_de, t_ie)
  tvar <- sum(t_vc)
  vc <- c(s_vc, t_vc, st_ie)
  lambda <- (svar + tvar) / (svar + tvar + st_ie)
  alpha <- svar / (svar + tvar)
  n_s <- s_ie / svar
  n_t <- t_ie / tvar
  pparm <- c(lambda = lambda, alpha = alpha, n_s = n_s, n_t = n_t)
  if (estmethod == "sv") {
    pparm <- c(pparm, ov_var = (svar + tvar + st_ie))
  }
  return(pparm)
}

r2p_product <- function(v_s, v_t, st_de, estmethod){
  pparm <- c(v_s = v_s, v_t = v_t)
  if (estmethod == "sv") {
    pparm <- c(pparm, ov_var = st_de)
  }
  return(pparm)
}

