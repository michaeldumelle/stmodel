# r2p <- function(s_de, s_ie, t_de, t_ie, st_de, st_ie, v_s, v_t, estmethod, stcov){
#   if (!(estmethod %in% c("reml", "ml", "sv"))){
#     stop("choose a valid estimation method")
#   }
#   stcov <- match.arg(stcov)
#   switch(stcov,
#          "productsum" = r2p_productsum(s_de = s_de, s_ie = s_ie, t_de = t_de,
#                                        t_ie = t_ie, st_de = st_de, st_ie = st_ie,
#                                        estmethod))
# }

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

