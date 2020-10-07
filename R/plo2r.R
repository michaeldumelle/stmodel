plo2r <- function(par, data_object){
  UseMethod("plo2r", object = data_object)
}

plo2r.svwls <- function(par, data_object){
  UseMethod("plo2r.svwls", object = data_object)
}

plo2r.svwls.productsum <- function(par, data_object){
  invlogit <- exp(par) / (1 + exp(par))
  lambda <- invlogit[["lambda"]]
  alpha <- invlogit[["alpha"]]
  n_s <- invlogit[["n_s"]]
  n_t <- invlogit[["n_t"]]
  n_st <- invlogit[["n_st"]]

  s_de <- lambda * alpha * (1 - n_s)
  s_ie <- lambda * alpha * n_s
  t_de <- lambda * (1 - alpha) * (1 - n_t)
  t_ie <- lambda * (1 - alpha) * n_t
  st_de <- (1 - lambda) * (1 - n_st)
  st_ie <- (1 - lambda) * n_st

  ov_var <- data_object$max_v * invlogit[["var_prop"]]
  rparm <- c(ov_var * c(s_de = s_de, s_ie = s_ie, t_de = t_de,
                      t_ie = t_ie, st_de = st_de, st_ie = st_ie),
             s_range = data_object$max_s_range * invlogit[["srange_prop"]],
             t_range = data_object$max_t_range * invlogit[["trange_prop"]])
  return(rparm)
}


plo2r.svwls.sum_with_error <- function(par, data_object){
  invlogit <- exp(par) / (1 + exp(par))
  lambda <- invlogit[["lambda"]]
  alpha <- invlogit[["alpha"]]
  n_s <- invlogit[["n_s"]]
  n_t <- invlogit[["n_t"]]


  s_de <- lambda * alpha * (1 - n_s)
  s_ie <- lambda * alpha * n_s
  t_de <- lambda * (1 - alpha) * (1 - n_t)
  t_ie <- lambda * (1 - alpha) * n_t
  st_ie <- (1 - lambda)

  ov_var <- data_object$max_v * invlogit[["var_prop"]]
  rparm <- c(ov_var * c(s_de = s_de, s_ie = s_ie, t_de = t_de,
                        t_ie = t_ie, st_ie = st_ie),
             s_range = data_object$max_s_range * invlogit[["srange_prop"]],
             t_range = data_object$max_t_range * invlogit[["trange_prop"]])
  return(rparm)
}

plo2r.svwls.product <- function(par, data_object){
  invlogit <- exp(par) / (1 + exp(par))
  v_s <- invlogit[["v_s"]]
  v_t <- invlogit[["v_t"]]
  ov_var <- data_object$max_v * invlogit[["var_prop"]]
  rparm <- c(v_s = v_s, v_t = v_t, st_de = ov_var,
             s_range = data_object$max_s_range * invlogit[["srange_prop"]],
             t_range = data_object$max_t_range * invlogit[["trange_prop"]])
  return(rparm)
}


# plo2r_productsum <- function(lambda, alpha, n_s, n_t, n_st, ov_var){
#   s_de <- lambda * alpha * (1 - n_s)
#   s_ie <- lambda * alpha * n_s
#   t_de <- lambda * (1 - alpha) * (1 - n_t)
#   t_ie <- lambda * (1 - alpha) * n_t
#   st_de <- (1 - lambda) * (1 - n_st)
#   st_ie <- (1 - lambda) * n_st
#   rparm <- ov_var * c(s_de = s_de, s_ie = s_ie, t_de = t_de,
#              t_ie = t_ie, st_de = st_de, st_ie = st_ie)
#   return(rparm)
# }
#
# p2r_sum_with_error <- function(lambda, alpha, n_s, n_t, ov_var){
#   s_de <- lambda * alpha * (1 - n_s)
#   s_ie <- lambda * alpha * n_s
#   t_de <- lambda * (1 - alpha) * (1 - n_t)
#   t_ie <- lambda * (1 - alpha) * n_t
#   st_ie <- (1 - lambda)
#   rparm <- ov_var * c(s_de = s_de, s_ie = s_ie, t_de = t_de,
#                       t_ie = t_ie, st_ie = st_ie)
#   return(rparm)
# }
#
# pr2_product <- function(v_s, v_t, ov_var) {
#   rparm(v_s = v_s, v_t = v_t, ov_var = ov_var)
# }
