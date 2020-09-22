p2r_productsum <- function(lambda, alpha, n_s, n_t, n_st, ov_var){
  s_de <- lambda * alpha * (1 - n_s)
  s_ie <- lambda * alpha * n_s
  t_de <- lambda * (1 - alpha) * (1 - n_t)
  t_ie <- lambda * (1 - alpha) * n_t
  st_de <- (1 - lambda) * (1 - n_st)
  st_ie <- (1 - lambda) * n_st
  rparm <- ov_var * c(s_de = s_de, s_ie = s_ie, t_de = t_de,
             t_ie = t_ie, st_de = st_de, st_ie = st_ie)
  return(rparm)
}

p2r_sum_with_error <- function(lambda, alpha, n_s, n_t, ov_var){
  s_de <- lambda * alpha * (1 - n_s)
  s_ie <- lambda * alpha * n_s
  t_de <- lambda * (1 - alpha) * (1 - n_t)
  t_ie <- lambda * (1 - alpha) * n_t
  st_ie <- (1 - lambda)
  rparm <- ov_var * c(s_de = s_de, s_ie = s_ie, t_de = t_de,
                      t_ie = t_ie, st_ie = st_ie)
  return(rparm)
}

pr2_product <- function(v_s, v_t, ov_var) {
  rparm(v_s = v_s, v_t = v_t, ov_var = ov_var)
}
