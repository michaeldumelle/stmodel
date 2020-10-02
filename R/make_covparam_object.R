make_covparam_object <- function(s_de, s_ie, t_de,
                                 t_ie, st_de, st_ie,
                                 v_s, v_t, srange,
                                 trange,
                                 stcov){
  covparam_object <- switch(stcov,
                            "productsum" = covparam_object_productsum(s_de = s_de, s_ie = s_ie,
                                                                      t_de = t_de, t_ie = t_ie,
                                                                      st_de = st_de, st_ie = st_ie,
                                                                      srange = srange, trange = trange),
                            "sum_with_error" = covparam_object_sum_with_error(s_de = s_de, s_ie = s_ie, t_de = t_de,
                                                                              t_ie = t_ie, st_ie = st_ie,
                                                                              srange = srange, trange = trange),
                            "product" = covparam_object_product(st_de = st_de, v_s = v_s, v_t = v_t,
                                                                srange = srange, trange = trange))
  covparam_object <- structure(covparam_object, class = stcov)
  return(covparam_object)
}

covparam_object_productsum <- function(s_de, s_ie, t_de,
                                t_ie, st_de, st_ie,
                                srange, trange){
  return(c(s_de = s_de, s_ie = s_ie, t_de = t_de,
              t_ie = t_ie, st_de = st_de, st_ie = st_ie,
              srange = srange, trange = trange))
}

covparam_object_sum_with_error <- function(s_de, s_ie, t_de,
                                       t_ie, st_ie,
                                       srange, trange){
  return(c(s_de = s_de, s_ie = s_ie, t_de = t_de,
           t_ie = t_ie, st_ie = st_ie,
           srange = srange, trange = trange))
}

covparam_object_product <- function(st_de, v_s, v_t,
                                       srange, trange){
  return(c(st_de = st_de, v_s = v_s, v_t = v_t,
           srange = srange, trange = trange))
}




# should just make this a named vector, not a list - works better with optim
# use this class to check the argument of stcov to stlm
