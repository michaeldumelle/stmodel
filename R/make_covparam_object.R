make_covparam_object <- function(s_de, s_ie, t_de,
                                 t_ie, st_de, st_ie,
                                 v_s, v_t, s_range,
                                 t_range,
                                 estmethod, stcov){
  # setting up flexible framework for other estimation methods
  if (estmethod %in% c("reml", "svwls")){
    # conditional call for type of spatio-temporal covariance
    covparam_object <- switch(stcov,
                              "productsum" = covparam_object_productsum(s_de = s_de, s_ie = s_ie,
                                                                        t_de = t_de, t_ie = t_ie,
                                                                        st_de = st_de, st_ie = st_ie,
                                                                        s_range = s_range, t_range = t_range),
                              "sum_with_error" = covparam_object_sum_with_error(s_de = s_de, s_ie = s_ie, t_de = t_de,
                                                                                t_ie = t_ie, st_ie = st_ie,
                                                                                s_range = s_range, t_range = t_range),
                              "product" = covparam_object_product(st_de = st_de, v_s = v_s, v_t = v_t,
                                                                  s_range = s_range, t_range = t_range),
                              stop("Use a valid error structure"))
    # giving the object the appropriate estimation method class
    covparam_object <- structure(covparam_object, class = c(estmethod, stcov))
  } else {
    # giving the object the appropriate estimation method class
    covparam_object <- structure(list(), class = c(estmethod))
  }
  return(covparam_object)
}

# covariance parameter productsum
covparam_object_productsum <- function(s_de, s_ie, t_de,
                                t_ie, st_de, st_ie,
                                s_range, t_range){
  return(c(s_de = s_de, s_ie = s_ie, t_de = t_de,
              t_ie = t_ie, st_de = st_de, st_ie = st_ie,
              s_range = s_range, t_range = t_range))
}

# covariance parameter sum with error
covparam_object_sum_with_error <- function(s_de, s_ie, t_de,
                                       t_ie, st_ie,
                                       s_range, t_range){
  return(c(s_de = s_de, s_ie = s_ie, t_de = t_de,
           t_ie = t_ie, st_ie = st_ie,
           s_range = s_range, t_range = t_range))
}

# covariance parameter product
covparam_object_product <- function(st_de, v_s, v_t,
                                       s_range, t_range){
  return(c(st_de = st_de, v_s = v_s, v_t = v_t,
           s_range = s_range, t_range = t_range))
}



