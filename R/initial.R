initial <- function(s_de, s_ie, t_de,
                       t_ie, st_de, st_ie,
                       v_s, v_t, s_range,
                       t_range,
                       estmethod, stcov){
  initialcov <- make_covparam_object(s_de = s_de, s_ie = s_ie,
                                     t_de = t_de, t_ie = t_ie,
                                     st_de = st_de, st_ie = st_ie,
                                     v_s = v_s, v_t = v_t,
                                     s_range = s_range, t_range = t_range,
                                     estmethod = estmethod, stcov = stcov)
  return(initialcov)
}
