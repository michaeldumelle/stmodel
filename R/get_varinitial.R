get_varinitial <- function(stempsv, stcov){
  inf_inf <- max(stempsv$gammahat)
  zeroplus_inf <- stempsv$gammahat[stempsv$h_s_avg == min(stempsv$h_s_avg[stempsv$h_s_avg > 0]) &
                                  stempsv$h_t_avg == max(stempsv$h_t_avg)]
  zero_inf <- stempsv$gammahat[stempsv$h_s_avg == 0 & stempsv$h_t_avg == max(stempsv$h_t_avg)]
  inf_zeroplus <- stempsv$gammahat[stempsv$h_t_avg == min(stempsv$h_t_avg[stempsv$h_t_avg > 0]) &
                                     stempsv$h_s_avg == max(stempsv$h_s_avg)]
  inf_zero <- stempsv$gammahat[stempsv$h_t_avg == 0 & stempsv$h_s_avg == max(stempsv$h_s_avg)]
  zeroplus_zeroplus <- stempsv$gammahat[stempsv$h_s_avg == min(stempsv$h_s_avg[stempsv$h_s_avg > 0]) &
                                         stempsv$h_t_avg == min(stempsv$h_t_avg[stempsv$h_t_avg > 0])]
  s_de <- inf_inf - zeroplus_inf
  s_ie <- zeroplus_inf - zero_inf
  t_de <- inf_inf - inf_zeroplus
  t_ie <- inf_zeroplus - inf_zero
  st_de <- inf_inf - zeroplus_zeroplus - s_de - t_de
  st_ie <- inf_inf - (s_de + s_ie + t_de + t_ie + st_de)
  varparams <- pmax(1e-10, c(s_de, s_ie, t_de, t_ie, st_de, st_ie))
  numzero <- sum(varparams == 1e-10)
  varparams <- (inf_inf - numzero * 1e-2) * varparams / (sum(varparams))
  varparams <- pmax(1e-2, varparams)
  names(varparams) <- c("s_de", "s_ie", "t_de", "t_ie", "st_de", "st_ie")
  varparams <- rescale_varinitial(varparams = varparams, stcov = stcov)
  return(varparams)
}
