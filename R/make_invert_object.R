make_invert_object <- function(stcov,
                               chol, co = NULL,
                               covparams, diag_tol,
                               f_s, f_t,
                               h_s, h_t,
                               logdet, m_index,
                               o_index, sp_cor,
                               t_cor, xo,
                               yo, ...){
  xyc_o <- cbind(xo, yo, co)
  if (chol){
    r_s <- NULL
    r_t <- NULL
    rf_s <-  r_make(h = f_s, range = covparams[["s_range"]], structure = sp_cor)
    rf_t <- r_make(h = f_t, range = covparams[["t_range"]], structure = t_cor)
  } else {
    r_s <-  r_make(h = h_s, range = covparams[["s_range"]], structure = sp_cor)
    r_t <- r_make(h = h_t, range = covparams[["t_range"]], structure = t_cor)
    rf_s <- NULL
    rf_t <- NULL
  }
  invert_object <- structure(list(chol = chol,
                                  covparams = covparams, diag_tol = diag_tol,
                                  logdet = logdet, m_index = m_index,
                                  o_index = o_index, r_s = r_s, r_t = r_t,
                                  rf_s = rf_s, rf_t = rf_t,
                                  sp_cor = sp_cor, t_cor = t_cor,
                                  xyc_o = xyc_o),
                             class = stcov)
  return(invert_object)
}
