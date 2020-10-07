make_invert_object <- function(stcov,
                               chol, co = NULL,
                               covparams, diag_tol,
                               h_s_large, h_t_large,
                               h_s_small, h_t_small,
                               n_s, n_t,
                               logdet, m_index,
                               o_index, sp_cor,
                               t_cor, xo,
                               yo, ...){
  xyc_o <- cbind(xo, yo, co)
  if (chol){
    r_s_small <- NULL
    r_t_small <- NULL
    r_s_large <-  make_r(h = h_s_large, range = covparams[["s_range"]], structure = sp_cor)
    r_t_large <- make_r(h = h_t_large, range = covparams[["t_range"]], structure = t_cor)
  } else {
    r_s_small <-  make_r(h = h_s_small, range = covparams[["s_range"]], structure = sp_cor)
    r_t_small <- make_r(h = h_t_small, range = covparams[["t_range"]], structure = t_cor)
    r_s_large <- NULL
    r_t_large <- NULL
  }
  invert_object <- structure(list(chol = chol,
                                  covparams = covparams, diag_tol = diag_tol,
                                  logdet = logdet, m_index = m_index,
                                  n_s = n_s, n_t = n_t,
                                  o_index = o_index, r_s_small = r_s_small, r_t_small = r_t_small,
                                  r_s_large = r_s_large, r_t_large = r_t_large,
                                  sp_cor = sp_cor, t_cor = t_cor,
                                  xyc_o = xyc_o),
                             class = stcov)
  return(invert_object)
}
