make_covest_object <- function(initial, ...){
  UseMethod("make_covest_object", object = initial)
}

make_covest_object.svwls <- function(initial, chol = NULL, diag_tol = NULL,
                                     f_s = NULL, f_t = NULL,
                                     h_s = NULL, h_t = NULL,
                                     logdet = NULL,
                                     m_index = NULL, max_srange,
                                     max_trange, max_v,
                                     o_index = NULL, sp_cor,
                                     sv = NULL, t_cor,
                                     weights = NULL, xo = NULL, yo = NULL, ...){
  xyc_o <- cbind(xo, yo)
  covest_object <- structure(list(chol = chol, diag_tol = diag_tol,
                                  f_s = f_s, f_t = f_t,
                                  h_s = h_s, h_t = h_t,
                                  initial = initial,
                                  logdet = logdet,
                                  m_index = m_index, max_srange = max_srange,
                                  max_trange = max_trange, max_v = max_v,
                                  o_index = o_index, sp_cor = sp_cor,
                                  sv = sv, t_cor = t_cor,
                                  weights = weights, xo = xo,
                                  xyc_o = xyc_o, yo = yo),
                             class = class(initial))
  covest_object$initial_plo <-  r2plo.svwls(covest_object)
  return(covest_object)
}


make_covest_object <- function(stcov,
                               chol = NULL, diag_tol = NULL,
                               f_s = NULL, f_t = NULL,
                               h_s = NULL, h_t = NULL,
                               initial, logdet = NULL,
                               m_index = NULL, max_srange,
                               max_trange, max_v,
                               o_index = NULL, sp_cor,
                               sv = NULL, t_cor,
                               weights = NULL, xo = NULL, yo = NULL, ...){
  xyc_o <- cbind(xo, yo)
  covest_object <- structure(list(chol = chol, diag_tol = diag_tol,
                                  f_s = f_s, f_t = f_t,
                                  h_s = h_s, h_t = h_t,
                                  initial = initial,
                                  logdet = logdet,
                                  m_index = m_index, max_srange = max_srange,
                                  max_trange = max_trange, max_v = max_v,
                                  o_index = o_index, sp_cor = sp_cor,
                                  sv = sv, t_cor = t_cor,
                                  weights = weights, xo = xo,
                                  xyc_o = xyc_o, yo = yo),
                             class = stcov)
  covest_object$initial_plo <-  r2plo_sv(covest_object)
  return(covest_object)
}


