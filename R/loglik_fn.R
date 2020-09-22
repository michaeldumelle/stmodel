loglik_fn_productsum <- function(par, xo, dim_xo, yo, dim_yo, n_s, n_t,
                      h_s, h_t, sp_cor, t_cor, max_srange, max_trange, o_index, m_index, chol = FALSE,
                      diag_tol, estmethod){
  invlogit <- exp(par) / (1 + exp(par))
  p2r <- p2r_productsum(lambda = invlogit[["lambda"]],
                        alpha = invlogit[["alpha"]],
                        n_s = invlogit[["n_s"]],
                        n_t = invlogit[["n_t"]],
                        n_st = invlogit[["n_st"]],
                        ov_var = 1)
  rparam <- c(p2r, s_range = max_srange * invlogit[["s_range"]],
              t_range = max_trange * invlogit[["t_range"]])
  r_s <- r_make(h = h_s, range = rparam[["s_range"]], structure = sp_cor)
  r_t <- r_make(h = h_t, range = rparam[["t_range"]], structure = t_cor)
  if (!chol) {
  ps_output <- invert_productsum(r_s = r_s, r_t = r_t, s_de = rparam["s_de"],
                    s_ie = rparam["s_ie"], t_de = rparam["t_de"],
                    t_ie = rparam["t_ie"], st_de = rparam["st_de"],
                    st_ie = rparam["st_ie"], o_index = o_index,
                    m_index = m_index, xyc_o = cbind(xo, yo), diag_tol = diag_tol,
                    logdet = TRUE)
  } else {
    ps_output <- invert_chol_productsum(r_s = r_s, r_t = r_t, s_de = rparam["s_de"],
                           s_ie = rparam["s_ie"], t_de = rparam["t_de"],
                           t_ie = rparam["t_ie"], st_de = rparam["st_de"],
                           st_ie = rparam["st_ie"], xyc_o = cbind(xo, yo),
                           diag_tol = diag_tol, logdet = TRUE)
  }
  m2ll <- minus2loglik(dim_xo = dim_xo, dim_yo = dim_yo, siginv_xyc_o = ps_output$siginv_o,
                       logdet_o = ps_output$logdet, diag_tol = diag_tol, estmethod = estmethod)
  return(m2ll)
}
