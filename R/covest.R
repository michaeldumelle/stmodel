covest_productsum <- function(rinit_vc, rinit_srange, max_srange, rinit_trange, max_trange,
                              xo, dim_xo, yo, dim_yo, n_s, n_t, h_s, h_t, sp_cor, t_cor, chol = FALSE,
                              diag_tol, estmethod){



  pinit_vc <- r2p_productsum(rinit_vc[1], rinit_vc[2], rinit_vc[3],
                             rinit_vc[4], rinit_vc[5], rinit_vc[6], estmethod)
  pinit_rc <- c(s_range = rinit_srange / max_srange, t_range = rinit_trange / max_trange)
  lopinit <- log(c(pinit_vc, pinit_rc) / (1 - c(pinit_vc, pinit_rc)))



  optim(par = lopinit, fn = loglik_fn_productsum, xo = xo, dim_xo = dim_xo,
        yo = yo, dim_yo = dim_yo, n_s = n_s, n_t = n_t,
        h_s = h_s, h_t = h_t, sp_cor = sp_cor, t_cor = t_cor, max_srange = max_srange,
        max_trange = max_trange, o_index = o_index, m_index = m_index, chol = chol,
        diag_tol = diag_tol, estmethod = estmethod, method = "Nelder-Mead")
}
