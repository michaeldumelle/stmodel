covest_productsum <- function(iv, s_range_max, t_range_max, v_max,
                              sp_cor, t_cor, estmethod, xo, yo, h_s, h_t, n_s, n_t, chol = FALSE,
                              diag_tol, sv){
  # initial values must be a named vector


  prof_iv <- r2p_productsum(s_de = iv[["s_de"]], s_ie = iv[["s_ie"]],
                             t_de = iv[["t_de"]], t_ie = iv[["t_ie"]],
                             st_de = iv[["st_de"]], st_ie = iv[["st_ie"]],
                             estmethod = estmethod)
  if (estmethod == "sv") {
    pinit_vc[["ov_var"]] <- pinit_vc[["ov_var"]] / max_v
  }
  pinit_rc <- c(s_range = rinit_srange / max_srange, t_range = rinit_trange / max_trange)
  lopinit <- log(c(pinit_vc, pinit_rc) / (1 - c(pinit_vc, pinit_rc)))


  # likelihod function
  if (estmethod %in% c("reml", "ml")){
    optim(par = lopinit, fn = loglik_fn_productsum, xo = xo, dim_xo = dim_xo,
          yo = yo, dim_yo = dim_yo, n_s = n_s, n_t = n_t,
          h_s = h_s, h_t = h_t, sp_cor = sp_cor, t_cor = t_cor, max_srange = max_srange,
          max_trange = max_trange, o_index = o_index, m_index = m_index, chol = chol,
          diag_tol = diag_tol, estmethod = estmethod, method = "Nelder-Mead")
  }
  # semivariogram function
  if (estmethod == "sv"){
    optim(par = lopinit, fn = sv_fn_productsum, sv = sv, sp_cor = sp_cor,
          t_cor = t_cor, max_srange = max_srange, max_trange = max_trange, max_v = max_v,
          method = "Nelder-Mead")
  }
}






















covest_productsum <- function(rinit_vc, rinit_srange, max_srange, rinit_trange, max_trange,
                              xo, dim_xo, yo, dim_yo, n_s, n_t, h_s, h_t, sp_cor, t_cor, max_v, chol = FALSE,
                              diag_tol, sv, estmethod){



  pinit_vc <- r2p_productsum(rinit_vc[1], rinit_vc[2], rinit_vc[3],
                             rinit_vc[4], rinit_vc[5], rinit_vc[6], estmethod)
  if (estmethod == "sv") {
    pinit_vc[["ov_var"]] <- pinit_vc[["ov_var"]] / max_v
  }
  pinit_rc <- c(s_range = rinit_srange / max_srange, t_range = rinit_trange / max_trange)
  lopinit <- log(c(pinit_vc, pinit_rc) / (1 - c(pinit_vc, pinit_rc)))


  # likelihod function
  if (estmethod %in% c("reml", "ml")){
  optim(par = lopinit, fn = loglik_fn_productsum, xo = xo, dim_xo = dim_xo,
        yo = yo, dim_yo = dim_yo, n_s = n_s, n_t = n_t,
        h_s = h_s, h_t = h_t, sp_cor = sp_cor, t_cor = t_cor, max_srange = max_srange,
        max_trange = max_trange, o_index = o_index, m_index = m_index, chol = chol,
        diag_tol = diag_tol, estmethod = estmethod, method = "Nelder-Mead")
  }
  # semivariogram function
  if (estmethod == "sv"){
    optim(par = lopinit, fn = sv_fn_productsum, sv = sv, sp_cor = sp_cor,
          t_cor = t_cor, max_srange = max_srange, max_trange = max_trange, max_v = max_v,
          method = "Nelder-Mead")
  }
}
