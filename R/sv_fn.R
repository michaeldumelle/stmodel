sv_fn_productsum <- function(par, sv, sp_cor, t_cor, max_srange, max_trange, max_v){
  # multiply overall var by exponentiated
  invlogit <- exp(par) / (1 + exp(par))
  p2r <- p2r_productsum(lambda = invlogit[["lambda"]],
                        alpha = invlogit[["alpha"]],
                        n_s = invlogit[["n_s"]],
                        n_t = invlogit[["n_t"]],
                        n_st = invlogit[["n_st"]],
                        ov_var = max_v * invlogit[["ov_var"]])
  rparam <- c(p2r, s_range = max_srange * invlogit[["s_range"]],
              t_range = max_trange * invlogit[["t_range"]])
  r_s <- r_make(h = sv$avg_hsp, range = rparam[["s_range"]], structure = sp_cor)
  r_t <- r_make(h = sv$avg_tsp, range = rparam[["t_range"]], structure = t_cor)

  sigma_s <- rparam[["s_de"]] * r_s + rparam[["s_ie"]] * (r_s == 1)
  sigma_t <- rparam[["t_de"]] * r_t + rparam[["t_ie"]] * (r_t == 1)
  sigma_st <- rparam[["st_de"]] * r_s * r_t + rparam[["st_ie"]] * (r_s == 1) * (r_t == 1)
  sigma <- sigma_s + sigma_t + sigma_st
  vparams <- c(rparam[["s_de"]], rparam[["s_ie"]], rparam[["t_de"]],
               rparam[["t_ie"]], rparam[["st_de"]], rparam[["st_ie"]])
  theo_sv <- sum(vparams) - sigma
  wts <- sv$n / theo_sv^2
  sumsq <- (sv$mean_sqdifs - theo_sv)^2
  return(sum(wts * sumsq))
}


weight_cressie <- function(sv, theo_sv) {
  wts <- sv$n / theo_sv^2
  return(wts)
}
