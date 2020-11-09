make_invert_object <- function(covparam_object,
                               chol, co = NULL,
                               condition,
                               h_s_large = NULL, h_t_large = NULL,
                               h_s_small = NULL, h_t_small = NULL,
                               logdet, m_index = NULL,
                               o_index = NULL, sp_cor,
                               t_cor, xo,
                               yo){
  xyc_o <- cbind(xo, yo, co)
  if (!chol & (is.null(h_s_small) | is.null(h_t_small))){
    stop("If not using Cholesky decomposition, h_s_small and h_t_small must be provided")
  }
  if (chol & (is.null(h_s_large) | is.null(h_t_large))){
    stop("If using Cholesky decomposition, h_s_large and h_t_large must be provided")
  }
    # nrow(NULL) = NULL
    n_s <- nrow(h_s_small)
    n_t <- nrow(h_t_small)
    invert_object <- structure(list(covparams = covparam_object, chol = chol,
                                    condition = condition, logdet = logdet,
                                    h_s_small = h_s_small, h_t_small = h_t_small,
                                    h_s_large = h_s_large, h_t_large = h_t_large,
                                    o_index = o_index, m_index = m_index,
                                    n_s = n_s, n_t = n_t,
                                    sp_cor = sp_cor, t_cor = t_cor,
                                    xo = xo, yo = yo, co = co,
                                    xyc_o = cbind(xo, yo, co)), class = class(covparam_object))
  return(invert_object)
}
