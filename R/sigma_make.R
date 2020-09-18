sigma_make <- function(de, r_mx, ie, v_ie, e, scale = FALSE) {
  if (!scale) {
    sigma <- de * r_mx
    diag(sigma) <- diag(sigma) + ie
    return(sigma)
  }
  else if (scale) {
    v_de <- 1 - v_ie
    sigma <- e * v_de * r_mx
    diag(sigma) <- diag(sigma) + e * v_ie
    return(sigma)
  } else {
    stop("scale must be TRUE or FALSE")
  }
}
