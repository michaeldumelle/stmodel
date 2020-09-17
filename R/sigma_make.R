sigma_make <- function(de = 0, r_mx, ie = 0) {
  sigma <- de * r_mx
  diag(sigma) <- diag(sigma) + ie
  return(sigma)
}
