sigma_make <- function(de, r_mx, ie, v_ie, e = 1, scale = FALSE) {
  if (!scale) {
    sigma <- de * r_mx + ie * (r_mx == 1)
    return(sigma)
  }
  else if (scale) {
    v_de <- 1 - v_ie
    sigma <- e * v_de * r_mx + e * v_ie * (r_mx == 1)
    return(sigma)
  } else {
    stop("scale must be TRUE or FALSE")
  }
}
