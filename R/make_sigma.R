make_sigma <- function(de, r_mx, ie, v_ie, e = 1, scale = FALSE) {
  if (any(r_mx > 1)) {
    stop("Bug - Diagonal or correlation matrix different from one")
  }
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

# could make this a generic based on class of r_mx
