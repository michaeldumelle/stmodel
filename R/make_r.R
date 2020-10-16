make_r <- function(h, range, structure = c("exponential", "spherical", "gaussian", "tent")) {
  switch(structure,
         exponential = r_exp(h, range),
         spherical = r_sph(h, range),
         gaussian = r_gau(h, range),
         tent = r_tent(h, range),
         stop("Choose a valid covariance structure"))
}


r_exp <- function(h, range) {
  return(exp(- (3 * (h / range))))
}

r_sph <- function(h, range) {
  return((1 - (3 / 2) * (h / range) + (1 / 2) * (h / range)^3) * (h <= range))
}

r_gau <- function(h, range) {
  return(exp(- (3 * (h / range)^2)))
}

r_tent <- function(h, range) {
  return((1 - h/range) * (h <= range))
}

