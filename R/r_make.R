## Make correlation matrix (r)
##
## @param h distance matrix
## @param range range (correlation decay) parameter defined in terms of effective range
## @param structure covariance structure for the correlation matrix
##
## @return a correlation matrix
## @export
##
## @examples
## h <- h_make(1:3)
## r <- r_make(h, 3, "exponential")
r_make <- function(h, range, structure = c("exponential", "spherical", "gaussian", "tent")) {
  structure <- match.arg(structure)
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

