## Return a distance matrix
##
## @param coord1 first coordinate
## @param coord2 second coordinate (if applicable)
## @param distmetric distance metric
##
## @return a distance matrix based on a distance metric
## @export
##
## @examples
## h_make(1:3)
## h_make(1:3, 1:3)
h_make <- function(coord1, coord2 = NULL, distmetric = "euclidean") {
  distmetric <- match.arg(distmetric)
  switch(distmetric,
         euclidean = euc_dist(coord1, coord2),
         stop("invalid distance metric"))
}




## helper functions for each distance metric
### euclidean
euc_dist <- function(coord1, coord2 = NULL) {

  eucdist_1d <- function(x) {
    abs(outer(x, x, "-"))
  }

  sqr_dif <- function(a, b) {
    (a - b)^2
  }

  eucdist_2d <- function(x, y) {
    sqrt(outer(x, x, sqr_dif) + outer(y, y, sqr_dif))
  }

  if (is.null(coord2)) {
    return(eucdist_1d(coord1))
  } else {
    return(eucdist_2d(coord1, coord2))
  }
}
