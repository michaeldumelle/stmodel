## Return a random effects design matrix
##
##
## @param index the appropriate index
##
## @return return a random effects design matrix
## @export
##
## @examples
## z_make(rep(1:4, each = 10))
z_make <- function(index) {
  # if index is not a character vector, store it as one
  if (!is.character(index)) {
    index <- as.character(index)
  }
  # make the fixed effect spatial design matrix
  z <- stats::model.matrix(~ index - 1)
  # remove attributes - useful to not carry over
  # various attributes from matrix sums
  attributes(z) <- list(dim = dim(z))
  # return the matrix
  return(z)
}
