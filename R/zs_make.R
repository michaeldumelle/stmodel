## Return spatial random effects design matrix
##
## Functionally the same as zs_make, but makes the
## temporal assignment explicit
##
## @param sp_index spatial index
##
## @return return the spatial random effects design matrix
## @export
##
## @examples
## zs_make(rep(1:4, each = 10))
zs_make <- function(sp_index) {
  # if spatial index is not a character vector, store it as one
  if (!is.character(sp_index)) {
    sp_index <- as.character(sp_index)
  }
  # make the fixed effect spatial design matrix
  zs <- stats::model.matrix(~ sp_index - 1)
  # remove attributes - useful to not carry over
  # various attributes from matrix sums
  attributes(zs) <- list(dim = dim(zs))
  # return the matrix
  return(zs)
}
