## Return temporal random effects design matrix
##
## Functionally the same as zs_make, but makes the
## temporal assignment explicit
##
## @param temp_index temporal index
##
## @return return the temporal random effects design matrix
##
## @export
##
## @examples
## zt_make(rep(1:4, times = 10))
zt_make <- function(temp_index) {
  # if spatial index is not a character vector, store it as one
  if (!is.character(temp_index)) {
    temp_index <- as.character(temp_index)
  }
  # make the fixed effect spatial design matrix
  zt <- stats::model.matrix(~ temp_index - 1)
  # remove attributes - useful to not carry over
  # various attributes from matrix sums
  attributes(zt) <- list(dim = dim(zt))
  # return the matrix
  return(zt)
}
