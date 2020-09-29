stmlist_make <- function(formula, estmethod, xcoord, ycoord, tcoord, spcor, tcor, data, ...){

}

stmlist_make_lik <- function(formula, xcoord, ycoord, tcoord, spcor, tcor, data, logdet = TRUE,
                             diagtol = 1e-4, ...)

stmlist_make_svwls <- function(formula, estmethod, xcoord, ycoord, tcoord, spcor, tcor, data, logdet = TRUE,
                               diagtol = 1e-4, ...)



test1 <- function(x, ...) {
  test2(x, ...)
}

test2 <- function(x, logdet = TRUE, ...) {
  return(list(x = x, logdet = logdet))
}

test1(1:4)
test1(1:4, logdet = FALSE)
