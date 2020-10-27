make_newdata_h <- function(newdata_coord1, data_coord1, newdata_coord2 = NULL, data_coord2 = NULL, distmetric = "euclidean"){
  distmetric <- match.arg(distmetric)
  switch(distmetric,
         euclidean = eucdist_newdata(newdata_coord1, data_coord1, newdata_coord2, data_coord2),
         stop("invalid distance metric"))
}

## helper functions for each distance metric
### euclidean
eucdist_newdata <- function(newdata_coord1, data_coord1, newdata_coord2, data_coord2) {

  if (is.null(newdata_coord2) & is.null(data_coord2)) {
    return(abs(outer(newdata_coord1, data_coord1, "-")))
  } else {
    return(sqrt(outer(newdata_coord1, data_coord1, sqr_dif) +
                  outer(newdata_coord2, data_coord2, sqr_dif)))
  }
}
