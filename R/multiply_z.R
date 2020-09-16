multiply_z <- function(mx, z_type, n_s) {
  # multiplication by an arbitrary dimension matrixon the left
  # we first store the total number of observations in the dense data rectangle
  n_st <- ncol(mx)
  n_row <- nrow(mx)
  if (z_type == "spatial"){
    # we take the appropriate entries in the matrix multiplied on the left and sum them
    return(vapply(1:n_s, function(a) rowSums(mx[, seq(a, n_st, by = n_s)]), double(n_row)))
  } else if (z_type == "temporal") {
    # can recover n_t with the matrix and n_s
    n_t <- n_st / n_s
    # we take the appropriate entries in the matrix multiplied on the left and sum them
    return(vapply(1:n_t, function(a) rowSums(mx[, seq(n_s * (a - 1) + 1, n_s * a, by = 1)]), double(n_row)))
  }
}


# multiply_z <- function(mx, z_type, n_s, density = "full", o_index) {
#   switch(density,
#          full = multiply_z_full(mx = mx, z_type = z_type, n_s = n_s),
#          observed = multiply_z_observed(mx = mx, z_type = z_type, n_s = n_s, o_index = o_index))
#   # taking observed vs unobserved means subsetting the "rows" in the switch functions
#   # this function is iterim, a final version will separate observed and uno
# }

# multiply_z_full <- function(mx, z_type, n_s) {
#   # we first store the total number of observations in the dense data rectangle
#   n_st <- ncol(mx)
#   n_row <- nrow(mx)
#   if (z_type == "spatial"){
#     # we take the appropriate entries in the matrix multiplied on the left and sum them
#     return(vapply(1:n_s, function(a) rowSums(mx[, seq(a, n_st, by = n_s)]), double(n_row)))
#   } else if (z_type == "temporal") {
#     # can recover n_t with the matrix and n_s
#     n_t <- n_st / n_s
#     # we take the appropriate entries in the matrix multiplied on the left and sum them
#     return(vapply(1:n_t, function(a) rowSums(mx[, seq(n_s * (a - 1) + 1, n_s * a, by = 1)]), double(n_row)))
#   }
# }

# multiply_z_observed <- function(mx, z_type, n_s, o_index) {
#   # we first store the total number of observations in the dense data rectangle
#   n_st <- ncol(mx)
#   if (z_type == "spatial"){
#     # we take the appropriate entries in the matrix multiplied on the left and sum them
#     return(vapply(1:n_s, function(a) rowSums(mx[o_index, seq(a, n_st, by = n_s)]), double(length(o_index))))
#   } else if (z_type == "temporal") {
#     # can recover n_t with the matrix and n_s
#     n_t <- n_st / n_s
#     # we take the appropriate entries in the matrix multiplied on the left and sum them
#     return(vapply(1:n_t, function(a) rowSums(mx[o_index, seq(n_s * (a - 1) + 1, n_s * a, by = 1)]), double(length(o_index))))
#   }
# }
