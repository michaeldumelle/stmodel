multiply_z <- function(mx, z_type, n_s, n_t, side = c("right", "left", "p_right", "p_left", "pz_z", "z_pz")) {
  # multiplication of an arbitrary dimension matrix by z on the right
  # we first store the total number of observations in the dense data rectangle
  n_st <- n_s * n_t
  side <- match.arg(side)
  switch(side,
         right = multiply_z_r(mx, z_type, n_s, n_t),
         left = multiply_z_l(mx, z_type, n_s, n_t),
         p_right = multiply_zp_r(mx, z_type, n_s, n_t),
         p_left = multiply_zp_l(mx, z_type, n_s, n_t),
         pz_z = multiply_zp_z(z_type, n_s, n_t),
         z_pz = multiply_z_zp(z_type, n_s, n_t),
         stop("Please choose an approporiate multiplication function for Z"))
}

multiply_z_r <- function(mx, z_type, n_s, n_t) {
  n_row <- nrow(mx)
  if (z_type == "spatial"){
    # we take the appropriate entries in the matrix multiplied on the left and sum them
    return(vapply(1:n_s, function(a) rowSums(mx[, seq(a, n_st, by = n_s)]), double(n_row)))
  } else if (z_type == "temporal") {
    # can recover n_t with the matrix and n_s
    # we take the appropriate entries in the matrix multiplied on the left and sum them
    return(vapply(1:n_t, function(a) rowSums(mx[, seq(n_s * (a - 1) + 1, n_s * a, by = 1), drop = FALSE]), double(n_row)))
  } else {
    stop("inappropriate z type")
  }
}

multiply_zp_l <- function(mx, z_type, n_s, n_t){
  n_col <- ncol(mx)
  if (z_type == "spatial"){
    return(t(vapply(1:n_s, function(a) colSums(mx[seq(a, n_st, by = n_s), , drop = FALSE]), double(n_col))))
  } else if (z_type == "temporal"){
    return(t(vapply(1:n_t, function(a) colSums(mx[seq(n_s * (a - 1) + 1, n_s * a, by = 1), , drop = FALSE]), double(n_col))))
  } else {
    stop("innapropriate Z type")
  }
}

multiply_z_l <- function(mx, z_type, n_s, n_t){
  if (z_type == "spatial"){
    return(mx[rep(seq(1, n_t), times = n_s), , drop = FALSE])
  } else if (z_type == "temporal"){
    return(mx[rep(seq(1, n_t), each = n_s), , drop = FALSE])
  } else {
    stop("innapropriate Z type")
  }
}

multiply_zp_r <- function(mx, z_type, n_s, n_t){
  if (z_type == "spatial"){
    return(mx[ , rep(seq(1, n_t), times = n_s), drop = FALSE])
  } else if (z_type == "temporal"){
    return(mx[ , rep(seq(1, n_t), each = n_s), drop = FALSE])
  } else {
    stop("innapropriate Z type")
  }
}

multiply_zp_z <- function(z_type, n_s, n_t){
  if (z_type == "spatial"){
    return(diag(n_t, nrow = n_s, ncol = n_s))
  } else if (z_type == "temporal") {
    return(diag(n_s, nrow = n_t, ncol = n_t))
  } else {
    stop("inappropriate z type")
  }
}

multiply_z_zp <- function(z_type, n_s, n_t){
  if (z_type == "spatial"){
    return(kronecker(matrix(1, nrow = n_t, ncol = n_t), diag(1, nrow = n_s, ncol = n_s)))
  } else if (z_type == "temporal") {
    return(kronecker(diag(1, nrow = n_t, ncol = n_t), matrix(1, nrow = n_s, ncol = n_s)))
  } else {
    stop("inappropriate z type")
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
