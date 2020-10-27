rescale_varinitial <- function(varparams, stcov){
  varparams <- switch(stcov,
                      "productsum" = rescale_varinitial_productsum(varparams = varparams),
                      "sum_with_error" = rescale_varinitial_sum_with_error(varparams = varparams),
                      "product" = rescale_varinitial_product(varparams = varparams),
                      stop("Not a valid stcov type for varinitial"))
  return(varparams)
}

rescale_varinitial_productsum <- function(varparams){
  varparams <- c(varparams, v_s = NA, v_t = NA)
  return(varparams)
}

rescale_varinitial_sum_with_error <- function(varparams){
  varsum <- sum(varparams)
  varparams["st_de"] <- 0
  varparams <- varsum * varparams / sum(varparams)
  varparams <- c(varparams, v_s = NA, v_t = NA)
  return(varparams)
}

rescale_varinitial_product <- function(varparams){
  v_s <- varparams[["s_ie"]] / (varparams[["s_ie"]] + varparams[["s_de"]])
  v_t <- varparams[["t_ie"]] / (varparams[["t_ie"]] + varparams[["t_de"]])
  varsum <- sum(varparams)
  varparams[which(names(varparams) != "st_de")] <- 0
  varparams <- varsum * varparams / sum(varparams)
  varparams <- c(varparams, v_s = v_s, v_t = v_t)
  return(varparams)
}
