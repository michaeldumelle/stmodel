residuals.stlmm <- function(stlmm_object, type = "raw"){
  if (type == "raw"){
  residuals <- stlmm_object$model$Response - stlmm_object$model$FixedDesignMatrix %*% stlmm_object$Coefficients
  }
  attr(residuals, "type") <- type
  return(residuals)
}
