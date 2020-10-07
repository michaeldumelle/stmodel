residuals.stlmm <- function(stlmm_object, resid_type = "raw"){
  if (resid_type == "raw"){
  residuals <- stlmm_object$model$Response - stlmm_object$model$FixedDesignMatrix %*% stlmm_object$Coefficients
  }
  attr(residuals, "resid_type") <- resid_type
  return(residuals)
}
