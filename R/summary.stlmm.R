summary.stlmm <- function(stlmm_object, ...) {

  call <- stlmm_object$formula

  NAvec <- rep(NA, times = length(stlmm_object$NamesCoefficients))

  regcoefs <- as.vector(stlmm_object$Coefficients)
  regvar <- stlmm_object$CovCoefficients
  p <- ncol(stlmm_object$model$FixedDesignMatrix)
  n <- nrow(stlmm_object$model$FixedDesignMatrix)

  sereg <- sqrt(diag(as.matrix(regvar)))

  tvec <- NAvec
  tvec <- regcoefs / sereg
  pvec <- NAvec
  pvec <- round(100000 * (1 - pt(abs(regcoefs / sereg),
                                 df = n - p)) * 2) / 100000

  fixed.effect.estimates <- data.frame(Estimate = regcoefs,
    Std.Error = sereg, t.value = tvec, prob.t = pvec)
  rownames(fixed.effect.estimates) <- stlmm_object$NamesCoefficients

  covmodels <- as.list(stlmm_object$CovarianceParameters)
  covmodelout <- data.frame(covmodels, stringsAsFactors = FALSE)


  covinfo <- as.list(stlmm_object$CovarianceForms)
  covinfoout <- data.frame(covinfo, stringsAsFactors = FALSE)

  resid_vec <- stlmm_object$Residuals

  objective <- stlmm_object$ObjectiveFn

  output <- structure(list(Call = call,
                FixedEffects = fixed.effect.estimates,
                CovarianceParameters = covmodelout,
                CovarianceForms = covinfoout,
                Residuals = resid_vec, ObjectiveFn = objective), class = "summary.stlmm")

  return(output)
}
