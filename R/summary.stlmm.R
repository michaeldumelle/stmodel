summary.stlmm <- function(object) {

  call <- object$formula

  NAvec <- rep(NA, times = length(object$NamesCoefficients))

  regcoefs <- as.vector(object$Coefficients)
  regvar <- object$CovCoefficients
  p <- ncol(object$model$FixedDesignMatrix)
  n <- nrow(object$model$FixedDesignMatrix)

  sereg <- sqrt(diag(as.matrix(regvar)))

  tvec <- NAvec
  tvec <- regcoefs / sereg
  pvec <- NAvec
  pvec <- round(100000 * (1 - pt(abs(regcoefs / sereg),
                                 df = n - p)) * 2) / 100000

  fixed.effect.estimates <- data.frame(Estimate = regcoefs,
    Std.Error = sereg, t.value = tvec, prob.t = pvec)
  rownames(fixed.effect.estimates) <- object$NamesCoefficients

  covmodels <- as.list(object$CovarianceParameters)
  covmodelout <- data.frame(covmodels, stringsAsFactors = FALSE)


  covinfo <- as.list(object$CovarianceForms)
  covinfoout <- data.frame(covinfo, stringsAsFactors = FALSE)

  resid_vec <- object$Residuals

  objective <- object$ObjectiveFn

  output <- structure(list(Call = call,
                FixedEffects = fixed.effect.estimates,
                CovarianceParameters = covmodelout,
                CovarianceForms = covinfoout,
                Residuals = resid_vec, ObjectiveFn = objective), class = "summary.stlmm")

  return(output)
}
