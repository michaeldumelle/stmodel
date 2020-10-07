print.summary.slmm <- function(x,
                                 digits = max(3L, getOption("digits") - 3L),
                                 signif.stars = getOption("show.signif.stars"), ...) {


  cat("\nCall:\n", paste(deparse(x$Call),
                         sep = "\n", collapse = "\n"),
      "\n", sep = "")

  cat("\nObjective Function:\n")
  print(x$ObjectiveFn)

  cat("\nResiduals:\n")
  resQ = c(min(x$Residuals), quantile(x$Residuals,
                                      p = c(0.25, 0.5, 0.75),
                                      na.rm = TRUE), max(x$Residuals))
  names(resQ) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(resQ, digits = digits)

  cat("\nCoefficients:\n")
  coefs = x$FixedEffects
  colnames(coefs) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
               na.print = "NA", ...)

  cat("\nCovariance Parameters:\n")
  print(x$CovarianceParameters)

  cat("\nCovariance Forms:\n")
  print(x$CovarianceForms)

}

print.slmm <- function(x,...) {
  print(summary(x,...))
}
