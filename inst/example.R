library(stmodel) # 6/10/2021 @ 3:05 pm
set.seed(2)

# generate some locations
n_s <- 10
n_t <- 10
x <- rep(rnorm(n_s), times = n_t)
y <- rep(rnorm(n_s), times = n_t)
t <- rep(1:n_t, each = n_s)
data <- data.frame(x, y, t, response)

# generate an iid response
data$response <- rnorm(n_s * n_t)

# estimate model (reml)
test <- stmodel::stlmm(data, response ~ 1, "x", "y", "t", stcov = "productsum")
summary(test)

# using strnorm to simulate correlated response -- first make covparam vector
covparam_object <- make_covparam_object(
  s_de = 1, s_ie = 1, t_de = 1, t_ie = 1, st_de = 1, st_ie = 1, s_range = 1, t_range = 1, stcov = "productsum"
)
# then simulate the response
data$response <- as.vector(strnorm(covparam_object, mu = 0, size = 1, error = "normal", xcoord = "x", ycoord = "y", tcoord = "t",
                                   data = data, s_cor = "exponential", t_cor = "exponential"))

# estimate model (reml)
test2 <- stmodel::stlmm(data, response ~ 1, "x", "y", "t", stcov = "productsum")
summary(test2)


# estimate model (svwls)
test3 <- stmodel::stlmm(data, response ~ 1, "x", "y", "t", stcov = "productsum", estmethod = "svwls")
summary(test3)

# plot st-semivariogram
ggplot(test3$stempsv, mapping = aes(x = h_s_avg, y = gammahat, col = as.factor(h_t_avg))) +
  geom_point() +
  geom_line() +
  ylim(0, NA)



