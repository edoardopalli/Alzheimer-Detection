set.seed(123456)
logistic <- function(x) 1 / (1 + exp(-x))
p <- function(x) logistic(1 - 3 * sin(x))
X <- runif(n = 200, -3, 3)
Y <- rbinom(n = 200, size = 1, prob = 0.4)

h <- 1
x <- seq(-3, 3, l = 501)

fitGlm <- sapply(x, function(x) {
  K <- dnorm(x = x, mean = X, sd = h)
  glm.fit(x = cbind(1, X - x), y = Y, weights = K,
          family = binomial())$coefficients
})

p.hat.test <- NULL
for (i in 1:length(x)) {
  K <- dnorm(x = X, mean = x[i], sd = 1)
  mdl <- glm(Y ~ X, family = binomial, weights = K)
  p.hat.test <- c(p.hat.test, predict(mdl, newdata = data.frame(X = x), type='response')[i])
}

x11()
plot(x, p.hat.test)

fit_locfit <- locfit(Y ~ lp(X, deg = 1, nn = h),
                     family = "binomial", kern = "gauss")
plot(fit_locfit)

