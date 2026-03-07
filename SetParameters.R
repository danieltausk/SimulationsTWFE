ngroups = 20
ntimes = 100

theta_AR = 0.8
sd_AR = 1
sd_error = 0.2

slopes = runif(ngroups,1/2,1)
intercepts = runif(ngroups,-20,20)

# true model generating the data:
# y_{it} = c_{it} + beta x_{it} + epsilon_{it},
# i = 1, ..., ngroups, t = 1, ..., ntimes
# epsilon_{it} independent from each other and independent from all x_{it}
# beta = 0 (so in fact y_{it} = c_{it} + epsilon_{it})
# epsilon_{it} is normal with mean zero and standard deviation sd_error

# in the additive model (assumptions of two-way fixed effects model satisfied),
# c_{it} = intercepts[i] + t,
# for each i, x_{it} is the sum of the time trend c_{it} with a stationary AR(1) time series with coefficient theta_AR and residual error sd_AR

# in the multiplicative model (assumptions of two-way fixed effects model violated),
# c_{it} = slopes[i]*t and x_{it} as above
