library(rstan)

mu = 13
sigma2 = 3
T = 300

# AR1 sim function
AR1.sim = function(mu, sigma2, T, phi) {
  x = numeric()
  x[1] = mu
  for(i in 2:T) {
    x[i] = mu + phi * (x[i-1] - mu) + rnorm(1, mean = 0, sd = sqrt(sigma2))
  }
  return(x)
}

# a)
# Simulating draws of x-values with different phis
phis = seq(-1,1,0.05)
nPhis = length(phis)
xVals = matrix(nrow=nPhis, ncol=T, 0)
for(i in 1:nPhis) {
  xVals[i,] = AR1.sim(mu, sigma2, T, phis[i])
}
# results
par(mfrow=c(2,3))
plot(xVals[1,], xlab="t", ylab="x_t", main="phi=-1")
plot(xVals[9,], xlab="t", ylab="x_t", main="phi=-0.6")
plot(xVals[15,], xlab="t", ylab="x_t", main="phi=-0.3")
plot(xVals[27,], xlab="t", ylab="x_t", main="phi=0.3")
plot(xVals[33,], xlab="t", ylab="x_t", main="phi=0.6")
plot(xVals[40,], xlab="t", ylab="x_t", main="phi=0.95")


# b)
x.phi_3 = AR1.sim(mu, sigma2, T, phi=0.3)
y.phi_95 = AR1.sim(mu, sigma2, T, phi=0.95)

ARStanModel = 'data {
  int<lower=0> N;
  vector[N] x;
}
parameters {
  real mu;
  real phi;
  real<lower=0> sigma;
}
model {
  for (n in 2:N)
    x[n] ~ normal(mu + phi * (x[n-1] - mu), sigma);
}'

# MCMC
x.fit = stan(model_code=ARStanModel, data=list(x=x.phi_3, N=T))
y.fit = stan(model_code=ARStanModel, data=list(x=y.phi_95, N=T))

# Summaries
x.summary = summary(x.fit)
y.summary = summary(y.fit)

# Get n_eff
x.nEff = x.summary$summary[,"n_eff"]
y.nEff = y.summary$summary[,"n_eff"]

# N effective samples
x.nEff.mu = x.nEff["mu"] 
# 4125

x.nEff.phi = x.nEff["phi"]
# 3741

x.nEff.sigma = x.nEff["sigma"] 
# 3894

y.nEff.mu = y.nEff["mu"] 
# 1305

y.nEff.phi = y.nEff["phi"] 
# 1322

y.nEff.sigma = y.nEff["sigma"] 
# 1600

# Get params
x.params = extract(x.fit)
y.params = extract(y.fit)

# CI:s
probs = c(0.025,0.975)

x.CI.mu <- apply(as.matrix(x.params$mu), 2, quantile, probs=probs) 
# (12.9, 13.5)

x.CI.phi <- apply(as.matrix(x.params$phi), 2, quantile, probs=probs) 
# (0.0896, 0.3148)

x.CI.sigma <- apply(as.matrix(x.params$sigma), 2, quantile, probs=probs) 
# (1.62, 1.9)

y.CI.mu <- apply(as.matrix(y.params$mu), 2, quantile, probs=probs) 
# (6.29, 14.91)

y.CI.phi <- apply(as.matrix(y.params$phi), 2, quantile, probs=probs) 
# (0.89, 0.983)

y.CI.sigma <- apply(as.matrix(y.params$sigma), 2, quantile, probs=probs) 
# (1.54, 1.82)


# Post means
x.postMean = get_posterior_mean(x.fit)
y.postMean = get_posterior_mean(y.fit)


# x post parameters
x.postMu = x.postMean[1,5] 
# 13.21

x.postPhi = x.postMean[2,5] 
# 0.20

x.postSigma = x.postMean[3,5] 
# 1.75

# y post parameters
y.postMu = y.postMean[1,5] 
# 10.86

y.postPhi = y.postMean[2,5] 
# 0.93

y.postSigma = y.postMean[3,5] 
# 1.67

#ii)
# plot density
par(mfrow=c(1,2))
plot(x=x.params$mu, y=x.params$phi,
     xlab="mu", ylab="phi", main="Posterior density, AR with phi=0.3")
plot(x=y.params$mu, y=y.params$phi,
     xlab="mu", ylab="phi", main="Posterior density, AR with phi=0.95")
