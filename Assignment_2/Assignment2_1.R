# --- Lab 2 / Assignment 1 --

# --- General ---
library(mvtnorm)

data = read.table("TempLambohov.txt", header=TRUE)

getBetas = function(my, sigma2, omega0Inv) {
  return(rmvnorm(1, mean=my, sigma=sigma2*omega0Inv))
}

getReg = function(betaVals, time) {
  return(betaVals[1] + betaVals[2]*time + betaVals[3]*(time^2))
}

par(mfrow=c(2,3))

# --- a ---
nDraws = 100
my0 = c(-10,100,-90) # -100 --> -90
v0 = 365
sigma02 = 1 # 2 --> 1
omega0 = 0.2 * diag(3)
omega0Inv = solve(omega0)
draw = rchisq(nDraws, v0)
sigma2 = v0 * sigma02 / draw

plot.new()
plot.window(xlim=c(0,1), ylim=c(-30, 40))
axis(side=1)
axis(side=2)

betas = matrix(0, nDraws, 3)
for (i in 1:nDraws) {
  betas[i,]=getBetas(my0, sigma2[i], omega0Inv)
  lines(data$time, getReg(betas[i,], data$time))
  lines(i, 1)
}
title(main="a - joint prior regression curves ", xlab="Time", ylab="Temp")


# --- b ---
X = as.matrix(data.frame(const=1, t=data$time, t2=data$time^2))
y = data$temp
n = length(data$time)
nV = v0 + n
nOmega = t(X)%*%X+omega0
nOmegaInv = solve(nOmega)
betaHat = solve(t(X)%*%X)%*%t(X)%*%y
nMy = solve(t(X)%*%X+omega0)%*%(t(X)%*%X%*%betaHat+omega0%*%my0)
nSigma2 = as.vector((v0*sigma02 + (t(y)%*%y+t(my0)%*%omega0%*%my0-t(nMy)%*%nOmega%*%nMy))/nV)

postSigma2 = nV * nSigma2 / draw
postBetas = matrix(0, nDraws, 3)
for (i in 1:nDraws) {
  postBetas[i,] = getBetas(nMy, postSigma2[i], nOmegaInv)
}

# --- b.i ---
# Marginal posteriors for betas
hist(postBetas[,1], breaks=50, main="b.i - B0 Marginal Posterior")
hist(postBetas[,2], breaks=50, main="b.i - B1 Marginal Posterior")
hist(postBetas[,3], breaks=50, main="b.i - B2 Marginal Posterior")
hist(postSigma2, breaks=50, main="b.i - sigma2 Marginal Posterior")

# --- b.ii ---
# calculate posterior mean temp and plot it as overlay over sample data
postTemp = matrix(0, nDraws, n)
for (i in 1:nDraws) {
  postTemp[i,] = getReg(postBetas[i,], data$time)
}
plot(data$time, y, main="b.ii - Observations, CI95, mean, maxTempTime", xlab="Time", ylab="Temp")
lines(data$time, colMeans(postTemp))

# calculate upper/lower bound of CI and overlay those as well
CI = matrix(0, n, 2)
for (i in 1:n) {
  CI[i,] = quantile(postTemp[,i], probs=c(0.025, 0.975))
}

lines(data$time, CI[, 1])
lines(data$time, CI[, 2])


# --- c ---
# due to negative quadratic polynomial function we have unique max where f'=0.
# f'=b2+2*B3*time=0 -> time=-B2/(2*B3)
maxTempTime = -postBetas[,2]/(2*postBetas[,3])
abline(v=mean(maxTempTime))


# --- d ---
'To avoid overfitting, my0 and omega0 can be set in the following ways:
 - my0 = 0: this parameter can be set to 0 to make the parameters
   stay closer to 0.
 - omega0 = lambda*I: the parameter omega0 can be set this way
   to give a smaller spread in the posterior distrubution. This
   way the parameters stay somewhat close to my0'

