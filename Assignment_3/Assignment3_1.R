library(mvtnorm)

# --- Lab 3 / Assignment 1 ---
nDraws = 1000

# --- a ---
data = readRDS("Precipitation.rds")
n = length(data)
logData = log(data)
logDataMean = mean(logData)

my0 = 1
sigma0 = 1
sigma = 1 # initial value updated later
v0 = 1
tao0 = 1
gibbsSamples = matrix(0, nDraws, 2)

for (i in 1:nDraws) {
  # draw next my
  w = (n/sigma) / ((n/sigma) + (1/tao0))
  tao = 1 / ((n/sigma0) + (1/sigma))
  myN = w * logDataMean + (1 - w)*my0
  my = rnorm(1, myN, tao)
  gibbsSamples[i, 1] = my
  
  # draw next sigma
  chiDraw = rchisq(1, df=(v0+n)) # v=v0+n
  sigmaVar = (v0*sigma0 + sum((logData-my)^2)) / (n + v0)
  sigma = n * sigmaVar / chiDraw
  gibbsSamples[i, 2] = sigma
}

# Inefficiency factor, IF=1+2*sum(rho(k)), where rho(k) is the sample
# autocorrelation at lag k calculated from sample values
acorrMyLagK = acf(gibbsSamples[,1], lag=nDraws, pl=FALSE)
acorrSigmaLagK = acf(gibbsSamples[,2], lag=nDraws, pl=FALSE)
MyIF = 1 + 2*sum(acorrMyLagK$acf[-1][1:30])
SigmaIF = 1 + 2*sum(acorrSigmaLagK$acf[-1][1:30])

# plot trajectories of sampled Markov Chains
par(mfrow=c(2,2))
plot(acorrMyLagK$acf[-1], type="l", main="1a - Trajectory my",
     xlab="Interation", ylab="value")
plot(acorrSigmaLagK$acf[-1], type="l", main="1a - Trajectory sigma",
     xlab="Interation", ylab="value")


# --- b ---
# plot histogram of daily precipitation & simulated draws density

simDraws = c()
for (i in 1:nDraws) {
  simDraws[i] = exp(rnorm(1, gibbsSamples[,1][i], gibbsSamples[,2][i]))
}
hist(data, 50, main="1b - Daily precipitation - observations",
     xlab="Precipitation")
plot(density(simDraws), main="1b - Daily precipipation - simulations", 
     xlim=c(0,90), xlab="Precipitation")

