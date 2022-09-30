# --- Lab 3 / Assignment 2 ---
library(mvtnorm)
data = read.table("eBayNumberOfBidderData.dat", header=TRUE)
data = data[,-2]



# --- a ---
eBayModel = glm(data$nBids~., data=data, poisson)
betas = eBayModel$coefficients



# --- b ---
Y = data[,1] # nBids
X = as.matrix(data[-1]) # features (rest)

# Zellner's g-prior
my = rep(0, dim(X)[2])
sigma = 100*solve(t(X)%*%X)


logPost = function(beta, X, Y, my, sigma) {
  # log likelihood
  tmp = beta%*%t(X)
  logLH = sum(Y*tmp - exp(tmp))
  
  # log prior
  logP = dmvnorm(x=beta, mean=my, sigma=sigma, log=TRUE)
  return(logLH + logP)
}

# use optim() for beta and hessian values
init = rep(0, dim(X)[2])
optimRes = optim(init, logPost, gr=NULL,  X, Y, my, sigma, method=c("BFGS"), 
               control=list(fnscale=-1), hessian=TRUE)

# optValue = optimRes$value
optBetas = optimRes$par
optHessians = -solve(optimRes$hessian)



# --- c ---
runMetropolis = function(nRuns, c, beta, hessian, postFunc, ...) {
  
  acceptedDraws = matrix(0, nrow=nRuns, ncol=length(beta))
  prevProposal = beta
  
  for (i in 1:nRuns) {
    # draw from proposal density
    proposal = rmvnorm(1, prevProposal, hessian*c)
    
    # draw from supplied posterior density
    drawRatio = min(1, exp(postFunc(proposal, ...) - postFunc(prevProposal, ...)))
    
    # draw a random ratio
    randomRatio = runif(1, 0, 1)
    
    # save proposal if OK, else stay on same
    if (drawRatio >= randomRatio) {
      prevProposal = proposal
    }
    acceptedDraws[i,] = prevProposal
  }
  return(acceptedDraws)
}

# use function to sample from earlier poisson regression for eBay data
nRuns = 1000
metroSample = runMetropolis(nRuns=nRuns, c=1, beta=my, hessian=optHessians, 
                       postFunc=logPost, X, Y, my, sigma)

# Assess MCMC convergence by graphical methods
par(mfrow=c(3,3))
for(i in 1:8){
  plot(1:nRuns, type='l', metroSample[,i], 
       main=paste(c("2c - Trajectory for covariate", i), collapse= " "),
       xlab="Iteration", ylab="Value")
}



# --- d ---
newData = c(1, 0, 1, 0, 1, 0, 1.2, 0.8)
newBetas = metroSample[nRuns,]
newLambda = exp(newData%*%newBetas)

steps = 0:5
newDist = dpois(steps, newLambda)
plot(steps, type="h", newDist, main="Probability for #bidders", 
     xlab="#Bidders", ylab="Probability")
cat(newDist[1]) # About 75% chance of 0 bidders
