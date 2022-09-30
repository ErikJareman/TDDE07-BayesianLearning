library("mvtnorm")

women = read.table("WomenAtWork.dat", header=TRUE)

# a)
nVar = nrow(women)
y = women$Work
X = as.matrix(women[,2:8])
nAtr <- ncol(X)
covNames=names(women[,2:8])

mu0 <- as.matrix(rep(0, nAtr))
tao = 5
sigma_pr = tao^2*diag(nAtr)

LogPostLogistic <- function(betas,y,X,mu,Sigma){
  linPred <- X%*%betas;
  logLik <- sum( linPred*y - log(1 + exp(linPred)) );
  #if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear 
  #the optimizer away from here!
  logPrior <- dmvnorm(betas, mu, Sigma, log=TRUE);
  
  return(logLik + logPrior)
}

# Initial values for beta
initVals <- rnorm(nAtr)

optiRes <- optim(initVals,LogPostLogistic,gr=NULL,y,X,mu0,sigma_pr,method=c("BFGS")
                 ,control=list(fnscale=-1),hessian=TRUE)

mode = optiRes$par
cov = -solve(optiRes$hessian)
names(mode) = covNames
print("The posterior mode is:")
print(mode)
print("Cov:")
print(cov)

# Marginal distribution for nSmallChild
NSC_mode = as.numeric(mode["NSmallChild"])
NSC_std = as.numeric(approx_PostStd["NSmallChild"])
credInterval = qnorm(p=c(0.025, 0.975), mean=NSC_mode, sd=NSC_std)
print(paste("Lower bound",
            round(credInterval[1], 6), "and the upper bound", 
            round(credInterval[2], 6)))

# Control calculations
glmModel = glm(Work ~ 0+., data = women, family = binomial)
print("Comparison:")
print(glmModel$coefficients)
print(mode)



# b)
makePred = function(data, nDraws, mean, sigma) {
  betaPred = rmvnorm(nDraws, mean=mean, sigma=sigma)
  pred = exp(t(data)%*%t(betaPred)) / (1 + exp(t(data)%*%t(betaPred)))
  return(pred)
}

woman=c(1, 20, 12, 8, 43, 0, 2)
nDraws = 10000
workPred = makePred(woman, nDraws, mode, cov)
plot(density(workPred))


# c)
makeMultiPred = function(data, nDraws, mean, sigma, n) {
  multiplePred=c()
  for (i in 1:nDraws) {
    betaPred = makePred(data, 1, mean, sigma)
    multiplePred=c(multiplePred, rbinom(1, n, betaPred))
  }
  barplot(table(multiplePred), main=paste("Predictive distribution for the number of women
  out of 11 that are working"), xlab="No. of women")
}

makeMultiPred(woman, 10000, mode, cov, 11)
