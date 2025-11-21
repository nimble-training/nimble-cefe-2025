## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                      cache = TRUE)
library(nimble)
library(nimbleMacros)
recalculate<-TRUE


## -----------------------------------------------------------------------------------------------------------------------------------------------
library(nimbleMacros)
data(mtcars) # a built-in dataset with base R
summary(mtcars) # am = 0 for automatic, 1 for manual transmission. mpg=miles per gallon
const <- mtcars[c("am", "mpg")]
pr <- setPriors(intercept = quote(dunif(-10,10)), 
                coefficient = quote(dnorm(0, 0.1)))
code <- nimbleCode({
  LM(am ~ scale(mpg), family = binomial, priors = pr) #LM is a model macro
})
mod <- nimbleModel(code = code, constants = const)
mod$getCode()


## -----------------------------------------------------------------------------------------------------------------------------------------------
set.seed(1)
x <- rnorm(10)
y <- rnorm(10, mean = 0.2 * x)
embed_dnorm <- nimbleFunction(
  setup=function(x, y) {},
  methods=list(
    dcustom=function(x=double(), alpha = double(), 
                     beta=double(), sd=double(), log=logical(0, default=FALSE)) {
      logProb <- sum(dnorm(y, mean = alpha+beta*x, sd = sd, log=TRUE))
      if(log) return(exp(logProb))
      return(logProb)
      returnType(double())
    }
  )
)
myDist <- embed_dnorm(x, y)
modelCode <- nimbleCode({
  beta ~ dnorm(0, sd=100)
  alpha ~ dnorm(0, sd = 100)
  sd ~ dunif(0, 20)
  dummy ~ myDist$dcustom(alpha, beta, sd)
})
model <- nimbleModel(modelCode,
                     data = list(dummy=1),
                     inits=list(alpha=0, beta=0, sd=1))
model$calculate()
# cmodel <- compileNimble(model) # This works too but is not run here.
# cmodel$calculate()


## ----eval=TRUE----------------------------------------------------------------------------------------------------------------------------------
weights <- c(1, .5, .25, .125)
weights <- weights / sum(weights)
centers <- c(-1, 0, 1, 2)
sds <- c(0.5, 0.6, 0.7, 0.8)
x <- seq(-3, 4, by = 0.05)
pdf <- x |> lapply(\(x) sum(weights * dnorm(x, centers, sds))) |> unlist()
{
  plot(x, pdf, type ='l', col = "red", lwd = 1.5,
       main = "Example of normal mixture (unimodal)",
       ylim = c(0, max(c(weights, pdf))))
  for(i in seq_along(weights)) {
    points(x, weights[i] * dnorm(x, centers[i], sds[i]), type = 'l' )
    points(c(centers[i], centers[i]), c(0, weights[i]), type = 'l', col = 'blue')
  }
}


## ----eval=TRUE----------------------------------------------------------------------------------------------------------------------------------
weights <- c(1, .7)
weights <- weights / sum(weights)
centers <- c(-1, 1)
sds <- c(0.5, 0.7)
x <- seq(-3, 4, by = 0.05)
pdf <- x |> lapply(\(x) sum(weights * dnorm(x, centers, sds))) |> unlist()
{
  plot(x, pdf, type ='l', col = "red", lwd = 1.5,
       main = "Example of normal mixture (bimodal)",
       ylim = c(0, max(c(weights, pdf))))
  for(i in seq_along(weights)) {
    points(x, weights[i] * dnorm(x, centers[i], sds[i]), type = 'l' )
    points(c(centers[i], centers[i]), c(0, weights[i]), type = 'l', col = 'blue')
  }
}












## -----------------------------------------------------------------------------------------------------------------------------------------------
library(nimble)
lmCode <- nimbleCode({
  psi ~ dunif(0,1)   # prior on inclusion probability
  sigma ~ dunif(0, 20)
  for(i in 1:numVars) {
    z[i] ~ dbern(psi) # indicator variable
    beta[i] ~ dnorm(0, sd = 100)
    zbeta[i] <- z[i] * beta[i]  # indicator * beta
  }
  for(i in 1:n) {
    pred.y[i] <- inprod(X[i, 1:numVars], zbeta[1:numVars])
    y[i] ~ dnorm(pred.y[i], sd = sigma)
  }
})
set.seed(1)
X <- matrix(rnorm(100*15), nrow = 100, ncol = 15)
lmConstants <- list(numVars = 15, n = 100, X = X)
lmModel <- nimbleModel(lmCode, constants = lmConstants)


## -----------------------------------------------------------------------------------------------------------------------------------------------
true_betas <- c(c(0.1, 0.2, 0.3, 0.4, 0.5),
                rep(0, 10))
lmModel$beta <- true_betas
lmModel$sigma <- 1
lmModel$z <- rep(1, 15)
lmModel$psi <- 0.5
lmModel$calculate()
set.seed(0) ## Make this reproducible
lmModel$simulate('y')
lmModel$y
lmModel$calculate() 
lmModel$setData('y')
lmData = list(y = lmModel$y)


## -----------------------------------------------------------------------------------------------------------------------------------------------
summary(lm(lmModel$y ~ lmModel$X))


## -----------------------------------------------------------------------------------------------------------------------------------------------
MCMCconf <- configureMCMC(lmModel)
MCMCconf$addMonitors('z')
MCMC <- buildMCMC(MCMCconf)
ClmModel <- compileNimble(lmModel)
CMCMC <- compileNimble(MCMC, project = lmModel, resetFunctions = TRUE)
set.seed(100)
system.time(samples_nimble <- runMCMC(CMCMC, niter = 100000, nburnin = 10000))


## -----------------------------------------------------------------------------------------------------------------------------------------------
inds <- seq(50000, 60000, by=10) ## Look at arbitrary 1001 iterations thinned
plot(samples_nimble[inds,'beta[1]'])
plot(samples_nimble[inds,'z[1]'])


## -----------------------------------------------------------------------------------------------------------------------------------------------
plot(samples_nimble[inds,'beta[4]'], ylim = c(-1, 1))
plot(samples_nimble[inds,'z[4]'])


## -----------------------------------------------------------------------------------------------------------------------------------------------
plot(samples_nimble[inds,'beta[5]'])
plot(samples_nimble[inds,'z[5]'])


## -----------------------------------------------------------------------------------------------------------------------------------------------
zCols <- grep("z\\[", colnames(samples_nimble))
posterior_inclusion_prob_nimble <- colMeans(samples_nimble[,zCols])
plot(true_betas, posterior_inclusion_prob_nimble)


## -----------------------------------------------------------------------------------------------------------------------------------------------
plot(density(samples_nimble[,'psi']))


## -----------------------------------------------------------------------------------------------------------------------------------------------
# make a new copy of the model to be totally independent
lmModel2 <- lmModel$newModel(replicate = TRUE)
MCMCconfRJ <- configureMCMC(lmModel2)
MCMCconfRJ$addMonitors('z')
configureRJ(MCMCconfRJ,
            targetNodes = 'beta',
            indicatorNodes = 'z',
            control = list(mean = 0, scale = .2))
MCMCRJ <- buildMCMC(MCMCconfRJ)


## -----------------------------------------------------------------------------------------------------------------------------------------------
ClmModel2 <- compileNimble(lmModel2)
CMCMCRJ <- compileNimble(MCMCRJ, project = lmModel2)
set.seed(100)
system.time(samples_nimble_RJ <- runMCMC(CMCMCRJ, niter = 100000, nburnin = 10000))


## -----------------------------------------------------------------------------------------------------------------------------------------------
inds <- seq(50000, 60000, by = 10)
plot(samples_nimble_RJ[inds,'beta[1]'])
plot(samples_nimble_RJ[inds,'z[1]'])


## -----------------------------------------------------------------------------------------------------------------------------------------------
plot(samples_nimble_RJ[inds,'beta[4]'])
plot(samples_nimble_RJ[inds,'z[4]'])


## -----------------------------------------------------------------------------------------------------------------------------------------------
plot(samples_nimble_RJ[inds,'beta[5]'])
plot(samples_nimble_RJ[inds,'z[5]'])


## -----------------------------------------------------------------------------------------------------------------------------------------------
zCols <- grep("z\\[", colnames(samples_nimble_RJ))
posterior_inclusion_prob_nimble_RJ <- colMeans(samples_nimble_RJ[,zCols])
plot(true_betas, posterior_inclusion_prob_nimble_RJ)


## -----------------------------------------------------------------------------------------------------------------------------------------------
plot(density(samples_nimble_RJ[,'psi']))

