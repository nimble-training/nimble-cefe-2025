## ----setup, include=FALSE-----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                      cache = TRUE)
library(nimble)
has_nimbleEcology <- require(nimbleEcology)
has_compareMCMCs <- require(compareMCMCs)
if(!has_nimbleEcology)
  message("This module will use nimbleEcology, which you don't have installed.")
if(!has_compareMCMCs)
  message("This module will use compareMCMCs, which you don't have installed.")


## -----------------------------------------------------------------------------------------------------
Section10p4_code <- nimbleCode({
  # Priors
  mean.p ~ dunif(0, 1)         # Detection intercept on prob. scale
  alpha0 <- logit(mean.p)      # Detection intercept
  alpha1 ~ dunif(-20, 20)      # Detection slope on wind
  mean.psi ~ dunif(0, 1)       # Occupancy intercept on prob. scale
  beta0 <- logit(mean.psi)     # Occupancy intercept
  beta1 ~ dunif(-20, 20)       # Occupancy slope on vegHt
  
  # Likelihood
  for (i in 1:M) {
    # True state model for the partially observed true state
    z[i] ~ dbern(psi[i])      # True occupancy z at site i
    logit(psi[i]) <- beta0 + beta1 * vegHt[i]
    for (j in 1:J) {
      # Observation model for the actual observations
      y[i,j] ~ dbern(p.eff[i,j])    # Detection-nondetection at i and j
      p.eff[i,j] <- z[i] * p[i,j]   # 'straw man' for WinBUGS
      logit(p[i,j]) <- alpha0 + alpha1 * wind[i,j]
    }
  }
  # Derived quantities are removed.
}
)

if(!exists("DO_PLOT"))
  DO_PLOT <- FALSE
# Choose sample sizes and prepare obs. data array y
set.seed(1)                   # So we all get same data set
M <- 100                      # Number of sites
J <- 3                        # Number of presence/absence measurements
y <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data

# Create a covariate called vegHt
vegHt <- sort(runif(M, -1, 1)) # sort for graphical convenience

# Choose parameter values for occupancy model and compute occupancy
beta0 <- 0                    # Logit-scale intercept
beta1 <- 3                    # Logit-scale slope for vegHt
psi <- plogis(beta0 + beta1 * vegHt) # Occupancy probability
# plot(vegHt, psi, ylim = c(0,1), type = "l", lwd = 3) # Plot psi relationship

# Now visit each site and observe presence/absence perfectly
z <- rbinom(M, 1, psi)        # True presence/absence

# Look at data so far
table(z)

# Plot the true system state
if(DO_PLOT) {
  par(mfrow = c(1, 3), mar = c(5,5,2,2), cex.axis = 1.5, cex.lab = 1.5)
  plot(vegHt, z, xlab="Vegetation height", ylab="True presence/absence (z)", frame = F, cex = 1.5)
  plot(function(x) plogis(beta0 + beta1*x), -1, 1, add=T, lwd=3, col = "red")
}

# Create a covariate called wind
wind <- array(runif(M * J, -1, 1), dim = c(M, J))

# Choose parameter values for measurement error model and compute detectability
alpha0 <- -2                        # Logit-scale intercept
alpha1 <- -3                        # Logit-scale slope for wind
p <- plogis(alpha0 + alpha1 * wind) # Detection probability
# plot(p ~ wind, ylim = c(0,1))     # Look at relationship

# Take J = 3 presence/absence measurements at each site
for(j in 1:J) {
  y[,j] <- rbinom(M, z, p[,j])
}
sum(apply(y, 1, max))               # Number of sites with observed presences

# Plot observed data and true effect of wind on detection probability
if(DO_PLOT) {
  plot(wind, y, xlab="Wind", ylab="Observed det./nondetection data (y)", frame = F, cex = 1.5)
  plot(function(x) plogis(alpha0 + alpha1*x), -1, 1, add=T, lwd=3, col = "red")
}
# Look at the data: occupancy, true presence/absence (z), and measurements (y)
cbind(psi=round(psi,2), z=z, y1=y[,1], y2=y[,2], y3=y[,3])

# Create factors
time <- matrix(rep(as.character(1:J), M), ncol = J, byrow = TRUE)
hab <- c(rep("A", 33), rep("B", 33), rep("C", 34))  # Must have M = 100

# Bundle and summarize data set
str( occupancy_data <- list(y = y, 
                            vegHt = vegHt,
                            wind = wind,
                            M = nrow(y),
                            J = ncol(y)))

# Initial values: must give for same quantities as priors given !
zst <- apply(y, 1, max)        # Avoid data/model/inits conflict
occupancy_inits <- function(){
  list(z = zst, 
       mean.p = runif(1), 
       alpha1 = runif(1), 
       mean.psi = runif(1), 
       beta1 = runif(1))
}


## -----------------------------------------------------------------------------------------------------
# Separate constants from data for better workflow.
occupancy_constants <- occupancy_data[c('M','J')]
occupancy_data2 <- occupancy_data[c('y','vegHt','wind')]
# Build the model
occ_model <- nimbleModel(Section10p4_code,
                         constants = occupancy_constants,
                         data = occupancy_data2,
                         inits = occupancy_inits())


## -----------------------------------------------------------------------------------------------------
# Build the MCMC configuration
occ_MCMCconf <- configureMCMC(occ_model)
# look at the samplers
occ_MCMCconf$printSamplers()


## -----------------------------------------------------------------------------------------------------
# Change the random walk samplers to slice samplers
paramNodes <- occ_model$getNodeNames(topOnly = TRUE)
paramNodes
occ_MCMCconf$removeSampler(paramNodes)
for(node in paramNodes)
  occ_MCMCconf$addSampler(target=node, type="slice")
# Alternative: occ_MCMCconf$addSampler(target=paramNodes, type="slice", targetByNode=TRUE)
occ_MCMCconf$printSamplers()


## -----------------------------------------------------------------------------------------------------
# Build the MCMC, compile and run
occ_MCMC <- buildMCMC(occ_MCMCconf)
Cocc_model <- compileNimble(occ_model)
Cocc_MCMC <- compileNimble(occ_MCMC, project=Cocc_model)
# Alternative: compiled <- compileNimble(occ_model, occ_MCMC)
samples <- runMCMC(Cocc_MCMC, niter = 10000, nburnin=1000)
summary(samples)


## -----------------------------------------------------------------------------------------------------
library(nimbleHMC)
occ_model <- nimbleModel(Section10p4_code,
                         constants = occupancy_constants,
                         data = occupancy_data2,
                         inits = occupancy_inits(),
                         buildDerivs = TRUE) # WE NEED AUTOMATIC DERIVATES FOR HMC
occ_MCMCconf <- configureMCMC(occ_model)
paramNodes <- occ_model$getNodeNames(topOnly = TRUE)
occ_MCMCconf$removeSampler(paramNodes)
occ_MCMCconf$addSampler(target=paramNodes, type="NUTS")
# Alternative: addHMC(occ_MCMCconf, paramNodes)
occ_MCMCconf$printSamplers()


## -----------------------------------------------------------------------------------------------------
occ_MCMC <- buildMCMC(occ_MCMCconf)
Cocc_model <- compileNimble(occ_model)
Cocc_MCMC <- compileNimble(occ_MCMC, project=Cocc_model) # much longer due to AD
HMCsamples <- runMCMC(Cocc_MCMC, niter=2000, nburnin=1000) # burn-in ("warmup") is typically half of niter (default)
summary(HMCsamples)


## -----------------------------------------------------------------------------------------------------
confSlice <- function(model) {
  conf <- configureMCMC(model)
  paramNodes <- model$getNodeNames(topOnly=TRUE)
  conf$replaceSamplers(paramNodes, type="slice", targetByNode=TRUE)
  conf
}
confHMC <- function(model) {
  conf <- configureMCMC(model)
  paramNodes <- model$getNodeNames(topOnly=TRUE)
  conf$removeSampler(paramNodes)
  conf$addSampler(paramNodes, type="NUTS")
  conf
}
confBarker <- function(model) {
  conf <- configureMCMC(model)
  paramNodes <- model$getNodeNames(topOnly=TRUE)
  conf$removeSampler(paramNodes)
  conf$addSampler(paramNodes, type="barker")
  conf
}


## -----------------------------------------------------------------------------------------------------
occ_model <- nimbleModel(Section10p4_code,
                         constants = occupancy_constants,
                         data = occupancy_data2,
                         inits = occupancy_inits(),
                         buildDerivs = TRUE) 
comparisons1 <- compareMCMCs(
  list(model=occ_model),
  nimbleMCMCdefs=list(myslice='confSlice',Barker='confBarker'),
  MCMCs=c('nimble','myslice','Barker'),
  MCMCcontrol = list(niter=20000, burnin=1000)
)
comparisons2 <- compareMCMCs(
  list(model=occ_model),
  nimbleMCMCdefs=list(HMC='confHMC'),
  MCMCs=c('HMC'),
  MCMCcontrol = list(niter=2000, burnin=1000)
)


## ----echo=FALSE---------------------------------------------------------------------------------------
comparisons <- c(comparisons1, comparisons2)
make_MCMC_comparison_pages(comparisons, dir="MCMC_comparisons_occ_orig")


## ----eval=FALSE---------------------------------------------------------------------------------------
# comparisons <- c(comparisons1, comparisons2)
# make_MCMC_comparison_pages(comparisons, dir="MCMC_comparisons_occ")


## -----------------------------------------------------------------------------------------------------
library(nimbleEcology)
Section10p4_code_dOcc <- nimbleCode({
  # Priors
  mean.p ~ dunif(0, 1)         # Detection intercept on prob. scale
  alpha0 <- logit(mean.p)      # Detection intercept
  alpha1 ~ dunif(-20, 20)      # Detection slope on wind
  mean.psi ~ dunif(0, 1)       # Occupancy intercept on prob. scale
  beta0 <- logit(mean.psi)     # Occupancy intercept
  beta1 ~ dunif(-20, 20)       # Occupancy slope on vegHt
  
  # Likelihood
  for (i in 1:M) {
    # True state model for the partially observed true state
    logit(psi[i]) <- beta0 + beta1 * vegHt[i]
    for (j in 1:J) {
      # Observation model for the actual observations
      logit(p[i,j]) <- alpha0 + alpha1 * wind[i,j]
    }
    y[i, 1:J] ~ dOcc_v(psi[i], p[i,1:J], len=J)
  }
  # Derived quantities are removed.
}
)


## -----------------------------------------------------------------------------------------------------
# Separate constants from data for better workflow.
occupancy_constants <- occupancy_data[c('M','J')]
occupancy_data2 <- occupancy_data[c('y','vegHt','wind')]
occ_model_dOcc <- nimbleModel(Section10p4_code_dOcc,
                         data = occupancy_data2,
                         constants = occupancy_constants,
                         inits = occupancy_inits(),
                         buildDerivs=TRUE)


## -----------------------------------------------------------------------------------------------------
library(compareMCMCs)
comparisons1_dOcc <- compareMCMCs(
  list(model=occ_model_dOcc),
  nimbleMCMCdefs=list(myslice='confSlice',Barker='confBarker'),
  MCMCs=c('nimble','myslice','Barker'),
  MCMCcontrol = list(niter=20000, burnin=1000)
)
comparisons2_dOcc <- compareMCMCs(
  list(model=occ_model_dOcc),
  nimbleMCMCdefs=list(HMC='confHMC'),
  MCMCs=c('HMC'),
  MCMCcontrol = list(niter=2000, burnin=1000)
)
comparisons_dOcc <- c(comparisons1_dOcc, comparisons2_dOcc)


## ----echo=FALSE---------------------------------------------------------------------------------------
comparisons_dOcc <- c(comparisons1_dOcc, comparisons2_dOcc)
make_MCMC_comparison_pages(comparisons_dOcc, dir="MCMC_comparisons_occ_dOcc_orig")


## ----eval=FALSE---------------------------------------------------------------------------------------
# comparisons_dOcc <- c(comparisons1_dOcc, comparisons2_dOcc)
# make_MCMC_comparison_pages(comparisons_dOcc, dir="MCMC_comparisons_occ_dOcc")

