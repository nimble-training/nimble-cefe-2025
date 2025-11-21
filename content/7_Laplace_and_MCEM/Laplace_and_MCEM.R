## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache = FALSE) # cache=TRUE can cause problems
library(nimble)
has_ggplot2 <- require(ggplot2)


## -----------------------------------------------------------------------------------------------------------------------------------------------
DeerEcervi <- read.table(file.path('..', 'examples', 'DeerEcervi', 'DeerEcervi.txt'), header = TRUE)
summary(DeerEcervi)

## Create presence/absence data from counts.
DeerEcervi$Ecervi_01 <- DeerEcervi$Ecervi
DeerEcervi$Ecervi_01[DeerEcervi$Ecervi>0] <- 1
## Set up naming convention for centered and uncentered lengths for exercises later
DeerEcervi$unctrLength <- DeerEcervi$Length
## Center Length for better interpretation
DeerEcervi$ctrLength <- DeerEcervi$Length - mean(DeerEcervi$Length)
## Make a factor version of Sex for plotting
DeerEcervi$fSex <- factor(DeerEcervi$Sex)
## Make a factor and id version of Farm
DeerEcervi$fFarm <- factor(DeerEcervi$Farm)
DeerEcervi$farm_ids <- as.numeric(DeerEcervi$fFarm)




## -----------------------------------------------------------------------------------------------------------------------------------------------
DEcode <- nimbleCode({
  for(i in 1:2) {
    # Priors for intercepts and length coefficients for sex = 1 (male), 2 (female)
    sex_int[i] ~ dnorm(0, sd = 1000)
    length_coef[i] ~ dnorm(0, sd = 1000)
  }

  # Priors for farm random effects and their standard deviation.
  farm_sd ~ dunif(0, 20)
  for(i in 1:num_farms) {
    farm_effect[i] ~ dnorm(0, sd = farm_sd)
  }

  # logit link and Bernoulli data probabilities
  for(i in 1:num_animals) {
    logit(disease_probability[i]) <-
      sex_int[ sex[i] ] +
      length_coef[ sex[i] ]*length[i] +
      farm_effect[ farm_ids[i] ]
    Ecervi_01[i] ~ dbern(disease_probability[i])
  }
})


## -----------------------------------------------------------------------------------------------------------------------------------------------
DEconstants <- list(num_farms = 24,
                    num_animals = 826,
                    length = DeerEcervi$ctrLength,
                    sex = DeerEcervi$Sex,
                    farm_ids = DeerEcervi$farm_ids)

DEmodel <- nimbleModel(DEcode,
                       constants = DEconstants,
                       buildDerivs=TRUE)
# We can set data AFTER buliding the model (but before build algorithms)
DEmodel$setData(list(Ecervi_01 = DeerEcervi$Ecervi_01))


## -----------------------------------------------------------------------------------------------------------------------------------------------
DElaplace <- buildLaplace(DEmodel)


## -----------------------------------------------------------------------------------------------------------------------------------------------
# knitr needs the model redefined here.
DEmodel <- nimbleModel(DEcode,
                       constants = DEconstants,
                       buildDerivs=TRUE)
DEmodel$setData(list(Ecervi_01 = DeerEcervi$Ecervi_01))

LAnodes <- setupMargNodes(DEmodel)
LAnodes$paramNodes # topological "top" nodes
LAnodes$randomEffectsNodes # topological "latent" nodes
LAnodes$calcNodes[c(1:3, 25:27, 851:853)] # most of the model
LAnodes$randomEffectsSets # conditionally independent random effects


## -----------------------------------------------------------------------------------------------------------------------------------------------
DEmodel <- nimbleModel(DEcode,
                       constants = DEconstants,
                       buildDerivs=TRUE)
DEmodel$setData(list(Ecervi_01 = DeerEcervi$Ecervi_01))
DElaplace <- buildLaplace(DEmodel)
comp <- compileNimble(DEmodel, DElaplace) # Slow due to AD
DEmle <- comp$DElaplace$findMLE()


## -----------------------------------------------------------------------------------------------------------------------------------------------
DEmle


## -----------------------------------------------------------------------------------------------------------------------------------------------
DEsummary <- comp$DElaplace$summary(DEmle)
DEsummary


## -----------------------------------------------------------------------------------------------------------------------------------------------
library(glmmTMB)
TMBfit <- glmmTMB(Ecervi_01 ~ ctrLength*fSex + (1|fFarm), data = DeerEcervi, family="binomial")
summary(TMBfit)
TMBcoefs <- fixef(TMBfit)$cond
## Arrange parameters from TMB into nimble's order to we can compare easily
NIMcoefs <- c(TMBcoefs[1], TMBcoefs[1]+TMBcoefs[3], TMBcoefs[2], TMBcoefs[2] + TMBcoefs[4],
              sqrt(summary(TMBfit)$varcor$cond$fFarm[1,1]))
NIMcoefs


## ----eval = TRUE--------------------------------------------------------------------------------------------------------------------------------
comp$DElaplace$calcLogLik(DEmle$par)     # Laplace approximation

comp$DElaplace$updateSettings(nQuad = 3) # Increase from 1 to 3 quadrature nodes
comp$DElaplace$calcLogLik(DEmle$par)     # Small change

comp$DElaplace$updateSettings(nQuad = 5) # 5 nodes
comp$DElaplace$calcLogLik(DEmle$par)     # Smaller change

comp$DElaplace$updateSettings(nQuad = 11) # 11 nodes
comp$DElaplace$calcLogLik(DEmle$par)     # Very small change

comp$DElaplace$updateSettings(nQuad = 15) # 5 nodes
comp$DElaplace$calcLogLik(DEmle$par)     # Tiny change


## -----------------------------------------------------------------------------------------------------------------------------------------------
DEcode2 <- nimbleCode({
  for(i in 1:2) {
    # Priors for intercepts and length coefficients for sex = 1 (male), 2 (female)
    sex_int[i] ~ dnorm(0, sd = 1000)
    length_coef[i] ~ dnorm(0, sd = 1000)
  }

  # Priors for farm random effects and their standard deviation.
  farm_sd ~ dunif(0, 1000) #gives a problem for MCEM
  for(i in 1:num_farms) {
    farm_effect[i] ~ dnorm(0, sd = farm_sd)
  }

  # logit link and Bernoulli data probabilities
  for(i in 1:num_animals) {
    disease_probability[i] <- 1e-6 + (0.99999*expit(
      sex_int[ sex[i] ] +
      length_coef[ sex[i] ]*length[i] +
      farm_effect[ farm_ids[i] ]))
    Ecervi_01[i] ~ dbern(disease_probability[i])
  }
})
if(!require(mcmcse))
  stop("You'll need mcmcse to run MCEM.")
DEmodel <- nimbleModel(DEcode2,
                       constants = DEconstants,
                       inits = list(farm_effect=rnorm(24,0,0.1),  # "friendly" inits (small sd)
                                    farm_sd=0.1, sex_int=c(0,0), length_coef=c(0,0)),
                       buildDerivs=TRUE)
DEmodel$setData(list(Ecervi_01 = DeerEcervi$Ecervi_01))
DEmcem <- buildMCEM(DEmodel) # The warning is not a problem
comp <- compileNimble(DEmodel, DEmcem) # Slow due to AD
#comp$DEmodel$farm_effect <- rnorm(24, 0, .1)
DEmle <- comp$DEmcem$findMLE(c(0,0,0,0,0.1), initM=2000, maxIter = 30)

