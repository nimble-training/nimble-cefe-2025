## ----setup, include=FALSE-----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache = TRUE)
library(nimble)
library(coda)


## -----------------------------------------------------------------------------------------------------
library(nimble)


## ----echo = TRUE--------------------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------------------------------
DO_PLOT <- TRUE # Comment-out this line if you don't want the plots
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
                            J = ncol(y),
                            XvegHt = seq(-1, 1, length.out=100),
                            Xwind = seq(-1, 1, length.out=100)) )

# Initial values: must give for same quantities as priors given !
zst <- apply(y, 1, max)        # Avoid data/model/inits conflict
occupancy_inits <- function(){
  list(z = zst, 
       mean.p = runif(1), 
       alpha1 = runif(1), 
       mean.psi = runif(1), 
       beta1 = runif(1))
}


## ----nimbleMCMC---------------------------------------------------------------------------------------
results <- nimbleMCMC(Section10p4_code,
                       constants = occupancy_data,
                       inits = occupancy_inits,
                       niter = 10000,
                       nburnin = 1000,
                       nchains = 2,
                       samplesAsCodaMCMC = TRUE,
                       WAIC = TRUE)
summary(results$samples) ## from coda
results$WAIC


## ----eval = FALSE-------------------------------------------------------------------------------------
# library(coda)
# {
#   pdf("occupancy_samples_mcmcplot.pdf")
#   plot(results$samples)
#   dev.off()
# }


## ----echo=FALSE---------------------------------------------------------------------------------------
library(coda)
{
  pdf("orig_occupancy_samples_mcmcplot.pdf")
  plot(results$samples)
  dev.off()
}


## ----eval = TRUE--------------------------------------------------------------------------------------
code <- nimbleCode({
  mu ~ dnorm(0, sd = 100)
  sigma ~ dhalfflat()
  for(i in 1:10) {
    x[i] ~ dnorm(mu, sd = sigma)
  }
  for(i in 1:8) {
    y[1:5, i] ~ ddirch(alpha[1:5])
  }
  # mean_y and cov_y assumed to be provided
})
model <- nimbleModel(code, calculate=FALSE)
model$getVarNames()
model$getNodeNames()


## ----eval = FALSE-------------------------------------------------------------------------------------
# code <- nimbleCode({
#   # Code snippet only
#   mu ~ dnorm(0, sd = 100)
#   sigma ~ dhalfflat()
#   for(i in 1:10) {
#     x[i] ~ dnorm(mu, sd = sigma)
#   }


## ----eval = FALSE-------------------------------------------------------------------------------------
# nimbleCode({
#   # Code snippet only:
#   for(i in 1:10) {
#     x[i] ~ dnorm(mu, sd = sigma)
#   }
#   sigma ~ dhalfflat()
#   mu ~ dnorm(0, sd = 100)


## ----eval=TRUE----------------------------------------------------------------------------------------
code <- nimbleCode({
  tau <- 1e-6
  x ~ dnorm(0, tau)
})
model <- nimbleModel(code, calculate=FALSE)
model$getNodeNames()


## ----eval=FALSE---------------------------------------------------------------------------------------
# nimbleCode({
#   tau <- 1e-6
#   lifted_d1_over_sqrt_oPtau_cP <- 1/sqrt(tau) # a lifted node
#   mu ~ dnorm(0, sd = lifted_d1_over_sqrt_oPtau_cP)
# })


## ----eval=FALSE---------------------------------------------------------------------------------------
# code <- nimbleCode({
#   for(i in 1:3) y[i] ~ dnorm(a + b*x[i], sd = sigma)
# })
# model <- nimbleModel(code, calculate=FALSE)
# model$getNodeNames()


## ----eval=FALSE---------------------------------------------------------------------------------------
# nimbleCode({
#   for(i in 1:3) {
#     lifted_a_plus_b_times_x_oBi_cB_L2[i] <- a + b*x[i] # lifted nodes
#     y ~ dnorm(lifted_a_plus_b_times_x_oBi_cB_L2[i], sd = sigma)
#   }})


## -----------------------------------------------------------------------------------------------------
# Build occupancy model
occ_model <- nimbleModel(Section10p4_code,
                       constants = occupancy_data,
                       inits = occupancy_inits())
# Get and set values
occ_model$calculate() # sum of all log probabilities
occ_model$z[6:10]
occ_model$calculate("z[6:10]") # sum of log probabilities of z[6]...z[10]
occ_model$mean.psi
occ_model$mean.psi <- 0.5
occ_model$calculate()
occ_model$getDependencies("p[1, 1]") # What depends on p[1, 1]


## ----eval=FALSE---------------------------------------------------------------------------------------
# code <- nimbleCode({
#   LM(weight ~ Time + (1|Chick))
# })


## ----eval=FALSE---------------------------------------------------------------------------------------
# # only showing relevant part of model code
# intercept ~ dnorm(0, sd = 100)
# for(i in 1:k)
#   group_effect[i] ~ dnorm(0, sd = sigma_beta)
# for(i in 1:N) {
#   predicted_y[i] <- intercept + group_effect[group[i] ] + slope*x[i]


## ----eval=FALSE---------------------------------------------------------------------------------------
# # only showing relevant part of model code
# intercept ~ dnorm(0, sd = 100)
# for(i in 1:k)
#   group_effect[i] ~ dnorm(intercept, sd = sigma_beta)
# for(i in 1:N) {
#   predicted_y[i] <- group_effect[group[i] ] + slope*x[i]


## ----eval=FALSE---------------------------------------------------------------------------------------
# for(i in 1:n)
#   predicted_y[i] <- intercept + slope*x[i]


## ----eval=FALSE---------------------------------------------------------------------------------------
# predicted_y[1:n] <- intercept + slope*x[1:n]

