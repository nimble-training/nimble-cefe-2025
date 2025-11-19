## ----setup, include=FALSE-----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache = TRUE)
library(nimble)
#library(compareMCMCs)
recalculate <- TRUE


## ----eval=FALSE---------------------------------------------------------------------------------------
# nimbleCode({
#   # ...model code snippet...
#   for (i in 1:M)
#     y[i, 1:J] ~ dOcc_v(probOcc = psi[i], probDetect = p[i, 1:J], len = J)
#   # ...


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
dOcc_R <- function(x, probOcc, probDetect, len, log=FALSE) {
  if (length(x) != length(probDetect))
    stop("Length of data does not match length of detection vector.")
  logProb_x_given_occupied <- sum(dbinom(x,
                                         prob = probDetect, 
                                         size = 1,
                                         log = TRUE))
  prob_x_given_unoccupied <- sum(x) == 0
  prob_x <- exp(logProb_x_given_occupied) * probOcc + 
    prob_x_given_unoccupied * (1 - probOcc)
  if (log)
    return(log(prob_x))
  return(prob_x)
}


## -----------------------------------------------------------------------------------------------------
y[9,] # A good example detection history
dOcc_R(y[9,], probOcc = 0.7, probDetect = c(0.5, 0.4, 0.3), log = TRUE)
# check the answer manually
log(0.7 * prod(dbinom(y[9,], prob = c(0.5, 0.4, 0.3), size = 1)))


## -----------------------------------------------------------------------------------------------------
dOcc <- nimbleFunction(
  run = function(x = double(1), # argument type declarations
                 probOcc = double(0),
                 probDetect = double(1),
                 len = integer(0, default = 0),
                 log = logical(0, default = 0)) {
    if (len != 0) 
      if (len != length(x))
        stop("Argument 'len' must match length of data, or be 0.")
    if (length(x) != length(probDetect))
      stop("Length of data does not match length of detection vector.")
    returnType(double(0)) # return type declaration (can be anywhere)
    logProb_x_given_occupied <- sum(dbinom(x,
                                           prob = probDetect, 
                                           size = 1,
                                           log = TRUE))
    prob_x_given_unoccupied <- sum(x) == 0
    prob_x <- exp(logProb_x_given_occupied) * probOcc + 
      prob_x_given_unoccupied * (1 - probOcc)
    if (log)
      return(log(prob_x))
    return(prob_x)
  }
)


## -----------------------------------------------------------------------------------------------------
dOcc(y[9,], probOcc = 0.7, probDetect = c(0.5, 0.4, 0.3), log = TRUE)


## ----eval=FALSE---------------------------------------------------------------------------------------
# debugonce(dOcc)
# dOcc(y[9,], probOcc = 0.7, probDetect = c(0.5, 0.4, 0.3), log = TRUE)


## -----------------------------------------------------------------------------------------------------
dOcc <- nimbleFunction(
  run = function(x = double(1), # argument type declarations
                 probOcc = double(0),
                 probDetect = double(1),
                 len = integer(0, default = 0),
                 log = logical(0, default = 0)) {
    if (len != 0) 
      if (len != length(x))
        stop("Argument 'len' must match length of data, or be 0.")
    if (length(x) != length(probDetect))
      stop("Length of data does not match length of detection vector.")
    returnType(double(0)) # return type declaration (can be anywhere)
    logProb_x_given_occupied <- sum(dbinom(x,
                                           prob = probDetect, 
                                           size = 1,
                                           log = TRUE))
    prob_x_given_unoccupied <- sum(x) == 0
    prob_x <- exp(logProb_x_given_occupied) * probOcc + 
      prob_x_given_unoccupied * (1 - probOcc)
    if (log)
      return(log(prob_x))
    return(prob_x)
  }
)


## -----------------------------------------------------------------------------------------------------
occCode2 <- nimbleCode({
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
    y[i, 1:J] ~ dOcc(probOcc = psi[i], probDetect = p[i, 1:J], len = J)
    for (j in 1:J) {
      logit(p[i,j]) <- alpha0 + alpha1 * wind[i,j]
    }
  }
})


## -----------------------------------------------------------------------------------------------------
occModel2 <- nimbleModel(occCode2,
                         constants = occupancy_data,
                         inits = occupancy_inits())


## -----------------------------------------------------------------------------------------------------
occModel2$calculate()
occModel2$calculate("y[9,]")
occModel2$psi[1:10]
occModel2$p[9,]
occModel2$psi[9] <- 0.7
occModel2$p[9,] <- c(.5,.4,.3)
occModel2$calculate("y[9,]") # same as above


## ----eval=FALSE---------------------------------------------------------------------------------------
# debug(dOcc)
# occModel2$calculate("y[9,]")
# undebug(dOcc)


## -----------------------------------------------------------------------------------------------------
occMCMC2 <- buildMCMC(occModel2)
occMCMC2$run(niter=5) # could use runMCMC(occMCMC2, niter=5, nchains=1)
as.matrix(occMCMC2$mvSamples)




## -----------------------------------------------------------------------------------------------------
C_dOcc(y[9,], probOcc = 0.7, probDetect = c(0.5, 0.4, 0.3), log = TRUE)










## -----------------------------------------------------------------------------------------------------
dist_code <- nimbleCode({
  for(i in 1:num_animals) {
    for(j in 1:num_detectors) {
      dist2[i, j] <- (sxy[i,1] - detector_xy[j,1])^2 + (sxy[i,2] - detector_xy[j,2])^2
    } # sxy are individual activity centers. detector_xy and detector locations.
  }
})

