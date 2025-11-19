## ----setup, include=FALSE-----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                      cache = TRUE)
library(nimble)
library(coda)
has_nimbleEcology <- require(nimbleEcology)
has_compareMCMCs <- require(compareMCMCs)
if(!has_nimbleEcology)
  message("This module will use nimbleEcology, which you don't have installed.")
doDHMMexample <- FALSE


## ----echo=FALSE---------------------------------------------------------------------------------------
if(!has_compareMCMCs) message("To run the code as given, you'll need to install compareMCMCs (see workshop installation information).  Otherwise you can run code via normal nimble workflows.")


## ----eval=FALSE---------------------------------------------------------------------------------------
# for (i in 1:M) {
#   z[i] ~ dbern(psi[i])      # True occupancy z at site i
#   logit(psi[i]) <- beta0 + beta1 * vegHt[i]
#   for (j in 1:J) {
#     y[i,j] ~ dbern(z[i] * p[i,j]) # Detection-nondetection at i and j
#     logit(p[i,j]) <- alpha0 + alpha1 * wind[i,j]
#   }
# }


## ----eval=FALSE---------------------------------------------------------------------------------------
# for (i in 1:M) {
#   logit(psi[i]) <- beta0 + beta1 * vegHt[i]
#   for (j in 1:J) {
#     logit(p[i,j]) <- alpha0 + alpha1 * wind[i,j]
#   }
#   y[i,1:J] ~ dOcc(psi[i], p[i,1:J])
# }


## ----eval=FALSE---------------------------------------------------------------------------------------
# for (i in 1:M) {
#   z[i] ~ dbern(psi[i])      # True occupancy z at site i
#   logit(psi[i]) <- beta0 + beta1 * vegHt[i]
#   for (j in 1:J) {
#     y[i,j] ~ dbern(z[i] * p[i,j]) # Detection-nondetection at i and j
#     logit(p[i,j]) <- alpha0 + alpha1 * wind[i,j]
#   }
# }


## ----eval=FALSE---------------------------------------------------------------------------------------
# y[i, 1:T] ~ dOcc_v(probOcc = psi[i],
#                    probDetect = p[i, 1:T], len = T)


## ----eval=FALSE---------------------------------------------------------------------------------------
# ## Not indexing transition and observation probabilities by individual, for simplicity of example
# for(t in 1:T) {
#   for(i in 1:num_states) {
#     for(j in 1:num_states) {
#       T[i, j, t] <- some_calculation() # Involving survival and transition probabilities
#                                        # that may depend on explanatory variables or random effects
#     }
#   }
# }
# ## Some similar definition of O, observation probabilities
# for(k in 1:num_individuals) {
#   for(t in 2:T) {
#     z[k, t] ~ dcat(T[z[k, t-1], 1:num_states, t])
#   }
#   for(t in 1:T) {
#     y[k, t] ~ dcat(O[z[k, t], 1:num_states, t])
#   }
# }




















## ----eval=FALSE---------------------------------------------------------------------------------------
# # Code for your own results:
# make_MCMC_comparison_pages(c(orchids_result, orchids_result_faster, orchids_result_DHMM),
#                            dir = "your_orchid_results",
#                            modelName = "orchids")

