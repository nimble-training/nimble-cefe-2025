## ----setup, include=FALSE-----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache = TRUE)
library(nimble)
#library(compareMCMCs)
recalculate <- TRUE


## ----eval=FALSE---------------------------------------------------------------------------------------
# x[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N],
#                      tau, c, zero_mean)




## ----eval=FALSE---------------------------------------------------------------------------------------
# num <- c(2, 2, 2, 2)  # each location has two neighbors
# adj <- c(2, 3,        # neighbors of loc'n 1
# 	      1, 4,         # neighbors of loc'n 2
# 	      1, 4,         # etc.
# 	      2, 3)

