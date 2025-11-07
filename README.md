# nimble-cefe-2025
This repository contains materials for a one-day hybrid `nimble` workshop (short course) hosted at the Centre d'Ecologie Fonctionnelle et Evolutive (CEFE), Centre National de la Recherche Scientifique (CNRS),  Montpellier, France, 21 November 2025.

This is scheduled as a one-day workshop from 9:00 - 16:30. The morning modules will be more introductory, while the afternoon modules will be more advanced. The workshop if free but requires registration. See [here](https://docs.google.com/forms/d/e/1FAIpQLSf1Y8hXdysE_VRYPA6HUObItGgTA-3JAOjEudJ8KzaxTshblg/viewform?usp=sharing&ouid=113227867823390804287). Attendees (in person or on Zoom) are welcome to attend only the portions of interest. The Zoom link will be provided to people who register.

To prepare for the workshop:

 - Install NIMBLE (see below)
 - Install additional packages (see below)
 - Download the materials provided in this repository (when ready)

All materials for the workshop will be in this GitHub repository. If you're familiar with Git/GitHub, you already know how to get all the materials on your computer. If you're not, go [here](https://github.com/nimble-training/nimble-cefe-2025) (where you probably are if you are reading this) click the (green) "Code" button, and choose the "Download ZIP" option.

## Background for the workshop

This workshop will focus on the `nimble` R package, not on statistical methodology per se.  The materials assume attendees have basic knowledge of statistical models and some ecological statistical models.

## Tentative Schedule

Friday, November 21st

1. (9:00 - 9:50) Introduction to NIMBLE: Basic concepts and workflows
2. (10:00 - 10:50) Extending NIMBLE: Introduction to nimbleEcology
3. (11:00 - 11:50) Extending NIMBLE: Writing new functions and distributions to use in models
4. (12:00 - 13:30) Break for lunch
5. (13:30 - 13:45) Brief overview of advanced features in NIMBLE: Bayesian nonparametric distributions; reversible jump sampling for variable selection; particle filtering; parameter transformations; automatic differentiation; advanced nimbleFunction programming concepts.
6. (13:40 - 14:20) Customizing MCMC configurations and using Hamiltonian Monte Carlo (nimbleHMC)
7. (14:30-15:20) Spatial models: conditional autoregressive (CAR) and Gaussian process models
8. (15:30-16:20) Maximum likelihood: Laplace approximation and Monte Carlo expectation maximization (MCEM)


## Help with NIMBLE

The NIMBLE web site is [here](https://r-nimble.org).

The NIMBLE user manual is [here](https://r-nimble.org/manual/cha-welcome-nimble.html).

A NIMBLE "cheatsheet" is available [here](https://r-nimble.org/documentation), where you will also find links to other resources.

## Installing NIMBLE

NIMBLE is an R package available on CRAN, so in general it will be straightforward to install as with any R package. However, NIMBLE does require a compiler and related tools installed on your system.

The steps to install NIMBLE are:

1. Install compiler tools on your system. [https://r-nimble.org/download](https://r-nimble.org/download) will point you to more details on how to install *Rtools* on Windows and how to install the command line tools of *Xcode* on a Mac. Note that if you have packages requiring a compiler (e.g., *Rcpp*) on your computer, you should already have the compiler tools installed. Linux users typically already have compiler tools on their systems.

2. Install the *nimble* package from CRAN in the usual fashion of installing an R package (e.g. `install.packages("nimble")`). More details (including troubleshooting tips) can also be found in Section 4 of the [NIMBLE manual](https://r-nimble.org/manual/cha-installing-nimble.html)

3) Test that the installation is working, by running the following code in R:

```
library(nimble)
code <- nimbleCode({ x ~ dnorm(0, 1) })
model <- nimbleModel(code)
cModel <- compileNimble(model)
```

If the above code runs without error, you're all set. If not, please see the troubleshooting tips.  The most common problems have to do with proper installation of the compiler tools.  On Windows, the `PATH` system variable must be set (see link to Rtools installation details from our download linked above).  On Mac OSX, command line tools must be installed as part of Xcode.  If you still have problems, please email the [nimble-users group](https://r-nimble.org/groups-and-issues.html) for help.

In general, we encourage you to update to the most recent version of NIMBLE (version 1.3.0).

## Installing additional packages

Prior to the workshop, you should also install the following R packages (beyond those automatically installed with `nimble`), which can be installed as follows:

```
install.packages(c("nimbleHMC", "mcmcplots", "coda", "nimbleEcology", "compareMCMCs"))
```

We will use some other packages to set up various examples. To be able to run everything in the workshop materials, you will also want:

```
install.packages(c("CARBayesdata","sp","spdep","classInt", "glmmTMB"))
```
