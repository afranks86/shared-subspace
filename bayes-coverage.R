rm(list=ls())

## Need these for Rscript
library("methods")
library("utils")
library(rstiefel)
library(sp)

source("fit-subspace.R")
source("generateData.R")
source("helper.R")

########################################################
############# Data Generation #########################
########################################################

n <- 50
S <- 2
R <- 2
P <- 200
ngroups <- 5

niters <- 1000
nwarmup <- niters/2

LambdaList <- Olist <- list()

LambdaList[[1]] <- c(100, 10)
LambdaList[[2]] <- c(100, 10)
LambdaList[[3]] <- c(100, 100/3)
LambdaList[[4]] <- c(100, 100/3)
LambdaList[[5]] <- c(100, 100)

Olist[[1]] <- matrix(c(sqrt(2)/2, sqrt(2)/2, -sqrt(2)/2, sqrt(2)/2), nrow=2)
Olist[[2]] <- matrix(c(-sqrt(2)/2, sqrt(2)/2, sqrt(2)/2, sqrt(2)/2), nrow=2)
Olist[[3]] <- matrix(c(-sqrt(2)/2, sqrt(2)/2, sqrt(2)/2, sqrt(2)/2), nrow=2)
Olist[[4]] <- diag(2)
Olist[[5]] <- diag(2)

nsim <- 100

for(i in 1:nsim) {

    dat <- generateData(P=P, S=S, R=2, ngroups=ngroups,
                        LambdaList=LambdaList, Olist=Olist,
                        nvec=rep(n, ngroups))

    dat$genType <- "Shared subspace"

    Ypooled <- c()
    for(i in 1:length(dat$Ylist)){
        Ypooled <- rbind(Ypooled, dat$Ylist[[i]])
    }
    dat$Ypooled <- Ypooled


########################################################
############# Data Fit #########################
########################################################

####### FIT DATA USING SHARED SUBSPACE ###############

    ## Initialize sampler
    ngroups <- dat$ngroups
    nvec <- dat$nvec
    Slist <- dat$Slist

    Vinit <- svd(do.call(cbind, lapply(1:ngroups, function(k) svd(t(dat$Ylist[[k]]))$u[, 1:R])))$u[, 1:S]
    EMFit <- subspaceEM(dat$Slist, P=P, S=S, R=R, nvec=dat$nvec,
                        PrecVec=rep(10, ngroups),
                        rho1=0.1, rho2=0.9, 
                        Vstart=Vinit, verbose=TRUE, maxIters=100)
    Vinit <- EMFit$V
    
    OmegaList <- Ulist <- list()
    for(k in 1:ngroups) {
        OmegaList[[k]] <- rep(1/2, R)
        
        Ok <- matrix(0, nrow=S, ncol=R)
        Ok[1:R, 1:R] <- diag(R)
        
        Ulist[[k]] <- Vinit %*% Ok    
    }
    s2vec <- rexp(ngroups)

    initSS <- list(V=Vinit, Ulist=Ulist, OmegaList=OmegaList, s2vec=s2vec)
    samples <- fitSubspace(dat$P, S, R, dat$Slist,
                           dat$nvec, dat$ngroups, init=initSS,
                           niters=niters, sigmaTruthList=dat$SigmaList,
                           draw=c(V=FALSE))

    sv <- svd(t(samples$V) %*% dat$V)
    Rot <- sv$v %*% t(sv$u)

    attach(samples)
    rotOsamps <- array(dim=dim(Osamps))
    for(i in 1:1000) {
        for(k in 1:ngroups) {
            rotOsamps[1:2, 1:2, k, i] <- t(Rot) %*% Osamps[1:2, 1:2, k, i]
        }
    }

    posteriorPlot(rotOsamps, omegaSamps[1:2, , ],
                  s2samps, nsamps=500, groups=1:ngroups, ymax=20,
                  probRegion=0.95, hline=1, plotPoints=FALSE)

    covered <- c()
    for(k in 1:ngroups) {

        eigRot <- eigen(t(Rot) %*%
                        t(samples$V) %*% dat$SigmaList[[k]] %*% samples$V %*% Rot)

        evals <- eigRot$values
        evecs <- eigRot$vectors

        truth <- c(atan(evecs[2]/evecs[1]), evals[1]/evals[2])

        hp <- getHullPoints(1000, omegaSamps[, k, ], rotOsamps[, , k, ])
        pts <- hp$pts
        hullPoints <- hp$hullPoints
        covered <- c(covered, point.in.polygon(truth[1], truth[2], pts[1, hullPoints], pts[2, hullPoints]))


    }

    write.table(matrix(covered, nrow=1), file="coverage-file.csv", append=TRUE,
                row.names=FALSE,
                col.names=FALSE)
    
    detach(samples)


}
