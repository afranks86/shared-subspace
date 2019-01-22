rm(list=ls())

## Need these for Rscript
library("methods")
library("utils")

library(rstiefel)

source("fit-subspace.R")
source("generateData.R")
source("helper.R")

########################################################
############# Data Generation #########################
########################################################

n <- 500
S <- 2
R <- 2
P <- 20000
ngroups <- 5

niters <- 1000
nwarmup <- niters/2

LambdaList <- Olist <- list()

LambdaList[[1]] <- c(1000, 100)
LambdaList[[2]] <- c(1000, 100)
LambdaList[[3]] <- c(1000, 1000/3)
LambdaList[[4]] <- c(1000, 1000/3)
LambdaList[[5]] <- c(1000, 1000)

Olist[[1]] <- matrix(c(sqrt(2)/2, sqrt(2)/2, -sqrt(2)/2, sqrt(2)/2), nrow=2)
Olist[[2]] <- matrix(c(-sqrt(2)/2, sqrt(2)/2, sqrt(2)/2, sqrt(2)/2), nrow=2)
Olist[[3]] <- matrix(c(-sqrt(2)/2, sqrt(2)/2, sqrt(2)/2, sqrt(2)/2), nrow=2)
Olist[[4]] <- diag(2)
Olist[[5]] <- diag(2)
    
dat <- generateData(P=P, S=S, R=2, ngroups=ngroups,
                    LambdaList=LambdaList, Olist=Olist,
                    nvec=rep(n, ngroups), saveSigmaList=FALSE, saveSList = FALSE)

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

P <- dat$P
S <- getRank(dat$Ypooled)
R <- max(sapply(dat$Ylist, getRank))
if(R > S)
    S <- R

ngroups <- dat$ngroups
nvec <- dat$nvec

## Vinit <- svd(do.call(cbind, lapply(1:ngroups, function(k) svd(t(dat$Ylist[[k]]))$u[, 1:R])))$u[, 1:S]
Vinit <- rustiefel(P, S)
EMFit <- subspaceEM(dat$Ylist, P=P, S=S, R=R, nvec=dat$nvec,
                    Vstart=Vinit, verbose=TRUE, maxIters=100, stiefelAlgo=2)
Vinit <- EMFit$V

tr((t(dat$V) %*%  EMFit$V) %*% (t(EMFit$V) %*% dat$V)) / S

OmegaList <- Ulist <- list()
for(k in 1:ngroups) {
    OmegaList[[k]] <- rep(1/2, R)
    
    Ok <- matrix(0, nrow=S, ncol=R)
    Ok[1:R, 1:R] <- diag(R)
    
    Ulist[[k]] <- Vinit %*% Ok    
}
s2vec <- rexp(ngroups)

initSS <- list(V=Vinit, Ulist=Ulist, OmegaList=OmegaList, s2vec=s2vec)
samples <- fitSubspace(dat$P, S, R, Q=S-R, dat$Ylist,
                   dat$nvec, dat$ngroups, init=initSS,
                   niters=niters)

save(samples, dat, file="bayes-simulations-1-19.RData")

pdf("paper/Figs/posteriorRegions-chull.pdf", width=5, height=5, font="Times")

sv <- svd(t(EMFit$V) %*% dat$V)
Rot <- sv$v %*% t(sv$u)

attach(samples)
rotOsamps <- array(dim=dim(Osamps))
for(i in 1:1000) {
    for(k in 1:ngroups) {
        rotOsamps[1:2, 1:2, k, i] <- Rot %*% Osamps[1:2, 1:2, k, i]
    }
}

posteriorPlot(rotOsamps, omegaSamps[1:2, , ],
              s2samps, nsamps=500, groups=1:ngroups, ymax=log(200),
              probRegion=0.95, hline=NULL, logRatio=TRUE, plotPoints=FALSE)

for(k in 1:ngroups) {
    ## eigRot <- eigen(Rot %*%
    ##                 (t(dat$Ylist[[k]] %*% EMFit$V) %*% (dat$Ylist[[k]] %*% EMFit$V) %*% t(Rot)))

    evals <- eigRot$values
    evecs <- eigRot$vectors

    ## points(atan(evecs[2]/evecs[1]), evals[1]/evals[2],
    ##       pch=19, col="white", cex=2)
    ## points(atan(evecs[2]/evecs[1]), evals[1]/evals[2],
    ##        pch="+:", col=k, cex=2)
    ## points(atan(evecs[2]/evecs[1]), evals[1]/evals[2],
    ##        pch=1, col="black", cex=2)

    points(atan(evecs[2]/evecs[1]), log(evals[1]/evals[2]),
           pch=19, col=k, cex=1.5)
}



detach(samples)
dev.off()

evalRatios <- sapply(1:ngroups, function(k) {
    numer <- sum(eigen(t(samples$V[, 1:R]) %*% dat$Slist[[k]] %*% samples$V[, 1:R])$values)/dat$nvec[k]
    denom <- sum(svd(dat$Ylist[[k]])$d[1:R]^2)/dat$nvec[k]
    correction <-  (denom - 1/mean(samples$s2samps[k, ]) * R * P / dat$nvec[k])  / denom
    (numer/denom) / correction
})

pdf(sprintf("paper/Figs/simRatio-%s.pdf", format(Sys.Date(), "%m-%d")))
par(mar=c(7.1, 4.1, 4.1, 2.1))
barplot(evalRatios, ylim=c(0, 1), xlab="", main="", space=0.2, cex.axis=2, col="#646464", col.axis="#646464", col.main="#646464", border=NA,
        names.arg=c("Black", "Red", "Green", "Blue", "Cyan"), cex.names=1.5)
dev.off()
