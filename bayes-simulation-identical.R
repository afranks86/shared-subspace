rm(list=ls())

## Need these for Rscript
library(microbenchmark)
library(rstiefel)

source("fit-subspace.R")
source("generateData.R")
source("subspace-functions.R")

########################################################
############# Data Generation #########################
########################################################

n <- 200
S <- 10
R <- 2
P <- 20000
ngroups <- 5

niters <- 1000
nwarmup <- niters/2

LambdaList <- Olist <- list()
Opooled <- matrix(0, nrow=S, ncol=S)
if(R  < S) { 
    Opooled[1:S, 1:S] <- rustiefel(S, S)
    LambdaPooled <- sort(rexp(S, 1/500), decreasing=TRUE)
    for(k in 1:ngroups) {
        Olist[[k]] <- Opooled
        LambdaList[[k]] <- LambdaPooled
    }
}

dat <- generateData(P=P, S=S, R=R, ngroups=ngroups,
                    V = rustiefel(P, S),
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

ngroups <- dat$ngroups
nvec <- dat$nvec

Vinit <- svd(do.call(cbind, lapply(1:ngroups, function(k) svd(t(dat$Ylist[[k]]))$u[, 1:S])))$u[, 1:S]
Vinit <- Vinit %*% rustiefel(S, S)

microbenchmark(EMFit <- subspaceEM(dat$Ylist, P=P, S=S, R=R, nvec=dat$nvec,
                    Vstart=Vinit, verbose=TRUE, EM_iters=10, M_iters=300, stiefelAlgo=2), times=1)

tr((t(dat$V) %*%  EMFit$V) %*% (t(EMFit$V) %*% dat$V)) / S
tr((t(dat$V) %*%  Vinit) %*% (t(Vinit) %*% dat$V)) / S
tr((t(EMFit$V) %*%  Vinit) %*% (t(Vinit) %*% EMFit$V)) / S

tr((t(dat$V[, 1:2]) %*%  EMFit$V[, 1:R]) %*% (t(EMFit$V[, 1:R]) %*% dat$V[, 1:2])) / R



OmegaList <- Ulist <- list()
for(k in 1:ngroups) {
    OmegaList[[k]] <- rep(1/2, S)
    
    Ok <- matrix(0, nrow=S, ncol=S)
    Ok <- diag(S) ## [1:R, 1:R] <- diag(R)
    
    Ulist[[k]] <- EMFit$V %*% Ok    
}
s2vec <- rexp(ngroups)

initSS <- list(V=EMFit$V, Ulist=Ulist, OmegaList=OmegaList, s2vec=s2vec)
samples <- fitSubspace(dat$P, S, R, Q=S-R, dat$Ylist,
                   dat$nvec, dat$ngroups, init=initSS,
                   niters=niters)


save(samples, dat, file="bayes-simulations-identical-2-19.RData")

pdf("paper/Figs/posteriorRegions-chull-identical.pdf", width=5, height=5, font="Times")

attach(samples)

indices <- c(6:7)

sv <- svd(t(EMFit$V[, indices]) %*% dat$V[, indices])
Rot <- sv$v %*% t(sv$u)


rotOsamps <- array(dim=c(length(indices), length(indices), ngroups, niters))
for(i in 1:1000) {
    for(k in 1:ngroups) {
        rotOsamps[1:length(indices), 1:length(indices), k, i] <- Rot %*% Osamps[indices, indices, k, i]
    }
}

posteriorPlot(rotOsamps, omegaSamps[indices, , ],
              s2samps, nsamps=500, groups=1:ngroups, ymax=log(200),
              probRegion=0.95, hline=NULL, logRatio=TRUE, plotPoints=FALSE)

for(k in 1:ngroups) {
    ## eigRot <- eigen(Rot %*%
    ##                 (t(dat$Ylist[[k]] %*% EMFit$V) %*% (dat$Ylist[[k]] %*% EMFit$V) %*% t(Rot)))

    evals <- dat$LambdaList[[k]]
    evecs <- dat$Olist[[k]][, 1]

    points(atan(evecs[2]/evecs[1]), log2(evals[1]/evals[2]),
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
