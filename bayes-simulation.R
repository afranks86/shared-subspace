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

LambdaList[[1]] <- c(1000, 100)
LambdaList[[2]] <- c(1000, 100)
LambdaList[[3]] <- c(1000, 1000/3)
LambdaList[[4]] <- c(1000, 1000/3)
LambdaList[[5]] <- c(1000, 1000)

Olist[[1]] <- matrix(c(sqrt(2)/2, sqrt(2)/2, -sqrt(2)/2, sqrt(2)/2), nrow=2)
Olist[[2]] <- matrix(c(-sqrt(2)/2, sqrt(2)/2, sqrt(2)/2, sqrt(2)/2), nrow=2)
Olist[[3]] <- diag(2)## matrix(c(-sqrt(2)/2, sqrt(2)/2, sqrt(2)/2, sqrt(2)/2), nrow=2)
Olist[[4]] <- diag(2)
Olist[[5]] <- diag(2)

Opooled <- matrix(0, nrow=S, ncol=S)
Opooled[1:2, 1:2] <- diag(2)
if(R  < S) { 
    Opooled[(R+1):S, (R+1):S] <- rustiefel(S-R, S-R)
    LambdaPooled <- rexp(S - R, 1/50)
    for(k in 1:ngroups) {
        Ok <- matrix(0, nrow=S, ncol=S)
        Ok[1:R, 1:R] <- Olist[[k]]
        Ok[(R+1):S, (R+1):S] <- Opooled[(R+1):S, (R+1):S]
        Olist[[k]] <- Ok
        LambdaList[[k]] <- c(LambdaList[[k]], LambdaPooled)
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

Vinit <- svd(do.call(cbind, lapply(1:ngroups, function(k) svd(t(dat$Ylist[[k]]))$u[, 1:R])))$u[, 1:S]

## ulst <- lapply(1:ngroups, function(k) svd(t(dat$Ylist[[k]]))$u[, 1:S])

## sapply(1:ngroups, function(n) max(abs(sapply(1:10, function(s) t(ulst[[1]][, 10]) %*% ulst[[n]][, s]))))

Vinit <- Vinit %*% rustiefel(S, S)


microbenchmark(EMFit <- subspaceEM(dat$Ylist, P=P, S=S, R=R, nvec=dat$nvec,
                                   Vstart=Vinit, verbose=TRUE, EM_iters=100, M_iters=100, stiefelAlgo=2), times=1)

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

## start from truth
OmegaList <- Ulist <- list()
for(k in 1:ngroups) {
    OmegaList[[k]] <- rep(1/2, S)
    
    Ok <- matrix(0, nrow=S, ncol=S)
    Ok <- diag(S) ## [1:R, 1:R] <- diag(R)
    
    Ulist[[k]] <- dat$V %*% Ok    
}


s2vec <- rexp(ngroups)

initSS <- list(V=EMFit$V, Ulist=Ulist, OmegaList=OmegaList, s2vec=s2vec)
initSS <- list(V=dat$V, Ulist=Ulist, OmegaList=OmegaList, s2vec=s2vec)
samples <- fitSubspace(dat$P, S, R, Q=S-R, dat$Ylist,
                   dat$nvec, dat$ngroups, init=initSS,
                   niters=niters)

save(samples, dat, file="bayes-simulations-2-4.RData")

pdf("paper/Figs/posteriorRegions-chull.pdf", width=5, height=5, font="Times")

attach(samples)

indices <- c(1:2)

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


    if(k == 3) {
        points(atan(evecs[2]/evecs[1]), log2(evals[1]/evals[2])+.05,
               pch=24, col=k, bg=k, cex=1.1)
    } else if (k == 4) {
        points(atan(evecs[2]/evecs[1]), log2(evals[1]/evals[2])-.05,
               pch=25, col=k, bg=k, cex=1.1)
    }
    else {
        points(atan(evecs[2]/evecs[1]), log2(evals[1]/evals[2]),
               pch=19, col=k, cex=1.5)
    }
    
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
