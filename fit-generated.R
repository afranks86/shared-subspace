rm(list=ls())
source("fit-subspace.R")
source("generateData.R")
source("helper.R")

## Subspace Shrinkage: U_k = VO_k
dat1 <- generateData(S=10, R=10, ngroups=10, nvec=rep(20, 10))
save(dat1, file="simDat1.RData")


## Complete Pooling U_k = VO
ngroups <- dat1$ngroups
S <- dat1$S
R <- dat1$R
Ok <- matrix(0, nrow=S, ncol=R)
Ok[1:R, 1:R] <- diag(R)

Ok <- rustiefel(S, R)
lam <- sort(rexp(R, 1/10), decreasing=TRUE)
omega <- lam/(1+lam)
Olist <- OmegaList <- list()
for( k in 1:ngroups) {
    Olist[[k]] <- Ok
    OmegaList[[k]] <- diag(omega)
}

dat2 <- generateData(S=S, R=R, Olist=Olist, ngroups=10, nvec=rep(20, 10))
save(dat2,file="simDat2.RData")

## No Shrinkage: U_k = V_KO_k (e.g. U_k is completely free)

## TODO

####### FIT DATA USING SHARED SUBSPACE ###############

## Initialize sampler

P <- dat1$P
S <- dat1$S
R <- dat1$R
ngroups <- dat1$ngroups
nvec <- dat1$nvec
Slist <- dat1$Slist

Vinit <- matrix(0,nrow=P,ncol=S)
Vinit[1:S,1:S] <- diag(S)

OmegaList <- Ulist <- list()
for(k in 1:ngroups) {
    OmegaList[[k]] <- rep(1/2, R)

    Ok <- matrix(0, nrow=S, ncol=R)
    Ok[1:R, 1:R] <- diag(R)
    
    Ulist[[k]] <- Vinit %*% Ok    
}
s2vec <- rexp(ngroups)

initSS <- list(V=Vinit, Ulist=Ulist, OmegaList=OmegaList, s2vec=s2vec)

## fit dat1
resDat1SS <- fitSubspace(dat1$P, dat1$S, dat1$R, dat1$Slist,
                         dat1$nvec, dat1$ngroups, init=initSS, niters=200)
save(resDat1SS, file="results_dat1_ss.RData")

## fit dat2 
resDat2SS <- fitSubspace(dat2$P, dat2$S, dat2$R, dat2$Slist,
                         dat2$nvec, dat2$ngroups, init=initSS, niters=200)
save(resDat2SS, file="results_dat2_ss.RData")

####### FIT DATA USING COMPLETE POOLING ###############

## fit dat1
Ylist <- dat1$Ylist
Ypooled <- c()
for(i in 1:length(Ylist)){
    Ypooled <- rbind(Ypooled, Ylist[[i]])
}

Slist <- Ulist <- OmegaList <- list()
Slist[[1]] <- Ypooled %*% t(Ypooled)
OmegaList[[1]] <- rep(1/2, S)
V <- initSS$V
Ulist[[1]] <- initSS$Ulist[[1]]
nvec <- sum(nvec)
s2vec <- rexp(1)

initCP <- list(V=Vinit, Ulist=Ulist, OmegaList=OmegaList, s2vec=s2vec)

## fit dat 1
resDat1CP <- fitSubspace(dat1$P, dat1$S, dat1$R,
                         Slist, nvec, ngroups=1, init=initCP, niters=200)
save(resDat1CP, file="results_dat1_cp.RData")

## fit dat2
Ylist <- dat2$Ylist
Ypooled <- c()
for(i in 1:length(Ylist)){
    Ypooled <- rbind(Ypooled, Ylist[[i]])
}

Slist <- list()
Slist[[1]] <- Ypooled %*% t(Ypooled)

resDat2CP <- fitSubspace(dat2$P, dat2$S, dat2$R,
                         Slist, nvec, ngroups=1, init=initCP, niters=200)
save(resDat2CP, file="results_dat2_cp.RData")

############################################################################

load("results_dat1_ss.RData")

for(k in 1:dat1$ngroups) {
    Uk <- resDat1SS$Usamps[, , k, 200]
    Omega <- resDat1SS$omegaSamps[, k, 200]
    s2 <- resDat1SS$s2samps[k, 200]
    SigmaInvHat <- getSigmaInv(P, Uk, Omega, s2)
    print(steinsLoss(dat1$SigmaList[[k]], SigmaInvHat))
}

load("results_dat1_cp.RData")

Uk <- resDat1CP$Usamps[, , 1, 200]
Omega <- resDat1CP$omegaSamps[, 1, 200]
s2 <- resDat1SS$s2samps[1, 200]
SigmaInvHat <- getSigmaInv(P, Uk, Omega, s2)

for(k in 1:dat1$ngroups) {
    print(steinsLoss(dat1$SigmaList[[k]], SigmaInvHat))
}


load("results_dat2_ss.RData")

sl <- 0
for(k in 1:dat2$ngroups) {
    Uk <- resDat2SS$Usamps[, , k, 200]
    Omega <- resDat2SS$omegaSamps[, k, 200]
    s2 <- resDat2SS$s2samps[k, 200]
    SigmaInvHat <- getSigmaInv(P, Uk, Omega, s2)
    sl <- sl+steinsLoss(dat2$SigmaList[[k]], SigmaInvHat)
}
print(sl/dat2$ngroups)

load("results_dat2_cp.RData")

Uk <- resDat2CP$Usamps[, , 1, 200]
Omega <- resDat2CP$omegaSamps[, 1, 200]
s2 <- resDat2SS$s2samps[1, 200]
SigmaInvHat <- getSigmaInv(P, Uk, Omega, s2)

sl <- 0
for(k in 1:dat2$ngroups) {
    sl <- sl + steinsLoss(dat2$SigmaList[[k]], SigmaInvHat)
}
print(sl / dat2$ngroups)
