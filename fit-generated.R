#!/usr/bin/env Rscript

#SBATCH -n 4  #Number of cores
#SBATCH -N 1  #Number of cores 
#SBATCH -t 15000  #Runtime in minutes
#SBATCH -J test_heatshock
#SBATCH -o outfiles/new_heatshock_%a.out #Standard output
#SBATCH -e outfiles/new_heatshock_%a.err #Standard error
#SBATCH -p airoldi,stats #Partition to submit to 
#SBATCH --mem=10000  #Memory per node in MB (see also --mem-per-cpu)
#SBATCH -a 1-6

rm(list=ls())
source("fit-subspace.R")
source("generateData.R")
source("helper.R")

idx <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

datType <- (idx-1) %% 3 + 1

########################################################
############# Data Generation #########################
########################################################
evals <- c(250, 125, 50, 30, 30 ,30, 20, 20, 20, 20, 20)
## For all models Sigma_k = s2(Psi_k + diag(P))

## Subspace Shrinkage: Psi_k = VO_kLam_kO_k^TV^T
dat1 <- generateData(S=10, R=10, ngroups=10, nvec=rep(20, 10))
save(dat1, file="simDat1.RData")

## Complete Pooling: Psi_k = Psi
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

## No pooling: Psi_k = Psi_k
P <- dat1$P
S <- P
R <- dat1$R
ngroups <- dat1$ngroups
Olist <- OmegaList <- list()
for( k in 1:ngroups) {
    Olist[[k]] <- rustiefel(200, 10)
    lamk <- diag(rexp(R, 1/10))
    OmegaList[[k]] <- lamk/(1+lamk)
}
dat3 <- generateData(P=P, S=S, R=R, V=diag(S), Olist=Olist,
                     ngroups=ngroups, nvec=rep(20, ngroups))

########################################################
############# Data Fit #########################
########################################################


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
## initSS <- list(V=dat1$V, Ulist=dat1$Ulist, OmegaList=dat1$OmegaList, s2vec=dat1$s2vec)
## fit dat1
resDat1SS <- fitSubspace(dat1$P, dat1$S, dat1$R, dat1$Slist,
                         dat1$nvec, dat1$ngroups, init=initSS, niters=200,
                         sigmaTruthList=dat1$SigmaList)
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
Slist[[1]] <- t(Ypooled) %*% Ypooled
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
Slist[[1]] <- t(Ypooled) %*% Ypooled

resDat2CP <- fitSubspace(dat2$P, dat2$S, dat2$R,
                         Slist, nvec, ngroups=1, init=initCP, niters=200)
save(resDat2CP, file="results_dat2_cp.RData")


#################  Fit Data Using No Pooling #############################

## fit dat1
nsamps <- niters <- 200
resDat1NP <- list(S=dat1$S, R=dat1$R, ngroups=dat1$ngroups)

Usamps <- array(dim=c(P, R, ngroups, nsamps))
omegaSamps <- array(dim=c(R, ngroups, ncol=nsamps))
s2samps <- matrix(nrow=ngroups, ncol=nsamps)
    
for(k in 1:dat1$ngroups) {
    
    resDat1NPk <- fitBayesianSpike(dat1$P, dat1$S, dat1$R, dat1$Slist[[k]], dat1$nvec[k], niters=niters)
    
    Usamps[, , k, ] <- resDat1NPk$Usamps
    omegaSamps[, k, ] <- resDat1NPk$omegaSamps
    s2samps[k, ] <- resDat1NPk$s2samps
}
resDat1NP$Usamps <- Usamps
resDat1NP$omegaSamps <- omegaSamps
resDat1NP$s2samps <- s2samps
save(resDat1NP, file="results_dat1_np.RData")

## fit dat2
nsamps <- niters <- 200
resDat2NP <- list(S=dat2$S, R=dat2$R, ngroups=dat2$ngroups)

Usamps <- array(dim=c(P, R, ngroups, nsamps))
omegaSamps <- array(dim=c(R, ngroups, ncol=nsamps))
s2samps <- matrix(nrow=ngroups, ncol=nsamps)
    
for(k in 1:dat2$ngroups) {

    resDat2NPk <- fitBayesianSpike(dat2$P, dat2$S, dat2$R, dat2$Slist[[k]], dat2$nvec[k], niters=niters)

    Usamps[, , k, ] <- resDat2NPk$Usamps
    omegaSamps[, k, ] <- resDat2NPk$omegaSamps
    s2samps[k, ] <- resDat2NPk$s2samps
}
resDat2NP$Usamps <- Usamps
resDat2NP$omegaSamps <- omegaSamps
resDat2NP$s2samps <- s2samps
save(resDat2NP, file="results_dat2_np.RData")


steinsLoss(dat1$SigmaList[[1]], getSigmaInv(dat1$P, resDat1NP$Usamps[, , 200], resDat1NP$omegaSamps[, 200], resDat1NP$s2samps[200]))

#################### Summarize Results ################################

load("results_dat1_ss.RData")

sl <- 0
for(k in 1:dat1$ngroups) {
    UkList <- resDat1SS$Usamps[, , k, 101:200]
    OmegaList <- resDat1SS$omegaSamps[, k, 101:200]
    s2list <- resDat1SS$s2samps[k, 101:200]
    SigmaInvPM <- getPostMeanSigmaInv(P, UkList, OmegaList, s2list, 100)
    sl <- sl + steinsLoss(dat1$SigmaList[[k]], SigmaInvPM)
}
print(sl/dat1$ngroups)

load("results_dat1_np.RData")

sl <- 0
for(k in 1:dat1$ngroups) {
    UkList <- resDat1NP$Usamps[, , k, 101:200]
    OmegaList <- resDat1NP$omegaSamps[, k, 101:200]
    s2list <- resDat1NP$s2samps[k, 101:200]
    SigmaInvPM <- getPostMeanSigmaInv(P, UkList, OmegaList, s2list, 100)
    sl <- sl+steinsLoss(dat1$SigmaList[[k]], SigmaInvPM)
}
print(sl/dat1$ngroups)

load("results_dat1_cp.RData")

UkList <- resDat1CP$Usamps[, , 1, 101:200]
OmegaList <- resDat1CP$omegaSamps[, 1, 101:200]
s2list <- resDat1CP$s2samps[1, 101:200]
SigmaInvPM <- getPostMeanSigmaInv(P, UkList, OmegaList, s2List, 100)

sl <- 0
for(k in 1:dat1$ngroups) {
    sl <- sl + steinsLoss(dat1$SigmaList[[k]], SigmaInvPM)
}
print(sl/dat1$ngroups) 

###############################################################

load("results_dat2_ss.RData")

sl <- 0
for(k in 1:dat2$ngroups) {
    UkList <- resDat2SS$Usamps[, , k, 101:200]
    OmegaList <- resDat2SS$omegaSamps[, k, 101:200]
    s2list <- resDat2SS$s2samps[k, 101:200]
    SigmaInvPM <- getPostMeanSigmaInv(P, UkList, OmegaList, s2list, 100)
    sl <- sl+steinsLoss(dat2$SigmaList[[k]], SigmaInvPM)
}
print(sl/dat2$ngroups)


load("results_dat2_np.RData")

sl <- 0
for(k in 1:dat2$ngroups) {
    UkList <- resDat2NP$Usamps[, , k, 101:200]
    OmegaList <- resDat2NP$omegaSamps[, k, 101:200]
    s2list <- resDat2NP$s2samps[k, 101:200]
    SigmaInvPM <- getPostMeanSigmaInv(P, UkList, OmegaList, s2list, 100)
    sl <- sl+steinsLoss(dat2$SigmaList[[k]], SigmaInvPM)
}
print(sl/dat2$ngroups)

load("results_dat2_cp.RData")

UkList <- resDat2CP$Usamps[, , 1, 101:200]
OmegaList <- resDat2CP$omegaSamps[, 1, 101:200]
s2list <- resDat2CP$s2samps[1, 101:200]
SigmaInvPM <- getPostMeanSigmaInv(P, UkList, OmegaList, s2list, 100)

sl <- 0
for(k in 1:dat2$ngroups) {
    sl <- sl + steinsLoss(dat2$SigmaList[[k]], SigmaInvPM)
}
print(sl / dat2$ngroups)
