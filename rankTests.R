#!/usr/bin/env Rscript

#SBATCH -n 4  #Number of cores
#SBATCH -N 1  #Number of cores 
#SBATCH -t 15000  #Runtime in minutes
#SBATCH -J ss3
#SBATCH -o outfiles/ssg_%a.out #Standard output
#SBATCH -e outfiles/ssg_%a.err #Standard error
#SBATCH -p airoldi,stats #Partition to submit to 
#SBATCH --mem=10000  #Memory per node in MB (see also --mem-per-cpu)
#SBATCH -a 1-100

rm(list=ls())
idx <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## Need these for Rscript
library("methods")
library("utils")

library(rstiefel)

source("fit-subspace.R")
source("generateData.R")
source("helper.R")

datType <- ((idx-1) %% 2) + 1


########################################################
############# Data Generation #########################
########################################################
n <- 50
P <- 200
ngroups <- 10

if(datType==1) {
  S <- 2
  R <- 2
} else if(datType==2) {
  S <- 10
  R <- 2
}

evals <- c(250, 25)
niters <- 1000
nwarmup <- niters/2

## Subspace Shrinkage: Psi_k = VO_kLam_kO_k^TV^T
LambdaList <- list()
lam <- evals
for( k in 1:ngroups) {
  LambdaList[[k]] <- lam
}
    
dat <- generateData(P=P, S=S, R=R, ngroups=ngroups,
                    LambdaList=LambdaList, nvec=rep(n, ngroups))

dat$genType <- "Shared subspace"

Ypooled <- c()
for(i in 1:length(dat$Ylist)){
    Ypooled <- rbind(Ypooled, dat$Ylist[[i]])
}
dat$Ypooled <- Ypooled

print("Saving data...")
save(dat, file=sprintf("/n/airoldifs2/lab/afranks/simDat%i.RData", datType))


########################################################
############# Data Fit #########################
########################################################

resList <- list()
lossVec <- numeric(3)

####### FIT DATA USING SHARED SUBSPACE ###############

## Initialize sampler

P <- dat$P
S <- 2
R <- 2

ngroups <- dat$ngroups
nvec <- dat$nvec
Slist <- dat$Slist

Vinit <- svd(do.call(cbind, lapply(1:ngroups, function(k) svd(t(dat$Ylist[[k]]))$u[, 1:R])))$u[, 1:S]
EMFit <- subspaceEM(dat$Slist, P=P, S=S, R=R, nvec=dat$nvec,
                    rho1=0.1, rho2=0.9, 
                    Vstart=Vinit, verbose=FALSE, maxIters=100)
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
res <- fitSubspace(dat$P, S, R, dat$Slist,
                   dat$nvec, dat$ngroups, init=initSS,
                   niters=niters, sigmaTruthList=dat$SigmaList,
                   draw=c(V=FALSE))

res$predictions <- makePrediction(res$Usamps, res$omegaSamps,
                                  res$s2samps, dat$SigmaList, n=100)
fullLoss <- ssLoss <- evalRatio <- c()
for(k in 1:dat$ngroups) {
  
  UkSamps <- res$Usamps[, , k, (nwarmup+1):niters]
  OmegaSamps <- res$omegaSamps[, k, (nwarmup+1):niters]
  s2list <- res$s2samps[k, (nwarmup+1):niters]

  ## How well do we recover the full covariance matrix
  SigmaInvK <- getSigmaInv(dat$P, dat$Ulist[[k]],
                           dat$OmegaList[[k]], dat$s2vec[k])
  SigmaHatInvPM <- getPostMeanSigmaInv(P, UkSamps, OmegaSamps,
                                       s2list, niters-nwarmup)
  fullLoss <- c(fullLoss, steinsLoss(solve(SigmaHatInvPM), SigmaInvK))

  ## How well do we recover the projected covariance matrix
  ProjSigmaInvK <- solve(t(Vinit) %*% dat$SigmaList[[k]] %*% Vinit)
  ProjSigmaHatInvPM <- getPostMeanSigmaProjInv(S, Vinit, UkSamps, OmegaSamps,
                                               s2list, niters-nwarmup)
  ssLoss <- c(ssLoss, steinsLoss(solve(ProjSigmaHatInvPM), ProjSigmaInvK))
  
  ## How much of the variation Sigma does projected Sigma account for
  evalRatio <- c(evalRatio,
                 sum(rowMeans(OmegaSamps / (1-OmegaSamps))) /
                 sum((eigen(dat$Slist[[k]])$values/nvec[k])[1:2]))
}

lossVec[1] <- sum(fullLoss / ngroups)
res$fullLoss <- fullLoss
res$ssLoss <- ssLoss
res$evalRatio <- evalRatio
resList[[1]] <- res

save(resList, dat, lossVec, datType, file=sprintf("/n/airoldifs2/lab/afranks/shared_subspace/rankTest-%i-%i.RData", datType, idx))


