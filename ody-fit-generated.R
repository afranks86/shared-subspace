#!/usr/bin/env Rscript

#SBATCH -n 4  #Number of cores
#SBATCH -N 1  #Number of cores 
#SBATCH -t 15000  #Runtime in minutes
#SBATCH -J test_ss
#SBATCH -o outfiles/ssg_%a.out #Standard output
#SBATCH -e outfiles/ssg_%a.err #Standard error
#SBATCH -p airoldi,stats #Partition to submit to 
#SBATCH --mem=10000  #Memory per node in MB (see also --mem-per-cpu)
#SBATCH -a 1-300

rm(list=ls())
idx <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## Need these for Rscript
library("methods")
library("utils")

library(rstiefel)

source("fit-subspace.R")
source("generateData.R")
source("helper.R")

datType <- ((idx-1) %% 3) + 1


########################################################
############# Data Generation #########################
########################################################

n <- 20
S <- 10
R <- 10
P <- 200
ngroups <- 10
evals <- c(250, 125, 50, 30, 30 ,30, 20, 20, 20, 20)

niters <- 1000
nwarmup <- niters/2

## For all models Sigma_k = s2(Psi_k + diag(P))
if(datType==1) {

  ## Subspace Shrinkage: Psi_k = VO_kLam_kO_k^TV^T
  LambdaList <- list()
  lam <- evals
  for( k in 1:ngroups) {
    LambdaList[[k]] <- lam
  }
  
  dat <- generateData(P=P, S=S, R=R, ngroups=ngroups,
                      LambdaList=LambdaList, nvec=rep(n, ngroups))
  
} else if(datType==2) {

  ## Complete Pooling: Psi_k = Psi
  Ok <- rustiefel(S, R)
  lam <- evals

  Olist <- LambdaList <- list()
  for( k in 1:ngroups) {
    Olist[[k]] <- Ok
    LambdaList[[k]] <- lam
  }

  dat <- generateData(P=P, S=S, R=R, Olist=Olist, ngroups=ngroups,
                      LambdaList=LambdaList, nvec=rep(n, ngroups))
  
} else if (datType==3) {

  ## No pooling: Psi_k = Psi_k
  Olist <- LambdaList <- list()
  for( k in 1:ngroups) {
    Olist[[k]] <- rustiefel(P, R)
    lamk <- evals 
    LambdaList[[k]] <- lamk
  }
  dat <- generateData(P=P, S=S, R=R, V=diag(P), Olist=Olist,
                      ngroups=ngroups, nvec=rep(n, ngroups),
                      LambdaList=LambdaList)
  
}

print("Saving data...")
save(dat, file=sprintf("/n/airoldifs2/lab/afranks/simDat%i.RData", datType))


########################################################
############# Data Fit #########################
########################################################


resList <- list()
lossVec <- numeric(3)

for(fitType in 1:3) {

####### FIT DATA USING SHARED SUBSPACE ###############

  if(fitType == 1) {
    ## Initialize sampler

    P <- dat$P
    S <- dat$S
    R <- dat$R
    ngroups <- dat$ngroups
    nvec <- dat$nvec
    Slist <- dat$Slist

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
    res <- fitSubspace(dat$P, dat$S, dat$R, dat$Slist,
                       dat$nvec, dat$ngroups, init=initSS,
                       niters=niters, sigmaTruthList=dat$SigmaList)
    loss <- 0
    for(k in 1:dat$ngroups) {
      
      UkSamps <- res$Usamps[, , k, (nwarmup+1):niters]
      OmegaSamps <- res$omegaSamps[, k, (nwarmup+1):niters]
      s2list <- res$s2samps[k, (nwarmup+1):niters]
      SigmaInvPM <- getPostMeanSigmaInv(P, UkSamps, OmegaSamps,
                                        s2list, (niters-nwarmup))
      loss <- loss + steinsLoss(dat$SigmaList[[k]], SigmaInvPM)
    }
    loss <- loss/dat$ngroups

    lossVec[1] <- loss
    resList[[1]] <- res
    
  } else if (fitType == 2) {

####### FIT DATA USING COMPLETE POOLING ###############

    Ylist <- dat$Ylist
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

    ## fit dat 
    res <- fitSubspace(dat$P, dat$S, dat$R, Slist, nvec,
                       ngroups=1, init=initCP, niters=niters,
                       sigmaTruthList=dat$SigmaList)

    loss <- 0
    UkSamps <- res$Usamps[, , 1, (nwarmup+1):niters]
    OmegaSamps <- res$omegaSamps[, 1, (nwarmup+1):niters]
    s2list <- res$s2samps[1, (nwarmup+1):niters]
    SigmaInvPM <- getPostMeanSigmaInv(P, UkSamps, OmegaSamps,
                                      s2list, (niters-nwarmup))
    loss <- steinsLoss(dat$SigmaList[[k]], SigmaInvPM)

    lossVec[2] <- loss
    resList[[2]] <- res
    
  } else if (fitType == 3) {

#################  Fit Data Using No Pooling ########################

    nsamps <- niters
    res <- list(S=dat$S, R=dat$R, ngroups=dat$ngroups)

    Usamps <- array(dim=c(P, R, ngroups, nsamps))
    omegaSamps <- array(dim=c(R, ngroups, ncol=nsamps))
    s2samps <- matrix(nrow=ngroups, ncol=nsamps)
    
    for(k in 1:dat$ngroups) {
      
      resk <- fitBayesianSpike(dat$P, dat$S, dat$R, dat$Slist[[k]],
                               dat$nvec[k], niters=niters)
      
      Usamps[, , k, ] <- resk$Usamps
      omegaSamps[, k, ] <- resk$omegaSamps
      s2samps[k, ] <- resk$s2samps

    }

    res$Usamps <- Usamps
    res$omegaSamps <- omegaSamps
    res$s2samps <- s2samps

    loss <- 0
    for(k in 1:dat$ngroups) {
      
      UkSamps <- res$Usamps[, , k, (nwarmup+1):niters]
      OmegaSamps <- res$omegaSamps[, k, (nwarmup+1):niters]
      s2list <- res$s2samps[k, (nwarmup+1):niters]
      SigmaInvPM <- getPostMeanSigmaInv(P, UkSamps, OmegaSamps,
                                        s2list, niters-nwarmup)
      loss <- loss + steinsLoss(dat$SigmaList[[k]], SigmaInvPM)
    }
    loss <- loss/dat$ngroups

    lossVec[3] <- loss
    resList[[3]] <- res

  }

  print(sprintf("Finished fitting %i", fitType))
  save(resList, dat, lossVec, datType,
       file=sprintf("/n/airoldifs2/lab/afranks/shared_subspace/results-%i-%i.RData", datType, idx))
}


