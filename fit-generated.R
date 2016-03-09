rm(list=ls())
## Need these for Rscript
library("methods")
library("utils")

library(rstiefel)

source("fit-subspace.R")
source("generateData.R")
source("helper.R")

datType <- 1

########################################################
############# Data Generation #########################
########################################################

n <- 20
S <- 10
R <- 2
P <- 20
ngroups <- 3

evals <- c(250, 20)

niters <- 100
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

  dat <- generateData(P=P, S=P, R=R, V=diag(P), Olist=Olist,
                      ngroups=ngroups, nvec=rep(n, ngroups),
                      LambdaList=LambdaList)
  
}

Ypooled <- c()
for(i in 1:length(dat$Ylist)){
  Ypooled <- rbind(Ypooled, dat$Ylist[[i]])
}
dat$Ypooled <- Ypooled

########################################################
############# Data Fit #########################
########################################################

resList <- list()
lossVec <- numeric(3)

for(fitType in 1:3) {

####### FIT DATA USING SHARED SUBSPACE ###############

  if(fitType == 1) {
    ## Initialize sampler

    fitP <- dat$P
    fitS <- fitR <- getRank(dat$Ypooled)
    ## S <- dat$S
    ## R <- dat$R
    ngroups <- dat$ngroups
    nvec <- dat$nvec
    Slist <- dat$Slist

    ##Vinit <- matrix(0,nrow=P,ncol=S)
    ##Vinit[(P-S+1):P, 1:S] <- diag(S)
    Vinit <- rustiefel(P, fitS)
#      Vinit <- dat$V %*% dat$Olist[[1]]
#      Vinit <- dat$V
    OmegaList <- Ulist <- list()
    for(k in 1:ngroups) {
      OmegaList[[k]] <- rep(0.996, fitR)

      Ok <- matrix(0, nrow=fitS, ncol=fitR)
      Ok[1:fitR, 1:fitR] <- diag(fitR)

      ##Ok <- rustiefel(S, R)
      #Ok <- dat$Olist[[k]]
      
      Ulist[[k]] <- Vinit %*% Ok    
    }
    s2vec <- rexp(ngroups)

    initSS <- list(V=Vinit, Ulist=Ulist, OmegaList=OmegaList, s2vec=s2vec)
    res <- fitSubspace(dat$P, fitS, fitR, dat$Slist,
                       dat$nvec, dat$ngroups, init=initSS,
                       niters=niters, sigmaTruthList=dat$SigmaList,
                       draw=c(V=FALSE, s2=TRUE, omega=TRUE, O=TRUE), Vmode="gibbs")
    
    res$predictions <- makePrediction(res$Usamps, res$omegaSamps,
                                  res$s2samps, dat$SigmaList, n=100)
    loss <- 0
    for(k in 1:dat$ngroups) {
      
      UkSamps <- res$Usamps[, , k, (nwarmup+1):niters]
      OmegaSamps <- res$omegaSamps[, k, (nwarmup+1):niters]
      s2list <- res$s2samps[k, (nwarmup+1):niters]

      SigmaInvK <- getSigmaInv(dat$P, dat$Ulist[[k]],
                               dat$OmegaList[[k]], dat$s2vec[k])
      SigmaHatInvPM <- getPostMeanSigmaInv(P, UkSamps, OmegaSamps,
                                        s2list, niters-nwarmup)
      loss <- loss + steinsLoss(solve(SigmaHatInvPM), SigmaInvK)
      
    }
    loss <- loss/ngroups
    print(loss)

    lossVec[1] <- loss
    resList[[1]] <- res
    
  } else if (fitType == 2) {

####### FIT DATA USING COMPLETE POOLING ###############

    P <- dat$P
    S <- R <- getRank(dat$Ypooled)
    ## S <- dat$S
    ## R <- dat$R

    Slist <- Ulist <- OmegaList <- list()
    Slist[[1]] <- t(Ypooled) %*% Ypooled
    OmegaList[[1]] <- rep(1/2, S)
    
    Vinit <- matrix(0,nrow=P,ncol=S)
    Vinit[1:S, 1:S] <- diag(S)
    Oinit <- rustiefel(S, R)
    Uinit <- Vinit %*% Oinit
    
    Ulist[[1]] <- Uinit
    nvec <- sum(dat$nvec)
    s2vec <- rexp(1)

    initCP <- list(V=Vinit, Ulist=Ulist, OmegaList=OmegaList, s2vec=s2vec)
    resk <- fitBayesianSpike(dat$P, S, R, Slist[[1]],
                             nvec, niters=niters,
                             SigmaTruth=dat$SigmaList[[1]])

    Usamps <- array(dim=c(dat$P, R, dat$ngroups, niters))
    omegaSamps <- array(dim=c(R, dat$ngroups, ncol=niters))
    s2samps <- matrix(nrow=dat$ngroups, ncol=niters)

    for(k in 1:dat$ngroups) {
      Usamps[, , k, ] <- resk$Usamps
      omegaSamps[, k, ] <- resk$omegaSamps
      s2samps[k, ] <- resk$s2samps
    }

    res <- list(S=S, R=R, ngroups=dat$ngroups)
    res$Usamps <- Usamps
    res$omegaSamps <- omegaSamps
    res$s2samps <- s2samps

    res$predictions <- makePrediction(res$Usamps, res$omegaSamps,
                                      res$s2samps, dat$SigmaList, n=100)

    
    SigmaHatInvPM <- getPostMeanSigmaInv(P, resk$Usamps,
                                         resk$omegaSamps,
                                         resk$s2samps, (niters-nwarmup))
    
    loss <- 0
    for(k in 1:dat$ngroups) {
      SigmaInvK <- getSigmaInv(dat$P, dat$Ulist[[k]],
                               dat$OmegaList[[k]], dat$s2vec[k])
      loss <- loss + steinsLoss(solve(SigmaHatInvPM), SigmaInvK)
    }
      loss <- loss/dat$ngroups
      print(loss)
      
    lossVec[2] <- loss
    resList[[2]] <- res
    
  } else if (fitType == 3) {

#################  Fit Data Using No Pooling ########################

    P <- dat$P

    ## All groups assumed to have same rank, take largest
    SRvec <- c()
    for(k  in 1:dat$ngroups) {
      SRvec <- c(SRvec, getRank(dat$Ylist[[k]]))
    }
    R <- max(SRvec)
    
    ## S <- dat$S
    ## R <- dat$R
    
    nsamps <- niters
    res <- list(S=S, R=R, ngroups=dat$ngroups)

    Usamps <- array(dim=c(P, R, ngroups, nsamps))
    omegaSamps <- array(dim=c(R, ngroups, ncol=nsamps))
    s2samps <- matrix(nrow=ngroups, ncol=nsamps)
    
    for(k in 1:dat$ngroups) {
      
      resk <- fitBayesianSpike(dat$P, R, dat$Slist[[k]],
                               dat$nvec[k], niters=niters,
                               SigmaTruth=dat$SigmaTruthList[[k]],
                               ngroups=ngroups)
      
      Usamps[, , k, ] <- resk$Usamps
      omegaSamps[, k, ] <- resk$omegaSamps
      s2samps[k, ] <- resk$s2samps

    }

    res$Usamps <- Usamps
    res$omegaSamps <- omegaSamps
    res$s2samps <- s2samps

    res$predictions <- makePrediction(res$Usamps, res$omegaSamps,
                                      res$s2samps, dat$SigmaList, n=100)

    loss <- 0
    for(k in 1:dat$ngroups) {

      UkSamps <- res$Usamps[, , k, (nwarmup+1):niters]
      OmegaSamps <- res$omegaSamps[, k, (nwarmup+1):niters]
      s2list <- res$s2samps[k, (nwarmup+1):niters]
      
      SigmaInvK <- getSigmaInv(dat$P, dat$Ulist[[k]],
                               dat$OmegaList[[k]], dat$s2vec[k])
      SigmaHatInvPM <- getPostMeanSigmaInv(P, UkSamps, OmegaSamps,
                                        s2list, niters-nwarmup)
      loss <- loss + steinsLoss(solve(SigmaHatInvPM), SigmaInvK)
    }
    loss <- loss/dat$ngroups
      print(loss)
      
    lossVec[3] <- loss
    resList[[3]] <- res

  }

  print(sprintf("Finished fitting %i with rank = %i", fitType, R))

}
