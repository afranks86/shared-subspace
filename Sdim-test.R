library(rstiefel)

source("fit-subspace.R")
source("generateData.R")
source("subspace-functions.R")

########################################################
############# Data Generation #########################
########################################################

n <- 50
S <- 2
R <- 2
P <- 200
ngroups <- 10

evals <- c(250, 25)

niters <- 1000
nwarmup <- niters/2



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
dat$genType <- "No pooling"
  



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

S <- idx
R <- 2


ngroups <- dat$ngroups
nvec <- dat$nvec
Slist <- dat$Slist

Vinit <- svd(do.call(cbind, lapply(1:ngroups, function(k) svd(t(dat$Ylist[[k]]))$u[, 1:R])))$u[, 1:min(S, R*ngroups)]
if(S > R*ngroups) {
    NC <- NullC(Vinit)
    Vinit <- cbind(Vinit, NC[, sample(1:ncol(NC), size=(S - R*ngroups))])
}
Vinit <- subspaceEM(dat$Slist, P=P, S=S, R=R, Q=0, dat$nvec,
                    Vstart=Vinit, verbose=FALSE,
                    rho1=0.1, rho2=0.9, maxIters=100)$V

##cancor(svd(do.call(cbind, dat$Olist))$u, Vinit)$cor

OmegaList <- Ulist <- list()
for(k in 1:ngroups) {
    OmegaList[[k]] <- rep(1/2, R)

    Ok <- matrix(0, nrow=S, ncol=R)
    Ok[1:R, 1:R] <- diag(R)
    
    Ulist[[k]] <- Vinit %*% Ok    
}
s2vec <- rexp(ngroups)

initSS <- list(V=Vinit, Ulist=Ulist, OmegaList=OmegaList, s2vec=s2vec)
res <- fitSubspace(dat$P, S, R, Q=0, Slist=dat$Slist,
                   nvec=dat$nvec, ngroups=dat$ngroups, init=initSS,
                   niters=niters, sigmaTruthList=dat$SigmaList,
                   draw=c(V=FALSE))

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

lossVec[1] <- loss
resList[[1]] <- res
    
print(sprintf("Finished fitting %i with rank = %i", fitType, R))
save(resList, dat, lossVec, datType,
     file=sprintf("data/sdimtest-%i-%i.RData", datType, idx))
