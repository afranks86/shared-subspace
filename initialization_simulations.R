rm(list=ls())

## Need these for Rscript
library(microbenchmark)
library(rstiefel)
library(tidyverse)

source("fit-subspace.R")
source("generateData.R")
source("helper.R")

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
Olist[[3]] <- matrix(c(-sqrt(2)/2, sqrt(2)/2, sqrt(2)/2, sqrt(2)/2), nrow=2)
Olist[[4]] <- diag(2)
Olist[[5]] <- diag(2)

Opooled <- matrix(0, nrow=S, ncol=S)
Opooled[1:2, 1:2] <- diag(2)
Opooled[, (R+1):S] <- NullC(Opooled[, 1:2]) %*% rustiefel(P, S - R)
LambdaPooled <- rexp(S -R, 1/500)
for(k in 1:ngroups) {
    Ok <- matrix(0, nrow=S, ncol=S)
    Ok[1:R, 1:R] <- Olist[[k]]
    Ok[, (R+1):S] <- Opooled[, (R+1):S]
    Olist[[k]] <- Ok
    LambdaList[[k]] <- c(LambdaList[[k]], LambdaPooled)
}
    
dat <- generateData(P=P, S=S, R=R, ngroups=ngroups,
                    V = rustiefel(P, S),
                    LambdaList=LambdaList, Olist=Olist,
                    nvec=rep(n, ngroups))

ngroups <- dat$ngroups
nvec <- dat$nvec

Vinit <- svd(do.call(cbind, lapply(1:ngroups, function(k) svd(t(dat$Ylist[[k]]))$u[, 1:S])))$u[, 1:S]
EMFit <- subspaceEM(dat$Ylist, P=P, S=S, R=R, nvec=dat$nvec,
                    Vstart=Vinit, verbose=TRUE, max_EM_iters=10, stiefelAlgo=2)
svd_init <- tr((t(dat$V) %*%  EMFit$V) %*% (t(EMFit$V) %*% dat$V)) / S

nruns <- 100
similarity <- numeric(runs)
for(i in nruns) {

    print(i)
    Vinit <- rustiefel(P, S)
    EMFit <- subspaceEM(dat$Ylist, P=P, S=S, R=R, nvec=dat$nvec,
                        Vstart=Vinit, verbose=TRUE, max_EM_iters=10, stiefelAlgo=2)
    random_init <- tr((t(dat$V) %*%  EMFit$V) %*% (t(EMFit$V) %*% dat$V)) / S
    similarity[i] <- random_init
}

tibble(score=similarity) %>% ggplot() + geom_density(aes(x=score), fill=alpha("red", 0.5)) + geom_vline(xintercept=svd_init, size=2, col="blue")  +xlim(c(0, 1)) + xlab("Subspace Similarity") + ylab("") + theme_grey(base_size=16)



