rm(list=ls())

## Need these for Rscript
library(microbenchmark)
library(rstiefel)
library(tidyverse)
library(cowplot)
library(colorspace)

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
                    Vstart=Vinit, verbose=TRUE, EM_iters=100, M_iters=100, stiefelAlgo=2)

svd_init <- tr((t(dat$V) %*%  EMFit$V) %*% (t(EMFit$V) %*% dat$V)) / S

nruns <- 100
similarity <- numeric(runs)

Vinit_base <- svd(do.call(cbind, lapply(1:ngroups, function(k) svd(t(dat$Ylist[[k]]))$u[, 1:n])))$u[, 1:(n*ngroups)]

for(i in nruns) {

    print(i)
    Vinit <- Vinit_base %*% rustiefel(n*ngroups, S)
    EMFit <- subspaceEM(dat$Ylist, P=P, S=S, R=R, nvec=dat$nvec,
                    Vstart=Vinit, verbose=TRUE, EM_iters=100, M_iters=1000, stiefelAlgo=2)
    random_init <- tr((t(dat$V) %*%  EMFit$V) %*% (t(EMFit$V) %*% dat$V)) / S
    similarity[i] <- random_init
}

load("results/similarity_scores.RData")

pdf("paper/Figs/initialization_fig.pdf")
cols <- qualitative_hcl(2, palette="Set 2")
tibble(score=similarity) %>% ggplot() + geom_density(aes(x=score), fill=cols[1]) + geom_vline(xintercept=svd_init, size=2, col=cols[2], linetype="dashed")  +xlim(c(0, 1)) + xlab("Subspace Similarity") + ylab("") + theme_classic(base_size=16) + ggtitle("Random initializations vs Eigen Initialization")
dev.off()


