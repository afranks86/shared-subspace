rm(list=ls())

## Need these for Rscript
library(methods)
library(utils)
library(microbenchmark)
library(rstiefel)
library(tidyverse)
library(tidyverse)
library(cowplot)
library(colorspace)

source("fit-subspace.R")
source("generateData.R")
source("helper.R")

########################################################
############# Data Generation #########################
########################################################



########################################################
############# Data Fit #########################
########################################################

nbench_times <- 10
Svec <- c(2, 10)
Pvec <- seq(from=1000, to=3000, by=1000)


## Initialize sampler
ngroups <- dat$ngroups
nvec <- dat$nvec

benchmark_array <- array(dim=c(length(Svec), length(Pvec), nbench_times))
for(sindex in 1:length(Svec)) {
    for(pindex in 1:length(Pvec)) {
        
        n <- 50
        S <- Svec[sindex]
        R <- S
        P <- Pvec[pindex]

        print(sprintf("P = %i, S = %i", P, S))
        
        ngroups <- 5

        niters <- 1000
        nwarmup <- niters/2

        dat <- generateData(P=P, S=S, R=S, ngroups=ngroups,
                            nvec=rep(n, ngroups), saveSigmaList=FALSE, saveSList = FALSE)

        dat$genType <- "Shared subspace"

        Ypooled <- c()
        for(i in 1:length(dat$Ylist)){
            Ypooled <- rbind(Ypooled, dat$Ylist[[i]])
        }
        dat$Ypooled <- Ypooled
        
        Vinit <- svd(do.call(cbind, lapply(1:ngroups, function(k) svd(t(dat$Ylist[[k]]))$u[, 1:R])))$u[, 1:S]
##        Vinit <- rustiefel(P, S)
        benchmarked <- microbenchmark(EMFit <- subspaceEM(dat$Ylist, P=P, S=S, R=R, nvec=dat$nvec, Vstart=Vinit, verbose=TRUE, maxIters=1000, stiefelAlgo=2), times=nbench_times)
        Vinit <- EMFit$V

        benchmark_array[sindex, pindex, ] <- benchmarked$time/1e9
    }
}

dimnames(benchmark_array) = list(S = Svec, P = Pvec, Iter =1:nbench_times)
as_tibble(reshape2::melt(benchmark_array)) %>% ggplot(aes(x=as.factor(P), y=value)) + geom_boxplot(aes(fill=factor(S))) + scale_fill_discrete_qualitative(palette="Set 2") +
    ylab("Seconds")
