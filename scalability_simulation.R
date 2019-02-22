rm(list=ls())

## Need these for Rscript
library(microbenchmark)
library(rstiefel)
library(tidyverse)
library(tidyverse)
library(cowplot)
library(colorspace)

source("fit-subspace.R")
source("generateData.R")
source("subspace-functions.R")


nbench_times <- 10
Svec <- c(2, 10, 25, 50)
Pvec <- seq(from=1000, to=10000, by=1000)

## Initialize sampler
ngroups <- dat$ngroups
nvec <- dat$nvec

benchmark_array <- array(dim=c(length(Svec), length(Pvec), nbench_times))
for(sindex in 1:length(Svec)) {
    for(pindex in 1:length(Pvec)) {
        for(nt in 1:nbench_times) {

            n <- 50
            S <- Svec[sindex]
            R <- S
            P <- Pvec[pindex]
            
            print("------------------------------")
            print(sprintf("P = %i, S = %i", P, S))
            print("------------------------------")
            
            ngroups <- 5

            niters <- 1000
            nwarmup <- niters/2

            dat <- generateData(P=P, S=S, R=S, ngroups=ngroups,
                                nvec=rep(n, ngroups))

            dat$genType <- "Shared subspace"

            Ypooled <- c()
            for(i in 1:length(dat$Ylist)){
                Ypooled <- rbind(Ypooled, dat$Ylist[[i]])
            }
            dat$Ypooled <- Ypooled

            Vinit <- rustiefel(P, S)
            benchmarked <- microbenchmark(EMFit <- subspaceEM(dat$Ylist, P=P, S=S, R=R, nvec=dat$nvec, Vstart=Vinit, verbose=FALSE, EM_iters=100, M_iters=50, stiefelAlgo=2), times=1)
            Vinit <- EMFit$V

            benchmark_array[sindex, pindex, nt] <- benchmarked$time/1e9
            
        }
    }
}

save(benchmark_array, file="results/benchmark_times.RData")

pdf("paper/Figs/run_times.pdf", height=4, width=8)
dimnames(benchmark_array) = list(S = Svec, P = Pvec, Iter = 1:nbench_times)
as_tibble(reshape2::melt(benchmark_array)) %>% ggplot(aes(x=as.factor(P), y=value)) +
    geom_jitter(aes(col=factor(S))) + scale_fill_discrete_qualitative(palette="Set 2") +
    ylab("Seconds") + ylim(c(0, 300)) +
    xlab("P") + scale_colour_discrete(name  = "S")

dev.off()

