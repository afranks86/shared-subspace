rm(list=ls())

## Need these for Rscript
library("methods")
library("utils")
library(foreach)
library(rstiefel)
library(doMC)
registerDoMC(cores=4)

source("fit-subspace.R")
source("generateData.R")
source("helper.R")

########################################################
############# Data Generation #########################
########################################################

n <- 50
S <- 2
R <- S
P <- 200

maxK <- 10
nrepl <- 30

## runs with three different sets of eigenvalues for psi_k
resMatList <- foreach(e=1:3) %dopar% {

    resMat <- matrix(0, nrow=4, ncol=maxK)
    ## 10 replicates per generated data
    for(repl in 1:nrepl) {
        for(ngroups in 1:maxK) {
            ## P <- ngroups * n

            if(e==1)
                evals <- c(5, 2)
            else if(e==2)
                evals <- c(50, 2)
            else
                evals <- c(5, 5)
            ## evals <- c(6, 5, 4, 3, 2)

            ## Subspace Shrinkage: Psi_k = VO_kLam_kO_k^TV^T
            LambdaList <- list()
            lam <- evals
            for( k in 1:ngroups) {
                LambdaList[[k]] <- lam
            }
            
            dat <- generateData(P=P, S=S, R=R, ngroups=ngroups,
                                LambdaList=LambdaList, nvec=rep(n, ngroups))
            dat$genType <- "Shared subspace"
            nvec <- dat$nvec

            Ypooled <- c()
            for(i in 1:length(dat$Ylist)){
                Ypooled <- rbind(Ypooled, dat$Ylist[[i]])
            }
            dat$Ypooled <- Ypooled

            Slist <- dat$Slist

            ## Vinit <- rustiefel(P, S)
            Vinit <- svd(do.call(cbind, lapply(1:ngroups, function(k) svd(t(dat$Ylist[[k]]))$u[, 1:S])))$u[, 1:S]
            EMFit <- subspaceEM(dat$Slist, P=P, S=S, R=R, nvec=dat$nvec,
                                PrecVec=rep(10, ngroups),
                                rho1=0.1, rho2=0.9, 
                                Vstart=Vinit, verbose=FALSE, maxIters=100)
            Vinit <- EMFit$V

            ## tr(VhatVhat^TVV^T)
            resMat[1, ngroups] <- resMat[1, ngroups] +
                tr(Vinit %*% t(Vinit) %*% dat$V %*% t(dat$V))/S
            ## lower bound
            cond <- evals <= sqrt(P/sum(nvec))
            resMat[2, ngroups] <- resMat[2, ngroups] +
                sum(ifelse(cond, 0, 
                           (1-P/sum(nvec)/evals[!cond]^2) / (1+P/sum(nvec)/evals[!cond])))/S
            ## crude upper bound
            resMat[3, ngroups] <- resMat[3, ngroups] +
                sum( (1-P/sum(nvec)/rep(max(evals), S)^2) / (1+P/sum(nvec)/rep(max(evals), S)) )/S
            ## better upper bound?
            resMat[4, ngroups] <- resMat[4, ngroups] +
                sum( (1-P/sum(nvec)/rep(mean(evals), S)^2) / (1+P/sum(nvec)/rep(mean(evals), S)) )/S

            print(sprintf("-------------- %i --------------", ngroups))

        }
        print(sprintf("-------------- Repl = %i --------------", repl))
    }
    resMat <- resMat/nrepl

    resMat
}
##save(resMatList, file="resMatList_K10.RData")

load("resMatList_K10.RData")

pdf("paper/Figs/asymptotics.pdf", width=14)
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(resMatList[[1]][1, 1:maxK], type="l", lwd=3, xlim=c(1, maxK), ylim=c(0, 1),
     lty=1, col="red",
     xlab="Number of Groups", ylab=expression("Subspace Accuracy"),
     cex.axis=2, cex.lab=2)
lines(resMatList[[1]][2, 1:maxK], lwd=3, col="red", lty=2)

lines(resMatList[[2]][1, 1:maxK], lwd=3, col="blue")
lines(resMatList[[2]][2, 1:maxK], lwd=3, col="blue", lty=2)

lines(resMatList[[3]][1, 1:maxK], lwd=3, col="dark green")
lines(resMatList[[3]][2, 1:maxK], lwd=3, col="dark green", lty=2)
grid()
legend("bottomright", legend=c("(5, 5)", "(5, 2)", "(50, 2)"), lwd=3, col=c("dark green", "red", "blue"), cex=1.5, bg="white", title=expression("("~lambda[1]~","~lambda[2]~")"))
dev.off()

par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(resMatList[[1]][1, ], type="l", lwd=3, ylim=c(0, 1), lty=1, col="red",
     xlab="Number of Groups", ylab=expression("Subspace Similarity"),
     cex.axis=2, cex.lab=2)
lines(resMatList[[1]][4, ], lwd=3, col="red", lty=2)

lines(resMatList[[2]][1, ], lwd=3, col="blue")
lines(resMatList[[2]][4, ], lwd=3, col="blue", lty=2)

lines(resMatList[[3]][1, ], lwd=3, col="dark green")
lines(resMatList[[3]][4, ], lwd=3, col="dark green", lty=2)
grid()
legend("bottomright", legend=c("(5, 5)", "(5, 2)", "(50, 2)"), lwd=3, col=c("dark green", "red", "blue"), cex=1.5, bg="white", title=expression("("~lambda[1]~","~lambda[2]~")"))
