rm(list=ls())

library(rstiefel)
source("helper.R")
source("fit-subspace.R")
##load("ae_XY_m527.RData")
load("C18_all.RData")

## number of eigenvectors if pooled space
S <- 12
## number of eigenvectors in group subspace
R <- 2
P <- ncol(Y)

table(X$Line)
table(X$Sex, X$Line)
table(X$Age, X$Line)


## New age groups
X$Sex <- as.factor(X$Sex)
X$Age <- as.factor(X$Age)
levels(X$Age) <- c("Young", "Young", "Middle", "Middle", "Old", "Old", "Old")

table(X$Age, X$Sex)

XDM <- model.matrix(~ -1 + Age:Sex, data=X )
ngroups <- ncol(XDM)
## XDM_All <- model.matrix(~ -1 + Age + Sex + Line, data=X )
## XDM_Age <- model.matrix(~ -1 + Age, data=X )

## fit <- lm(Y~-1+XDM_All)
fit <- lm(Y~-1+XDM)

## ngroups <- ncol(XDM_Age) 

nvec <- rep(0, ngroups)
residualList <- Slist <- list()
for( i in 1:ngroups ) {
    ##indices <- which(XDM_Age[, i]==1)
    indices <- which(XDM[, i]==1)
    res <- fit$residuals[indices, ]
    getRank(res)
    residualList[[i]] <- res
    Slist[[i]] <-     cur <- t(res) %*% res
    nvec[i] <- nrow(res)
    
}
getRank(fit$residuals)

## Vinit <- svd(do.call(cbind, lapply(1:ngroups, function(k) svd(t(residualList[[k]]))$u[, 1:R])))$u[, 1:min(ngroups*R, S)]
Vinit <- svd(do.call(cbind, lapply(1:ngroups, function(k) svd(t(residualList[[k]]))$u[, 1:S])))$u[, 1:S]
PrecVec = rep(5, ngroups)
EMFit <- subspaceEM(Slist, P=P, S=S, R=R, nvec=nvec, Vstart=Vinit, PrecVec=PrecVec, verbose=TRUE, rho1=1e-5, rho2=1-1e-5)
Vinit <- EMFit$V

## In paper using, 04-28
save(EMFit, Vinit, S, R, P, nvec, residualList,
     file=sprintf("dmelanEM-%s.RData", format(Sys.Date(), "%m-%d")))

evalRatios <- sapply(1:ngroups, function(k) {
    numer <- sum(eigen(t(Vinit[, 1:S]) %*% Slist[[k]] %*% Vinit[, 1:S])$values)/nvec[k]
    denom <- sum(svd(residualList[[k]])$d[1:S]^2)/nvec[k]
    correction <-  (denom - 1/EMFit$PrecVec[k] * S * P / nvec[k])  / denom
    (numer/denom) / correction
})

pdf(sprintf("paper/Figs/dmelanRatio-%s.pdf", format(Sys.Date(), "%m-%d")), font="Times")
par(mar=c(6.1, 4.1, 4.1, 2.1))
barplot(evalRatios, ylim=c(0, 1), xlab="", main="", space=0.2, cex.axis=1.5, col="#646464", col.axis="black", col.main="black", border=NA)
text(seq(1, 7, length.out=ngroups) , par("usr")[1], labels=paste(rep(c("Young", "Middle", "Old"), 2), rep(c("F", "M"), each=ngroups/2), sep="/"),
     srt=45, xpd=TRUE, cex=1.5, col="black", adj=c(1, NA))
dev.off()

OmegaList <- Ulist <- list()
for(k in 1:ngroups) {

    eigK <- eigen(diag(S) - EMFit$PhiList[[k]] / EMFit$PrecVec[k])

    OmegaList[[k]] <- eigK$values
    Ulist[[k]] <- Vinit %*% eigK$vectors

}
s2vec <- 1/EMFit$PrecVec

initSS <- list(V=Vinit, Ulist=Ulist, OmegaList=OmegaList, s2vec=s2vec)
samples <- fitSubspace(P, S, R, Q=S-R, Slist,
                   nvec, ngroups, init=initSS,
                   niters=200, draw=c(V=FALSE))
save(samples, initSS, file=sprintf("dmelanBayes-%s.RData",
                                   format(Sys.Date(), "%m-%d")))


## Posterior plots
pdf(sprintf("paper/Figs/dmelanPosterior-%s.pdf", format(Sys.Date(), "%m-%d")), font="Times")
with(samples,
     posteriorPlot(Osamps[1:2, 1:2, , ], omegaSamps[1:2, , ],
                s2samps, nsamps=100, groups=1:ngroups,
                probRegion=0.95, hline=NULL,
                col=rep(1:3, 2), pch=rep(c(2, 3), each=3), lty=rep(c(1, 2), each=3), ymax=log(8), logRatio=TRUE, plotPoints=FALSE))
## legend("topright", legend=paste(rep(c("Young", "Middle", "Old"), 2), rep(c("F", "M"), each=3), sep="/"), col=rep(1:3, 2), pch=rep(c(2, 3), each=3), bty="n")
legend("topright", legend=paste(rep(c("Young", "Middle", "Old"), 2), rep(c("F", "M"), each=3), sep="/"), col=rep(1:3, 2), lty=rep(c(1, 2), each=3), bty="n", lwd=3)
dev.off()

for(k in 1:ngroups) {

    eigK <- eigen(solve(EMFit$PhiList[[k]][1:2, 1:2])*EMFit$PrecVec[k]-diag(2))
    lambda <- eigK$values
    evecs <- eigK$vectors

    maxIndex <- which.max(lambda)
    lamRatio <- lambda[maxIndex]/lambda[-maxIndex]
    angle <- atan(evecs[2, maxIndex]/evecs[1, maxIndex])

    points(angle, lamRatio, col=k, pch=19, cex=2)
    print(sprintf("%f, %f", angle, lamRatio))
}

with(samples,
     vectorPlot(Osamps[c(3,7), c(3,7), , ], omegaSamps[c(3,6), , ],
                s2samps, nsamps=100, groups=1:ngroups,
                probRegion=0.95, type=2, hline=1))

mag <- apply(EMFit$V[, 1:2], 1, function(x) sqrt(sum(x^2)))
indices <- which(mag > quantile(mag, 0.95))
plot(EMFit$V[indices, 1], EMFit$V[indices, 2], xlab="V1", ylab="V2", cex=0.5,
     xlim=c(-0.1, 0.1), ylim=c(-0.1, 0.1))
text(EMFit$V[indices, 1], EMFit$V[indices, 2], labels=rownames(EMFit$V)[indices],
     xlab="V1", ylab="V2", cex=0.5, pos=4)
abline(h=0, v=0)

## Check magnitude of eigenvalues vs P/n
eigVals <- sapply(1:ngroups, function(g) with(samples, head(rowMeans(omegaSamps[, g, ]/(1-omegaSamps[, g, ])), n=S)))

rowMeans(sapply(1:R, function(i) (1 - sqrt(P/nvec)/eigVals[i, ]^2) / (1 + sqrt(P/nvec)/eigVals[i, ])))
sapply(1:R, function(i) (1 - sqrt(P/(ngroups*min(nvec)))/min(eigVals[i, ])^2) / (1 + sqrt(P/(ngroups*min(nvec)))/min(eigVals[i, ])))
