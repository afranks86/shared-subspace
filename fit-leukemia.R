rm(list=ls())

library(car)
library(rstiefel)
library(xtable)
source("helper.R")
source("fit-subspace.R")
source("https://bioconductor.org/biocLite.R")
## biocLite("hgu95av2.db")
## library(hgu95av2.db)
load("LeukemiaData/leukemia.RData")

## number of eigenvectors if pooled space
S <- 14
## number of eigenvectors in group subspace
R <- 2
P <- ncol(Ypooled)

residualList <- Slist <- list()
for( k in 1:ngroups ) {
    residual <- lm(Ylist[[k]] ~ 1)$residuals
    residual <- t(lm(t(residual) ~ 1)$residuals)
    residualList[[k]] <- residual
    Slist[[k]] <- t(residual) %*% residual
}
pooledResiduals <- do.call(rbind, residualList)
plot(svd(residualList[[ngroups]]/sqrt(nvec[ngroups]))$d^2, pch=19, cex=0.75)
for(k in 1:6) {
    points(svd(residualList[[k]]/sqrt(nvec[k]))$d^2, col=k+1, pch=19, cex=0.75)
}
svd(residualList[[ngroups]]/sqrt(nvec[ngroups]))$d^2

Vinit <- svd(do.call(cbind, lapply(1:ngroups, function(k) svd(t(residualList[[k]]))$u[, 1:R])))$u[, 1:min(ngroups*R, S)]
## Vinit <- Vinit[, sample(1:ncol(Vinit))]
## Vinit <- rustiefel(P, S)
PrecVec = rep(100, ngroups)
EMFit <- subspaceEM(Slist, P=P, S=S, R=R, nvec=nvec, Vstart=Vinit, PrecVec=PrecVec, verbose=TRUE, rho1=1e-5, rho2=1-1e-5)
Vinit <- EMFit$V

save(EMFit, Vinit, S, R, P, nvec, residualList,
     file=sprintf("leukemiaEM-%s.RData", format(Sys.Date(), "%m-%d")))

sapply(1:ngroups, function(k) {
    sum(eigen(t(Vinit) %*% Slist[[k]] %*% Vinit)$values) / sum(svd(residualList[[k]])$d[1:S]^2)
})

evalRatios <- sapply(1:ngroups, function(k) {
    numer <- sum(eigen(t(Vinit[, 1:S]) %*% Slist[[k]] %*% Vinit[, 1:S])$values)/nvec[k]
    denom <- sum(svd(residualList[[k]])$d[1:S]^2)/nvec[k] - 1/EMFit$PrecVec[k] * S * P / nvec[k]
    (numer/denom)
})


evalRatiosQuad <- sapply(1:ngroups, function(k) {
    numer <- sum(eigen(t(Vinit[, 1:S]) %*% Slist[[k]] %*% Vinit[, 1:S])$values)/nvec[k]

    evals <- svd(residualList[[k]])$d[1:S]^2/nvec[k]
    b <- (1/EMFit$PrecVec[k]*P/nvec[k] - evals - 1)

    quadSol <- (-b + sqrt(b^2 - 4*evals))/2
    quadSol[is.nan(quadSol)] <- (evals - (1/EMFit$PrecVec[k]*P/nvec[k]))[is.nan(quadSol)]
    denom <- sum(quadSol)
    

    (numer/denom)
})



pdf(sprintf("paper/Figs/leukemiaRatio-%s.pdf", format(Sys.Date(), "%m-%d")), font="Times")
par(mar=c(6.1, 4.1, 4.1, 2.1))
barplot(evalRatios, ylim=c(0, 1), xlab="", main="", space=0.2, cex.axis=1.5, col="#646464", col.axis="black", col.main="black", border=NA)
text(seq(.8, 8.8, by=1.2) , par("usr")[1]+0.1, labels=unique(typeVec), srt=45, xpd=TRUE, cex=1.5, col="black", adj=c(1, NA))
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
                   niters=1000, draw=c(V=FALSE))
samples$Usamps <- NULL
save(samples, initSS, file=sprintf("leukemiaBayes-%s.RData", format(Sys.Date(), "%m-%d")))

## Posterior plots
pdf(sprintf("paper/Figs/leukemiaPosterior-%s.pdf", format(Sys.Date(), "%m-%d")), font="Times")
with(samples,
     posteriorPlot(Osamps[1:R, 1:R, 2:7, ], omegaSamps[1:R, 2:7, ],
                s2samps[2:7, ], nsamps=1000, groups=1:(ngroups-1),
                probRegion=0.95, hline=NULL, logRatio=TRUE,
                plotPoints=FALSE, col=2:7, ymax=log(80)))
group1pts <- with(samples, getHullPoints(nsamps=100, omegaSamps[1:2, 1, ], Osamps[1:2, 1:2, 1, ], logRatio=TRUE)$allPts)

split1 <- which(group1pts[1, ] < 0)
split2 <- which(group1pts[1, ] >= 0)
hull1 <- chull(group1pts[1, split1], group1pts[2, split1])
hull2 <- chull(group1pts[1, split2], group1pts[2, split2])

pts1 <- group1pts[, split1[hull1]]
pts1[2, 4] <- min(pts1[2, ])
pts2 <- group1pts[, split2[hull2]]
pts2 <- pts2[, 3:ncol(pts2)]
pts2 <- cbind(c(pi/2, min(pts1[2, ])), pts2)
polygon(pts1[1, ], pts1[2, ], lwd=0.01, border="black", col=alpha("black", 1/4), lty=1)
polygon(pts2[1, ], pts2[2, ], lwd=0.01, border="black", col=alpha("black", 1/4), lty=1)
lines(pts1[1, c(7:9, 1:4)], pts1[2, c(7:9, 1:4)], col="black", lwd=3)
lines(pts2[1, 1:5], pts2[2, 1:5], col="black", lwd=3)
legend("topright", legend=unique(typeVec), col=1:7, pch=19, bty="n")
dev.off()

pdf(sprintf("paper/Figs/leukemia-biplot-%s.pdf", format(Sys.Date(), "%m-%d")), font="Times")
mag <- apply(EMFit$V[, 1:2], 1, function(x) sqrt(sum(x^2)))
indices <- which(mag > quantile(mag, 0.99))

Vsub <- EMFit$V[indices, ]
topright <- which(Vsub[, 1] > 0 & Vsub[, 2] > 0)
topleft <- which(Vsub[, 1] < 0 & Vsub[, 2] > 0)
bottom <- which(Vsub[, 2] < 0)

plot(EMFit$V[, 1:2], col="light grey", pch=19, cex=0.5,
     xlim=c(-0.1, 0.1), ylim=c(-0.1, 0.1), xlab="V1", ylab="V2")
points(Vsub[topright, 1], Vsub[topright, 2], xlab="V1", ylab="V2", cex=1.5,
     xlim=c(-0.1, 0.1), ylim=c(-0.1, 0.1), pch=15, col="black")
points(Vsub[topleft, 1], Vsub[topleft, 2], xlab="V1", ylab="V2", cex=1.5,
       xlim=c(-0.1, 0.1), ylim=c(-0.1, 0.1), pch=17, col="black")
points(Vsub[bottom, 1], Vsub[bottom, 2], xlab="V1", ylab="V2", cex=1.5,
     xlim=c(-0.1, 0.1), ylim=c(-0.1, 0.1), pch=16, col="black")
#text(EMFit$V[indices, 1], EMFit$V[indices, 2], labels=rownames(EMFit$V)[indices],
                                        #     xlab="V1", ylab="V2", cex=0.5, pos=4)

abline(h=0, v=0)
legend("topright", legend=unique(typeVec)[1:3], lty=1, lwd=2, col=1:3, cex=1, bty='n', ncol=1)


probeMap <- as.list(hgu95av2SYMBOL)
res <- unlist(probeMap[rownames(EMFit$V)])[indices]


trnames <- sort(res[topright])
tlnames <- sort(res[topleft])
bnames <- sort(res[bottom])
length(trnames) <- length(tlnames) <- length(bnames) <-
    max(length(trnames), length(tlnames), length(bnames))

geneList <- cbind(trnames, tlnames, bnames)
rownames(geneList) <- rep("", nrow(geneList))
colnames(geneList) <- c("Top Right ($\\blacksquare$)", "Top Left ($\\blacktriangle$)", "Bottom (\\tikz\\draw[black,fill=black] (0,0) circle (.7ex);)")

print(xtable(geneList, na="", quote=FALSE),
      include.rownames=FALSE, include.colnames=TRUE,
      floating=FALSE, sanitize.text.function=identity,
      file="paper/Figs/biplot.tab")

for(k in 1:3) {

    eigK <- eigen(solve(EMFit$PhiList[[k]][1:2, 1:2])*EMFit$PrecVec[k]-diag(2))
    lambda <- eigK$values
    print(lambda)
    evecs <- eigK$vectors

    maxIndex <- which.max(lambda)
    lamRatio <- lambda[maxIndex]/lambda[-maxIndex]
    angle <- atan(evecs[2, maxIndex]/evecs[1, maxIndex])
    print(angle)
    ellipse(center=c(0,0), shape=evecs[, 1:2] %*% diag(lambda[1:2]) %*% t(evecs[, 1:2])/lambda[maxIndex], radius=0.01, col=k, center.cex=000)
    ellipse(center=c(0,0), shape=evecs[, 1:2] %*% diag(lambda[1:2]) %*% t(evecs[, 1:2])/lambda[maxIndex], radius=0.025, col=k, center.cex=000)
    ellipse(center=c(0,0), shape=evecs[, 1:2] %*% diag(lambda[1:2]) %*% t(evecs[, 1:2])/lambda[maxIndex], radius=0.04, col=k, center.cex=000)
    ## points(angle, lamRatio, col=k, pch=19, cex=2)

}
dev.off()

plot(residualList[[1]] %*% EMFit$V[, 1:2], xlim=c(-10, 10), ylim=c(-10, 10), pch=19)
for(k in 2:ngroups) {
    points(residualList[[k]] %*% EMFit$V[, 1:2], xlim=c(-10, 10), ylim=c(-10, 10),
         col=k, pch=19)
}

## Check magnitude of eigenvalues vs P/n
eigVals <- sapply(1:ngroups, function(g) with(samples, head(rowMeans(omegaSamps[, g, ]/(1-omegaSamps[, g, ])), n=2)))

rowMeans(sapply(1:R, function(i) (1 - sqrt(P/nvec)/eigVals[i, ]^2) / (1 + sqrt(P/nvec)/eigVals[i, ])))
sapply(1:R, function(i) (1 - sqrt(P/(ngroups*min(nvec)))/min(eigVals[i, ])^2) / (1 + sqrt(P/(ngroups*min(nvec)))/min(eigVals[i, ])))
