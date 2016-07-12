load("SdimLoss.Rdata")
pdf("paper/Figs/LossVsDimension.pdf", font="Times")
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(5:50, lossLst, type="l", lwd=3, ylim=c(0, 100), xlim=c(0, 50), ylab="Stein's Risk", xlab=expression(hat(s)), cex.axis=2, cex.lab=2)
grid()
lines(5:50, lossLst, type="l", lwd=3, ylim=c(0, 120))
abline(h=3, lty=2, col="blue", lwd=3)
abline(h=0.89, lty=2, col="red", lwd=3)

legend("topright", lty=2, lwd=3, col=c("blue", "red"),
       legend=c(expression(VV^T==I[p]),
                expression(paste(VV^T, "= span(", U[1], ",..., ", U[K], ")"))),
                bty='n', cex=2)
dev.off()

pdf(sprintf("paper/Figs/simRatio-s5.pdf", format(Sys.Date(), "%m-%d")), font="Times")
load("sdimtest-5-1.RData")
evalRatios <- with(resList[[1]],
    sapply(1:ngroups, function(k) {
    numer <- sum(eigen(t(V[, 1:S]) %*% dat$Slist[[k]] %*% V[, 1:S])$values)/dat$nvec[k]
    denom <- sum(svd(dat$Ylist[[k]])$d[1:S]^2)/dat$nvec[k]
    correction <-  (denom - 1/mean(s2samps[k, ]) * S * P / dat$nvec[k])  / denom
    (numer/denom) / correction
    }))

par(mar=c(5.1, 4.5, 4.1, 2.1))
barplot(evalRatios, ylim=c(0, 1), xlab="Group", ylab="Goodness of Fit", main="",
        space=0.2, cex.axis=2, col="#646464", col.axis="black",
        col.main="#646464", border=NA,
        names.arg=1:10, cex.names=2, cex.main=1, cex.lab=2, col.lab="black")
dev.off()

pdf(sprintf("paper/Figs/simRatio-s20.pdf", format(Sys.Date(), "%m-%d")), font="Times")
load("sdimtest-20-1.RData")
evalRatios <- with(resList[[1]],
    sapply(1:ngroups, function(k) {
    numer <- sum(eigen(t(V[, 1:S]) %*% dat$Slist[[k]] %*% V[, 1:S])$values)/dat$nvec[k]
    denom <- sum(svd(dat$Ylist[[k]])$d[1:S]^2)/dat$nvec[k]
    correction <-  (denom - 1/mean(s2samps[k, ]) * S * P / dat$nvec[k])  / denom
    (numer/denom) / correction
    }))

par(mar=c(5.1, 4.5, 4.1, 2.1))
barplot(evalRatios, ylim=c(0, 1), xlab="Group", ylab="Goodness of Fit", space=0.2, cex.axis=2, col="#646464", col.axis="black", col.lab="black", border=NA,
        names.arg=1:10, cex.names=2, cex.lab=2)
dev.off()


## gen: S=R=2, fit: S=R=2
load("rankResults1.RDat_a")
ssMat1 <- ssMat
fullMat1 <- fullMat
ratioMat1 <- ratioMat

## gen: S=10, R=2, fit:S=R=2
load("rankResults2.RData")
ssMat2 <- ssMat
fullMat2 <- fullMat
ratioMat2 <- ratioMat

pdf("lossBoxplot.pdf")
boxplot(as.numeric(ssMat1), as.numeric(ssMat2), pch=19, cex=0.5, ylab="Stein's Loss", outline=FALSE, names=c("S=2", "S=10"))
dev.off()

pdf("paper/Figs/evalRatios.pdf", width=14)
par(mfrow=c(1,2))
barplot(ratioMat1[1, ], ylim=c(0, 1), xlab="Group", main="S=2", space=0.2, cex.axis=2, col="#646464", col.axis="#646464", col.main="#646464", border=NA, cex.main=2,
        names.arg=1:10)
barplot(ratioMat2[1, ], ylim=c(0, 1), xlab="Group", main="S=10", space=0.2, cex.axis=2, col="#646464", col.axis="#646464", col.main="#646464", border=NA, cex.main=2,
        names.arg=1:10)
dev.off()


load("rankTest-1-1.RData")
res <- resList[[1]]
pdf("paper/Figs/posteriorRegions.pdf")
## with(res, vectorPlot(Osamps, omegaSamps, s2samps, nsamps=1000, groups=1:10, probRegion=0.95, type=1))
with(res, posteriorPlot(Osamps, omegaSamps, s2samps, nsamps=1000, groups=1:10, probRegion=0.95, plotPoints=FALSE))
sv <- svd(t(res$V) %*% dat$V)
Rot <- sv$v %*% t(sv$u)
## abline(v=sapply(dat$Olist, function(O) {
##     O <- t(Rot) %*% O
##     atan(O[, 1][2]/O[, 1][1])}),
##     col=1:10, lwd=2, lty=2)

segments(x0=sapply(dat$Olist, function(O) {
    O <- t(Rot) %*% O
    atan(O[, 1][2]/O[, 1][1])}),
    y0=0, y1=2,
    col=1:10, lwd=3, lty=1)
dev.off()
