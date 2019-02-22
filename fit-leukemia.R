rm(list=ls())

library(car)
library(rstiefel)
library(xtable)
library(made4)
library(coda)
library(colorspace)
library(microbenchmark)
library(cowplot)
library(colorspace)
library(ggdendro)
library(tidyverse)

source("subspace-functions.R")
source("fit-subspace.R")
source("leukemia_functions.R")

load("LeukemiaData/leukemia.RData")

## Look at PCA for pooled data to evaluate differences in means
YpooledDemeaned <- sweep(Ypooled, 2, colMeans(Ypooled), '-')
indices <- which(typeVec %in% c("TEL-AML1", "T-ALL", "MLL"))
plotarrays(YpooledDemeaned[indices, ], axis1=8, axis2=9, classvec=typeVec[indices],
           star=FALSE, arraycol=c("red", "blue", "dark green"))

## typeVec[typeVec=="TEL-AML1"] <- 1
## typeVec[typeVec=="T-ALL"] <- 2
## typeVec[typeVec=="MLL"] <- 3
## plotarrays(YpooledDemeaned[indices, ], axins2=3, classvec=typeVec[indices],
##            star=FALSE, arraycol=qualitative_hcl(5, palette="Dark 3")[1:3], ellipse=2, labelsize=2)

## indices <- which(typeVec %in% c("1", "2", "Hyperdip50"))
## plotarrays(YpooledDemeaned[indices, ], classvec=typeVec[indices],
##            star=FALSE, labelsize=1.5, arraycol=c("red", "blue", "dark green"))

## plotarrays(YpooledDemeaned, classvec=typeVec,
##            star=FALSE)

## number of eigenvectors if pooled space

S <- 45
## number of eigenvectors in group subspace
R <- 10
P <- ncol(Ypooled)

residualList <- Slist <- list()
for( k in 1:ngroups ) {
    residual <- lm(Ylist[[k]] ~ 1)$residuals
    residual <- t(lm(t(residual) ~ 1)$residuals)
    residualList[[k]] <- residual
    Slist[[k]] <- t(residual) %*% residual
}
pooledResiduals <- do.call(rbind, residualList)

getRank(pooledResiduals)

plot(svd(residualList[[ngroups]]/sqrt(nvec[ngroups]))$d^2, pch=19, cex=0.75, ylim=c(0, 50))
for(k in 1:6) {
    points(svd(residualList[[k]]/sqrt(nvec[k]))$d^2, col=k+1, pch=19, cex=0.75)
}
svd(residualList[[ngroups]]/sqrt(nvec[ngroups]))$d^2

Vinit <- svd(do.call(cbind, lapply(1:ngroups, function(k) svd(t(residualList[[k]]))$u[, 1:min(S, nvec[k])])))$u[, 1:S]

isoVar <- sapply(1:ngroups, function(i) median(apply(residualList[[i]], 2, var)))
weightsList <- sapply(1:ngroups, function(i) {
    evals <- svd(residualList[[i]]/sqrt(nvec[i]))$d^2
    dp <- ifelse(evals/isoVar[i] <= (1+sqrt(P/nvec[i])), 0, sqrt((1-P/nvec[i]/(evals/isoVar[i]-1)^2) / (1+P/nvec[i]/(evals/isoVar[i]-1))))
    weights <- 1/(1-dp) - 1
    weights[1:min(nvec[i], S)]
})

## Weighted initialization
Vinit <- svd(do.call(cbind, lapply(1:ngroups, function(k) {
    svdk <- svd(t(residualList[[k]]))$u[, 1:min(nvec[k], S)] %*%  diag(weightsList[[k]][1:min(nvec[k], S)])
    })))$u[, 1:S]

R <- 10
Vinit <- Vinit %*% rustiefel(S, S)
microbenchmark(EMFit <- subspaceEM(residualList, P=P, S=S, R=R, nvec=nvec, Vstart=Vinit, verbose=TRUE, stiefelAlgo=2), times=1)
Vinit <- EMFit$V

save(EMFit, Vinit, S, R, P, nvec, residualList,
     file=sprintf("results/leukemiaEM-%s.RData", format(Sys.Date(), "%m-%d")))

evalRatiosQuad <- compute_variance_explained(cbind(Vstar, V2), residualList, nvec, 1/EMFit$PrecVec)


## GF PLOT
pdf(sprintf("paper/Figs/leukemia_cumulative_ratio-%s.pdf", format(Sys.Date(), "%m-%d")))
rownames(evalRatiosQuad) <- unique(typeVec)
colnames(evalRatiosQuad) <- 1:S
tib <- as_tibble(evalRatiosQuad)
tib$Type <- unique(typeVec)
tib %>% gather(key=S, value=GF, -Type) %>%
    ggplot(aes(x=as.numeric(S), y=GF, col=Type)) + geom_line(size=1.5) + scale_color_discrete_qualitative(alpha=0.75) +
    xlab("Dimension") + ylab("Goodness of Fit") +
    geom_hline(yintercept=1, linetype="dashed")
dev.off()

## Bar plot
pdf(sprintf("paper/Figs/leukemiaRatio-%s.pdf", format(Sys.Date(), "%m-%d")), font="Times")
par(mar=c(6.1, 4.1, 4.1, 2.1))
barplot(evalRatiosQuad[, S], ylim=c(0, 1), xlab="", main="", space=0.2, cex.axis=1.5, col="#646464", col.axis="black", col.main="black", border=NA, names.arg="")
abline(h=1, lty="dashed", lwd=2)
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
samples <- fitSubspace(P, S, R, Q=S-R, residualList,
                   nvec, ngroups, init=initSS,
                   niters=10000, nskip=10)

save(samples, initSS, file=sprintf("results/leukemiaBayes-%s.RData", format(Sys.Date(), "%m-%d")))

#########################
## Posterior plots
#########################

regionColors <- c(4, 2:3, 1, 5:7)

g1 <- 7
g2 <- 2
pmPsi6 <- getPostMeanPsi(samples$Osamps[, , g1 , ],
                            samples$omegaSamps[, g1, ],
                         samples$s2samps[g1, ], 100)
pmPsi7 <- getPostMeanPsi(samples$Osamps[, , g2 , ],
                            samples$omegaSamps[, g2, ],
                         samples$s2samps[g2, ], 100)

O <- svd(pmPsi6[1:R, 1:R] - pmPsi7[1:R, 1:R])$u[, 3:4]


## view data projected onto a 2D subspace
unique(typeVec)
g <- 1


O <- samples$Osamps[1:R, 1:R, g, 1000]
O2 <- samples$Osamps[(R+1):S, (R+1):S, 7, 1000]

Vstar <- Vinit[, 1:R] %*% O
V2 <- Vinit[, (R+1):S] %*% O2
Osamps_proj <- array(dim = c(2, 2, ngroups, 1000))
omegaSamps_proj <- array(dim = c(2, ngroups, 1000))
view_indices <- c(1, 2)
for(i in 100:1000) {
    for(k in 1:ngroups) {
        ## Osamps2[1:R, 1:R, k, i ] <- t(O) %*% samples$Osamps[1:R, 1:R, k, i]

        tmp <- t(Vstar[, view_indices]) %*% Vstar %*% t(O) %*% samples$Osamps[1:R, 1:R, k, i]
        omega <- samples$omegaSamps[1:R, k, i]
        eig <- eigen(tmp %*% diag(omega / (1-omega)) %*% t(tmp))
        Osamps_proj[, , k, i] <- eig$vectors
        lambda <- eig$values
        omegaSamps_proj[, k, i] <- lambda/(lambda+1)
    }
}


pdf(sprintf("paper/Figs/leukemiaPosterior-%-i-%s.pdf", g, format(Sys.Date(), "%m-%d")), font="Times")
with(samples,
     posteriorPlot(Osamps_proj[, , 1:7, 901:1000], omegaSamps_proj[, 1:7, 901:1000],
                s2samps[1:7, 901:1000], nsamps=100, groups=1:ngroups,
                probRegion=0.95, hline=NULL, type="mag",
                plotPoints=TRUE, col=regionColors[1:7], ymax=500, cex.pts=1))

unique(typeVec)
legend("topright", legend=unique(typeVec), col=regionColors, pch=19, bty="n", title="Leukemia type", cex=1.2)

omegaSamps_proj[, , 1000]/(1-omegaSamps_proj[, , 1000])

dev.off()


samples$omegaSamps[, , 1000]/(1-samples$omegaSamps[, , 1000])
samples$omegaSamps[, , 1000]/(1-samples$omegaSamps[, , 1000])

pdf(sprintf("paper/Figs/leukemia-biplot-%i-%s.pdf", g, format(Sys.Date(), "%m-%d")), font="Times")

grps <- c(1, 4, 6)
## grps <- c(1, 2, 3, 7)

symbol_names <- AnnotationDbi::select(hgu95av2.db, rownames(Vstar), c("SYMBOL"))

symbol_names %>% dplyr::filter(!duplicated(PROBEID))

mag <- apply(Vstar[, 1, drop=FALSE], 1, function(x) sqrt(sum(x^2)))
pos_mag <- apply(Vstar[, 1, drop=FALSE], 1, function(x) sqrt(sum(x^2)))

pos_indices <- order(Vstar[, 1], decreasing=TRUE)[1:10]
neg_indices <- order(Vstar[, 1], decreasing=FALSE)[1:10]

Vsub <- Vstar[c(pos_indices, neg_indices), ]

rownames(Vsub)

print(
    xtable(
        tibble(letters1 = toupper(letters[1:10]), symbols1 = symbol_names[pos_indices, "SYMBOL"],
               letters2= toupper(letters[11:20]), symobls2 = symbol_names[neg_indices, "SYMBOL"]),
    na="", quote=FALSE),
    include.rownames=FALSE, include.colnames=TRUE,
    floating=FALSE, sanitize.text.function=identity,
    file=sprintf("paper/Figs/genes-%i.txt", g))

plot(Vstar[, 1:2], col="light grey", pch=19, cex=0.5,
     xlim=c(-1.1, 1.1)*max(abs(Vstar[, 1:2])), ylim=c(-1.1, 1.1)*max(abs(Vstar[, 1:2])), xlab="V1", ylab="V2")

text(Vstar[pos_indices, 1], Vstar[pos_indices, 2], labels=toupper(letters)[1:10], col="black")
text(Vstar[neg_indices, 1], Vstar[neg_indices, 2], labels=toupper(letters)[11:20], col="black")


abline(h=0, v=0)
legend("topleft", legend=unique(typeVec)[grps], lty=1, lwd=2, col=regionColors[grps],
       cex=1.3, bty='o', ncol=1, title="Leukemia Type", bg="white", box.col="white")


for(k in grps) {

    pmPsi <- getPostMeanPsi(Osamps_proj[, , k , ],
                            omegaSamps_proj[, k, ],
                            samples$s2samps[k, ], 100)

    eigK <- eigen(pmPsi)
    lambda <- eigK$values
    print(lambda)
    evecs <- eigK$vectors

    maxIndex <- which.max(lambda)
    lamRatio <- lambda[maxIndex]/lambda[-maxIndex]
    angle <- atan(evecs[2, maxIndex]/evecs[1, maxIndex])
    print(angle)

    ellipse(center=c(0,0), shape=evecs[, 1:2] %*% diag(lambda[1:2]) %*% t(evecs[, 1:2])/25, radius=0.01, col=regionColors[k], center.cex=000)
    ellipse(center=c(0,0), shape=evecs[, 1:2] %*% diag(lambda[1:2]) %*% t(evecs[, 1:2])/25, radius=0.025, col=regionColors[k], center.cex=000)
    ellipse(center=c(0,0), shape=evecs[, 1:2] %*% diag(lambda[1:2]) %*% t(evecs[, 1:2])/25, radius=0.04, col=regionColors[k], center.cex=000)
    points(angle, lamRatio, col=k, pch=19, cex=2)

}


dev.off()

#############################################################
################### GO ANALYSIS #######################
#############################################################


symbol_names <- sapply(rownames(EMFit$V), function(x) {
    na.omit(select(hgu95av2.db, x, c("SYMBOL"))$SYMBOL)[1]
})


v1 <- Vstar[, view_indices[1]]
v2 <- Vstar[, view_indices[2]]
mag <- rowSums(Vstar[, view_indices]^2)
names(mag) <- rownames(EMFit$V)

gp_v1 <- find_go_groups(abs(v1), 0.05, go_name_length=60)
gp_v1

print(xtable(gp_v1[, c("Name", "Q-value", "Numer of Genes")], na="", quote=FALSE),
      include.rownames=TRUE, include.colnames=TRUE,
      floating=FALSE, sanitize.text.function=identity,
      file=sprintf("paper/Figs/go-%i.txt", g))




gp_v2 <- find_go_groups(abs(v2), 0.1, go_name_length=60)
gp_v2

gp_mag <- find_go_groups(mag, 0.1, go_name_length=60)
xtable(gp_mag)

mag_all <- rowSums(Vstar^2)
names(mag_all) <- rownames(EMFit$V)
gp_all <- find_go_groups(abs(mag_all), 0.05, go_name_length=60)
gp_all



print(xtable(gp_all[, c("Name", "Q-value", "Numer of Genes")], na="", quote=FALSE),
      include.rownames=TRUE, include.colnames=TRUE,
      floating=FALSE, sanitize.text.function=identity,
      file=sprintf("paper/Figs/go-all.txt"))


#############################################################
################ Check Covariance Similarity #################
#############################################################

similarity_matrix <- matrix(0, nrow=ngroups, ncol=ngroups)
for(s in 500:1000){
    for(i in 1:ngroups) {
        for(j in 1:ngroups) {
            
            O1 <- samples$Osamps[1:R, 1:R, i, s]
            O2 <- samples$Osamps[1:R, 1:R, j, s]

            omega1 <-  samples$omegaSamps[1:R, i, s]
            omega2 <-  samples$omegaSamps[1:R, j, s]

            l1 <- omega1/(1-omega1)
            l2 <- omega2/(1-omega2)

            sig1 <- O1 %*% diag(l1) %*% t(O1)
            sig2 <- O2 %*% diag(l2) %*% t(O2)

            d <- sig1 - sig2

            similarity_matrix[i, j] <- similarity_matrix[i, j] + tr(t(d) %*% d)
        }
    }
}

as.dist(similarity_matrix)
sort(similarity_matrix[1, ])
sort(similarity_matrix[2, ])
sort(similarity_matrix[7, ])

diag(similarity_matrix) <- NA
library(superheat)
rownames(similarity_matrix) <- colnames(similarity_matrix) <- unique(typeVec)


pdf("paper/Figs/leukemia_dendro.pdf")
hc <- hclust(as.dist(similarity_matrix), "complete")
ggdendrogram(hc, aes(size=100))
dev.off()

pdf("paper/Figs/leukemia_heat.pdf")
ord <- hc$order
superheat(similarity_matrix[ord, ord]/500, heat.pal=sequential_hcl(10, palette="Blue"), heat.na.col="black")
dev.off()

#############################################################
################ Check Predictions  ########################
#############################################################

Y <- pooledResiduals
preds <- makePrediction(Y, EMFit$V, samples$Osamps, samples$omegaSamps, samples$s2samps,
               nvec=nvec, ngroups=length(nvec), nsamps=1000,
               numToAvg=500)

## For how many 
sapply(unique(typeVec), function(type) {
       mean((as.integer((as.factor(typeVec))) == apply(preds_med, 1, which.max))[typeVec==type])})
mean((as.integer((as.factor(typeVec))) == apply(preds_med, 1, which.max)))



sort(preds_med[1, ], index.return=TRUE, decreasing=TRUE)$ix



order(preds_med[1, ], decreasing=FALSE)
preds_med[1, ]

apply(preds_med, 1, which.max)
sort(preds_med,partial=6)[6]



preds_med <- apply(preds, c(1, 3), median)
colnames(preds) <- unique(typeVec)
tib <- as_tibble(preds)
tib$sample <- 1:nrow(Y)

tv <- unique(typeVec)
tv[1] <- "BCR"
tv[2] <- "E2A"
labelDat <- tibble(Group = tv, x = c(0, cumsum(nvec))[1:7] + nvec/2, y=rep(1.05, 7))

pdf("paper/Figs/membership_probs.pdf", width=12, height=4)
tib %>% gather(key=group, value=prob, -sample) %>%
    ggplot() + geom_bar(aes(x=sample, y=prob, fill=group), stat="identity",width=1) +
    geom_vline(xintercept=c(0, cumsum(nvec))+1/2) +
    scale_fill_discrete_qualitative() +
    xlab("Samples") + ylab("Estimated Probabilities") +
    geom_text(data=labelDat, aes(x=x, y=y, label=Group), angle=0, size=5) +
    theme(legend.position="top", legend.title=element_blank())
dev.off()


apply(apply(preds, c(1, 3), median), 1, which.max) %>% table

plot(residualList[[7]] %*% Vstar[, 1:2], ylim=c(-15, 15), xlim=c(-15,15))
points(Ylist[[1]] %*% Vstar[, 1:2], ylim=c(-15, 15), xlim=c(-15,15), col="blue")

dim(Ylist[[1]])
unique(typeVec)
nvec

library(Rtsne)
rtsne_pooled <- Rtsne(pooledResiduals)

plot(rtsne_pooled$Y, col=as.factor(typeVec), pch=19)

Z <- pooledResiduals %*% Vinit[, 1:10]

res <- sapply(1:7, function(g) {
    psi <- getPostMeanPsi(samples$Osamps[, , g, ], samples$omegaSamps[, g, ], samples$s2samps[, g], 200)
    psi <- mean(samples$s2samps[, g]) + psi
    dmat <- sqrt(mahalanobis(Z, center=FALSE, cov=psi[1:10, 1:10]))
    plot(dmat, col=as.factor(typeVec), ylim=c(0, 25))
    dmat
})


dist_clust <- Rtsne(scale(res, center=FALSE), perplexity=50, dims=3)

plot(dist_clust$Y[, c(1, 2)], col=as.factor(typeVec), pch=19)

plot3D::scatter3D(x=dist_clust$Y[,1], y=dist_clust$Y[,2], z=dist_clust$Y[,3],
                  bty = "g", pch = 18, 
                  col.var = as.integer(typeVec), 
                  col = 1:7,
                  pch = 18)







                  col.var=as.factor(typeVec), col=1:7)












rownames(Vstar)

Vstar2 <- Vstar
rownames(Vstar2) <- symbol_names[!duplicated(symbol_names[, 1]), 2]
head(Vstar2)
Vstar["STAT1", ]
