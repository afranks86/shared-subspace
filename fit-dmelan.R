rm(list=ls())

library(tidyverse)
library(rstiefel)
library(colorspace)
source("subspace-functions.R")
source("fit-subspace.R")
##load("ae_XY_m527.RData")
load("results/C18_all.RData")


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

group <- apply(XDM, 1, function(x) which(x==1))

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

S <- getRank(fit$residuals)
R <- S
P <- ncol(Y)

Vinit <- svd(do.call(cbind, lapply(1:ngroups, function(k) svd(t(residualList[[k]]))$u[, 1:S])))$u[, 1:S]
PrecVec = rep(5, ngroups)
EMFit <- subspaceEM(residualList, P=P, S=S, R=R, nvec=nvec, Vstart=Vinit, verbose=TRUE, M_iters=1000)
Vinit <- EMFit$V

evalRatiosQuad <- compute_variance_explained(Vinit, residualList, nvec, 1/EMFit$PrecVec)
rownames(evalRatiosQuad) <- colnames(XDM)
colnames(evalRatiosQuad) <- 1:S
tib <- as_tibble(evalRatiosQuad)
tib$Type <- colnames(XDM)
tib %>% gather(key=S, value=GF, -Type) %>%
    ggplot(aes(x=as.numeric(S), y=GF, col=Type)) + geom_line(size=1.5) + scale_color_discrete_qualitative(alpha=0.75) +
    xlab("Dimension") + ylab("Goodness of Fit") +
    geom_hline(yintercept=1, linetype="dashed")

save(EMFit, Vinit, S, R, P, nvec, residualList,
     file=sprintf("dmelanEM-%s.RData", format(Sys.Date(), "%m-%d")))

pdf(sprintf("paper/Figs/dmelanRatio-%s.pdf", format(Sys.Date(), "%m-%d")), font="Times")
par(mar=c(6.1, 4.1, 4.1, 2.1))
barplot(evalRatiosQuad[, S], ylim=c(0, 1), xlab="", main="", space=0.2, cex.axis=1.5, col="#646464", col.axis="black", col.main="black", border=NA)
text(seq(1, 7, length.out=ngroups) , par("usr")[1], labels=paste(rep(c("Young", "Middle", "Old"), 2), rep(c("F", "M"), each=ngroups/2), sep="/"),
     srt=45, xpd=TRUE, cex=1.5, col="black", adj=c(1, NA))
dev.off()



OmegaList <- Ulist <- list()
for(k in 1:ngroups) {

    eigK <- eigen( Matrix::nearPD(diag(S) - EMFit$PhiList[[k]] / EMFit$PrecVec[k])$mat)

    OmegaList[[k]] <- eigK$values
    Ulist[[k]] <- Vinit %*% eigK$vectors

}
s2vec <- 1/EMFit$PrecVec

initSS <- list(V=Vinit, Ulist=Ulist, OmegaList=OmegaList, s2vec=s2vec)
samples <- fitSubspace(P, S, R, Q=S-R, residualList,
                   nvec, ngroups, init=initSS,
                   niters=10000, nskip=10)



save(samples, initSS, file=sprintf("dmelanBayes-%s.RData",
                                   format(Sys.Date(), "%m-%d")))

regionColors <- rep(sequential_hcl(palette="Heat", 3), 2)
pch <- rep(c(19, 15), each=3)

g1 <- 1
g2 <- 3

pmList <- lapply(1:ngroups, function(g) {
    getPostMeanPsi(samples$Osamps[, , g , ],
                   samples$omegaSamps[, g, ],
                   samples$s2samps[g, ], 100)})

view_indices <- c(5, 6)

## age difference
O <- svd(pmList[[1]][1:R, 1:R] + pmList[[4]] -
         pmList[[3]][1:R, 1:R] - pmList[[6]])$u[, view_indices]

## gender difference
O <- svd(pmList[[1]][1:R, 1:R] + pmList[[2]] + pmList[[3]] -
         pmList[[4]][1:R, 1:R] - pmList[[5]] - pmList[[6]])$u[, view_indices]


Vstar <- Vinit[, 1:R] %*% O
Osamps_proj <- array(dim = c(2, 2, ngroups, 1000))
omegaSamps_proj <- array(dim = c(2, ngroups, 1000))

for(i in 1:1000) {
    for(k in 1:ngroups) {
        ## Osamps2[1:R, 1:R, k, i ] <- t(O) %*% samples$Osamps[1:R, 1:R, k, i]

        tmp <- t(O) %*% samples$Osamps[1:R, 1:R, k, i]
        omega <- samples$omegaSamps[1:R, k, i]
        eig <- eigen(tmp %*% diag(omega / (1-omega)) %*% t(tmp))
        Osamps_proj[, , k, i] <- eig$vectors
        lambda <- eig$values
        omegaSamps_proj[, k, i] <- lambda/(lambda+1)
    }
}


pdf(sprintf("paper/Figs/dmelanPosterior-%-i-%i-%s.pdf", g1, g2, format(Sys.Date(), "%m-%d")), font="Times")
toshow <- 1:6
with(samples,
     posteriorPlot(Osamps_proj[, , , 501:1000], omegaSamps_proj[, , 501:1000],
                s2samps[, 501:1000], nsamps=100, groups=toshow,
                probRegion=0.95, hline=NULL, type="mag",
                plotPoints=TRUE, col=regionColors, ymax=NULL, cex.pts=1, pch=pch))

legend("topleft", legend=colnames(XDM), col=regionColors, pch=pch, title="Leukemia type", cex=1.2, bty="o", box.col="white")

omegaSamps_proj[, , 1000]/(1-omegaSamps_proj[, , 1000])


    

dev.off()

preds <- makePrediction(fit$residuals, EMFit$V,
                        samples$Osamps,
                        samples$omegaSamps,
                        samples$s2samps,
               nvec=nvec, ngroups=length(nvec), nsamps=1000,
               numToAvg=500)

## For how many 
sapply(1:6, function(type) {
       mean(type == apply(preds, 1, which.max)[group == type])})
mean(group == apply(preds, 1, which.max))


t(apply(preds, 1, order, decreasing=TRUE))
sort(preds,partial=6)[6]



