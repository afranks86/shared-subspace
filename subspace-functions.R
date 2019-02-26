library(rstiefel)
library(mvtnorm)
library(Matrix)
library(magrittr)
library(scales)


## stanstart <- stan(file="vectorBMF.stan", data=list(M=9,lam=rep(1,10), gamma=rep(1, 10)), chains=1, iter=1)

tr <- function(X) { sum(diag(X)) }

## Compute frobenius norms
## H = Hierachical eigen pooling
## e.g. compare trace (U_k^T U_J)^2
## SS = Shares space
## e.g. compare || U_k^T U_J || (frobenius)

createNormMatrix <- function(eigenlist, vindices, type="H") {
    ngroups <- length(eigenlist)
    normsMat <- matrix(NA, nrow=ngroups, ncol=ngroups)

    for( i in 1:ngroups ) {
        for( j in i:ngroups ){
            Ui <- eigenlist[[i]]$vectors[, vindices]
            Uj <- eigenlist[[j]]$vectors[, vindices]
            if(type == "H") {
                normsMat[i, j] <- sum((t(Ui)%*%Uj)^2)/sqrt(ncol(Ui)) ## these are the same
            } else if (type == "SS") {
                
                normsMat[i, j] <- norm(t(Ui)%*%Uj, type="F")/sqrt(ncol(Ui))
            }
        }
    }
    diag(normsMat) <- NA
    normsMat[lower.tri(normsMat)] <- t(normsMat)[lower.tri(normsMat)]

    return(normsMat)
}

posteriorPlot <- function(Osamps, OmegaSamps, s2samps, nsamps, groups,
                          probRegion=0.95, hline=NULL, col=NULL,
                          pch=NULL, lty=NULL, ymax=30, type = "mag", 
                          plotPoints=TRUE, polar=FALSE, cex.axis=1.5, cex.pts=0.5,
                          splitGroups=c(), splitPts=rep(0, length(splitPts))) {

    ngroups <- length(groups)
    
    if(is.null(col)){
        col <- 1:ngroups
    }
    if(is.null(pch)){
        pch=rep(19, ngroups)
    }
    if(is.null(lty)){
        lty=rep(1, ngroups)
    }
    par(mar=c(5.1, 5.1, 4.1, 2.1))

    if(type=="mag") {
        ylab <- expression(lambda[1])
    }
    else if(type=="logmag") {
        ylab <- expression("log"[2]~"(" ~ lambda[1] ~ ")")
    } else if (type =="logratio") {
        ylab <- expression("log"[2]~"(" ~ lambda[1]/lambda[2] ~ ")")
    } else {
        ylab <- expression("(" ~ lambda[1]/lambda[2] ~ ")")
    }

    if(is.null(ymax))
        ymax <- 1.1*max(OmegaSamps/(1-OmegaSamps))
    
    plot(0, 0, xlim=c(-pi/2, pi/2),
         ylim=c(0, ymax), cex=0, xlab=expression("angle, acos("~U[1]^T*V[1]~")"), ylab=ylab, xaxt="n", cex.axis=cex.axis, cex.lab=1.5)
    axis(1, at=seq(-pi/2, pi/2, by=pi/4), labels=expression(-pi/2, -pi/4, 0, pi/4, pi/2), cex.axis=cex.axis, cex.lab=1.5)

    for(g in setdiff(groups, splitGroups)) {

        pmPsi <- getPostMeanPsi(Osamps[, , g , ], OmegaSamps[, g, ],
                                s2samps[g, ], nsamps)

        eigPsi <- eigen(pmPsi)
        pmValues <- eigPsi$values
        pmVectors <- eigPsi$vectors
        maxIndex <- which.max(pmValues)

        pmPoint <- pmValues[maxIndex]*pmVectors[, maxIndex]
        if(pmPoint[1] < 0)
            pmPoint <- pmPoint * c(-1, -1)

        hp <- getHullPoints(nsamps, OmegaSamps[, g, ], Osamps[, , g, ],
                            type=type, probRegion=probRegion)
        pts <- hp$pts
        hullPoints <- hp$hullPoints

        if(plotPoints) {
            points(pts[1, ], pts[2, ], col=alpha(col[g], 1/2), pch=pch[g], cex=cex.pts)
        } else {
            polygon(pts[1, hullPoints], pts[2, hullPoints], lwd=3, border=col[g], col=alpha(col[g], 1/4), lty=lty[g])
        }
    }

    if(!is.null(hline))
        abline(h=hline, lty=2)

    for(g in splitGroups) {
        browser()
        allPts <- getHullPoints(nsamps, OmegaSamps[1:2, g, ], Osamps[1:2, 1:2, g, ], type=type)$allPts
        split1 <- which(allPts[1, ] < 0)
        split2 <- which(allPts[1, ] >= 0)



        allPts[2, split1]
        
        hull1 <- chull(group1pts[1, split1], group1pts[2, split1])
        numPtsToRemove <- round(nsamps*(1-probRegion))
        while(numPtsToRemove > 0) {
            hullPoints1 <- chull(pts[1, ], pts[2, ])
            if(length(hullPoints) > numPtsToRemove) {
                hullPoints <- sample(hullPoints, numPtsToRemove)
                pts <- pts[, -hullPoints]
                numPtsToRemove <- 0
            } else{
                pts <- pts[, -hullPoints]
                numPtsToRemove <- numPtsToRemove - length(hullPoints)
            }
        }
        
        hull2 <- chull(group1pts[1, split2], group1pts[2, split2])
        pts1 <- group1pts[, split1[hull1]]
        pts1[2, 4] <- min(pts1[2, ])
        pts2 <- group1pts[, split2[hull2]]
        pts2 <- pts2[, 3:ncol(pts2)]
        pts2 <- cbind(c(pi/2, min(pts1[2, ])), pts2)
        polygon(pts1[1, ], pts1[2, ], lwd=0.01, border="black", col=alpha(regionColors[g], 1/4), lty=1)
        polygon(pts2[1, ], pts2[2, ], lwd=0.01, border="black", col=alpha(regionColors[g], 1/4), lty=1)
        lines(pts1[1, c(5:9, 1:4)], pts1[2, c(5:9, 1:4)], col=regionColors[g], lwd=3)
        lines(pts2[1, 1:10], pts2[2, 1:10], col=regionColors[g], lwd=3)
    }

    
}


eigenvalueDists <- function(OmegaSamps, nsamps, groups,
                            probRegion=0.95, hline=NULL, col=NULL,
                            pch=NULL, lty=NULL, ymax=30, logRatio=FALSE,
                            plotPoints=TRUE, polar=FALSE, cex.axis=1.5) {

    ngroups <- length(groups)
    
    if(is.null(col)){
        col <- 1:ngroups
    }
    if(is.null(pch)){
        pch=rep(19, ngroups)
    }
    if(is.null(lty)){
        lty=rep(1, ngroups)
    }
    par(mar=c(5.1, 5.1, 4.1, 2.1))

    plot(0, 0, xlim=c(-pi/2, pi/2),
         ylim=c(0, ymax), cex=0, xlab=expression("angle, acos("~U[1]^T*V[1]~")"), ylab=ylab, xaxt="n", cex.axis=cex.axis, cex.lab=1.5)
    axis(1, at=seq(-pi/2, pi/2, by=pi/4), labels=expression(-pi/2, -pi/4, 0, pi/4, pi/2), cex.axis=cex.axis, cex.lab=1.5)

    LambdaSamps <- OmegaSamps[1, , ] / (1 - OmegaSamps[1, , ])

    tib <- as_tibble(LambdaSamps)
    tib
    tib$type <- types
    tib %>% gather(key=Sample, value=Lambda, -type) %>% mutate(Lambda = as.numeric(Lambda)) %>% ggplot(aes(type, Lambda)) + geom_violin(aes(fill=type)) + ylim(c(0, 300))

}

posteriorPlotPolar <- function(Osamps, OmegaSamps, s2samps, nsamps, groups,
                          probRegion=0.95, hline=NULL, col=NULL,
                          pch=NULL, lty=NULL, ymax=NULL, logRatio=FALSE,
                          plotPoints=TRUE, polar=FALSE, cex.axis=1.5) {

    ngroups <- length(groups)
    
    if(is.null(col)){
        col <- 1:ngroups
    }
    if(is.null(pch)){
        pch=rep(19, ngroups)
    }
    if(is.null(lty)){
        lty=rep(1, ngroups)
    }

    radius_vec <- angle_vec <- type_vec <- c()
                              
    for(g in groups) {
        
        pmPsi <- getPostMeanPsi(Osamps[, , g , ], OmegaSamps[, g, ],
                                s2samps[g, ], nsamps)

        eigPsi <- eigen(pmPsi)
        pmValues <- eigPsi$values
        pmVectors <- eigPsi$vectors
        maxIndex <- which.max(pmValues)

        pmPoint <- pmValues[maxIndex]*pmVectors[, maxIndex]
        if(pmPoint[1] < 0)
            pmPoint <- pmPoint * c(-1, -1)

        hp <- getHullPoints(nsamps, OmegaSamps[, g, ], Osamps[, , g, ],
                            logRatio=logRatio)
        pts <- hp$pts

        radius_vec <- c(radius_vec, pts[2, ], pts[2, ])
        angle_vec <- c(angle_vec, pts[1, ], pts[1, ]-pi)
        type_vec <- c(type_vec, rep(g, 2*length(pts[1,])))
        
    }

    tibble(radius=radius_vec, angle=angle_vec, type=type_vec) %>% ggplot(aes(y=radius, x=angle, col=as.factor(type))) +
        geom_point() + coord_polar(theta="x", clip="off") + ylim(0, 2) + xlim(c(-pi, pi)) + theme_bw() + scale_color_discrete_qualitative(alpha=0.3)

}


## type is "mag", "ratio",
getHullPoints <- function(nsamps, OmegaSamps, Osamps, type="mag",
                          probRegion=0.95) {

    PointsList <- lapply(1:nsamps, function(i) {
        LambdaSamp <- OmegaSamps[, i]/(1-OmegaSamps[, i])
        maxIndex <- which.max(LambdaSamp)

        if(type == "mag")
            yval <- LambdaSamp[maxIndex]
        else{
            yval <- LambdaSamp[maxIndex]/LambdaSamp[-maxIndex]
        }
        
        O1 <- Osamps[, maxIndex, i]
        angle <- atan(O1[2]/O1[1])
        
        c(angle, yval)
        
    })

    pts <- simplify2array(PointsList)
    allPts <- pts
    if(type == "logratio") {
        allPts[2, ] <- log2(allPts[2, ])
        pts[2, ] <- log2(pts[2, ])
    }
    
    numPtsToRemove <- round(nsamps*(1-probRegion))
    while(numPtsToRemove > 0) {
        hullPoints <- chull(pts[1, ], pts[2, ])
        if(length(hullPoints) > numPtsToRemove) {
            hullPoints <- sample(hullPoints, numPtsToRemove)
            pts <- pts[, -hullPoints]
            numPtsToRemove <- 0
        } else{
            pts <- pts[, -hullPoints]
            numPtsToRemove <- numPtsToRemove - length(hullPoints)
        }
    }

    
    
    hullPoints <- chull(pts[1, ], pts[2, ])

    list(allPts=allPts, pts=pts, hullPoints=hullPoints)
    
}


getHullPointsProbabilityRegion <- function(pts, probRegion) {
    
}




compute_variance_explained <- function(V, Ylist, nvec, s2vec) {

    P <- nrow(V)
    S <- ncol(V)
    ngroups <- length(Ylist)
    
    ## Use qudratic fomr to compute goodness of fit
    evalRatiosQuad <- sapply(1:S, function(M) {
        sapply(1:ngroups, function(k) {

            YVt <- t(V[, 1:M]) %*% t(Ylist[[k]])
            numer <- tr(YVt %*% t(YVt))/nvec[k]
            
            evals <- svd(Ylist[[k]])$d[1:min(M, nvec[k])]^2/nvec[k]
            b <- (s2vec[k]*P/nvec[k] - evals - 1)

            quadSol <- ifelse(b^2 - 4*evals < 0, (evals - (s2vec[k]*P/nvec[k])), suppressWarnings((-b + sqrt(b^2 - 4*evals))/2))
            denom <- sum(quadSol)

            (numer/denom)
        })
    })

    evalRatiosQuad
}






##############################################################
### Functions for Sampling
#############################################################

## Sample from 1/sigma^2

sampleSigma2 <- function(Y, YV, O, omega, n, nu0=1, s20=1) {

    YU <- YV %*% O
    
    p <- ncol(Y)
    a <- ( nu0 + n*p)/2

    b <- ( nu0*s20 + sum(Y^2) - tr(t(YU) %*% YU %*% diag(omega, nrow=length(omega))))/2
    
    1/rgamma(1, a, b)
}

## trucnated gamma
tgamma <- function(a, b, scale) {

    logM <- -1 * pgamma(scale, a, b, log.p=TRUE)
    qgamma(exp(log(runif(length(a))) - logM), a, b)

}

sampleOmega <- function(YV, O, s2, n) {

    R <- ncol(O)
    
    Rtilde <- min(n, R)
    O <- O[, 1:Rtilde]
    
    YVO <- YV %*% O

    cvec <- diag(t(YVO) %*% YVO) /(2*s2)

    ## add divide by n for stability?
    g <- tgamma(rep(n/2 + 1, length(cvec)), cvec/n, n)
    g <- g/n

    omega <- 1 - g
    
    if(n < R)
        omega <- c(omega, rep(0, R - n))

    omega

}

## Shared subspace Sampling
## phi is the prior concentration
sampleO <- function(YV, O, s2, omega) {

    A <- t(YV) %*% YV / (2*s2)
    B <- diag(omega, nrow=length(omega))

    ord <- order(omega, decreasing=TRUE)
    revOrd <- order(ord)

    Btilde <- B[ord, ord]
    Otilde <- O[, ord]

    O <- rbing.matrix.gibbs(A, Btilde, Otilde)

    O[, revOrd]
}

## Draw the shared orthonormal matrix V given other terms
## Ok is an S x R matrix,  V is a P x S and U is P x R
## B_k = Ok*omega_k*Ok^T
## A_k is sample covariance matrix for group k

sampleV <- function(Slist, Ulist, s2vec, OmegaList, V, method="gibbs") {

    K <- length(Ulist)
    S <- ncol(V)
    P <- nrow(V)

    ## First save the Bk matrices for each k
    BkList <- OkList <- list()
    for( k in 1:K ){
         OkList[[k]] <- t(V) %*% Ulist[[k]]
         omegaK <- OmegaList[[k]]
         BkList[[k]] <- OkList[[1]] %*% (diag(omegaK, nrow=length(omegaK)) /
                                         (2*s2vec[k])) %*% t(OkList[[1]])
    }

    if(K == 1 & method=="bing") {

        A <- Slist[[1]]
        B <- BkList[[1]]

        Beigen <- eigen(B)
        Bvals <- Beigen$values
        Bvecs <- Beigen$vectors
        
        ord <- order(Bvals, decreasing=TRUE)
        revOrd <- order(ord)

        Btilde <- diag(Bvals)
        Vtilde <- V %*% Bvecs[, ord]

        Vtilde <- rbing.matrix.gibbs(A, Btilde, Vtilde)
        V <- Vtilde %*% t(Bvecs[, ord])

    } else {
    
    ## Randomly sample a column,  i
    for( i in sample(1:S) ) {

      N <- NullC(V[, -i])

        ## Get the linear term
        ## for j neq i
        C <- rep(0, P)
        for( j in setdiff(1:S, i)) {
            Ak_bij <- matrix(0, P, P)
            for(k in 1:K) {
                Ak <- Slist[[k]]
                bk_ij <- BkList[[k]][i, j]
                Ak_bij <- Ak_bij + Ak*bk_ij
            }
            C <- C+t(V[, j]) %*% Ak_bij
        }
        C <- 2*C

        ## Get the quadratic term
        Ak_bii <- matrix(0, P, P)
        for(k in 1:K) {
            Ak <- Slist[[k]]
            bk_ii <- BkList[[k]][i, i]
            Ak_bii <- Ak_bii+Ak*bk_ii
        }
        A <- Ak_bii
        
        Ctilde <- as.vector(C %*% N)
        Atilde <- t(N) %*% A %*% N

        NV <- as.vector(t(N) %*% V[, i])

        if( method == "gibbs" ) {
            V[, i] <- N %*% R.rbmf.vector.gibbs(Atilde, Ctilde, NV)
        } else if( method == "hmc" ) {
            V[, i] <- N %*% R.rbmf.vector.hmc(Atilde, Ctilde, NV)
        } else if(method == "mh" ){
            V[, i] <- N %*% R.rbmf.vector.mh(Atilde, Ctilde, NV)
        } else {
            stop("Unknown method")
        }
    }}

    ## Propose swaps to handle multi-modality
    if(ncol(V) > 1) {

        if(ncol(V)==2) {
            Vswap <- V[, 2:1]
        } else {
            perm <- sample(2:ncol(V))
            onePos <- sample(2:ncol(V), size=1)
            perm <- c(perm[1:(onePos-1)], onePos, perm[onePos:length(perm)])
            Vswap <- V[, sample(1:ncol(V))]
        }

        logdens <- function(Vcur) {
            sum(
                sapply(1:length(Slist), function(k) {
                    Ok <- OkList[[k]]
                    omegaK <- OmegaList[[k]]
                    B <- Ok %*% (diag(omegaK, nrow=length(omegaK)) /
                                 (2*s2vec[k])) %*% t(Ok)

                    tr(Slist[[k]] %*% Vcur %*%
                       B %*% t(Vcur))
                })
            )
        }

        logMFswap <- logdens(Vswap)
        logMFcur <- logdens(V)
        log.mh <- logMFswap - logMFcur
        u <- runif(1,0,1)
        if ( log(u) < log.mh ) {
            V <- Vswap
            print("SWAP")
        }
    }
    
    V
}

proposeBinaryO <- function(S, U, V, Sig, s2, omega, n, flipProb=0.1) {
    
    Ocur <- round(t(V) %*% U, digits=10)
    
    if(ncol(Ocur)!=S)
        stop("R=S for binary sampling")
    
    ## change nflip entries of Ok
    cur <- apply(Ocur, 1, sum)
    
    if(!any(cur==1))
        stop("O is not binary")

    flip <- which(runif(S) < 0.1)
    cur[flip] <- 1-cur[flip]

    Oprop <- matrix(0, nrow=S, ncol=S)
    for( i in 1:S ) {
        Oprop[i, i] <- cur[i]
    }
    
    ## k = t(O)t(V)SkVO
    omega_cur <- sampleOmega(Sig, V %*% Ocur, s2, n)
    omega_prop <- sampleOmega(Sig, V %*% Oprop, s2, n)

    H <-  t(V) %*% Sig %*% V / (2*s2)
    
    Kprop <- diag(t(Oprop) %*% H %*% Oprop)
    Kcur <-  diag(t(Ocur) %*% H  %*% Ocur)

    trans1 <- sum(log(sapply(1:length(omega_prop),
                             function(i) debeta(omega_prop[i], 1, 1+n/2, Kprop[i]))))
    trans2 <- sum(log(sapply(1:length(omega_cur),
               function(i) debeta(omega_cur[i], 1, 1+n/2, Kcur[i]))))

    
    ll1 <- tr(omega %*% t(Oprop) %*% H %*% Oprop %*% omega) +
        n/2*sum(log(omega_prop))
    ll2 <- tr(omega %*% t(Ocur) %*% H %*% Ocur %*% omega) +
        n/2*sum(log(omega_cur))
    
    if(-rexp(1) < ll1+trans2-ll2-trans1 ) {
        Ocur <- Oprop
        omega_cur <- omega_prop
    } 

    list(O=Ocur, omega=omega_cur)

}

steinsLoss <- function(C1, C2inv) {

    sum(diag(C1 %*% C2inv)) - log(det(C1 %*% C2inv)) - nrow(C1)

}




getSigmaInv <- function(P, U, Omega, s2) {
    
    1/s2*(diag(P) - U %*% diag(Omega) %*% t(U))
    
}

getPostMeanSigmaInv <- function(P, USamps, OmegaSamps, s2vec, nsamps) {
    SigmaInvList <- lapply(1:nsamps, function(i) {
        1/s2vec[i]*(diag(P) - USamps[, , i] %*% diag(OmegaSamps[, i]) %*% t(USamps[, , i]))
    })
    apply(simplify2array(SigmaInvList), c(1, 2), mean)
}

## posterior mean of (t(V) x Sig x V)^(-1)
getPostMeanSigmaProjInv <- function(S, V, USamps, OmegaSamps, s2vec, nsamps) {
    SigmaProjInvList <- lapply(1:nsamps, function(i) {
        
        1/s2vec[i]*(diag(S) - (t(V) %*% USamps[, , i]) %*%
                    diag(OmegaSamps[, i]) %*%
                    t(t(V) %*% USamps[, , i]))
    })
    apply(simplify2array(SigmaProjInvList), c(1, 2), mean)
}


getPostMeanPsi <- function(Osamps, OmegaSamps, s2vec, nsamps) {
    PsiList <- lapply(1:nsamps, function(i) {
        s2vec[i]*(Osamps[, , i] %*% diag(OmegaSamps[, i]/(1-OmegaSamps[, i])) %*% t(Osamps[, , i]))
    })

    apply(simplify2array(PsiList), c(1, 2), function(x) mean(x, na.rm=TRUE))
}

getPostMeanSigma <- function(P, Ulist, OmegaList, s2vec) {

    SigmaList <- sapply(1:nsamps, function(i) {
        s2vec[i]*(Ulist[, , i] %*% diag((1-omega)/omega) %*% t(Ulist[, , i]) + diag(P))
    })
    apply(simplify2array(SigmaList), c(1, 2), mean)
}

## Get mean predictions for each group relative to truth
makePrediction <- function(Y, V, Osamps, omegaSamps, s2samps, nvec,
                           ngroups=length(nvec),
                           nsamps=dim(Osamps)[4],
                           numToAvg=nsamps/2) {
    if(ngroups==1) {
        return(matrix(1, ncol=1, nrow=1))
    }

    ##predictionMat <- matrix(0, nrow=length(genGroups), ncol=ngroups)
    ## Y <- Y / sum(Y^2)
    YV <- Y %*% V
    SigmaHat <- list()
        
    ll_array <- array(0, dim=c(nrow(Y), numToAvg, ngroups))
    for(i in (nsamps-numToAvg+1):nsamps) {
            
        for(k  in 1:ngroups) {
            Ok <- Osamps[, , k, i]
            om <- omegaSamps[, k, i]
            Lam <- diag(om/(1-om))
            s2 <- s2samps[k, i]
            SigmaHat[[k]] <-  s2*(Ok %*% Lam %*% t(Ok) + diag(nrow(Ok)))
        }

        ll <- computeLL(YV, SigmaHat)
        ll_array[, nsamps - i, ] <- ll
    }
    median_ll <- apply(ll_array, c(1, 3), median)
    dim(median_ll)

    priorWeights <- rep(1/length(SigmaHat), length(SigmaHat))
    probs <- apply(median_ll, 1, function(probVec) {
        normProbs <- probVec - max(probVec)
        normProbs <- exp(normProbs) * priorWeights /
            sum(exp(normProbs) * priorWeights)
        
        normProbs
    })
      
    t(probs)

}

computeLL <- function(Y, SigmaList, priorWeights=rep(1/length(SigmaList),
                                                     length(SigmaList))) {

    priorWeights <- priorWeights / sum(priorWeights)
    unnormalizedProbs <- sapply(SigmaList, function(Sigma) {
        dmvnorm(Y, sigma=Sigma, log=TRUE)
    })

    unnormalizedProbs
}


computeMembershipProbabilities <-  function(Y, SigmaList,
                                            priorWeights=rep(1/length(SigmaList),
                                                             length(SigmaList))) {

    priorWeights <- priorWeights / sum(priorWeights)
    unnormalizedProbs <- sapply(SigmaList, function(Sigma) {
        dmvnorm(Y, sigma=Sigma, log=TRUE)
    })

    if( !class(unnormalizedProbs) == "matrix")
        unnormalizedProbs <- matrix(unnormalizedProbs, nrow=1, ncol=length(SigmaList))
    
    probs <- apply(unnormalizedProbs, 1, function(probVec) {
        normProbs <- probVec - max(probVec)
        normProbs <- exp(normProbs) * priorWeights /
            sum(exp(normProbs) * priorWeights)
        
        normProbs
    })

    t(probs)
}


## Optimal threshold from Gavish, Donoho 2014
getRank <- function(Y) {

  svals <- svd(Y)$d

  m <- max(nrow(Y), ncol(Y))
  n <- min(nrow(Y), ncol(Y))
  
  if(m==n) {
    rank <- sum(svals > 2.858*median(svals))
  } else {
    beta <- n/m
    omeg <- 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43
    rank <- sum(svals > omeg*median(svals))
  }

  rank
}


sampleV2 <- function(Slist, Ulist, s2vec, OmegaList, V, method="gibbs") {

    K <- length(Ulist)
    S <- ncol(V)
    P <- nrow(V)

    ## First save the Bk matrices for each k
    BkList <- OkList <- list()
    for( k in 1:K ){
        OkList[[k]] <- t(V) %*% Ulist[[k]]
         omegaK <- OmegaList[[k]]
         BkList[[k]] <- OkList[[k]] %*% (diag(omegaK, nrow=length(omegaK)) /
                                (2*s2vec[k])) %*% t(OkList[[k]])
    }
    
    ## Randomly sample a column,  i
    for( i in sample(1:S) ) {

      N <- NullC(V[, -i])

        ## Get the linear term
        ## for j neq i
        C <- rep(0, P)
        for( j in setdiff(1:S, i)) {
            Ak_bij <- matrix(0, P, P)
            for(k in 1:K) {
                Ak <- Slist[[k]]
                bk_ij <- BkList[[k]][i, j]
                Ak_bij <- Ak_bij + Ak*bk_ij
            }
            C <- C+t(V[, j]) %*% Ak_bij
        }
        C <- 2*C

        ## Get the quadratic term
        Ak_bii <- matrix(0, P, P)
        for(k in 1:K) {
            Ak <- Slist[[k]]
            bk_ii <- BkList[[k]][i, i]
            Ak_bii <- Ak_bii+Ak*bk_ii
        }
        A <- Ak_bii
        
        Ctilde <- as.vector(C %*% N)
        Atilde <- t(N) %*% A %*% N

        NV <- as.vector(t(N) %*% V[, i])

        V[, i] <- N %*% R.rbmf.vector.mises(Atilde, Ctilde, NV)

    }

    ## Propose swaps to handle multi-modality
    if(ncol(V) > 1) {

        if(ncol(V)==2) {
            Vswap <- V[, 2:1]
        } else {
            perm <- sample(2:ncol(V))
            onePos <- sample(2:ncol(V), size=1)
            perm <- c(perm[1:(onePos-1)], onePos, perm[onePos:length(perm)])
            Vswap <- V[, sample(1:ncol(V))]
        }
        
        logdens <- function(Vcur) {
            sum(
                sapply(1:length(Slist), function(k) {
                    Ok <- OkList[[k]]
                    omegaK <- OmegaList[[k]]
                    B <- Ok %*% (diag(omegaK, nrow=length(omegaK)) /
                                 (2*s2vec[k])) %*% t(Ok)

                    tr(Slist[[k]] %*% Vcur %*%
                       B %*% t(Vcur))
                })
            )
        }

        logMFswap <- logdens(Vswap)
        logMFcur <- logdens(V)
        log.mh <- logMFswap - logMFcur
        u <- runif(1,0,1)
        if ( log(u) < log.mh ) {
            V <- Vswap
            print("SWAP")
        }
    }

    
    V
}

##############################################################
##############  Optimization on the Stiefel Manifold #########
##############  This is the M-step in the EM algo ############
##############################################################


subspaceEM <- function(Ylist, P, S, R=S, Q=S-R, nvec, rho1=0.1, rho2=0.9,
                       Vstart=NULL, stiefelAlgo=1, lambda=0,
                       EM_iters=10,
                       M_iters=100,
                       verbose=FALSE) {


    ## Function to optimize
    F <- function(V, PhiList, PrecVec) {
        obj <- 0
        for(k in 1:length(PhiList)) {
            VY <- t(V) %*% t(Ylist[[k]])
            obj <- obj +
                1/2 * tr(  VY %*% t(VY) %*% PhiList[[k]] ) -
                1/2 * PrecVec[k] * tr(VY %*% t(VY))
        }
        obj + lambda*sum(abs(V))
    }

    ## dF(X)/dX
    dF <- function(V, PhiList, PrecVec) {
        G <- 0
        for(k in 1:length(PhiList)) {

            G <- G + t(Ylist[[k]]) %*% ((Ylist[[k]] %*% V) %*% PhiList[[k]]) -
                (PrecVec[k] * t(Ylist[[k]])) %*%  (Ylist[[k]] %*% V)

        }
        G + lambda*sign(V)
    }
    
    
    if(is.null(Vstart)) {
        Vstart = rustiefel(P, S)
    }

    Vnew <-  Vstart

    convCheck <- Inf
    iter <- 0
    while(convCheck > 1e-6 & iter < EM_iters ) {
        ## ---------- E-step -----------

        ## E[ 1/sigma^2 * (psi+I)^(-1) | V]
        PrecVec <- rep(1, length(Ylist))
        PhiList <- list()
        V1 <- Vnew[, 1:(S-Q)] 
        if(Q > 0) {
            V2 <- Vnew[, (S-Q+1):S]
            V2Yk <- lapply(1:length(Ylist), function(k) t(V2) %*% t(Ylist[[k]])*sqrt(PrecVec[k]))
            V2Ssum <- Reduce('+', lapply(V2Yk, function(x) x %*% t(x)))
            
            ## Diag Q  to ensure inversion is stable
            PhiShared <- solve(V2Ssum + 1e-5*diag(Q)) * (sum(nvec) + Q + 1)
            
        } else{
            V2 <- matrix(nrow=nrow(V1), ncol=0)
            PhiShared <- matrix(nrow=0, ncol=0)
        }
        
        for(k in 1:length(Ylist)) {

            ## Diag (S-Q)  to ensure inversion is stable
            V1Y <- t(V1) %*% t(Ylist[[k]])
            PsiK <- solve(V1Y %*% t(V1Y) + 1e-5*diag(S-Q)) * (nvec[k] + (S-Q) +1)
            PhiList[[k]] <- as.matrix(bdiag(PsiK, 1/PrecVec[k]*PhiShared))

        }
        
        ## E[ 1/sigma^2  | V]
        for(k in 1:length(Ylist)) {
            PrecVec[k] <- (nvec[k] * (P - S) + 2) /
                (sum(Ylist[[k]]^2) - tr((t(Vnew) %*% t(Ylist[[k]])) %*% (Ylist[[k]] %*% Vnew )))
        }

        ## Objective function to optimize in M-step
        F_t <- function(V) F(V, PhiList, PrecVec)
        dF_t <- function(V) dF(V, PhiList, PrecVec)
        
        ## ------- M-step -----------

        Vnew <- optStiefel(F_t, dF_t, Vinit=Vstart, method="bb",
                           maxIters= M_iters,
                           searchParams = NULL, verbose=verbose)

        ## ---- Check for convergence ------
        
        convCheck <- 1 - (norm(t(Vstart) %*% Vnew, type="F")/sqrt(S))
        Vstart <- Vnew
        iter <- iter + 1
        if(verbose) {
            print(PrecVec)
            print(convCheck)
        }

    }
    
    list(V=Vnew, PhiList=PhiList, PrecVec=PrecVec)
}


## for finding first p eigenvalues of M
sumFirstP <- function(M, n, p, rho1=0.1, rho2=0.9, tol=1-1e-6) {
    
    F <- function(X) { - sum(diag(t(X) %*% M %*% X)) }

    X <- rbind(diag(p), matrix(0, nrow=n-p, p))
    G <- -2*M %*% X

    Fprev <- 0
    while(Fprev/F(X) < tol ) {

        Xprev <- X
        Fprev <- F(Xprev)
        
        tau <- 1
        X <- lineSearch(n, p, X, G, F, rho1, rho2, tau)
p
        G <- -2*M %*% X

        print(sprintf("Objective is %f", F(X)))
    
    }

    X
}

##################################################3


