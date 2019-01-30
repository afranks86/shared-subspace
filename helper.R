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
                          pch=NULL, lty=NULL, ymax=30, logRatio=FALSE,
                          plotPoints=TRUE, cex.axis=1.5) {

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
    ylab <- ifelse(logRatio,
                   expression("log"[2]~"(" ~ lambda[1]/lambda[2] ~ ")"),
                   expression(lambda[1]/lambda[2]))
    plot(0, 0, xlim=c(-pi/2, pi/2),
         ylim=c(0, ymax), cex=0, xlab=expression("angle, acos("~U[1]^T*V[1]~")"), ylab=ylab, xaxt="n", cex.axis=cex.axis, cex.lab=1.5)
    axis(1, at=seq(-pi/2, pi/2, by=pi/4), labels=expression(-pi/2, -pi/4, 0, pi/4, pi/2), cex.axis=cex.axis, cex.lab=1.5)

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
        hullPoints <- hp$hullPoints



        if(plotPoints) {
            points(pts[1, ], pts[2, ], col=alpha(col[g], 1/2), pch=pch[g], cex=0.5)
        } else {
            polygon(pts[1, hullPoints], pts[2, hullPoints], lwd=3, border=col[g], col=alpha(col[g], 1/4), lty=lty[g])
        }
    }

    if(!is.null(hline))
        abline(h=hline, lty=2)

}

getHullPoints <- function(nsamps, OmegaSamps, Osamps, logRatio=FALSE,
                          probRegion=0.95) {

    PointsList <- lapply(1:nsamps, function(i) {
        LambdaSamp <- OmegaSamps[, i]/(1-OmegaSamps[, i])
        maxIndex <- which.max(LambdaSamp)
        LambdaRatio <- LambdaSamp[maxIndex]/LambdaSamp[-maxIndex]
        O1 <- Osamps[, maxIndex, i]
        angle <- atan(O1[2]/O1[1])
        
        c(angle, LambdaRatio)
        
    })

    pts <- simplify2array(PointsList)
    allPts <- pts
    if(logRatio == TRUE) {
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
tgamma <- function(a, b) {

    M <- 1/pgamma(1, a, b)
    qgamma(runif(length(a))/M, a, b)

}

sampleOmega <- function(YV, O, s2, n) {

    YVO <- YV %*% O
    R <- ncol(O)
    cvec <- diag( t(YVO) %*% YVO /(2*s2) )
    g <- tgamma(rep(n/2 + 1, length(cvec)), cvec/n)
    g <- g/n
    
    omega <- 1 - g
    
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


R.rbmf.vector.gibbs <- function (A, c, x) {

    evdA <- eigen(A, symmetric=TRUE)
    E <- evdA$vec
    l <- evdA$val
    y <- t(E) %*% x
    d <- t(E) %*% c
    x <- E %*% R.ry_bmf(y, l, d, length(y))
    x/sqrt(sum(x^2))

}

R.rbmf.vector.mh <- function (A, c, x) {

    evdA <- eigen(A, symmetric=TRUE)
    E <- evdA$vec
    l <- evdA$val

    cnorm <- c/sqrt(sum(c^2))
    c %*% cnorm

    Esigned <- E*sign(cnorm)*sign(E)
    
    mxVec <- which.max(c(c %*% cnorm + cnorm %*% A %*% cnorm, c %*% Esigned + l))
    if(mxVec==1) {
        x <- cnorm
    } else{
        x <- Esigned[, mxVec-1]
    }

    x
}

R.rbmf.vector.hmc <- function (A, c, x, iter=1000) {

    evdA <- eigen(A)
  
    E <- evdA$vec
    l <- evdA$val
    l[l < 0] <- 0
    y <- t(E) %*% x
    d <- as.numeric( t(E) %*% c )
  
    sink("/dev/null");
    results <- stan(fit=stanstart,
                    data=list(M=length(y)-1, lam=l, gamma=d),
                    chains=1, iter=iter)
    sink()
  
    samps <- extract(results, permuted=FALSE)
    samps <- samps[dim(samps)[1], 1, ]
    y <- samps[grep("Y", names(samps))]
    
    x <- E %*% y
    x/sqrt(sum(x^2))

}



R.ry_bmf <- function(y, l, d, n) {

    k <-  (n-1)/2

    for(i in sample(1:n)) {

        omyi <- 1/(1-y[i]^2)
        smyi <- sqrt(omyi)
        a <- l[i]+l[i]*y[i]^2*omyi
        b <- -1*y[i]*d[i]*smyi
        for(j in 1:n) {
            a <- a-l[j]*y[j]^2*omyi
            b <- b+y[j]*d[j]*smyi
        }
##        print(y)
##        print(i)
##        print(l)
        theta <- R.rtheta_bmf.mh(k, a, b, abs(d[i]))
        if(theta==1e-16) {

        }
        
        for(j in 1:n){
            y[j] <- y[j]*sqrt(1-theta)*smyi
        }
        y[i] <- sqrt(theta)*(-1^(rbinom(1,1,1/(1+exp(2*sqrt(theta)*d[i]) ))) )
    }
    y
    ##rstiefel::ry_bmf(y, l, d)
}


R.rtheta_bmf.mh <- function(k, a, b, cc, steps=50) {

    f <- function(x) {
        -1/2*log(x) + k*log(1-x) + a*x+b*sqrt(1-x) + sqrt(x)*cc +
            log(1.0 + exp(-2*sqrt(x)*cc))
    }

    mode <- optimize(function(x) - f(x),
                     lower=0, upper=1,
                     tol=.Machine$double.eps)$minimum

    if(f(1e-16) > f(mode)) {
        mode <- 1e-16
    }

    fdoubleprime <- function(x) {
        -b/(4*(1-x)^(3/2)) + cc^2/(x*(exp(cc*sqrt(x))+1)) - k/(1-x)^2 + 1/(2*x^2) - cc*tanh(cc*sqrt(x))/(4*x^(3/2) )
        
    }

    dd <- fdoubleprime(mode)
    if( dd >= 0 ) {
        p <- 1/2
        q <- 1## 1/(2*mode)-1/2
    } else {
        pplusq <- -1*mode*(1-mode)*dd-1
        p <- mode*(pplusq-2)+1
        if( p < 1 ) {
            browser()
            print("P < 1")
            p <- 1
            q <- (1-mode)/mode
        } else {
            q <- pplusq - p
        }
    }
    
    if(mode < 0.1) {
        ##browser()
    }

    lprop <- function(x) { dbeta(x, p, q, log=TRUE) }
    thCur <- mode
    reject <- 0
    for(i in 1:steps) {

        th <- rbeta(1, p, q)

        log.mh <- f(th) - f(thCur) +
            lprop(thCur) - lprop(th)
        
        if( log.mh > 0 ) {
            thCur <- th
        } else {
            u <- runif(1,0,1)
            if ( log(u) < log.mh ) {
                thCur <- th
            } else {
                reject <- reject+1
            }
        }
    }

    if( reject==steps ) {
        print("ALL REJECTS")
        ## dn <- function(x) dnorm(x, mode, sd=sqrt(-1/dd), log=TRUE)
        ## print(mode)
        ## curve(f(x), from=0, to=0.001)
        ## curve(lprop(x)-(lprop(mode)-f(mode)), col="red", add=TRUE)
        ## lprop2 <- function(x) { dnorm(x, mean=mode, sd=sqrt(-1/dd), log=TRUE) }
        ## curve(lprop2(x)-(lprop2(mode)-f(mode)), col="green", add=TRUE)
        ## curve(-1/2*log(x) - (-1/2*log(1e-5)-f(1e-5)), from=0, to=0.001, col="red", add=TRUE)
        thCur <- mode
        
    }

  ##print(reject/steps)
    
    thCur
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
    apply(simplify2array(PsiList), c(1, 2), mean)
}

getPostMeanSigma <- function(P, Ulist, OmegaList, s2vec) {

    SigmaList <- sapply(1:nsamps, function(i) {
        s2vec[i]*(Ulist[, , i] %*% diag((1-omega)/omega) %*% t(Ulist[, , i]) + diag(P))
    })
    apply(simplify2array(SigmaList), c(1, 2), mean)
}

## Get mean predictions for each group relative to truth
makePrediction <- function(Usamps, omegaSamps, s2samps, SigmaTrueList,
                           ngroups=dim(Usamps)[3],
                           genGroups=(1:ngroups), n=20,
                           nsamps=dim(Usamps)[4],
                           numToAvg=nsamps/2) {
    if(ngroups==1) {
        return(matrix(1, ncol=1, nrow=1))
    }
    
  predictionMat <- matrix(0, nrow=length(genGroups), ncol=ngroups)
  for(group in 1:length(genGroups)) {
    Yk <- rmvnorm(n, sigma=SigmaTrueList[[genGroups[group]]])
    SigmaHatList <- list()
    
    probs <- rep(0, ngroups)
    for(i in (nsamps-numToAvg+1):nsamps) {
      
      for(k  in 1:ngroups) {

        U <- Usamps[, , k, i]
        om <- omegaSamps[, k, i]
        Lam <- diag(om/(1-om))
        s2 <- s2samps[k, i]
        SigmaHatList[[k]] <- U %*% Lam %*% t(U) + s2*diag(nrow(Usamps))
      }
      probs <- probs +
        rowMeans(computeMembershipProbabilities(SigmaHatList, Yk))
    }
    probs <- probs / numToAvg

    predictionMat[group, ] <- probs
  }

  predictionMat
}

computeMembershipProbabilities <-

  function(SigmaList, Y,
           priorWeights=rep(1/length(SigmaList), length(SigmaList))) {

    priorWeights <- priorWeights / sum(priorWeights)
    unnormalizedProbs <- sapply(SigmaList, function(Sigma) {
      dmvnorm(Y, sigma=Sigma, log=TRUE)
    })

    if( !class(unnormalizedProbs) == "matrix")
      unnormalizedProbs <- matrix(unnormalizedProbs, nrow=1, ncol=10)
    
    probs <- apply(unnormalizedProbs, 1, function(probVec) {
      normProbs <- probVec - max(probVec)
      normProbs <- exp(normProbs) * priorWeights /
        sum(exp(normProbs) * priorWeights)
    
      normProbs
    })

    probs
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

R.rbmf.vector.mises <- function(Atilde, Ctilde, xinit) {

    xcur <- xinit
    reject <- 0
    for(i in 1:50) {

        xprop <- rmf.vector(kmu=xcur * 10000)
        logMFprop <- Ctilde %*% xprop + t(xprop) %*% Atilde %*% xprop
        logMFcur <- Ctilde %*% xcur + t(xcur) %*% Atilde %*% xcur
        log.mh <- logMFprop - logMFcur

        if( log.mh > 0 ) {
            xcur <- xprop
        } else {
            u <- runif(1,0,1)
            if ( log(u) < log.mh ) {
                xcur <- xprop
            } else {
                reject <- reject+1
            }
        }
    }

    xcur
    
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
        PhiList <- list()
        V1 <- Vnew[, 1:(S-Q)] 
        if(Q > 0) {
            V2 <- Vnew[, (S-Q+1):S]
            V2Yk <- lapply(Ylist, function(x) t(V2) %*% t(x))
            V2Ssum <- Reduce('+', lapply(V2Yk, function(x) x %*% t(x)))
            PhiShared <- solve(V2Ssum) * (sum(nvec) + Q + 1 )
        } else{
            V2 <- matrix(nrow=nrow(V1), ncol=0)
            PhiShared <- matrix(nrow=0, ncol=0)
        }
        
        for(k in 1:length(Ylist)) {

            V1Y <- t(V1) %*% t(Ylist[[k]])
            PsiK <- solve(V1Y %*% t(V1Y)) * (nvec[k] + (S-Q) +1)
            PhiList[[k]] <- as.matrix(bdiag(PsiK, PhiShared))

        }
        
        ## E[ 1/sigma^2  | V]
        PrecVec <- c()
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


optimV3 <- function(Ylist, P, S, nvec, Vinit=NULL, rho=0.1,
                    maxIters=50, verbose=FALSE, tauInit=2, subsetSize=0.3, lambda=0) {
    
    if(is.null(Vinit)) 
        V <- rustiefel(P, S)
    else
        V <- Vinit## G <- G +
            ##     nvec[k]/2 * (2*VSVinv - diag(diag(VSVinv))) %*% t(V) %*% Slist[[k]]
            ## + nvec[k] * (P - S)/2 * diag(-1/(tr(Slist[[k]]) - tr(VSV)), S) %*% t(V) %*% Slist[[k]]

    Slist <- lapply(Ylist, function(Y) t(Y) %*% Y)
    ## Function to optimize
    F <- function(V) {
        obj <- 0
        for(k in 1:length(Slist)) {
            VSV <- t(V) %*% Slist[[k]] %*% V
            obj <- obj + nvec[k]/2 * determinant(VSV)$modulus  +
                nvec[k]*(P-S)/2 * log(tr(Slist[[k]]) - tr(VSV))
        }
        1/sum(nvec)*obj + lambda*sum(abs(V))
    }

    ## dF(X)/dX
    dF <- function(V) {

        ## Subsample to get stochastic gradientn
        SlistSub <- lapply(Ylist, function(Y) {
            nk <- nrow(Y)
            indices <- sample(1:nk, round(nk*subsetSize))
            Ysub <- Y[indices, ]
            
            t(Ysub) %*% Ysub
        })
        
        G <- 0
        for(k in 1:length(SlistSub)) {

            VSV <- t(V) %*% SlistSub[[k]] %*% V
            VSVinv <- solve(VSV)
            
            G <- G +
                nvec[k]/2 * (2*VSVinv - diag(diag(VSVinv))) %*% t(V) %*% SlistSub[[k]]
            + nvec[k] * (P - S)/2 * diag(-1/(tr(SlistSub[[k]]) - tr(VSV)), S) %*% t(V) %*% SlistSub[[k]] 

        }
        -1/sum(nvec)*t(G) - lambda*sign(V)
    }

    
    eta <- 0.1
    iter <- 1
    Vprev <- rustiefel(P, S)
    Fprev <- Inf
    Fcur <- F(V)

    Ccur <- Fcur
    Gcur <- dF(V)
    Gprev <- Gcur

    Qcur <- 1
    tau <- tauInit
    ##while(Fcur/Fprev < 1-1e-15 & iter < maxIters) {
     while(iter < maxIters) {

        if(verbose) {
            print(sprintf("Iteration %i: Obj=%f, Tau=%f", iter, Fprev, tau))
        }
        
        Fprev <- Fcur
        ## newV <- linreSearchBB(P, S, V, Vprev, Gcur, Gprev, F, rho=0.01, Ccur)
        ##newV <- lineSearch(P, S, V, Gcur, F, rho1=0.01, rho2=0.9, tauStart=1)
        newV <- takeStep(S, V, Gcur, tau)
        if(iter %% 100 == 0)
        tau <- tauInit * 1/(iter/100+1)
         
        ## Update
        Vprev <- V
        V <- newV
        Fcur <- F(V)
        Gprev <- Gcur
        Gcur <- dF(V)

        Qprev <- Qcur
        Qcur <- eta*Qcur + 1
        
        Ccur <- (eta*Qprev*Ccur + Fcur) / Qcur
        
        iter <- iter + 1
    }

    V
}

takeStep <- function(p, X, G, tau) {

    A <- G %*% t(X) - X %*% t(G)
    U <- cbind(G, X)
    V <- cbind(X, -1*G)
    
    HV <- solve(diag(2*p) + tau/2 * t(V) %*% U) %*% t(V)
    Ytau <- X - tau * U %*% (HV %*% X)

}
