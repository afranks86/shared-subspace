library(rstan)

library(mvtnorm)

stanstart <- stan(file="vectorBMF.stan", data=list(M=9,lam=rep(1,10), gamma=rep(1, 10)), chains=1, iter=1)

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

getGenderStat <- function(normsMat ) {

    withinGender <- c(as.vector(normsMat[1:6, 1:6]), as.vector(normsMat[7:12, 7:12]))
    betweenGender <- normsMat[1:6, 7:12]

    median(withinGender, na.rm=TRUE)-median(betweenGender)

    
}

getAgeStat <- function(normsMat) {

    ages <- unique(X$Age)
    ageStats <- list()
    ageDiff <- c()
    
    UL <- normsMat[1:6, 1:6]
    LR <- normsMat[7:12, 7:12]
    for( i in 1:5 ) {

        ageStats[[i]] <- c(UL[row(UL)==(col(UL)-i)], 
                           (LR)[row(LR)==(col(LR)-i)])
        ageDiff <- c(ageDiff, rep(diff(ages, lag=i), 2))
    }

    Yage <- unlist(ageStats)
    Xage <- ageDiff
    
    ## return slope of norms against age
    slope <- coef(lm(formula = Yage ~ Xage))[2]

    list(slope=slope, ages=Xage, norms=Yage)

}

##############################################################
### Functions for Sampling
#############################################################

## Sample from 1/sigma^2
sampleSigma2 <- function(S, U, omega, n, nu0=1, s20=1) {

    p <- nrow(S)
    a <- ( nu0 + n*p)/2
    b <- ( nu0*s20 + tr(S%*%( diag(p)  - U %*% diag(omega, nrow=length(omega)) %*% t(U))) )/2
    1/rgamma(1, a, b)
}

## Unormalized eigenvalues
debeta_un <- function(w, a, b, c, log=FALSE) {
    x <- (a-1)*log(w) + (b-1)*log(1-w) + c*(w-1)
    if(!log){
        x <- exp(x)
    }
    x
}

## normalizing constant
debeta_nc <- function(a, b, c) {
    1/(integrate(debeta_un, 0, 1, a=a, b=b, c=c)$val)
}

## normalized
debeta <- function(w, a, b, c){
    debeta_un(w, a, b, c)*debeta_nc(a, b, c)
}

pebeta<-function(w, a, b, c) {
    dint<- function(w, a, b, c) {
        int <- 0
        if(w > 0) {
            int <- integrate(debeta_un, 0, w, a, b, c, 
                             rel.tol = .Machine$double.eps^0.75,
                             subdivisions=200L)$value
        }
        int
    }
    pmin(sapply(w, dint, a, b, c)/dint(1, a, b, c), 1)
}

qebeta<-function(p, a, b, c) {

    f <- function(w, p, a, b, c){ p - pebeta(w, a, b, c) }
    w <- p*0
    for(i in 1:length(p)) {
        w[i] <- uniroot(f, interval=c(1e-10, 1-1e-10), p[i], a, b, c)$root
    }
    w
}

rebeta<-function(n, a, b, c, interval=c(0, 1), MH=1000) {

    uinterval <- sort( pebeta(interval, a, b, c) )
    if(uinterval[1] < uinterval[2] ) {
        u <- runif(n,  uinterval[1],   uinterval[2])
        w <- qebeta(u, a, b, c)
    }
    if(uinterval[1] == uinterval[2] ) {
        w<-runif(n, interval[1], interval[2])
    }

    w[w>=interval[2]]<-interval[2]-1e-6
    w[w<=interval[1]]<-interval[1]+1e-6

    for(s in seq(1, MH, length=MH)) {
        wp <- pmin(1-1e-6, pmax(0+1e-6, w+runif(n, -1e-3, 1e-3) ))

        lhr <- debeta_un(wp, a, b, c, log=TRUE) -
            debeta_un(w, a, b, c, log=TRUE) -
            log(interval[1] < wp | wp < interval[2])
        lu <- log(runif(n))
        w[lu < lhr] <- wp[lu < lhr]
    }
    
    w
}

sampleOmega <- function(Sig, U, s2, n, a=1, b=1) {

    R <- ncol(U)
    cvec <- diag( t(U) %*% Sig %*% U/(2*s2) )
    omega <- numeric(R)

    for(r in 1:R) {
        if(cvec[r]==0) {
            omega[r] <- rbeta(1, a, n/2+b)
        } else {
            omega[r] <- rebeta(1, a, b+n/2, cvec[r])
        }
    }

    omega
}

## If we want to maintain order of eigenvalues
rsomega_gibbs <- function(S, U, s2, n, omega, a=1, b=1) {

    R<-ncol(U) ; c<-diag( t(U)%*%S%*%U/(2*s2) )
    for(r in 1:R) {
        lb<-omega[r+1] ; if (is.na(lb)){ lb<-0 }
        ub<-omega[r-1] ; if (length(ub)==0){ ub<-1 }
        if(any(is.na(c(lb, ub))))
            browser()
        omega[r]<-rebeta(1, a, b+n/2, c[r], c(lb, ub))
        
        if(is.na(omega[r])){
            omega[r]<-runif(1, lb, ub)
        }
    }
    omega
}

## Shared subspace Sampling
## phi is the prior concentration
sampleO <- function(SC, U, s2, omega, V, phi=0) {

    O <- t(V) %*% U

    ## if phi > 0, prior concentrates around first R eigenvectors of V
    prior.target <- matrix(0, nrow=ncol(V), ncol=ncol(V))
    prior.target[1:ncol(U), 1:ncol(U)] <- phi*diag(ncol(U))

    A <- t(V) %*% (SC/(2*s2))%*%V+prior.target
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
    BkList <- list()
    for( k in 1:K ){
         Ok <- t(V) %*% Ulist[[k]]
         omegaK <- OmegaList[[k]]
         BkList[[k]] <- Ok %*% (diag(omegaK, nrow=length(omegaK)) /
                                (2*s2vec[k])) %*% t(Ok)
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
        
        Ctilde <- as.vector(C%*%N)
        Atilde <- t(N)%*%A%*%N

        NV <- as.vector(t(N)%*%V[, i])
      
        if( method == "gibbs" ) {
            V[, i] <- N %*% R.rbmf.vector.gibbs(Atilde, Ctilde, NV)
        } else if( method == "hmc" ) {
            V[, i] <- N %*% R.rbmf.vector.hmc(Atilde, Ctilde, NV)
        } else {
            stop("Unknown method")
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

    evdA <- eigen(A)
    E <- evdA$vec
    l <- evdA$val
    y <- t(E) %*% x
    d <- t(E) %*% c
    x <- E %*% R.ry_bmf(y, l, d, length(y))
    x/sqrt(sum(x^2))

}

R.rbmf.vector.hmc <- function (A, c, x, iter=10) {

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

    ## Use the last sample
    y <- samps[nrow(samps), ]

    x <- E %*% y
    x/sqrt(sum(x^2))

}



R.ry_bmf <- function(y, l, d, n) {

    k <-  (n-1)/2

    for(i in 1:n) { 
        omyi <- 1/(1-y[i]^2)
        smyi <- sqrt(omyi)
        a <- l[i]+l[i]*y[i]^2*omyi
        b <- -1*y[i]*d[i]*smyi
        for(j in 1:n) {
            a <- a-l[j]*y[j]^2*omyi
            b <- b+y[j]*d[j]*smyi
        }

        theta <- R.rtheta_bmf.mh(k, a, b, abs(d[i]))
        for(j in 1:n){
            y[j] <- y[j]*sqrt(1-theta)*smyi
        }
        y[i] <- sqrt(theta)*(-1^(rbinom(1,1,1/(1+exp(2*sqrt(theta)*d[i]) ))) )
    }
    y
}


R.rtheta_bmf.mh <- function(k, a, b, c, steps=50) {

    w <- c
    u <- Inf
    g <- k

    if(a>0) {
        g <- max(1/(1+log(2+a)),k-a)
    }
    count <- 1
    
    f <- function(x) { 1/2*log(x)+k*log(1-x)+a*x+b*sqrt(1-x)+sqrt(x)*c+log(1.0+exp(-2*sqrt(x)*c)) }
    mode <- optimize(function(x) -f(x),lower=0,upper=1, tol=.Machine$double.eps)$minimum

    fprime <- function(x) {
        a - b/(2*sqrt(1-x)) - k/(1-x) + 1/(2*x) + c/(2*sqrt(x))
    }
    
    fdoubleprime <- function(x) {
        -b/(4*(1-x)^(3/2)) + c^2/(x*(exp(c*sqrt(x))+1)) - k/(1-x)^2 - 1/(2*x^2) - c*tanh(c*sqrt(x))/(4*x^(3/2) )
        
    }

    dd <- fdoubleprime(mode)
    if( dd >= 0 ) {
        p <- 1/2
        q <- 1/(2*mode)-1/2
    } else {
        pplusq <- -1*mode*(1-mode)*dd-1
        p <- mode*pplusq
        if( p < 1 ) {
            p <- 1
            q <- (1-mode)/mode
        } else {
            q <- pplusq - p
        }
    }

    lprop <- function(x) { dbeta(x,p,q, log=TRUE) }
    thCur <- mode
    reject <- 0
    for(k in 1:steps ) {

        u <- runif(1,0,1)
        th <- ifelse(runif(1,0,1)<1/2,rbeta(1,1/2,g),rbeta(1,p,q))
        th <- rbeta(1,p,q)

        log.mh <- f(th)-f(thCur)+lprop(thCur)-lprop(th)
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
        
        browser()
        print(mode)
        curve(f(x))
        curve(lprop(x)-optimize(function(x) lprop(x)-f(x),interval=c(0,1))$objective,col="red",add=TRUE)
        
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

getPostMeanSigma <- function(P, Ulist, OmegaList, s2vec) {

    SigmaList <- sapply(1:nsamps, function(i) {
        s2vec[i]*(Ulist[, , i] %*% diag((1-omega)/omega) %*% t(Ulist[, , i]) + diag(P))
    })
    apply(simplify2array(SigmaList), c(1, 2), mean)
}

## Get mean predictions for each group relative to truth
makePrediction <- function(Usamps, omegaSamps, s2samps, SigmaTrueList,
                           genGroup=1, n=20,
                           ngroups=dim(Usamps)[3],
                           nsamps=dim(Usamps)[4],
                           numToAvg=nsamps/2) {

  
  Yk <- rmvnorm(n, sigma=SigmaTrueList[[genGroup]])
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

  probs

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
