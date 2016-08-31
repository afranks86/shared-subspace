## P - number of features (cols) for each group of measurements
## S - Number of columns of V
## R - Number of columns of U = VO
## nvec = sample size for each group
## s2vec - measurement noise variance
## LambdaList - list of length ngroups, each a vec of evals length R
## Olist - 

library(rstiefel)
library(mvtnorm)
library(rstiefel)
source("helper.R")

generateData <- function(P=200, S=10, R=S, Q=S-R, ngroups=10,
                         nvec=rep(100, ngroups),
                         V=rbind(diag(S), matrix(0, nrow=P-S, ncol=S)),
                         s2vec=rep(1, ngroups), LambdaList=NULL,
                         Olist=NULL) {

    if(Q + R != S)
        stop("Q + R must be equal S")

    if( is.null(LambdaList) ) {
        LambdaList <- vector("list", ngroups)
    } else if (length(LambdaList) != ngroups) {
        stop(sprintf("LambdaList must be of length ngroups=%s"), ngroups)
    } 

    if( is.null(Olist) ){
        Olist <- vector("list", ngroups)
    } else if (length(Olist) != ngroups) {
        stop(sprintf("Olist must be of length ngroups=%s", ngroups))
    }
    
    Slist <- Ylist <- Ulist <- OmegaList <- SigmaList <- list()

    if(Q > 0) {
        Oshared <- rbind(matrix(0, nrow=(S-Q), ncol=Q), rustiefel(Q, Q))
    } else {
        Oshared <- matrix(ncol=0, nrow=S)
    }

    LambdaShared <- sort(rexp(Q, 1/10), decreasing=TRUE)
    for ( k in 1:ngroups ) {

        if( is.null(Olist[[k]]) ) {
            print(k)
            ## Generate Ok
            Odiff <- rbind(rustiefel(R, R), matrix(0, nrow=S-R, ncol=R))
            Ok <- cbind(Odiff, Oshared)
            Olist[[k]] <- Ok
        }


        if( is.null(LambdaList[[k]]) ){
            LambdaDiff <- sort(rexp(R, 1/10), decreasing=TRUE)
            LambdaList[[k]] <- c(LambdaDiff, LambdaShared)
        }

        Ok <- Olist[[k]]
        Uk <- V %*% Ok
        Ulist[[k]] <- Uk
        Lamk <- LambdaList[[k]]
        OmegaList[[k]] <- Lamk/(Lamk+1)

        LamMat <- diag(Lamk, nrow=length(Lamk), ncol=length(Lamk))
        SigmaList[[k]] <- s2vec[k]*(Uk %*% LamMat %*% t(Uk) + diag(P))

        ## Generate sample covariance matrix and save
        
        Y <- rmvnorm(n=nvec[k], sigma=SigmaList[[k]])
        ## Y <- sqrt(s2vec[k])*(matrix(rnorm(nvec[k]*P), nrow=nvec[k], ncol=P) %*% V %*% Ok %*% sqrt(LamMat) %*% t(Ok) %*% t(V) + matrix(rnorm(nvec[k]*P), nrow=nvec[k], ncol=P))
        Ylist[[k]] <- Y
        Slist[[k]] <- t(Y) %*% Y

    }

    list(V=V, Slist=Slist, Ylist=Ylist, Ulist=Ulist, Olist=Olist,
         OmegaList=OmegaList, LambdaList=LambdaList, SigmaList=SigmaList,
         s2vec=s2vec, S=S, R=R, Q=Q, P=P, ngroups=ngroups, nvec=nvec)
    
}



if( FALSE ) {


    genData <- generateData()
    ## Sample only Conditionals
    
    S <- genData$S
    R <- genData$R
    U <- genData$Ulist[[1]]
    V <- genData$V
    Sig <- genData$Slist[[1]]
    s2 <- genData$s2vec[1]
    omega <- genData$OmegaList[[1]]
    nvec <- genData$nvec

    count <- 0
    ## Ocur <- round(t(V) %*% U, digits=10)
    Ok <- matrix(0,nrow=S,ncol=S)
    rows <- sort(sample(1:S,R))
    for(i in rows) {
        Ok[i,i] <- 1
    }
    U <- V %*% Ok
    for( i in 1:1000 ) {
        Ok <- proposeBinaryO(S, U, V, Sig, s2, omega, nvec[1], flipProb=0.1)
        U <- V %*% Ok
        if(any(new!=Ocur)) {
            count <- count+1
        }
    }
    print(count/1000)
    
    k <- 1

    ## Sample s2vec
    s2vec.samps <- sapply(1:100, function(i) rs2_fc(Slist[[k]], Ulist[[k]], OmegaList[[k]], nvec[k]))
    hist(s2vec.samps,col="grey")
    abline(v=1,col="red",lwd=3)

    ## Sample omegas
    omegas <- lapply(1:100, function(i) romega_fc(Slist[[k]], Ulist[[k]], s2vec[k], nvec[k]))
    hist(sapply(omegas,function(x) x[1]),col="red",xlim=c(0.8,1))
    abline(v=OmegaList[[k]][1],col="red",lwd=3,lty=2)

    hist(sapply(omegas,function(x) x[R]),col="blue",add=TRUE)
    abline(v=OmegaList[[k]][R],col="blue",lwd=3,lty=2)

    ## Sample O's
    U <- Ulist[[k]]
    O.init <- t(V)%*%U
    O.samps <- list()
    for(i in 1:1000) {
        print(i)
        Ok <- rU_csm_gibbs(Slist[[k]], U, s2vec[k], OmegaList[[k]], V)
        U <- V%*%Ok
        O.samps[[i]] <- Ok
    }
    hist(sapply(1:1000, function(i) norm(t(O.init)%*%O.samps[[i]],type="f")),breaks=10,xlim=c(2,3.5),col="red")
    hist(sapply(1:1000, function(i) norm(t(t(V)%*%Ulist[[k+1]])%*%O.samps[[i]],type="f")),breaks=10,xlim=c(2.5,3.5),add=TRUE)
    hist(sapply(1:1000, function(i) norm(t(t(V)%*%Ulist[[k+2]])%*%O.samps[[i]],type="f")),breaks=10,xlim=c(2.5,3.5),add=TRUE)
    abline(v=norm(t(O.init)%*%(t(V)%*%Ulist[[k]]),type="F"),lwd=3,col="red")

    ## Sample V's succeeds...
    Vinit <- rustiefel(P, S)
    Vcur <- Vinit
    dp <- numeric(100)
    for(i in 1:100) {
        print(sprintf("---- Iteration %i------",i))
        Vcur <- sampleV(dat$Slist, dat$Ulist, dat$s2vec, dat$OmegaList, Vcur, method="gibbs")
        dp[i] <- norm(t(Vcur)%*%dat$V,type="F")
    }
    plot(c(norm(t(Vinit) %*% dat$V,type="F"),dp[1:100]),type="l",ylim=c(0,4.5))
    Vcur <- Vinit
    dp <- numeric(100)
    for(i in 1:100) {
        print(sprintf("---- Iteration %i------",i))
        Vcur <- sampleV(dat$Slist, dat$Ulist, dat$s2vec, dat$OmegaList, Vcur, method="hmc")
        dp[i] <- norm(t(Vcur) %*% dat$V,type="F")
    }
    lines(c(norm(t(Vinit) %*% V,type="F"), dp[1:100]), type="l",col="blue")

    Vinit <- matrix(0,nrow=nrow(V),ncol=ncol(V))
    Vinit[1:ncol(V),1:ncol(V)] <- diag(ncol(V))
    Vcur <- Vinit
    dp <- numeric(10)
    for(i in 1:10) {
        print(sprintf("---- Iteration %i------",i))
        Vcur <- sampleV(Slist, Ulist,s2vec, OmegaList, Vcur, method="gibbs")
        dp[i] <- norm(t(Vcur)%*%V,type="F")
    }
    lines(c(norm(t(Vinit)%*%V,type="F"),dp[1:10]),lty=2)
    Vcur <- Vinit
    dp <- numeric(10)
    for(i in 1:10) {
        print(sprintf("---- Iteration %i------",i))
        Vcur <- sampleV(Slist, Ulist,s2vec, OmegaList, Vcur, method="hmc")
        dp[i] <- norm(t(Vcur)%*%V,type="F")
    }
    lines(c(norm(t(Vinit)%*%V,type="F"),dp[1:10]),lty=2,col="blue")

    Vinit <- matrix(0,nrow=nrow(V),ncol=ncol(V))
    Vinit[1:ncol(V),1:ncol(V)] <- diag(ncol(V))
    dp <- numeric(10)
    for(i in 1:10) {
        print(sprintf("---- Iteration %i------",i))
        Vcur <- sampleV(lapply(SigmaList,function(s) s*500), Ulist,s2vec, OmegaList, Vcur)
        dp[i] <- norm(t(Vcur)%*%V,type="f")
    }
    plot(dp,type="l")
}
