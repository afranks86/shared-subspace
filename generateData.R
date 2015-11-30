 generateData <- function(P=100, S=10, R=5, ngroups=10) {

    library(rstiefel)
    library(mvtnorm)

    library(rstiefel)
    source("helper.R")


    load("ae_XY_m527.RData")
    XDM <- model.matrix(~ -1 + as.factor(Age):as.factor(Sex), data=X )
    fit<-lm(Y ~ -1+XDM)

    ## number of eigenvectors if pooled space
    S <- 10
    ## number of eigenvectors in group subspace
    R <- 5
    nsamps <- 10
    P <- ncol(Y)

    ngroups <- ncol(XDM)

    res <- fit$residuals
    eigen.all <- eigen(t(res) %*% res/(nrow(res)-1))
    V <- eigen.all$vectors[1:P, 1:S]
    V <- sweep(V,2,sqrt(colSums(V^2)),'/')

    ## Generate data to look like real data

    Slist <- Ulist <- Olist <- OmegaList <- SigmaList <- eigenlist <- list()
    nvec <- numeric(ngroups)
    s2vec <- numeric(ngroups)

    for ( k in 1:ngroups ) {

        print(k)

        indices <- which(XDM[, k]==1)
        res <- fit$residuals[indices, ]
        cur <- t(res) %*% res/(nrow(res)-1)
        ecur <- eigen(cur)
        eigenlist[[k]] <- ecur
        
        s2vec[k] <- 1
        
        ## Generate Ok
        Ok <- matrix(0,nrow=S,ncol=R)
        rows <- sample(1:S,R)
        for(r in 1:R) {
            Ok[rows[r],r] <- 1
        }
        Olist[[k]] <- Ok
        Uk <- V%*%Ok
        Ulist[[k]] <- Uk
        
        Lam <- diag(head(ecur$values,n=R))
        OmegaList[[k]] <- diag(Lam/(Lam+1))
        nvec[k] <- 20
        SigmaList[[k]] <- s2vec[k]*(Uk%*%Lam%*%t(Uk)+diag(P))

        res <- rmvnorm(n=nvec[k],sigma=SigmaList[[k]])
        Slist[[k]] <- t(res)%*%res

    }

    list(V=V, Slist=Slist, Ulist=Ulist, Olist=Olist, OmegaList=OmegaList,
         SigmaList=SigmaList, eigenlist=eigenlist, s2vec=s2vec,
         S=S, R=R, P=P, ngroups=ngroups, nvec=nvec)
    
}


if( FALSE ) {

    source("helper.R")
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
        Ok <- proposeBinaryO(S, R, U, V, Sig, s2, omega, nvec[1], nflip=2)
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
    Vinit <- V
    Vcur <- Vinit
    dp <- numeric(10)
    for(i in 1:10) {
        print(sprintf("---- Iteration %i------",i))
        Vcur <- sampleV(Slist, Ulist,s2vec, OmegaList, Vcur, method="gibbs")
        dp[i] <- norm(t(Vcur)%*%V,type="F")
    }
    plot(c(norm(t(Vinit)%*%V,type="F"),dp[1:10]),type="l",ylim=c(0,4.5))
    Vcur <- Vinit
    dp <- numeric(10)
    for(i in 1:10) {
        print(sprintf("---- Iteration %i------",i))
        Vcur <- sampleV(Slist, Ulist,s2vec, OmegaList, Vcur, method="hmc")
        dp[i] <- norm(t(Vcur)%*%V,type="F")
    }
    lines(c(norm(t(Vinit)%*%V,type="F"),dp[1:10]),type="l",col="blue")

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
