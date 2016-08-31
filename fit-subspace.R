## P=number of features, S=subspace dimension
## psi_k has dimension R + Q <= S, where Q eigenvectors are common
## and R eigenvectors vary across groups. 
##
## Slist: a list of p x p matrices corresponding to Y^TY
## init: a list containing initial V, U's, Omega's,
## nvec: vector of number of samples per group

fitSubspace <- function(P, S, R, Q=S-R, Slist, nvec, ngroups=length(Slist),
                        niters=100, nskip=1, init=NULL, binaryO=FALSE,
                        verbose=TRUE, sigmaTruthList=NULL, draw=c(),
                        saveU=FALSE,
                        printLoss=FALSE,
                        Vmode="hmc") {
    

    niters <- niters
    nskip <- nskip
    nsamps <- floor(niters/nskip)

    if(R + Q > S)
        stop("R + Q must be less than S")
    

    if(saveU)
        Usamps <- array(dim=c(P, RQ, ngroups, nsamps))
    else
        Usamps <- NULL
    Osamps <- array(dim=c(S, R + Q, ngroups, nsamps))
    Vsamps <- array(dim=c(P, S, nsamps))
    omegaSamps <- array(dim=c(R + Q, ngroups, ncol=nsamps))
    s2samps <- matrix(nrow=ngroups, ncol=nsamps)

    initItems <- c("V", "Ulist", "OmegaList")
    if( is.null(init) ) {

        V1 <- matrix(0, nrow=P, ncol=(S-Q))
        V1[1:(S-Q), 1:(S-Q)] <- diag(S-Q)

        V2 <- matrix(0, nrow=P, ncol=Q)
        V2[(S-Q+1):S, (S-Q+1):S] <- diag(Q)

        V <- cbind(V1, V2)
        
        OmegaList <- Olist <- Ulist <- list()

        SigInvInitList <- list()
        s2vec <- rexp(ngroups)
        for(k in 1:ngroups) {
            OmegaList[[k]] <- rep(1/2, (R+Q))

            Ok <- matrix(0, nrow=S, ncol=(R+Q))
            Ok[1:R, 1:R] <- diag(R)
            if(Q > 0)
                Ok[(S-Q+1):S, (R+1):(R+Q)] <- diag(Q)

            Olist[[k]] <- Ok
            Ulist[[k]] <- V %*% Ok

            SigInvInitList[[k]] <-
                getSigmaInv(P, Ulist[[k]], OmegaList[[k]], s2vec[k])
        }


        init <- list(V=Vinit, Ulist=Ulist, OmegaList=OmegaList, s2vec=s2vec)
        
    } else if (any(!(initItems %in% names(init)))) {
        stop("Some parameters unintialized!")
    }

    sigmaTruthInvList <- NULL
    if(!is.null(sigmaTruthList)) {
        sigmaTruthInvList <- list()
        for(k in 1:ngroups) {
            sigmaTruthInvList[[k]] <- solve(sigmaTruthList[[k]])
        }
    }
    
    V <- init$V
    Ulist <- init$Ulist
    Olist <- lapply(Ulist, function(u) t(V) %*% u)
    OmegaList <- init$OmegaList
    s2vec <- init$s2vec
    
    initSigmaList <- list()
    for(k  in 1:ngroups) {
        omega <- OmegaList[[k]]
        Lambda <- diag((omega/(1-omega)))
        initSigmaList[[k]] <-
            init$s2vec[k] * (Ulist[[k]] %*% Lambda %*% t(Ulist[[k]]) + diag(P))
    }


    draw <- c(draw, V=TRUE, O=TRUE, s2=TRUE, omega=TRUE)
    draw <- draw[unique(names(draw))]
    
    for ( i in 1:niters ) {

        if(draw["V"])
            V <- sampleV(Slist, Ulist, s2vec, OmegaList,
                         V, method=Vmode)
        
        Ssum <- lapply(1:ngroups, function(k) Slist[[k]]/s2vec[k]) %>%
            Reduce('+', .)

        ## Sample common omega and common O
        Uk <- Ulist[[1]]

        if(Q > 0) {
            ## common omega
            Omega2 <- sampleOmega(Ssum, Uk[, (R+1):(R+Q)],
                                  1, sum(nvec[k]))
        
            ## Sample common O
            O2 <- sampleO(Ssum, Uk[, (R+1):(R+Q)], 1,
                          OmegaList[[k]][(R+1):(R+Q)],
                          V[, (S-Q+1):S])
        } else {
            Omega2 <- c()
            O2 <- matrix(nrow=0, ncol=0)
        }
        
        for ( k in 1:ngroups ) {
            
            Uk <- V %*% Olist[[k]]
                
            
            ## Sample sigma^2_k
            if(draw["s2"])
                s2vec[k] <- sampleSigma2(Slist[[k]], Uk, OmegaList[[k]], nvec[k])

            
            ## Sample omegas_k, don't requrie ordered
            if(draw["omega"]) {

                Omega1 <- sampleOmega(Slist[[k]], Uk[, 1:R],
                                         s2vec[k], nvec[k])

                OmegaList[[k]] <- c(Omega1, Omega2)
                
            }
                ## Sample O_k
            if(draw["O"]) {
                
                O1 <- sampleO(Slist[[k]], Uk[, 1:R], s2vec[k],
                                  (OmegaList[[k]])[1:R], V[, 1:(S-Q)])
                                  
                Ok <- as.matrix(bdiag(O1, O2))
            }

            
            Olist[[k]] <- Ok
            Ulist[[k]] <- V %*% Ok
            
            ## save samples        
            if(i %% nskip==0) {
                if(saveU)
                    Usamps[, , k, i/nskip] <- Ulist[[k]]

                Osamps[, , k, i/nskip] <- t(V) %*% Ulist[[k]]
                
                omegaSamps[, k, i/nskip] <- OmegaList[[k]]
                s2samps[k, i/nskip] <- s2vec[k]
            }
            
        }

        if (i%%nskip==0) {
            Vsamps[, , i/nskip] <- V
            if(verbose & !printLoss)
                print(sprintf("Iteration %i", i))
        }
        
        if(i%%nskip==0 & verbose & printLoss) {
            
            sl <- 0
            if(is.null(sigmaTruthInvList)) {
                ## Print loss relative to starting point
                for(k in 1:ngroups) {

                    omegak <- OmegaList[[k]]
                    SigHat <- s2vec[k]*(Ulist[[k]] %*%
                                        diag(omegak/(1-omegak)) %*%
                                        t(Ulist[[k]])+diag(P))

                    SigmaStartInvK <- with(initSS,
                                           1/s2vec[k]*(diag(P) - Ulist[[k]] %*% diag(OmegaList[[k]]) %*% t(Ulist[[k]])))

                    sl <- sl + steinsLoss(SigHat, SigmaStartInvK)
                }
            } else {
                for(k in 1:ngroups) {
                    omegak <- OmegaList[[k]]
                    SigHat <- s2vec[k]*(Ulist[[k]] %*%
                                        diag(omegak/(1-omegak)) %*%
                                        t(Ulist[[k]])+diag(P))

                    sl <- sl + steinsLoss(SigHat, sigmaTruthInvList[[k]])

                }
            }
            print(sprintf("Iteration %i, Loss %f", i, sl/ngroups))
        }
    }

    list(S=S, R=R, Q=Q, V=V, init=init, ngroups=ngroups, Usamps=Usamps,
         Osamps=Osamps, Vsamps=Vsamps, omegaSamps=omegaSamps, s2samps=s2samps)
}

## bayesian single group estimation
fitBayesianSpike <- function(P, R, SC, n, niters=100, nskip=1,
                             init=NULL, SigmaTruth=NULL,
                             verbose=TRUE, ngroups=1) {

    niters <- niters
    nskip <- nskip
    nsamps <- floor(niters/nskip)

    Usamps <- array(dim=c(P, R, nsamps))
    omegaSamps <- matrix(nrow=R, ncol=nsamps)
    s2samps <- numeric(nsamps)

    initItems <- c("U", "omega", "s2")
    if( is.null(init) ) {
        U <- rustiefel(P, R)
        s2 <- rexp(1)
        omega <- rep(1/2, R)
        init <- list(U=U, omega=omega, s2=s2)
        
    } else if (any(!(initItems %in% names(init)))) {
        stop("Some parameters unintialized!")
    }

    SigmaTruthInv <- NULL
    if( !is.null(SigmaTruth) ) {
        SigmaTruthInv <- solve(SigmaTruth)
    }
    
    s2 <- init$s2
    U <- init$U
    omega <- init$omega

    SigInvInit <- getSigmaInv(P, init$U, init$omega, init$s2)
    
    for ( i in 1:niters ) {

        ## Sample sigma^2_k
        
        s2 <- sampleSigma2(SC, U, omega, n)

        A <- SC/(2*s2)
        B <- diag(omega, nrow=length(omega))

        ## order according to eigenvalues (matters when U is full rank)
        ord <- order(omega, decreasing=TRUE)
        revOrd <- order(ord)
        Btilde <- B[ord, ord]

        U <- rbing.matrix.gibbs(A, Btilde, U[, ord])
        U <- U[, revOrd]
        
        omega <- sampleOmega(SC, U, s2, n)

        if(i%%nskip==0) {

            Usamps[, , i/nskip] <- U
            omegaSamps[, i/nskip] <- omega
            s2samps[i/nskip] <- s2
            
            sl <- 0

            if(is.null(SigmaTruthInv)) {

                SigHat <- s2*(U %*% diag(omega/(1-omega)) %*% t(U)+diag(P))
                sl <- steinsLoss(SigHat, SigInvInit)
                
            } else {
                SigHat <- s2*(U %*% diag(omega/(1-omega)) %*% t(U)+diag(P))
                sl <- steinsLoss(SigHat, SigmaTruthInv)

            }
            if(verbose)
                print(sprintf("Iteration %i, Loss %f", i, sl/ngroups))
            
        }


    }

    list(R=R, init=init, Usamps=Usamps, omegaSamps=omegaSamps, s2samps=s2samps)

}



