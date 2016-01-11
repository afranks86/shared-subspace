## init is a list containing initial V, U's, Omega's,
## samples per group (nvec) and Sample covariance (Slist)

fitSubspace <- function(P, S, R, Slist, nvec, ngroups=length(Slist),
                        niters=100, nskip=1, init=NULL, binaryO=FALSE,
                        verbose=TRUE, sigmaTruthList=NULL) {
    

    niters <- niters
    nskip <- nskip
    nsamps <- floor(niters/nskip)

    Usamps <- array(dim=c(P, R, ngroups, nsamps))
    Osamps <- array(dim=c(S, R, ngroups, nsamps))
    Vsamps <- array(dim=c(P, S, nsamps))
    omegaSamps <- array(dim=c(R, ngroups, ncol=nsamps))
    s2samps <- matrix(nrow=ngroups, ncol=nsamps)

    initItems <- c("V", "Ulist", "OmegaList")
    if( is.null(init) ) {

        Vinit <- matrix(0,nrow=P, ncol=S)
        Vinit[1:S,1:S] <- diag(S)

        OmegaList <- Olist <- Ulist <- list()
        for(k in 1:ngroups) {
            OmegaList[[k]] <- rep(1/2, R)

            Ok <- matrix(0, nrow=S, ncol=R)
            Ok[1:R, 1:R] <- diag(R)

            Olist[[k]] <- Ok
            Ulist[[k]] <- V %*% Ok
        }
        s2vec <- rexp(ngroups)  

        init <- list(V=Vinit, Ulist=Ulist, OmegaList=OmegaList, s2vec=s2vec)
        
    } else if (any(!(initItems %in% names(init)))) {
        stop("Some parameters unintialized!")
        ## To Finish
    }

    V <- init$V
    Ulist <- init$Ulist
    Olist <- lapply(Ulist, function(u) t(V) %*% u)
    OmegaList <- init$OmegaList

    initSigmaList <- list()
    for(k  in 1:ngroups) {
        omega <- OmegaList[[k]]
        Lambda <- diag((omega/(1-omega)))
        initSigmaList[[k]] <-
            init$s2vec[k] * (Ulist[[k]] %*% Lambda %*% t(Ulist[[k]]) + diag(P))
    }
    
    for ( i in 1:niters ) {

        V <- sampleV(Slist, Ulist, s2vec, OmegaList, V, method="hmc")    

        for ( k in 1:ngroups ) {

            Uk <- V %*% Olist[[k]]
            
            ## Sample sigma^2_k
            s2vec[k] <- sampleSigma2(Slist[[k]], Uk, OmegaList[[k]], nvec[k])
            
            if(binaryO) {

                samp <- proposeBinaryO(S, Uk, V, Slist[[k]], s2vec[k],
                                       OmegaList[[k]], nvec[k], flipProb=0.1)

                Ok <- samp$O
                OmegaList[[k]] <- samp$omega
                
            } else {
                
                ## Sample omegas_k,  do not necessarily maintain order
                OmegaList[[k]] <- sampleOmega(Slist[[k]], Uk, s2vec[k], nvec[k])
                ## Sample O_k
                Ok <- sampleO(Slist[[k]], Uk, s2vec[k], OmegaList[[k]], V)
                
            }

            Olist[[k]] <- Ok
            Ulist[[k]] <- V %*% Ok
            
            ## save samples        
            if(i %% nskip==0) {
                Usamps[, , k, i/nskip] <- Ulist[[k]]
                Osamps[, , k, i/nskip] <- t(V) %*% Ulist[[k]]

                omegaSamps[, k, i/nskip] <- OmegaList[[k]]
                s2samps[k, i/nskip] <- s2vec[k]
            }
            
        }

        if (i%%nskip==0)
            Vsamps[, , i/nskip] <- V    
        
        
        if(i%%nskip==0 & verbose) {
            sl <- 0
            if(is.null(sigmaTruthList)) {
                for(k in 1:ngroups) {
                    SigInv <- getSigmaInv(P, Ulist[[k]], OmegaList[[k]], s2vec[[k]])
                    sl <- sl + steinsLoss(initSigmaList[[k]], SigInv)
                }
            } else {
                for(k in 1:ngroups) {
                    SigInv <- getSigmaInv(P, Ulist[[k]], OmegaList[[k]], s2vec[[k]])
                    sl <- sl + steinsLoss(sigmaTruthList[[k]], SigInv)
                }
            }
            print(sprintf("Iteration %i, Loss %f", i, sl/ngroups))
            
        }
    }

    list(S=S, R=R, V=V, init=init, ngroups=ngroups, Usamps=Usamps,
         Osamps=Osamps, Vsamps=Vsamps, omegaSamps=omegaSamps, s2samps=s2samps)
}

## bayesian single group estimation
fitBayesianSpike <- function(P, S, R, SC, n, niters=100, nskip=1, init=NULL,
                             verbose=FALSE) {

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

    s2 <- init$s2
    U <- init$U
    omega <- init$omega
    
    for ( i in 1:niters ) {

        ## Sample sigma^2_k
        
        s2 <- sampleSigma2(SC, U, omega, n)

        U <- rbing.matrix.gibbs(SC/(2*s2), diag(omega), U)

        omega <- sampleOmega(SC, U, s2, n)
        
        ## save samples        
        if(i %% nskip==0) {
            Usamps[, , i/nskip] <- U
            omegaSamps[, i/nskip] <- omega
            s2samps[i/nskip] <- s2

            if( verbose ) { print(i) }
        }
        
    }
    
    list(S=S, R=R, init=init, Usamps=Usamps, omegaSamps=omegaSamps, s2samps=s2samps)
}

