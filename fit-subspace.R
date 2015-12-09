## init is a list containing initial V, U's, Omega's,
## samples per group (nvec) and Sample covariance (Slist)

fitSubspace <- function(P, S, R, Slist, nvec, ngroups=length(Slist),
                        niters=100, nskip=1, init=NULL, binaryO=FALSE) {


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
    OmegaList <- init$OmegaList
    
    for ( i in 1:niters ) {

        V <- sampleV(Slist, Ulist, s2vec, OmegaList, V, method="hmc")    

        for ( k in 1:ngroups ) {

            Uk <- V %*% Olist[[k]]
            
            ## Sample sigma^2_k
            s2vec[k] <- sampleSigma2(Slist[[k]], Uk, OmegaList[[k]], nvec[k])

            if(binaryO) {
                samp <- proposeBinaryO(S, Uk, V, Sig, s2, Omega, n, flipProb=0.1)

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
        
        
        if(i%%nskip==0) {
            print(i)
        }
    }

    list(S=S, R=R, V=V, init=init, ngroups=ngroups, Usamps=Usamps,
         Osamps=Osamps, Vsamps=Vsamps, omegaSamps=omegaSamps, s2samps=s2samps)
}
