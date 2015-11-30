## TODO
fitSubspace <- function(init, niters=100, nskip=1, outfile="fit.RData") {

    niters <- niters
    nskip <- nskip
    nsamps <- floor(niters/nskip)

    Usamps <- array(dim=c(P, R, ngroups, nsamps))
    Osamps <- array(dim=c(S, R, ngroups, nsamps))
    Vsamps <- array(dim=c(P, S, nsamps))
    omegaSamps <- array(dim=c(R, ngroups, ncol=nsamps))
    s2samps <- matrix(nrow=ngroups, ncol=nsamps)

    defined <- c("V", "Ulist", "OmegaList", "nvec", "Slist") %in% names(init)
    names(defined) <- c("V", "Ulist", "OmegaList", "nvec", "Slist")
    if(any(!defined)) {
        stop(sprintf("%s not initialized!",
                     paste(names(defined)[!defined], collapse=", ")))
    }
    

    V <- init$V
    Ulist <- init$Ulist
    Slist <- init$Slist
    OmegaList <- init$OmegaList
    nvec <- init$nvec
    
    for ( i in 1:niters ) {

        V <- sampleV(Slist, Ulist, s2vec, OmegaList, V, method="hmc")        
        
        for ( k in 1:ngroups ) {
            
            ## Sample sigma^2_k
            s2vec[k] <- rs2_fc(Slist[[k]], Ulist[[k]], OmegaList[[k]], nvec[k])
            
            ## Sample omegas_k,  do not necessarily maintain order
            OmegaList[[k]] <- sampleOmega(Slist[[k]], Ulist[[k]], s2vec[k], nvec[k])
            
            ## Sample O_k
            Ok <- rU_csm_gibbs(Slist[[k]], Ulist[[k]], s2vec[k], OmegaList[[k]], V)
            
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

    save(S, R, V, init, ngroups, Usamps, Osamps, Vsamps, 
         omegaSamps, s2samps, file=outfile)
    
}
