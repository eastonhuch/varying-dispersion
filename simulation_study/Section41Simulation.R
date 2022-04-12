#######################
# Code for Section 4.1 of "Efficient estimation of the GP-P regression model"
# Simulation for different P, phi, n
# NOTE: Since all data are generated and models are fit sequentially, 
# this code will take several hours to run
#######################

myfilename <- "simulation41"

source('gpp.R')
library(MASS)
library(COUNT)
library(VGAM)

set.seed(1056223)

# Values for model-fitting functions
tol <- 1e-10
max.iter <- 100
max.stephalving <- 10

#maximizing function for "optim" and "nlm"
mygpp_for_optim <- function(myparvec, y, X, P, log=T, omit_constant=T){
	thismu <- exp(X%*%myparvec[1:3])
	-sum(dgpoisP(y, mu=thismu, phi=myparvec[4], P=P, log=T, omit_constant=omit_constant))
}

#different sample sizes to consider
sampsizes <- c(150, 500, 1000) 
nsamps <- length(sampsizes)

#different phis to consider
phis <- c(-0.1, 0.5, 1)
nphis <- length(phis)

#different values of P to consider
Ps <- c(0.26, 1, 2, 3.4)
nPs <- length(Ps)

#How many simulations? 
nreps <- 1000

#Model names
modnames <- c("GP1", "GP2", "GP-0.26", "GP-3.4", "GP-mom", "nlm-T", "optim-T", "vglm-1", "vglm-2", "Pois", "QPois", "NB-1", "NB-2")
nmods <- length(modnames)

#WHAT TO SAVE
comptime <- array(0, dim=c(3,nreps, nmods,nPs, nsamps, nphis))
betahat <- array(0, dim=c(nreps, nmods, nPs, nsamps, nphis))
betalow <- array(0, dim=c(nreps, nmods, nPs, nsamps, nphis))
betahigh <- array(0, dim=c(nreps, nmods, nPs, nsamps, nphis))
niters <- array(0, dim=c(nreps, nmods, nPs, nsamps, nphis)) 
phihat <- array(0, dim=c(nreps, nmods, nPs, nsamps, nphis))
philow <- array(0, dim=c(nreps, nmods, nPs, nsamps, nphis))
phihigh <- array(0, dim=c(nreps, nmods, nPs, nsamps, nphis))


for(thisP in 1:nPs){
	P <- Ps[thisP]
	
	for(thisphi in 1:nphis){	
	  phi <- phis[thisphi] 
		
	  for(thissampsize in 1:nsamps){
		  n <- sampsizes[thissampsize]
		  beta <- 1 #true covariate 
	    
		  for(iter in 1:nreps){
		    
		    #Generate Data
	      x <- runif(n) -0.5 #covariate #centering at 0 => mean of ~e^0 = 1
	      eta <- x*beta #+ 2.3 #(add beta0 to increase mean...doesn't seem to help)
	      mu <- exp(eta) #log(mu) = X%*%beta 
	      x2 <- runif(n)-.5 #nuissance covariate
	      X <- cbind(rep(1, n), x, x2)
	      y <- rgpoisP(n, mu, phi, P) #response variable

	      #starting values for the algorithms (start at the Poisson mle)
        betastart <- c(glm(y~x+x2, family="poisson")$coefficient)
	      phistart <- 0 
	
	      #GP-1 with gpp
	      tim1 <- proc.time()
	      GP1fit <- gpp(y, X, betastart=betastart, phistart=phistart, P=1, tol=tol, max_iter=max.iter, phi_method="joint", stephalving_max=10, verbose=F, penalty=NULL, regularize_intercept=F, offset=NULL)
	      tim2 <- proc.time()
	      comptime[,iter,1,thisP,thissampsize,thisphi] <- (tim2-tim1)[1:3]
	      betahat[iter,1,thisP,thissampsize,thisphi] <- GP1fit$beta[2]	
	      phihat[iter,1,thisP,thissampsize,thisphi] <- GP1fit$phi
	      covpar <- solve(GP2fit$J_star)
	      betalow[iter,1,thisP,thissampsize,thisphi] <- GP1fit$beta[2] + qnorm(.025)*sqrt(covpar[2,2])
	      betahigh[iter,1,thisP,thissampsize,thisphi] <- GP1fit$beta[2] + qnorm(.975)*sqrt(covpar[2,2])
	      philow[iter,1,thisP,thissampsize,thisphi] <- GP1fit$phi + qnorm(.025)*sqrt(covpar[4,4])
	      phihigh[iter,1,thisP,thissampsize,thisphi] <- GP1fit$phi + qnorm(.975)*sqrt(covpar[4,4])
	      niters[iter,1,thisP,thissampsize,thisphi] <- GP1fit$iters
	
	      #GP-2 with gpp
	      tim1 <- proc.time()
	      GP2fit <- gpp(y, X, betastart=betastart, phistart=phistart, P=2, tol=tol, max_iter=max.iter) 
	      tim2 <- proc.time()
	      comptime[,iter,2,thisP,thissampsize,thisphi] <- (tim2-tim1)[1:3]
	      betahat[iter,2,thisP,thissampsize,thisphi] <- GP2fit$beta[2]	
	      phihat[iter,2,thisP,thissampsize,thisphi] <- GP2fit$phi
	      covpar <- solve(GP2fit$J_star)
	      betalow[iter,2,thisP,thissampsize,thisphi] <- GP2fit$beta[2] + qnorm(.025)*sqrt(covpar[2,2])
	      betahigh[iter,2,thisP,thissampsize,thisphi] <- GP2fit$beta[2] + qnorm(.975)*sqrt(covpar[2,2])
	      philow[iter,2,thisP,thissampsize,thisphi] <- GP2fit$phi + qnorm(.025)*sqrt(covpar[4,4])
	      phihigh[iter,2,thisP,thissampsize,thisphi] <- GP2fit$phi + qnorm(.975)*sqrt(covpar[4,4])
	      niters[iter,2,thisP,thissampsize,thisphi] <- GP2fit$iters
	
	      #GP-0.26 with gpp
	      tim1 <- proc.time()
	      GP26fit <- gpp(y, X, betastart=betastart, phistart=phistart, P=0.26, tol=tol, max_iter=max.iter) 
	      tim2 <- proc.time()
	      comptime[,iter,3,thisP,thissampsize,thisphi] <- (tim2-tim1)[1:3]
	      betahat[iter,3,thisP,thissampsize,thisphi] <- GP26fit$beta[2]	
	      phihat[iter,3,thisP,thissampsize,thisphi] <- GP26fit$phi
	      covpar <- solve(GP26fit$J_star)
	      betalow[iter,3,thisP,thissampsize,thisphi] <- GP26fit$beta[2] + qnorm(.025)*sqrt(covpar[2,2])
	      betahigh[iter,3,thisP,thissampsize,thisphi] <- GP26fit$beta[2] + qnorm(.975)*sqrt(covpar[2,2])
	      philow[iter,3,thisP,thissampsize,thisphi] <- GP26fit$phi + qnorm(.025)*sqrt(covpar[4,4])
	      phihigh[iter,3,thisP,thissampsize,thisphi] <- GP26fit$phi + qnorm(.975)*sqrt(covpar[4,4])
	      niters[iter,3,thisP,thissampsize,thisphi] <- GP26fit$iters
	      
	      #GP-3.4 with gpp
	      tim1 <- proc.time()
	      GP34fit <- gpp(y, X, betastart=betastart, phistart=phistart, P=3.4, tol=tol, max_iter=max.iter) 
	      tim2 <- proc.time()
	      comptime[,iter,4,thisP,thissampsize,thisphi] <- (tim2-tim1)[1:3]
	      betahat[iter,4,thisP,thissampsize,thisphi] <- GP34fit$beta[2]	
	      phihat[iter,4,thisP,thissampsize,thisphi] <- GP34fit$phi
	      covpar <- solve(GP34fit$J_star)
	      betalow[iter,4,thisP,thissampsize,thisphi] <- GP34fit$beta[2] + qnorm(.025)*sqrt(covpar[2,2])
	      betahigh[iter,4,thisP,thissampsize,thisphi] <- GP34fit$beta[2] + qnorm(.975)*sqrt(covpar[2,2])
	      philow[iter,4,thisP,thissampsize,thisphi] <- GP34fit$phi + qnorm(.025)*sqrt(covpar[4,4])
	      phihigh[iter,4,thisP,thissampsize,thisphi] <- GP34fit$phi + qnorm(.975)*sqrt(covpar[4,4])
	      niters[iter,4,thisP,thissampsize,thisphi] <- GP34fit$iters
	
	      #GP-P (P at true value), estimate phi with moment estimator
	      #For some values of phi and P, this returned errors so used a tryCatch wrapper so the code would not fail
	      moveon <- F
	      tim1 <- proc.time()
	      tryCatch(GPmoment <- gpp(y, X, betastart=betastart, phistart=phi, P=P, tol=tol, max_iter=max.iter, phi_method="moment", stephalving_max=10, verbose=F, penalty=NULL, regularize_intercept=F, offset=NULL), error=function(ee){moveon<<-T})
	      tim2 <- proc.time()
	      if(moveon==F){
	        comptime[,iter,5,thisP,thissampsize,thisphi] <- (tim2-tim1)[1:3]
	        betahat[iter,5,thisP,thissampsize,thisphi] <- GPmoment$beta[2]	
	        phihat[iter,5,thisP,thissampsize,thisphi] <- GPmoment$phi
	        betalow[iter,5,thisP,thissampsize,thisphi] <- GPmoment$beta[2] + qnorm(.025)*sqrt(solve(GPmoment[[3]])[2,2])
	        betahigh[iter,5,thisP,thissampsize,thisphi] <- GPmoment$beta[2] + qnorm(.975)*sqrt(solve(GPmoment[[3]])[2,2])
	        niters[iter,5,thisP,thissampsize,thisphi] <- GPmoment$iters
	      }else{
		      comptime[,iter,5,thisP,thissampsize,thisphi] <- NA
		      betahat[iter,5,thisP,thissampsize,thisphi] <- NA
		      phihat[iter,5,thisP,thissampsize,thisphi] <- NA
		      betalow[iter,5,thisP,thissampsize,thisphi] <- NA
		      betahigh[iter,5,thisP,thissampsize,thisphi] <- NA
		      niters[iter,5,thisP,thissampsize,thisphi] <- NA
	      }
	      
	      #GP-P (P set to true value) by optimizing via nlm
	      #again had to use tryCatch wrapper so code would not fail
	      moveon <- F
	      tim1 <- proc.time()
	      nlm.fit <- tryCatch(nlm(mygpp_for_optim, p=c(betastart, phistart), hessian=T, gradtol=tol, iterlim=max.iter, y=y, X=X, P=P), error=function(ee){moveon<<-T})
	      tim2 <- proc.time()
	      if(moveon==F){
	        comptime[,iter,6,thisP,thissampsize,thisphi] <- (tim2-tim1)[1:3]
	        betahat[iter,6,thisP,thissampsize,thisphi] <- nlm.fit$estimate[2]
	        phihat[iter,6,thisP,thissampsize,thisphi] <- nlm.fit$estimate[4]
	        niters[iter,6,thisP,thissampsize,thisphi] <- nlm.fit$iterations
	        moveon2 <- F
	        covpar <- tryCatch(solve(nlm.fit$hessian), error=function(ee){moveon2<<-T})
		      if(moveon2==F){
			      betalow[iter,6,thisP,thissampsize,thisphi] <- nlm.fit$estimate[2] + qnorm(.025)*sqrt(covpar[2,2])
			      betahigh[iter,6,thisP,thissampsize,thisphi] <- nlm.fit$estimate[2] + qnorm(.975)*sqrt(covpar[2,2])
			      philow[iter,6,thisP,thissampsize,thisphi] <- nlm.fit$estimate[4] + qnorm(.025)*sqrt(solve(nlm.fit$hessian)[4,4])
			      phihigh[iter,6,thisP,thissampsize,thisphi] <- nlm.fit$estimate[4] + qnorm(.975)*sqrt(covpar[4,4])
			    }else{
			      betalow[iter,6,thisP,thissampsize,thisphi] <- NA
			      betahigh[iter,6,thisP,thissampsize,thisphi] <- NA
			      philow[iter,6,thisP,thissampsize,thisphi] <- NA
			      phihigh[iter,6,thisP,thissampsize,thisphi] <- NA
			    }
	      }else{
		      comptime[,iter,6,thisP,thissampsize,thisphi] <- NA
		      betahat[iter,6,thisP,thissampsize,thisphi] <- NA
		      phihat[iter,6,thisP,thissampsize,thisphi] <- NA
		      betalow[iter,6,thisP,thissampsize,thisphi] <- NA
		      betahigh[iter,6,thisP,thissampsize,thisphi] <- NA
		      philow[iter,6,thisP,thissampsize,thisphi] <- NA
		      phihigh[iter,6,thisP,thissampsize,thisphi] <- NA
		      niters[iter,6,thisP,thissampsize,thisphi] <- NA
	      }
	      
	      #GP-P model with optim (makes use of bounds on phi)
	      # also had to use tryCatch so code wouldn't fail
	      moveon <- F
	      if(P==1 | P==2){
		      tim1 <- proc.time()	
		      optim.fit <- tryCatch(optim(par=c(betastart, phistart), mygpp_for_optim, y=y, X=X, P=2, log=T, omit_constant=T, method="L-BFGS-B", hessian=T, lower=c(-Inf,-Inf,-Inf, -2^{-P})), error=function(ee){moveon<<-T})
		      tim2 <- proc.time()
	      }else{
		      tim1 <- proc.time()
		      optim.fit <- tryCatch(optim(par=c(betastart, phistart), mygpp_for_optim, y=y, X=X, P=2, log=T, omit_constant=T, method="BFGS", hessian=T), error=function(ee){moveon<<-T})
		      tim2 <- proc.time()
	      }
	      if(moveon==F){
		      comptime[,iter,7,thisP,thissampsize,thisphi] <- (tim2-tim1)[1:3]
		      betahat[iter,7,thisP,thissampsize,thisphi] <- optim.fit$par[2]
		      phihat[iter,7,thisP,thissampsize,thisphi] <- optim.fit$par[4]
		      niters[iter,7,thisP,thissampsize,thisphi] <- optim.fit$count[2] # number of times called gr...also a number that counts calls to the fn
		      moveon2 <- F
		      covpar <- tryCatch(solve(optim.fit$hessian), error=function(ee){moveon2<<-T})
		      if(moveon2==F){
			      betalow[iter,7,thisP,thissampsize,thisphi] <- optim.fit$par[2] + qnorm(.025)*sqrt(covpar[2,2])
			      betahigh[iter,7,thisP,thissampsize,thisphi] <- optim.fit$par[2] + qnorm(.975)*sqrt(covpar[2,2])
			      philow[iter,7,thisP,thissampsize,thisphi] <- optim.fit$par[4] + qnorm(.025)*sqrt(covpar[4,4])
			      phihigh[iter,7,thisP,thissampsize,thisphi] <- optim.fit$par[4] + qnorm(.975)*sqrt(covpar[4,4])
		      }else{
			      betalow[iter,7,thisP,thissampsize,thisphi] <- NA
			      betahigh[iter,7,thisP,thissampsize,thisphi] <- NA
			      philow[iter,7,thisP,thissampsize,thisphi] <- NA
			      phihigh[iter,7,thisP,thissampsize,thisphi] <- NA
		      }	
	      }else{
		      comptime[,iter,7,thisP,thissampsize,thisphi] <- NA
		      betahat[iter,7,thisP,thissampsize,thisphi] <- NA
		      phihat[iter,7,thisP,thissampsize,thisphi] <- NA
		      betalow[iter,7,thisP,thissampsize,thisphi] <- NA
		      betahigh[iter,7,thisP,thissampsize,thisphi] <- NA
		      philow[iter,7,thisP,thissampsize,thisphi] <- NA
		      phihigh[iter,7,thisP,thissampsize,thisphi] <- NA
		      niters[iter,7,thisP,thissampsize,thisphi] <- NA	
	      }
	
	      #GP-1 fit via vglm function
	      tim1 <- proc.time()
	      vglm1fit <- vglm(y~x+x2, family='genpoisson1')
	      tim2 <- proc.time()
	      comptime[,iter,8,thisP,thissampsize,thisphi] <- (tim2-tim1)[1:3]
	      betahat[iter,8,thisP,thissampsize,thisphi] <- coef(vglm1fit)[3]
	      betalow[iter,8,thisP,thissampsize,thisphi] <- coef(vglm1fit)[3] + qnorm(.025)*coef(summary(vglm1fit))[3,2]
	      betahigh[iter,8,thisP,thissampsize,thisphi] <- coef(vglm1fit)[3] + qnorm(.975)*coef(summary(vglm1fit))[3,2]
	      phihat[iter,8,thisP,thissampsize,thisphi] <- sqrt(exp(exp(coef(vglm1fit)[4])))-1 #https://rdrr.io/cran/VGAM/man/genpois1UC.html  dispind = (1+phi)^2 and they model the loglog link 
	      philow[iter,8,thisP,thissampsize,thisphi] <- sqrt(exp(exp(coef(vglm1fit)[4] + qnorm(.025)*coef(summary(vglm1fit))[2,2] )))-1
	      phihigh[iter,8,thisP,thissampsize,thisphi] <- sqrt(exp(exp(coef(vglm1fit)[4] + qnorm(.975)*coef(summary(vglm1fit))[2,2] )))-1
	      niters[iter,8,thisP,thissampsize,thisphi] <- vglm1fit@iter 
	      
	      #GP-2 fit via vglm function
	      tim1 <- proc.time()	
	      vglm2fit <- vglm(y~x+x2, family='genpoisson2')
	      tim2 <- proc.time()
	      comptime[,iter,9,thisP,thissampsize,thisphi] <- (tim2-tim1)[1:3]
	      betahat[iter,9,thisP,thissampsize,thisphi] <- coef(vglm2fit)[3]
	      betalow[iter,9,thisP,thissampsize,thisphi] <- coef(vglm2fit)[3] + qnorm(.025)*coef(summary(vglm2fit))[3,2]
	      betahigh[iter,9,thisP,thissampsize,thisphi] <- coef(vglm2fit)[3] + qnorm(.975)*coef(summary(vglm2fit))[3,2]
	      phihat[iter,9,thisP,thissampsize,thisphi] <- exp(coef(vglm2fit)[4]) #https://rdrr.io/cran/VGAM/man/genpois1UC.html  disppar = phi (they call it alpha) and they model the log link 
	      philow[iter,9,thisP,thissampsize,thisphi] <- exp(coef(vglm2fit)[4] + qnorm(.025)*coef(summary(vglm2fit))[2,2])
	      phihigh[iter,9,thisP,thissampsize,thisphi] <- exp(coef(vglm2fit)[4] + qnorm(.975)*coef(summary(vglm2fit))[2,2])
	      niters[iter,9,thisP,thissampsize,thisphi] <- vglm2fit@iter 
	
	      #Poisson model fit via glm
	      tim1 <- proc.time()
	      Poisfit <- glm(y~x+x2, family=poisson)
	      tim2 <- proc.time()
	      comptime[,iter,10,thisP,thissampsize,thisphi] <- (tim2-tim1)[1:3]
	      betahat[iter,10,thisP,thissampsize,thisphi] <- coef(Poisfit)[2]	
	      betalow[iter,10,thisP,thissampsize,thisphi] <- coef(Poisfit)[2]	+ qnorm(.025)*summary(Poisfit)$coefficients[2,2]
	      betahigh[iter,10,thisP,thissampsize,thisphi] <- coef(Poisfit)[2] + qnorm(.975)*summary(Poisfit)$coefficients[2,2]
	      niters[iter,10,thisP,thissampsize,thisphi] <- Poisfit$iter
	
	      #Quasi-Poisson fit via glm
	      tim1 <- proc.time()
	      QPoisfit <- glm(y~x+x2, family=quasipoisson)
	      tim2 <- proc.time()
	      comptime[,iter,11,thisP,thissampsize,thisphi] <- (tim2-tim1)[1:3]
	      betahat[iter,11,thisP,thissampsize,thisphi] <- coef(QPoisfit)[2]
	      betalow[iter,11,thisP,thissampsize,thisphi] <- coef(QPoisfit)[2]	+ qnorm(.025)*summary(QPoisfit)$coefficients[2,2]
	      betahigh[iter,11,thisP,thissampsize,thisphi] <- coef(QPoisfit)[2] + qnorm(.975)*summary(QPoisfit)$coefficients[2,2]	
	      niters[iter,11,thisP,thissampsize,thisphi] <- QPoisfit$iter
	      
	      #NB-1 fit via ml.nb1 
	      #Again needed to use tryCatch
	      moveon <- F
	      tim1 <- proc.time()
	      NB1fit <- tryCatch(ml.nb1(y~x+x2,data=data.frame(y, x, x2)), error=function(ee){moveon<<-T}) #library(COUNT)
	      tim2 <- proc.time()
	      if(moveon==F){
		      comptime[,iter,12,thisP,thissampsize,thisphi] <- (tim2-tim1)[1:3]
		      betahat[iter,12,thisP,thissampsize,thisphi] <- NB1fit[2,1]	
		      betalow[iter,12,thisP,thissampsize,thisphi] <- NB1fit[2,1]	+ qnorm(.025)*NB1fit[2,2]
		      betahigh[iter,12,thisP,thissampsize,thisphi] <- NB1fit[2,1]	+ qnorm(.975)*NB1fit[2,2]
	      }else{
		      comptime[,iter,12,thisP,thissampsize,thisphi] <- NA
		      betahat[iter,12,thisP,thissampsize,thisphi] <- NA	
		      betalow[iter,12,thisP,thissampsize,thisphi] <- NA
		      betahigh[iter,12,thisP,thissampsize,thisphi] <- NA
	      }
	
	      #NB-2 fit via glm.nb
	      tim1 <- proc.time()
	      NB2fit <- glm.nb(y~x+x2) #library(MASS)
	      tim2 <- proc.time()
	      comptime[,iter,13,thisP,thissampsize,thisphi] <- (tim2-tim1)[1:3]
	      betahat[iter,13,thisP,thissampsize,thisphi] <- coef(NB2fit)[2]	
	      betalow[iter,13,thisP,thissampsize,thisphi] <- coef(NB2fit)[2]	+ qnorm(.025)*summary(NB2fit)$coefficients[2,2]
	      betahigh[iter,13,thisP,thissampsize,thisphi] <- coef(NB2fit)[2] + qnorm(.975)*summary(NB2fit)$coefficients[2,2]
	      niters[iter,13,thisP,thissampsize,thisphi] <- NB2fit$iter

      }#end iters
      cat("Finished: n =", thissampsize, "; phi =", thisphi, "\n") #Print progress
    }#end thissampsize
  }#end thisphi
}#end thisP


save(comptime, niters,hess, betahat, phihat, betalow, betahigh, philow, phihigh, sampsizes, phis, P, modnames, file=paste(myfilename,".Rdata", sep=""))


