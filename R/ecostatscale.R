#' Simulate time series for a single species in a single patch
#'
#' Function for simulating dynamics from Eq.1 in the main text.
#' @param r per-capita growth rate (r in Eq.1)
#' @param f the waiting time (or average waiting time) between disturbance events (equal to 1/lambda in Eq.1)
#' @param d mean size of disturbance function (mu in Eq.1)
#' @param d_sd standard deviation of disturbance function (sigma in Eq.1)
#' @param sf waiting time between sampling events
#' @param tmax the time series length to be simulated
#' @param stochd a logical variable, indicating whether disturbance size should be stochastic - otherwise, all disturbances are of magnitude d - defaults to TRUE
#' @param stocht a logical variable, indicating whether waiting time between disturbance events should be stochastic - otherwise, waiting time is always f - defaults to TRUE
#' @param as.matrix indicates whether results should be returned as matrix (potentially faster for some applications) - defaults to FALSE
#' @param oscillate_dist a logical variable indicating whether the sign of the disturbance should oscillate between positive and negative - ignored if stochd==TRUE - defaults to FALSE
#' @return a matrix or data.frame with columns for sampling times, abundances, and number of disturbances for each time interval
#' @export
#' @import stats
#' @examples
#' # see xt2fun

symdyn<-function(r, f, d, d_sd, sf, tmax, stochd=TRUE, stocht=TRUE, as.matrix=FALSE, oscillate_dist=FALSE) {
  st<-seq(0, tmax, by=sf)
  nobs<-length(st)

  datout<-matrix(nrow=length(st), ncol=3)
  colnames(datout)<-c("time", "state", "disturbed")
  datout[,"time"]<-st
  datout[1,"state"]<-0
  datout[,"disturbed"]<-0

  x<-0 #standardized abundance
  tm<-0 #time
  n<-2 #sample position
  md<-1 #number of disturbances

  while(n <= nobs) {
    if(n==2) {
      if(stocht) {
        tdist<-tm+rexp(1, 1/f) #time to next disturbance
      } else {
        tdist<-tm+f
      }
      tsamp<-st[n]
    }

    while((n <= nobs) & (tsamp<tdist)) {
      x<-xt(tm, tsamp, x, r)
      datout[n,"state"]<-x
      tm<-tsamp
      n<-n+1
      if(n <= length(st)) {
        tsamp<-st[n]
      }
    }

    while((n <= nobs) & (tdist<=tsamp)) {
      x<-xt(tm, tdist, x, r)
      tm<-tdist
      if(stochd) {
        rd<-rnorm(1, d, d_sd)
      } else {
        if(!oscillate_dist | (md%%2 == 0)) {
          rd<-d
        } else {
          rd<-(-d)
        }
        md<-md+1
      }
      x<-x+rd
      datout[n,"disturbed"]<-datout[n,"disturbed"]+1

      if(stocht) {
        tdist<-tm+rexp(1, 1/f) #time to next disturbance
      } else {
        tdist<-tm+f
      }
    }
  }

  if(!as.matrix) {
    return(data.frame(datout))
  } else {
    return(datout)
  }
}


#' Simulate deterministic dynamics of a single species
#'
#' Helper function for symdyn, used to simulate dynamics between disturbance events.
#' @param t0 initial time step
#' @param t1 desired time step
#' @param B0 initial abundance
#' @param r relative growth rate
#' @return predicted value of x at time t1
#' @export
#' @examples
#' # see xt2fun

xt <- function(t0, t1, B0, r) {
  dt<-(t1-t0)
  return(exp(-r*dt)*B0)
}


#' Simulate deterministic dynamics of a N competing species
#'
#' Helper function for symdynN, used to simulate dynamics between disturbance events.
#' @param time time step - ignored, but required for consistency with ode function
#' @param state vector of current states
#' @param pars parameter list, including matrix A with full interaction matrix (and growth rates along the diagonal)
#' @return rate of change for each species
#' @export

df0<-function(time, state, pars) {
  list(c(t(as.matrix(state))%*%pars$A))
}

#' Simulate deterministic dynamics of N patches with dispersal
#'
#' Helper function for symdynN, used to simulate dynamics between disturbance events.
#' Note that K is set to 1 for all species.
#' @param time time step - ignored, but required for consistency with ode function
#' @param state vector of current states
#' @param pars parameter list, including matrix A with full interaction matrix (and growth rates along the diagonal),
#' and a value Ifrac, which is the dispersal rate (D in the main text), and Ksim, which is the carrying capacity
#' @return rate of change for each species
#' @export

df_col<-function(time, state, pars) {
  list(diag(pars$A)*state-pmax(pars$Ifrac*(state+pars$Ksim),0)+pars$Ifrac*mean(pmax(state+pars$Ksim, 0)))
}

#' Simulate deterministic dynamics of N patches with dispersal and loss
#'
#' Helper function for symdynN, used to simulate dynamics between disturbance events.
#' Note that K is set to 1 for all species.
#' @param time time step - ignored, but required for consistency with ode function
#' @param state vector of current states
#' @param pars parameter list, including matrix A with full interaction matrix (and growth rates along the diagonal),
#' and a value Ifrac, which is the dispersal rate (D in the main text), and Ksim, which is the carrying capacity, and Iloss
#' which is the loss rate.
#' @return rate of change for each species
#' @export

df_col_loss<-function(time, state, pars) {
  list(diag(pars$A)*state-pmax((pars$Ifrac+pars$Iloss)*(state+pars$Ksim),0)+mean(pars$Ifrac*pmax(state+pars$Ksim, 0)))
}


#' Simulate deterministic dynamics
#'
#' Helper function for symdynN. Simulate an ODE given parameters, starting value, and times.
#' @param t0 initial time step
#' @param t1 desired time step
#' @param B0 initial abundance
#' @param odepars parameter list
#' @param dffun function for calculating derivatives
#' @param nsteps number of time steps to return - defaults to 2
#' @return a matrix of species abundances
#' @export

xtN <- function(t0, t1, B0, odepars, dffun, nsteps=2) {
  out<-ode(y=B0, times=seq(t0, t1, length=nsteps), parms=odepars, func = dffun)
  out[-1,-1]
}


#' Simulate time series for N species or patches
#'
#' Function for simulating dynamics from Eq.2-3 in the main text.
#' @param r per-capita growth rate (r in Eq.2-3)
#' @param amu the mean interaction strength
#' @param asd standard deviation used to generate interaction strengths
#' @param f the waiting time (or average waiting time) between disturbance events (equal to 1/lambda in Eq.2-3)
#' @param d mean size of disturbance function (mu in Eq.2-3)
#' @param d_sd standard deviation of disturbance function (sigma in Eq.5)
#' @param d_cov the covariance for generating disturbances
#' @param N number of species or patches
#' @param sf waiting time between sampling events
#' @param tmax the time series length to be simulated
#' @param stochd a logical variable, indicating whether disturbance size should be stochastic - otherwise, all disturbances are of magnitude d - defaults to TRUE
#' @param stocht a logical variable, indicating whether waiting time between disturbance events should be stochastic - otherwise, waiting time is always f - defaults to TRUE
#' @param as.matrix indicates whether results should be returned as matrix (potentially faster for some applications) - defaults to FALSE
#' @param amax the maximum value allowed for interaction coefficients - defaults to zero
#' @param amin the minimum value allowed for interaction coefficients - defaults to -Inf
#' @param Ifrac dispersal rate (D in Eq. 2) - defaults to NULL (i.e. no dispersal)
#' @param Iloss loss rate from dispersal that falls outside of the patch - defaults to NULL
#' @param dffun the function handed to the ODE solver - should be df_col for spatial simulations, and df0 for multi-species simulations - defaults to df0
#' @param fullout a logical, determining whether the full output or just a summary is returned - defaults to fullout
#' @param xstart  optional vector of starting abundances - defaults to NULL (i.e. no values)
#' @param Ksim carrying capacities - defaults to 1
#' @return a matrix or data.frame with columns for sampling times, abundances, and number of disturbances for each time interval
#' @export
#' @import deSolve
#' @import stats
#' @examples
#' ### Example 1: 10 patches
#' r<-1 #rate of recovery
#' d<-(0) #mean size of disturbance (mu in text)
#' d_sd<-sqrt(0.1) #SD of disturbances (sigma in text)
#' f<-1 #average time between disturbances (1/lambda in text)
#' sf<-0.1 #sampling interval
#' tmax<-120 #maximum time for simulation
#' d_cov<-d_cov0<-(d_sd)^2/2 #covariance in disturbance size among patches
#'
#' xtNpatches<-symdynN(r = r, amu=0, asd=0, f=f, d=d,
#'               d_sd=d_sd, d_cov=d_cov, N=10,
#'               sf=sf, tmax=tmax, Ifrac=0, dffun = df_col)
#'
#'
#' ### Example 2: 30 species
#' r<-1 #rate of recovery
#' d<-(0) #mean size of disturbance (mu in text)
#' d_sd<-sqrt(0.1) #SD of disturbances (sigma in text)
#' f<-1 #average time between disturbances (1/lambda in text)
#' sf<-0.1 #sampling interval
#' tmax<-120 #maximum time for simulation
#' d_cov<-0 #covariance in disturbances among species
#' amu<-(-r/2) #average interaction coefficient
#' asd<-0.1 #standard deviation of interaction coefficient
#'
#' xtNsp<-symdynN(r = r, amu=amu, asd=asd, f=f, d=d,
#'              d_sd=d_sd, d_cov=d_cov, N=30,
#'              sf=sf, tmax=tmax)

symdynN<-function(r, amu, asd, f, d, d_sd, d_cov, N, sf, tmax, stochd=TRUE, stocht=TRUE, as.matrix=FALSE, amax=0, amin=-Inf, Ifrac=NULL, Iloss = NULL, dffun=df0, fullout=FALSE, xstart=NULL, Ksim = 1) {
  st<-seq(0, tmax, by=sf)
  nobs<-length(st)

  datout<-matrix(nrow=length(st), ncol=2+N)
  colnames(datout)<-c("time", "disturbed", paste("N", 1:N, sep="_"))
  datout[,"time"]<-st
  datout[,-1]<-0
  sppos<-1:N+2

  if(is.null(xstart)) {
    x<-rep(0, N) #standardized abundance
  } else {
    x<-xstart
    datout[1,-c(1:2)]<-x
  }
  tm<-0 #time
  n<-2 #sample position
  m<-1 #disturbance position

  #disturbance times
  if(stocht) {
    dtime<-cumsum(rexp(round(2*(tmax/f)), 1/f))
    mm<-3
    while(max(dtime)<tmax) {
      dtime<-cumsum(rexp(round(mm*(tmax/f)), 1/f))
      mm<-mm+1
    }
    dtime<-dtime[dtime<=tmax]
  } else {
    dtime<-seq(f, tmax, by=f)
  }

  #disturbance quantities
  if(stochd) {
    covmat<-diag(N)*d_sd^2
    covmat[row(covmat)!=col(covmat)]<-d_cov
    dquant<-rmvnorm(length(dtime), mean=rep(d, N), sigma = covmat)
  } else {
    dquant<-matrix(nrow=length(dtime), ncol=N, data=d)
  }

  #interaction matrix
  A<-(-diag(N)*r)
  ps<-which(row(A)!=col(A))
  A[ps]<-rnorm(N^2-N, amu, asd)
  while(any(A[ps]<amin) | any(A[ps]>amax)) {
    A[ps][A[ps]<amin]<-rnorm(sum(A[ps]<amin), amu, asd)
    A[ps][A[ps]>amax]<-rnorm(sum(A[ps]>amax), amu, asd)
  }

  odepars<-list(A=A, Ifrac=Ifrac, Iloss=Iloss, Ksim=Ksim)

  while(n <= nobs) {
    if(n==2) {
      tdist<-dtime[m] #time to next disturbance
      tsamp<-st[n]
    }

    while((n <= nobs) & (tsamp<tdist)) {
      x<-xtN(tm, tsamp, x, odepars, dffun)
      datout[n,sppos]<-x
      tm<-tsamp
      n<-n+1
      if(n <= length(st)) {
        tsamp<-st[n]
      }
    }

    while((n <= nobs) & (tdist<=tsamp)) {
      x<-xtN(tm, tdist, x, odepars, dffun)
      tm<-tdist

      x<-x+dquant[m,]
      datout[n,"disturbed"]<-datout[n,"disturbed"]+1

      m<-m+1
      if(m<=length(dtime)) {
        tdist<-dtime[m]
      } else {
        tdist<-tmax+1
      }
    }
  }

  if(fullout) {
    return(list(datout=datout, A=A, dquant=dquant))
  } else if(!as.matrix) {
    return(data.frame(datout))
  } else {
    return(datout)
  }
}







#' Unbiased stability paramter estimation
#'
#' Function for solving for stability paramter values from observed time series.
#' Equivalent to Eq.5 in the main text.
#' @param x0 value of x^2 at time t (x(t) in Eq.5)
#' @param r per-capita growth rate (r in Eq.5)
#' @param d mean size of disturbance function (mu in Eq.5)
#' @param d_sd standard deviation of disturbance function (sigma in Eq.5)
#' @param dt time step (i.e. time between x0 and x1) - can be a vector of the same length as x0, or a number if all time steps are of equal length
#' @param ndist number of disturbances in each time step (equivalent to p(t+tau) in Eq.5) - must be same length as x0
#' @return predicted value of x^2 at time t+dt
#' @export
#' @examples
#' # simulate dynamics, with r=1, d=0, and d_sd=0.1
#' xtout<-symdyn(r=1, f=1, d=0, d_sd=0.1, sf=0.1, tmax=100)
#'
#' # abundance in current time step
#' x0<-xtout$state[1:(nrow(xtout)-1)]
#' # abundance at t+1
#' x1<-xtout$state[2:nrow(xtout)]
#'
#' dt<-diff(xtout$time)
#' ndist<-xtout$disturbed[-1]
#'
#' # fit model - note square root transform of response variable,
#' # and log transform of parameter values
#'
#' mod<-nls(sqrt(x1^2)~sqrt(xt2fun(x0, r=exp(lr), d=0, d_sd=exp(ld_sd), dt, ndist)),
#'          start=c(lr=log(1), ld_sd=log(0.1)))
#' exp(coef(mod)) # model estimates

xt2fun<-function(x0, r, d, d_sd, dt, ndist) {
  if(length(dt)==1) {
    dt<-rep(dt, length(x0))
  }
  if(length(r)==1) {
    r<-rep(r, length(x0))
  }
  if(length(d_sd)==1) {
    d_sd<-rep(d_sd, length(x0))
  }

  tstep<-dt/(ndist+1)
  stmp1<-(x0*exp(-r*tstep))^2
  if(max(ndist)>0) {
    for(k in 1:max(ndist)) {
      ps<-which(ndist>=k)
      stmp1[ps]<-(stmp1[ps]+2*d*sqrt(stmp1[ps])+d^2+d_sd[ps]^2)*exp(-2*r[ps]*tstep[ps])
    }
  }

  return(stmp1)
}


#' Variance scaling function
#'
#' Extrapolate variance observed at spatial or ecological scale b
#' to a different scale, B. Equivalent to Eq.7a in the main text.
#' @param mvar_b Mean variance of abundance values observed at scale b.
#' @param murho_b Mean Pearson correlation coefficient of abundance values observed at scale b, calculated as mucov_b/mvar_b. If NULL, mucov_b is used instead.
#' @param mucov_b Mean covariance of abundances observed at scale b. Ignored if mrho_b is not NULL. Defaults to NULL.
#' @param b Size of observed scale. Defaults to 1.
#' @param B Size of desired scale for extrapolation.
#' @return Extrapolated variance at scale B.
#' @import mvtnorm
#' @export
#' @examples
#' # extrapolate from scale of 1 to 10 - e.g. from a 1m2 patch to a 10m2 patch
#' var_scale(mvar_b = 1, murho_b = 0.5, b = 1, B = 10)
#'
#' # example with 100 simulated species
#' nsp<-100 # number of species
#' var_b<-1 # species-level abundance variance
#' cov_b<-(-0.01) # between-specie abundance covariance
#' # note - if nsp is large, cov_b must be near zero
#' # this is because, e.g. many variables cannot all be
#' # simultaneously negatively correlated
#'
#' # make a covariance matrix based on var_b and cov_b
#' sigmamat<-diag(nsp)*var_b+(1-diag(nsp))*cov_b
#' # simulate 1000 observations of 100 species
#' sim_x<-mvtnorm::rmvnorm(n=1e3, mean = rep(0,100), sigma = sigmamat)
#'
#' # calculate mean variance, covariance, and correlation from sim_x
#' cvmat<-cov(sim_x)
#' mvar_b<-mean(diag(cvmat))
#' mucov_b<-mean(cvmat[row(cvmat)!=col(cvmat)])
#' murho_b<-mucov_b/mvar_b
#'
#' # test function vs. observation
#' # note - answers match exactly
#' var(rowSums(sim_x))
#' var_scale(mvar_b, murho_b = murho_b, b=1, B=100)

var_scale<-function(mvar_b, murho_b, mucov_b=NULL, b=1, B) {
  if(is.null(murho_b)) {
    murho_b<-mucov_b/mvar_b
  }

  return(mvar_b*(B/b)*(1 + murho_b*((B/b) - 1)))
}


#' Sigma scaling function
#'
#' Extrapolate disturbance standard deviation observed at spatial or ecological scale b
#' to a different scale, B (inversely related to resistance). Equivalent to Eq.7b in the main text.
#' @param msd_b Mean disturbance standard deviation observed at scale b.
#' @param murho_b Mean Pearson correlation coefficient of disturbances observed at scale b, calculated as mucov_b/mvar_b. If NULL, mucov_b is used instead.
#' @param mucov_b Mean covariane of disturbances observed at scale b. Ignored if mrho_b is not NULL. Defaults to NULL.
#' @param b Size of observed scale. Defaults to 1.
#' @param B Size of desired scale for extrapolation.
#' @return Extrapolated disturbance standard deviation at scale B.
#' @export
#' @examples
#' #extrapolate from scale of 1 to 10 - e.g. from a 1m2 patch to a 10m2 patch
#' sd_scale(msd_b = 1, murho_b = 0.5, b = 1, B = 10)

sd_scale<-function(msd_b, murho_b, mucov_b=NULL, b=1, B) {
  if(is.null(murho_b)) {
    murho_b<-mucov_b/(msd_b^2)
  }

  return(msd_b*sqrt((B/b)*(1 + murho_b*((B/b) - 1))))
}


#' Resilience scaling function
#'
#' Extrapolate resilience observed at the scale of a single spatial or ecological scale b (e.g. a patch or species)
#' to a larger scale, B (e.g. functional group or landscape).
#' @param mvar_b Mean abundance variance observed at scale b
#' @param murho_b_abundance Mean Pearson correlation coefficient of abundance values observed at scale b
#' @param mucov_b_abundance Mean covariance of abundance values observed at scale b. Ignored unless murho_b_abundance is NULL Defaults to NULL.
#' @param msd_b Mean disturbance standard deviation observed at scale b
#' @param murho_b_disturbance Mean Pearson correlation coefficient of disturbances observed at scale b
#' @param mucov_b_disturbance Mean covariance of disturbances observed at scale b. Ignored unless murho_b_abundance is NULL Defaults to NULL.
#' @param b Size of observed scale. Defaults to 1.
#' @param B Larger scale being extrapolated to (e.g. total number of species, or size of patch B relative to b)
#' @param lambda Mean disturbance frequency.
#' @return Extrapolated median resilience at scale of M species.
#' @export
#' @examples
#' # extrapolate from scale of 1 species to 10 species
#' res_scale(mvar_b = 0.25, murho_b_abundance = -0.034, msd_b = sqrt(0.1),
#'            murho_b_disturbance = 0, B = 30, lambda=1)
#'
#' # plot relationship for groups of 1 to 30 species
#' plot(1:30, res_scale(mvar_b = 0.25, murho_b_abundance = -0.034, msd_b = sqrt(0.1),
#'       murho_b_disturbance = 0, B = 1:30, lambda=1),
#'       xlab="ecological scale", ylab="resilience, r", type="b")

res_scale<-function(mvar_b, murho_b_abundance, mucov_b_abundance=NULL,
                    msd_b, murho_b_disturbance, mucov_b_disturbance=NULL,
                    b = 1, B, lambda) {
  B = B/b
  if(is.null(murho_b_abundance)) {
    murho_b_abundance<-mucov_b_abundance/mvar_b
  }
  if(is.null(murho_b_disturbance)) {
    murho_b_disturbance<-mucov_b_disturbance/(msd_b^2)
  }

  return((lambda/2)*(msd_b^2)*(1+murho_b_disturbance*(B-1))/
           (mvar_b*(1+murho_b_abundance*(B-1))))
}
