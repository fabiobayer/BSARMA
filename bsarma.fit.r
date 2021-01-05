# Implemented by Fabio M Bayer (bayer@ufsm.br) 
# Date: 06/2015

barma.fit<- function (y, ar, ma, AR, MA, link, 
                      names_phi,names_PHI,names_theta,names_THETA,diag,h1,resid,ljunbox_lag)
{
  #auxiliar functions
  operator<-function(phi,PHI,ar,AR)
  {
    parameters<- c(phi,PHI)
    index<- c(ar,AR)
    
    j1<-1
    for(j in ar)
    {
      J1<-1
      for(J in AR)
      {      
        parameters<- c(parameters, -phi[j1]*PHI[J1])
        index<- c(index, (j+J))
        
        J1<-J1+1
      }
      j1<-j1+1
    }
    
    z<-c()
    z$parameters<-parameters
    z$index<-index
    return(z)
    
  }
  Monti.test<-function (x, lag = 1, type = c("Box-Pierce", "Ljung-Box"), fitdf = 0) 
  {
    if (NCOL(x) > 1) 
      stop("x is not a vector or univariate time series")
    DNAME <- deparse(substitute(x))
    type <- match.arg(type)
    cor <- pacf(x, lag.max = lag, plot = FALSE, na.action = na.pass)
    n <- sum(!is.na(x))
    PARAMETER <- c(df = lag - fitdf)
    obs <- cor$acf[1:(lag )]
    if (type == "Box-Pierce") {
      METHOD <- "Box-Pierce test"
      STATISTIC <- n * sum(obs^2)
      PVAL <- 1 - pchisq(STATISTIC, lag - fitdf)
    }
    else {
      METHOD <- "Box-Ljung test"
      STATISTIC <- n * (n + 2) * sum(1/seq.int(n - 1, n - lag) * 
                                       obs^2)
      PVAL <- 1 - pchisq(STATISTIC, lag - fitdf)
    }
    names(STATISTIC) <- "X-squared"
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
                   p.value = PVAL, method = METHOD, data.name = DNAME), 
              class = "htest")
  }
  
  
  maxit1<-1000
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  mu.eta <-  link$mu.eta
  diflink <- link$diflink
  
  ynew = linkfun(y)
  ystar = log(y/(1-y))
  
  S<-frequency(y)
  p <- max(ar)
  q <- max(ma)
  P <- max(AR)
  Q <- max(MA)
  n <- length(y)
  m <- max(p,q,S*P,S*Q,S*P+p,S*Q+q,na.rm=T) 
  p1 <- length(ar)
  q1 <- length(ma)
  P1 <- length(AR)
  Q1 <- length(MA)
  
  
  y_prev <- c(rep(NA,(n+h1)))
  
  if(any(is.na(c(ar))==F)) 
  {
    A <- matrix(rep(NA,(n-m)*(p1)),ncol=(p1))
    
    for(i in 1:(n-m))
    {
      A[i,] <- ynew[i+m-ar]
    }
    
    Z <- cbind(rep(1,(n-m)),A)
  }else{
    Z <- as.matrix(rep(1,(n-m)))
  }
  
  if(any(is.na(c(AR))==F)) 
  {
    As <- matrix(rep(NA,(n-m)*(P1)),ncol=(P1))
    
    for(i in 1:(n-m))
    {
      As[i,] <- ynew[i+m-(S*AR)]
    }
    Z <- cbind(Z,As)
  }
  
  x <- as.matrix(Z)
  Y <- y[(m+1):n]
  Ynew = linkfun(Y)
  Ystar = log(Y/(1-Y))
  ajuste = lm.fit(x, Ynew)
  mqo = c(ajuste$coef)
  k = length(mqo)
  n1 = length(Y)
  mean = fitted(ajuste)
  mean = linkinv(mean)
  dlink = diflink(mean)
  er = residuals(ajuste)
  sigma2 = sum(er^2)/((n1 - k) * (dlink)^2)
  prec = 1/n1 * sum(mean * (1 - mean)/sigma2 - 1)
  
  
  #################################################################
  ######### BSARMA complete model
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(AR)==F) && any(is.na(MA)==F))
  { 
    if(diag>0)print("BSARMA model",quote=F)
    reg <- c(mqo, rep(0,(q1+Q1)), prec) # initializing the parameter values
    
    loglik <- function(z) 
    {
      alpha <- z[1]
      phi = z[2:(p1+1)] 
      PHI = z[(p1+2):(p1+P1+1)]
      theta = z[(p1+P1+2):(p1+P1+q1+1)]
      THETA = z[(p1+P1+q1+2):(p1+P1+q1+Q1+1)]
      prec <- z[(p1+P1+q1+Q1+2)] # precision parameter
      
      error<-rep(0,n) # E(error)=0 
      eta<-rep(NA,n)
      
      ar_par_index <- operator(phi,PHI,ar,S*AR)
      ma_par_index <- operator(theta,THETA,ma,S*MA)
      
      ar_par <- ar_par_index$parameters
      ar_ind <- ar_par_index$index
      
      ma_par <- ma_par_index$parameters
      ma_ind <- ma_par_index$index
      
      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + sum(ar_par*ynew[i-ar_ind]) - sum(ma_par*error[i-ma_ind])
        error[i] <- ynew[i]-eta[i] 
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings(dbeta(y1, mu * prec, (1 - mu) * prec, log = TRUE))
      sum(ll)
    } 
    
    escore <- function(z) 
    {
      
      alpha <- z[1]
      phi = z[2:(p1+1)] 
      PHI = z[(p1+2):(p1+P1+1)]
      theta = z[(p1+P1+2):(p1+P1+q1+1)]
      THETA = z[(p1+P1+q1+2):(p1+P1+q1+Q1+1)]
      prec <- z[(p1+P1+q1+Q1+2)] # precision parameter
      
      error<-rep(0,n) # E(error)=0 
      eta<-rep(NA,n)
      
      ar_par_index <- operator(phi,PHI,ar,S*AR)
      ma_par_index <- operator(theta,THETA,ma,S*MA)
      
      ar_par <- ar_par_index$parameters
      ar_ind <- ar_par_index$index
      
      ma_par <- ma_par_index$parameters
      ma_ind <- ma_par_index$index
      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + sum(ar_par*ynew[i-ar_ind]) - sum(ma_par*error[i-ma_ind])
        error[i]<- ynew[i]-eta[i] 
      }
      
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      ystar <- log(y1/(1-y1))
      mustar <- digamma(mu * prec) - digamma((1 - mu) * prec)
      
      mT <- diag(mu.eta(mu))
      
      R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
      for(i in 1:(n-m))
      {
        for(j in 1:q1)
        {
          R[i,j] <- -sum(error[i+m-c(ma[j],ma[j]+(S*MA))]*c(1,-THETA))
        }
      }
      Rs <- matrix(rep(NA,(n-m)*Q1),ncol=Q1)
      for(i in 1:(n-m))
      {
        for(j in 1:Q1)
        {
          Rs[i,j] <- -sum(error[i+m-c(S*MA[j],S*MA[j]+(ma))]*c(1,-theta))
        }
      }
      
      A <- matrix(rep(NA,(n-m)*(p1)),ncol=(p1))
      for(i in 1:(n-m))
      {
        for(j in 1:p1)
        {
          A[i,j] <- sum(ynew[i+m-c(ar[j],ar[j]+(S*AR))]*c(1,-PHI))
        }
      }    
      
      As <- matrix(rep(NA,(n-m)*(P1)),ncol=(P1))
      for(i in 1:(n-m))
      {
        for(j in 1:P1)
        {
          As[i,j] <- sum(ynew[i+m-c(S*AR[j],S*AR[j]+(ar))]*c(1,-phi))
        }
      }    
      
      deta.dalpha <- rep(0,n)
      deta.dphi <- matrix(0, ncol=p1,nrow=n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      deta.dPHI <- matrix(0, ncol=P1,nrow=n)
      deta.dTHETA <- matrix(0, ncol=Q1,nrow=n)
      
      for(i in (m+1):n)
      {
        deta.dalpha[i]<- 1 + theta%*%deta.dalpha[i-ma] + THETA%*%deta.dalpha[i-S*MA]- ma_par%*%deta.dalpha[i-ma_ind]
        
        deta.dphi[i,]<- A[(i-m),] + theta%*%deta.dphi[i-ma,]+ THETA%*%deta.dphi[i-S*MA,]- ma_par%*%deta.dphi[i-ma_ind,]
        
        deta.dPHI[i,]<- As[(i-m),] + theta%*%deta.dPHI[i-ma,]+ THETA%*%deta.dPHI[i-S*MA,]- ma_par%*%deta.dPHI[i-ma_ind,]
        
        deta.dtheta[i,]<- R[(i-m),] + theta%*%deta.dtheta[i-ma,]+ THETA%*%deta.dtheta[i-S*MA,]- ma_par%*%deta.dtheta[i-ma_ind,]
        
        deta.dTHETA[i,]<- Rs[(i-m),] + theta%*%deta.dTHETA[i-ma,]+ THETA%*%deta.dTHETA[i-S*MA,]- ma_par%*%deta.dTHETA[i-ma_ind,]
      }
      
      a <- matrix(deta.dalpha[(m+1):n],ncol=1)
      rP <- matrix(deta.dphi[(m+1):n,],ncol=p1)
      rPs <- matrix(deta.dPHI[(m+1):n,],ncol=P1)
      rR <- matrix(deta.dtheta[(m+1):n,],ncol=q1)
      rRs <- matrix(deta.dTHETA[(m+1):n,],ncol=Q1)
      
      ymstar <- matrix((ystar-mustar),ncol=1)
      
      Ualpha <- prec * t(a) %*% mT %*% ymstar
      Uphi <-   prec * t(rP) %*% mT %*% ymstar
      UPHI <-   prec * t(rPs) %*% mT %*% ymstar
      Utheta <- prec * t(rR) %*% mT %*% ymstar
      UTHETA <- prec * t(rRs) %*% mT %*% ymstar
      Uprec <-  sum(mu * (ystar-mustar) + log(1-y1)
                    - digamma((1 - mu) * prec) + digamma(prec) )
      
      rval <- c(Ualpha,Uphi,UPHI,Utheta,UTHETA,Uprec)
    }
    names_par <- c("alpha",names_phi,names_PHI,names_theta,names_THETA,"precision")
    
    opt <- optim(reg, loglik, escore, method = "BFGS", 
                 control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    
    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")
    
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+q1+P1+Q1+2)]
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <- coef[1]
    phi = coef[2:(p1+1)] 
    PHI = coef[(p1+2):(p1+P1+1)]
    theta = coef[(p1+P1+2):(p1+P1+q1+1)]
    THETA = coef[(p1+P1+q1+2):(p1+P1+q1+Q1+1)]
    prec <- coef[(p1+P1+q1+Q1+2)] # precision parameter
    
    z$alpha <- alpha
    z$phi <- phi
    z$PHI <- PHI
    z$theta <- theta
    z$THETA <- THETA
    z$prec <- prec
    
    errorhat<-rep(0,n) # E(error)=0 
    etahat<-rep(NA,n)
    
    ar_par_index <- operator(phi,PHI,ar,S*AR)
    ma_par_index <- operator(theta,THETA,ma,S*MA)
    
    ar_par <- ar_par_index$parameters
    ar_ind <- ar_par_index$index
    
    ma_par <- ma_par_index$parameters
    ma_ind <- ma_par_index$index
    
    for(i in (m+1):n)
    {
      etahat[i] <- alpha + (ar_par%*%ynew[i-ar_ind]) - (ma_par%*%errorhat[i-ma_ind])
      errorhat[i] <- ynew[i]-etahat[i] # predictor scale
      
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    z$mustarhat <- digamma(muhat * prec) - digamma((1 - muhat) * prec)
    
    R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m))
    {
      for(j in 1:q1)
      {
        R[i,j] <- -sum(errorhat[i+m-c(ma[j],ma[j]+(S*MA))]*c(1,-THETA))
      }
    }
    Rs <- matrix(rep(NA,(n-m)*Q1),ncol=Q1)
    for(i in 1:(n-m))
    {
      for(j in 1:Q1)
      {
        Rs[i,j] <- -sum(errorhat[i+m-c(S*MA[j],S*MA[j]+(ma))]*c(1,-theta))
      }
    }
    A <- matrix(rep(NA,(n-m)*(p1)),ncol=(p1))
    for(i in 1:(n-m))
    {
      for(j in 1:p1)
      {
        A[i,j] <- sum(ynew[i+m-c(ar[j],ar[j]+(S*AR))]*c(1,-PHI))
      }
    }    
    As <- matrix(rep(NA,(n-m)*(P1)),ncol=(P1))
    for(i in 1:(n-m))
    {
      for(j in 1:P1)
      {
        As[i,j] <- sum(ynew[i+m-c(S*AR[j],S*AR[j]+(ar))]*c(1,-phi))
      }
    }    
    
    deta.dalpha <- rep(0,n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    deta.dPHI <- matrix(0, ncol=P1,nrow=n)
    deta.dTHETA <- matrix(0, ncol=Q1,nrow=n)
    
    for(i in (m+1):n)
    {
      deta.dalpha[i]<- 1 + theta%*%deta.dalpha[i-ma] + THETA%*%deta.dalpha[i-S*MA]- ma_par%*%deta.dalpha[i-ma_ind]
      
      deta.dphi[i,]<- A[(i-m),] + theta%*%deta.dphi[i-ma,]+ THETA%*%deta.dphi[i-S*MA,]- ma_par%*%deta.dphi[i-ma_ind,]
      
      deta.dPHI[i,]<- As[(i-m),] + theta%*%deta.dPHI[i-ma,]+ THETA%*%deta.dPHI[i-S*MA,]- ma_par%*%deta.dPHI[i-ma_ind,]
      
      deta.dtheta[i,]<- R[(i-m),] + theta%*%deta.dtheta[i-ma,]+ THETA%*%deta.dtheta[i-S*MA,]- ma_par%*%deta.dtheta[i-ma_ind,]
      
      deta.dTHETA[i,]<- Rs[(i-m),] + theta%*%deta.dTHETA[i-ma,]+ THETA%*%deta.dTHETA[i-S*MA,]- ma_par%*%deta.dTHETA[i-ma_ind,]
    }
    
    a <- matrix(deta.dalpha[(m+1):n],ncol=1)
    rP <- matrix(deta.dphi[(m+1):n,],ncol=p1)
    rPs <- matrix(deta.dPHI[(m+1):n,],ncol=P1)
    rR <- matrix(deta.dtheta[(m+1):n,],ncol=q1)
    rRs <- matrix(deta.dTHETA[(m+1):n,],ncol=Q1)
    
    vI <- matrix(rep(1,n-m),ncol=1) 
    
    psi1 = trigamma(muhat * prec)
    psi2 = trigamma((1 - muhat) * prec)
    vc = matrix(prec * (psi1 * muhat - psi2 * (1 - muhat)),ncol=1)
    C = diag( as.vector(vc) )
    mT = diag(mu.eta(muhat))
    W = diag( c((prec^2)*( psi1 + psi2 )) )
    
    D = diag( c(psi1 * (muhat^2) + psi2 * (1 - muhat)^2 - trigamma(prec)) )
    T2<-mT^2
    
    Kaa <- t(a) %*% W %*% T2 %*% a 
    Kap <- t(a) %*% W %*% T2 %*% rP 
    Kpa <- t(Kap)
    KaP <- t(a) %*% W %*% T2 %*% rPs
    KPa <- t(KaP)
    Kat <- t(a) %*% W %*% T2 %*% rR
    Kta <- t(Kat)
    KaT <- t(a) %*% W %*% T2 %*% rRs
    KTa <- t(KaT)
    Kaprec <- Kpreca <- t(vI) %*% C %*% mT %*% a
    
    Kpp <- t(rP) %*% W %*% T2 %*% rP
    KpP <- t(rP) %*% W %*% T2 %*% rPs
    KPp <- t(KpP)
    Kpt <- t(rP) %*% W %*% T2 %*% rR
    Ktp <- t(Kpt)
    KpT <- t(rP) %*% W %*% T2 %*% rRs
    KTp <- t(KpT)
    Kpprec <- t(rP) %*% C %*% mT %*% vI
    Kprecp <- t(Kpprec)
    
    KPP <- t(rPs) %*% W %*% T2 %*% rPs
    KPt <- t(rPs) %*% W %*% T2 %*% rR
    KtP <- t(KPt)
    KPT <- t(rPs) %*% W %*% T2 %*% rRs
    KTP <- t(KPT)
    KPprec <- t(rPs) %*% C %*% mT %*% vI
    KprecP <- t(KPprec)
    
    Ktt <- t(rR) %*% W %*% T2 %*% rR
    KtT <- t(rR) %*% W %*% T2 %*% rRs
    KTt <- t(KtT)
    Ktprec <- t(rR) %*% C %*% mT %*% vI
    Kprect <- t(Ktprec)
    
    KTT <- t(rRs) %*% W %*% T2 %*% rRs
    KTprec <- t(rRs) %*% C %*% mT %*% vI
    KprecT <- t(KTprec)
    
    Kprecprec <- sum(diag(D))
    
    K <- rbind(
      cbind(Kaa,Kap,KaP,Kat,KaT,Kaprec),
      cbind(Kpa,Kpp,KpP,Kpt,KpT,Kpprec),
      cbind(KPa,KPp,KPP,KPt,KPT,KPprec),
      cbind(Kta,Ktp,KtP,Ktt,KtT,Ktprec),
      cbind(KTa,KTp,KTP,KTt,KTT,KTprec),
      cbind(Kpreca,Kprecp,KprecP,Kprect,KprecT,Kprecprec)
    )
    
    #### Seasonality test
    #vcov <- solve(K)
    vcov <- chol2inv(chol(K,pivot=T))
    z$vcov <- vcov
    
    z$wald <- as.numeric(t(coef[c((p1+2):(p1+P1+1),(p1+P1+q1+2):(p1+P1+q1+Q1+1)) ] ) %*% 
                           chol2inv(chol(z$vcov[c((p1+2):(p1+P1+1),(p1+P1+q1+2):(p1+P1+q1+Q1+1)),
                                                c((p1+2):(p1+P1+1),(p1+P1+q1+2):(p1+P1+q1+Q1+1))])) 
                         %*% coef[c((p1+2):(p1+P1+1),(p1+P1+q1+2):(p1+P1+q1+Q1+1))] )
    
    
    z$p_wald <- (1 - pchisq( z$wald, (Q1+P1) ) )
    
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha + (ar_par%*%ynew_prev[n+i-ar_ind]) - (ma_par%*%errorhat[n+i-ma_ind])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 # original scale
    }
    
    z$forecast <- y_prev[(n+1):(n+h1)]
  } 
  
  
  
  
  
  
  
  #################################################################
  ######### BSARMA sem ma
  if(any(is.na(ar)==F) && any(is.na(ma)==T) && any(is.na(AR)==F) && any(is.na(MA)==F))
  { 
    if(diag>0)print("BSARMA model without ma",quote=F)
    reg <- c(mqo, rep(0,(Q1)), prec) # initializing the parameter values
    
    loglik <- function(z) 
    {
      alpha <- z[1]
      phi = z[2:(p1+1)] 
      PHI = z[(p1+2):(p1+P1+1)]
      THETA = z[(p1+P1+2):(p1+P1+Q1+1)]
      prec <- z[(p1+P1+Q1+2)] # precision parameter
      
      error<-rep(0,n) # E(error)=0 
      eta<-rep(NA,n)
      
      ar_par_index <- operator(phi,PHI,ar,S*AR)
      
      ar_par <- ar_par_index$parameters
      ar_ind <- ar_par_index$index
      
      ma_par <- THETA
      ma_ind <- S*MA
      
      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + sum(ar_par*ynew[i-ar_ind]) - sum(ma_par*error[i-ma_ind])
        error[i] <- ynew[i]-eta[i] 
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings(dbeta(y1, mu * prec, (1 - mu) * prec, log = TRUE))
      sum(ll)
    } 
    
    escore <- function(z) 
    {
      
      alpha <- z[1]
      phi = z[2:(p1+1)] 
      PHI = z[(p1+2):(p1+P1+1)]
      THETA = z[(p1+P1+2):(p1+P1+Q1+1)]
      prec <- z[(p1+P1+Q1+2)] # precision parameter
      
      error<-rep(0,n) # E(error)=0 
      eta<-rep(NA,n)
      
      ar_par_index <- operator(phi,PHI,ar,S*AR)
      
      ar_par <- ar_par_index$parameters
      ar_ind <- ar_par_index$index
      
      ma_par <- THETA
      ma_ind <- S*MA
      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + sum(ar_par*ynew[i-ar_ind]) - sum(ma_par*error[i-ma_ind])
        error[i]<- ynew[i]-eta[i] 
      }
      
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      ystar <- log(y1/(1-y1))
      mustar <- digamma(mu * prec) - digamma((1 - mu) * prec)
    
      mT <- diag(mu.eta(mu))
      
      Rs <- matrix(rep(NA,(n-m)*Q1),ncol=Q1)
      for(i in 1:(n-m))
      {
        for(j in 1:Q1)
        {
          Rs[i,j] <- -sum(error[i+m-(S*MA[j])])
        }
      }
      
      A <- matrix(rep(NA,(n-m)*(p1)),ncol=(p1))
      for(i in 1:(n-m))
      {
        for(j in 1:p1)
        {
          A[i,j] <- sum(ynew[i+m-c(ar[j],ar[j]+(S*AR))]*c(1,-PHI))
        }
      }    
      
      As <- matrix(rep(NA,(n-m)*(P1)),ncol=(P1))
      for(i in 1:(n-m))
      {
        for(j in 1:P1)
        {
          As[i,j] <- sum(ynew[i+m-c(S*AR[j],S*AR[j]+(ar))]*c(1,-phi))
        }
      }    
      
      deta.dalpha <- rep(0,n)
      deta.dphi <- matrix(0, ncol=p1,nrow=n)
      deta.dPHI <- matrix(0, ncol=P1,nrow=n)
      deta.dTHETA <- matrix(0, ncol=Q1,nrow=n)
      
      for(i in (m+1):n)
      {
        deta.dalpha[i]<- 1  + THETA%*%deta.dalpha[i-S*MA]
        
        deta.dphi[i,]<- A[(i-m),] + THETA%*%deta.dphi[i-S*MA,]
        
        deta.dPHI[i,]<- As[(i-m),] + THETA%*%deta.dPHI[i-S*MA,]
        
        deta.dTHETA[i,]<- Rs[(i-m),] +  THETA%*%deta.dTHETA[i-S*MA,]
      }
      
      a <- matrix(deta.dalpha[(m+1):n],ncol=1)
      rP <- matrix(deta.dphi[(m+1):n,],ncol=p1)
      rPs <- matrix(deta.dPHI[(m+1):n,],ncol=P1)
      rRs <- matrix(deta.dTHETA[(m+1):n,],ncol=Q1)
      
      ymstar <- matrix((ystar-mustar),ncol=1)
      
      Ualpha <- prec * t(a) %*% mT %*% ymstar
      Uphi <-   prec * t(rP) %*% mT %*% ymstar
      UPHI <-   prec * t(rPs) %*% mT %*% ymstar
      UTHETA <- prec * t(rRs) %*% mT %*% ymstar
      Uprec <-  sum(mu * (ystar-mustar) + log(1-y1)
                    - digamma((1 - mu) * prec) + digamma(prec) )
      
    
      rval <- c(Ualpha,Uphi,UPHI,UTHETA,Uprec)
    }
    names_par <- c("alpha",names_phi,names_PHI,names_THETA,"precision")
    
    opt <- optim(reg, loglik, escore, method = "BFGS", 
                 control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    
    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")
    
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+P1+Q1+2)]
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <- coef[1]
    phi = coef[2:(p1+1)] 
    PHI = coef[(p1+2):(p1+P1+1)]
    THETA = coef[(p1+P1+2):(p1+P1+Q1+1)]
    prec <- coef[(p1+P1+Q1+2)] # precision parameter
    
    z$alpha <- alpha
    z$phi <- phi
    z$PHI <- PHI
    z$THETA <- THETA
    z$prec <- prec
    
    errorhat<-rep(0,n) # E(error)=0 
    etahat<-rep(NA,n)
    
    ar_par_index <- operator(phi,PHI,ar,S*AR)
    
    ar_par <- ar_par_index$parameters
    ar_ind <- ar_par_index$index
    
    ma_par <- THETA #ma_par_index$parameters
    ma_ind <- S*MA #ma_par_index$index
    
    for(i in (m+1):n)
    {
      etahat[i] <- alpha + (ar_par%*%ynew[i-ar_ind]) - (ma_par%*%errorhat[i-ma_ind])
      errorhat[i] <- ynew[i]-etahat[i] # predictor scale
      
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    z$mustarhat <- digamma(muhat * prec) - digamma((1 - muhat) * prec)
    
    Rs <- matrix(rep(NA,(n-m)*Q1),ncol=Q1)
    for(i in 1:(n-m))
    {
      for(j in 1:Q1)
      {
        Rs[i,j] <- -sum(errorhat[i+m-(S*MA[j])])
      }
    }
    A <- matrix(rep(NA,(n-m)*(p1)),ncol=(p1))
    for(i in 1:(n-m))
    {
      for(j in 1:p1)
      {
        A[i,j] <- sum(ynew[i+m-c(ar[j],ar[j]+(S*AR))]*c(1,-PHI))
      }
    }    
    As <- matrix(rep(NA,(n-m)*(P1)),ncol=(P1))
    for(i in 1:(n-m))
    {
      for(j in 1:P1)
      {
        As[i,j] <- sum(ynew[i+m-c(S*AR[j],S*AR[j]+(ar))]*c(1,-phi))
      }
    }    
    
    deta.dalpha <- rep(0,n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dPHI <- matrix(0, ncol=P1,nrow=n)
    deta.dTHETA <- matrix(0, ncol=Q1,nrow=n)
    
    for(i in (m+1):n)
    {
      deta.dalpha[i]<- 1 +  THETA%*%deta.dalpha[i-S*MA]
      
      deta.dphi[i,]<- A[(i-m),] + THETA%*%deta.dphi[i-S*MA,]
      
      deta.dPHI[i,]<- As[(i-m),] + THETA%*%deta.dPHI[i-S*MA,]
      
      deta.dTHETA[i,]<- Rs[(i-m),] +  THETA%*%deta.dTHETA[i-S*MA,]
    }
    
    a <- matrix(deta.dalpha[(m+1):n],ncol=1)
    rP <- matrix(deta.dphi[(m+1):n,],ncol=p1)
    rPs <- matrix(deta.dPHI[(m+1):n,],ncol=P1)
    rRs <- matrix(deta.dTHETA[(m+1):n,],ncol=Q1)
    
    vI <- matrix(rep(1,n-m),ncol=1) 
    
    psi1 = trigamma(muhat * prec)
    psi2 = trigamma((1 - muhat) * prec)
    vc = matrix(prec * (psi1 * muhat - psi2 * (1 - muhat)),ncol=1)
    C = diag( as.vector(vc) )
    mT = diag(mu.eta(muhat))
    W = diag( c((prec^2)*( psi1 + psi2 )) )
    
    D = diag( c(psi1 * (muhat^2) + psi2 * (1 - muhat)^2 - trigamma(prec)) )
    T2<-mT^2
    
    Kaa <- t(a) %*% W %*% T2 %*% a 
    Kap <- t(a) %*% W %*% T2 %*% rP 
    Kpa <- t(Kap)
    KaP <- t(a) %*% W %*% T2 %*% rPs
    KPa <- t(KaP)
    KaT <- t(a) %*% W %*% T2 %*% rRs
    KTa <- t(KaT)
    Kaprec <- Kpreca <- t(vI) %*% C %*% mT %*% a
    
    Kpp <- t(rP) %*% W %*% T2 %*% rP
    KpP <- t(rP) %*% W %*% T2 %*% rPs
    KPp <- t(KpP)
    KpT <- t(rP) %*% W %*% T2 %*% rRs
    KTp <- t(KpT)
    Kpprec <- t(rP) %*% C %*% mT %*% vI
    Kprecp <- t(Kpprec)
    
    KPP <- t(rPs) %*% W %*% T2 %*% rPs
    KPT <- t(rPs) %*% W %*% T2 %*% rRs
    KTP <- t(KPT)
    KPprec <- t(rPs) %*% C %*% mT %*% vI
    KprecP <- t(KPprec)
    
    KTT <- t(rRs) %*% W %*% T2 %*% rRs
    KTprec <- t(rRs) %*% C %*% mT %*% vI
    KprecT <- t(KTprec)
    
    Kprecprec <- sum(diag(D))
    
    K <- rbind(
      cbind(Kaa,Kap,KaP,KaT,Kaprec),
      cbind(Kpa,Kpp,KpP,KpT,Kpprec),
      cbind(KPa,KPp,KPP,KPT,KPprec),
      cbind(KTa,KTp,KTP,KTT,KTprec),
      cbind(Kpreca,Kprecp,KprecP,KprecT,Kprecprec)
    )
    
    #### Seasonality test
    #vcov <- solve(K)
    vcov <- chol2inv(chol(K,pivot=T))
    z$vcov <- vcov
    
    
    z$wald <- as.numeric(t(coef[c((p1+2):(p1+P1+1),(p1+P1+2):(p1+P1+Q1+1)) ] ) %*% 
                           chol2inv(chol(z$vcov[c((p1+2):(p1+P1+1),(p1+P1+2):(p1+P1+Q1+1)),
                                                c((p1+2):(p1+P1+1),(p1+P1+2):(p1+P1+Q1+1))])) 
                         %*% coef[c((p1+2):(p1+P1+1),(p1+P1+2):(p1+P1+Q1+1))] )
    
    
    z$p_wald <- (1 - pchisq( z$wald, (Q1+P1) ) )
    
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha + (ar_par%*%ynew_prev[n+i-ar_ind]) - (ma_par%*%errorhat[n+i-ma_ind])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 # original scale
    }
    
    z$forecast <- y_prev[(n+1):(n+h1)]
  } 
  
  
  
  
  ############################################################################
  ############################################################################
  
  z$serie <- y
  z$barma <- names_par
  
  
  # residuals
  res1 <- y-z$fitted
  res2 <- ynew-z$etahat
  res3 <- as.vector(ystar)[(m+1):n]-z$mustarhat
  z$resid1 <- (res1[(m+1):n])/sqrt(z$fitted[(m+1):n]*(1-z$fitted[(m+1):n])/(1+prec))
  z$resid2 <- (res2[(m+1):n])/sqrt((mu.eta(z$fitted[(m+1):n])^{-2})*(z$fitted[(m+1):n]*(1-z$fitted[(m+1):n])/(1+prec)) )
  z$resid3 <- (res3/sqrt( trigamma(z$fitted[(m+1):n]*prec)+trigamma((1-z$fitted[(m+1):n])*prec) ) )
  
  l_tilde <- (dbeta(y[(m+1):n], y[(m+1):n] * prec, (1 - y[(m+1):n]) * prec, log = TRUE))
  l_hat <- (dbeta(y[(m+1):n], z$fitted[(m+1):n] * prec, (1 - z$fitted[(m+1):n]) * prec, log = TRUE))
  
  dt <- (l_tilde-l_hat)
  dt[which(dt<0)]<-0
  
  if(resid==1) resid2 <- z$resid1
  if(resid==2) resid2 <- z$resid2
  if(resid==3) resid2 <- z$resid3
  if(resid==4) resid2 <- z$resid4 <- sign(y[(m+1):n]-z$fitted[(m+1):n])*sqrt(2*(dt))
  
  z$deviance <- 2*sum(dt)
  
  k<- sum(c(p1,q1,P1,Q1),na.rm=T) # numero de parametros
 
  z$m<- m
  z$rank <-k
  
  stderror <- sqrt(diag(z$vcov))
  z$stderror <- stderror
  
  z$zstat <- abs(z$coef/stderror)
  z$pvalues <- 2*(1 - pnorm(z$zstat) )
  
  #z$loglik <- opt$value
  z$loglik <- opt$value*(n/(n-m))
  z$counts <- as.numeric(opt$counts[1])

  z$aic <- -2*z$loglik+2*(p1+q1+P1+Q1+2)
  z$bic <- -2*z$loglik+log(n)*(p1+q1+P1+Q1+2)
  
  model_presentation <- cbind(round(z$coef,4),round(z$stderror,4),round(z$zstat,4),round(z$pvalues,4))
  colnames(model_presentation)<-c("Estimate","Std. Error","z value","Pr(>|z|)")
  
  z$model <- model_presentation
  z$link <- link
  
  if(is.na(ljunbox_lag)==T)
  {
    ljunbox_lag <- 10 # non-seasonal
    
    if(any(is.na(AR)==F) || any(is.na(AR)==F)){
      ljunbox_lag <- 2*S
    }
  } 
  
  ljungbox<- Box.test(resid2, lag = ljunbox_lag, type = "Ljung-Box", fitdf = k)
  monti<- Monti.test(resid2, lag = ljunbox_lag, type = "Ljung-Box", fitdf = k)
  
  z$ljungbox<-ljungbox$statistic
  z$p_ljungbox<-ljungbox$p.value
  
  z$monti<-monti$statistic
  z$p_monti<-monti$p.value
  
  ###################################################
  ######### GRAPHICS ################################
  
  if(diag>0)
  {
    
    print(model_presentation)
    print(" ",quote=F)
    print(c("Log-likelihood =",round(z$loglik,4)),quote=F)
    print(c("Number of iterations in BFGS optim =",z$counts),quote=F)
    print(c("MAIC =",round(z$aic,4),"MSIC =",round(z$bic,4)),quote=F)
    print(c("Deviance =",round(z$deviance,4)," DF:",n-m-k),quote=F)
    if(any(is.na(AR)==F) || any(is.na(MA)==F))
      print(paste("Seasonality test: W = ",round(z$wald,4),", p-value = ",round(z$p_wald,4),sep=""),quote=F)
    print("Residuals:",quote=F)
    print(summary(resid2))
    print(ljungbox)
    print(monti)
    
    t<-seq(-5,n+6,by=1)
    
    par(mfrow=c(1,1))
    par(mar=c(2.8, 2.7, 1.2, 1)) # margens c(baixo,esq,cima,direia)
    par(mgp=c(1.7, 0.45, 0))
    plot(resid2,main=" ",xlab="Index",ylab="Residuals", pch = "+",ylim=c(-4,4))
    lines(t,rep(-3,n+12),lty=2,col=1)
    lines(t,rep(3,n+12),lty=2,col=1)
    lines(t,rep(-2,n+12),lty=3,col=1)
    lines(t,rep(2,n+12),lty=3,col=1)
    
    
    max_y<- max(c(z$fitted,y),na.rm=T)
    min_y<- min(c(z$fitted,y),na.rm=T)
    plot(as.vector(z$fitted), as.vector(y), main=" ", pch = "+",
         xlab="Fitted values",ylab="Observed data",
         xlim=c(0.95*min_y,max_y*1.05),
         ylim=c(0.95*min_y,max_y*1.05))
    lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
    
    densidade<-density(resid2)
    plot(densidade,ylab="density",main=" ")
    lines(densidade$x,dnorm(densidade$x),lty=2)
    legend("topleft",c("Exact distribution of residuals","Normal approximation"),#pch=vpch,
           pt.bg="white", lty=c(1,2), bty="n")
    
    acf(resid2,ylab="ACF",xlab="Lag") # função de autocorrelação
    
    pacf(resid2,ylab="PACF",xlab="Lag") # função de autocorrelação parcial
    
    max_r<- max(resid2,na.rm=T)
    min_r<- min(resid2,na.rm=T)
    qqnorm(resid2, pch = "+",
           xlim=c(0.95*min_r,max_r*1.05),
           ylim=c(0.95*min_r,max_r*1.05),
           main="",xlab="Normal quantiles",ylab="Empirical quantiles")
    lines(c(-10,10),c(-10,10),lty=2)
    
    par(mfrow=c(1,1))
    plot(y,type="l",ylab="serie",xlab="tempo")
    lines(z$fitted,col="red")
    
    fim<-end(y)[1]+end(y)[2]/12
    
    y_prev <- ts(y_prev, start=start(y), frequency=frequency(y))
    par(mfrow=c(1,1))
    plot(y_prev,type="l",col="red", ylim=c(min(y),max(y)),ylab="Serie",xlab="Time")
    abline(v=fim,lty=2)
    lines(y)
    
    w1<-5
    h1<-4
    
    if(diag>1)
    {
      #pdf(file = "resid_v_ind.pdf",width = w1, height = h1,family = "Times")
      postscript(file = "resid_v_ind.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
        par(mgp=c(1.7, 0.45, 0))
        plot(resid2,main=" ",xlab="Index",ylab="Residuals", pch = "+",ylim=c(-4,4))
        lines(t,rep(-3,n+12),lty=2,col=1)
        lines(t,rep(3,n+12),lty=2,col=1)
        lines(t,rep(-2,n+12),lty=3,col=1)
        lines(t,rep(2,n+12),lty=3,col=1)
      }
      dev.off()
      
      #pdf(file = "obs_v_fit.pdf",width = w1, height = h1,family = "Times")
      postscript(file = "obs_v_fit.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
        par(mgp=c(1.7, 0.45, 0))
        plot(as.vector(z$fitted), as.vector(y), main=" ", pch = "+",
             xlab="Fitted values",ylab="Observed data",
             xlim=c(0.95*min_y,max_y*1.05),
             ylim=c(0.95*min_y,max_y*1.05))
        lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
      }
      dev.off()
      
      #pdf(file = "resid_density.pdf",width = w1, height = h1,family = "Times")
      postscript(file = "resid_density.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(1.5, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
        par(mgp=c(1.7, 0.45, 0))
        
        plot(densidade,ylab="Density",main=" ",xlab=" ",ylim=c(0,1.15*max(densidade$y)))
        lines(densidade$x,dnorm(densidade$x),lty=2)
        #legend("topleft",c("Densidade estimada","Normal padrão"),#pch=vpch,
        legend("topleft",c("Exact distribution of residuals","Normal approximation"),#pch=vpch,
               pt.bg="white", lty=c(1,2), bty="n")
      }
      dev.off()
      
      #pdf(file = "resid_FAC.pdf",width = w1, height = h1,family = "Times")
      postscript(file = "resid_FAC.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
        par(mgp=c(1.7, 0.45, 0))
        acf(resid2,ylab="ACF",xlab="Lag") # função de autocorrelação
      }
      dev.off()
      
      #pdf(file = "resid_FACP.pdf",width = w1, height = h1,family = "Times")
      postscript(file = "resid_FACP.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
        par(mgp=c(1.7, 0.45, 0))
        pacf(resid2,ylab="PACF",xlab="Lag") # função de autocorrelação parcial
      }
      dev.off()
      
      #pdf(file = "qq_plot.pdf",width = w1, height = h1,family = "Times")
      postscript(file = "qq_plot.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {  
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
        par(mgp=c(1.7, 0.45, 0))
        qqnorm(resid2, pch = "+",
               xlim=c(0.95*min_r,max_r*1.05),
               ylim=c(0.95*min_r,max_r*1.05),
               main="",xlab="Normal quantiles",ylab="Empirical quantiles")
        lines(c(-10,10),c(-10,10),lty=2)
      }
      dev.off()
      
      #pdf(file = "adjusted.pdf",width = 6, height = 4,family = "Times")
      postscript(file = "adjusted.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
        par(mgp=c(1.7, 0.45, 0))
        plot(y,type="l",ylab="RH",xlab="Time",ylim=c(0.58,0.9))
        lines(z$fitted,lty=2,col=1)
        
        legend("bottomleft",c("Observed data","Fitted values"),#pch=vpch,
               pt.bg="white", lty=c(1,2), bty="n",col=c(1,1))
        
      }
      dev.off()
      
      
      postscript(file = "forecast.eps",horizontal=F,paper="special",width = 6, height = 4.7,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
        par(mgp=c(1.7, 0.45, 0))
        plot(y_prev,type="l",lty=2,col="red", ylim=c(min(y),max(y)),ylab="RH",xlab="Times")
        abline(v=fim,lty=2)
        lines(y)
        legend("bottomleft",c("Observed data","Fitted and forecast values"),#pch=vpch,
               pt.bg="white", lty=c(1,2), bty="n",col=c(1,"red"))
      }
      dev.off()
      
    }    
  }  
  return(z)
}