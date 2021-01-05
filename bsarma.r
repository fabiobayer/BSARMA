# Implemented by Fabio M Bayer (bayer@ufsm.br) 
# Date: 06/2015
#
# Reference:
# BAYER, F. M.; CINTRA, R. J.; CRIBARI-NETO, F.. 
# Beta seasonal autoregressive moving average models. 
# Journal of Statistical Computation and Simulations, 
# v. 88, p. 2961–2981, 2018.
#
# Input arguments:
# ar: AR orders
# AR: seasonal AR orders
# ma: MA orders
# MA: seasonal MA orders
#
# h: forecast horizon
#
# diag: diagnostic graphs:
#   diag = 0 : do not plot any graph (useful for simulations)
#   diag = 1 : plot graphs 
#   diga = 2 : save graphs on ps files
#
# resid: there are three types of residuals: 
#   resid = 1 : standardized residual
#   resid = 2 : standardized residual 2
#   resid = 3 : standardized weighted residual
#   resid = 4 : deviance residual
#
# The input data needs to be a time series object (ts)
#
# Usage examples:
#
# 1) BARMA(2,3) with logit link and 6 out-of-sample forecasts
# fit <- barma(y,ar=c(1,2),ma=c(1,2,3))
# 
# 2) BSARMA(p,q)(P,Q) with cloglog link and 12 out-of-sample forecasts
# fit <- barma(y,ar=c(1,2,5),ma=c(1,3),AR=c(1,2),MA=1,link="cloglog",h=12)
# 

barma<- function(y, ar=NA, ma=NA, AR=NA, MA=NA, link = "logit",diag=1,h=6,resid=3,ljungbox=NA)
{  
  source("bsarma.fit.r")
  
  if( sum(is.na(y)) >0 )
  {
    fit1<-c()
    fit1$conv <- 1
    warning("OUT OF RANGE (0,1)!")
    return(fit1)
  }else{
    if (min(y) <= 0 || max(y) >= 1) 
    {
      fit1<-c()
      fit1$conv <- 1
      warning("OUT OF RANGE (0,1)!")
      return(fit1)
      #stop("OUT OF RANGE (0,1)!")
    }else{
      
      
      if(is.ts(y)==T)
      {
        S<-frequency(y)
      }else stop("data can be a time-series object")
      
      
      if(any(is.na(ar))==F)
      {
        names_phi<-c(paste("phi",ar,sep=""))
      }else{
        names_phi <- NA
      }
      
      if(any(is.na(ma))==F)
      {
        names_theta<-c(paste("theta",ma,sep=""))
      }else{
        names_theta <- NA
      }
      
      if(any(is.na(AR))==F)
      {
        names_PHI<-c(paste("PHI",AR,sep=""))
      }else{
        names_PHI <- NA
      }
      
      if(any(is.na(MA))==F)
      {
        names_THETA<-c(paste("THETA",MA,sep=""))
      }else{
        names_THETA <- NA
      }
      
      #   p <- max(ar)
      #   q <- max(ma)
      #   P <- max(AR)
      #   Q <- max(MA)
      #   n <- length(y)
      #   m <- max(p,q,S*P,S*Q,na.rm=T)
      #   p1 <- length(ar)
      #   q1 <- length(ma)
      #   P1 <- length(AR)
      #   Q1 <- length(MA)
      
      linktemp <- substitute(link)
      if (!is.character(linktemp))
      {
        linktemp <- deparse(linktemp)
        if (linktemp == "link")
          linktemp <- eval(link)
      }
      if (any(linktemp == c("logit", "probit", "cloglog")))
        stats <- make.link(linktemp)
      else stop(paste(linktemp, "link not available, available links are \"logit\", ",
                      "\"probit\" and \"cloglog\""))
      
      link1 <- structure(list(link = linktemp, 
                              linkfun = stats$linkfun,
                              linkinv = stats$linkinv, 
                              mu.eta = stats$mu.eta, 
                              diflink = function(t) 1/(stats$mu.eta(stats$linkfun(t)))
      )
      )
      
      
      fit1 <- barma.fit(y, ar, ma, AR, MA, link1, names_phi, names_PHI, names_theta, names_THETA, diag, h,resid,ljunbox_lag=ljungbox) # model estimation
      
      return(fit1)
    }
  }
}

