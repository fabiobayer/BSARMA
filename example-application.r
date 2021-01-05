# Application explored in 
# BAYER, F. M.; CINTRA, R. J.; CRIBARI-NETO, F.. 
# Beta seasonal autoregressive moving average models. 
# Journal of Statistical Computation and Simulations, 
# v. 88, p. 2961–2981, 2018.

source("bsarma.r") # read the code with barma function

ur<-scan(file="ur-sm-mensal.txt") # read the data
y<-ur/100 # transformation to (0,1)
y<-ur[13:180]/100 # Start in January 2003
Y<-ts(y,start=c(2003,1),frequency=12) # time series object

# descrptive analysis
mean(y) 
summary(y) 
hist(y) 
plot(Y) # plota a serie temporal

library(forecast) # useful package
seasonplot(Y) 
monthplot(Y) 
plot(decompose(Y)) 

# correlogram
acf(Y,lag.max = 24) 
pacf(Y,lag.max = 24) 

# zera as opcoes graficas
# op <- par(no.readonly = TRUE) # the whole list of settable par's.
# par(op)

## Fit Beta SARMA model
# Input arguments:
# ar: AR orders
# AR: seasonal AR orders
# ma: MA orders
# MA: seasonal MA orders
# h: forecast horizon
# diag: diagnostic graphs:
#   diag = 0 : do not plot any graph (useful for simulations)
#   diag = 1 : plot graphs 
#   diga = 2 : save graphs on ps files
# resid: there are three types of residuals: 
#   resid = 1 : standardized residual
#   resid = 2 : standardized residual 2
#   resid = 3 : standardized weighted residual

# fit the model
fit1<-barma(Y,ar=c(1),AR=c(1),MA=c(1),h=12,diag=2,resid=3)

