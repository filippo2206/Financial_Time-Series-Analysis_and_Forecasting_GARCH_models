################################################################################
###### ANALISI TIME SERIES FINANZIARIA SULLA SOCIETA' ORACLE CORPORATION    ######

rm(list = ls())

library(tseries)  
library(sandwich)
library(lmtest)
library(urca)   
library(rugarch)  
library(FinTS)   
library(car)
library(forecast) 
library(xts)     
library(quantmod) 

source("C:/Users/alida/Desktop/UNI/3 anno/Statistica Economica/TSA-Predict-Student-Functions.R")
source("C:/Users/alida/Desktop/UNI/3 anno/Statistica Economica/TSA-Finance-Functions.R")


file.data <- "C:/Users/alida/Desktop/ORCL ts finanziaria.csv"

data <- read.table(file = file.data, header = TRUE, sep = ",", 
                   check.names = FALSE, comment.char = "", na.strings = ".")

colnames(data)[6] <- "Adjusted" 
data <- data.frame(data, cc.ret = c(NA, diff(log(data$Adjusted))),    
            gkVol = .garmanklass(data = data, sd = TRUE),    
            check.names = TRUE)
################################################################################

###ESTRAGGO PERIODO
ind   <- as.Date(x = "2015-01-03") <= as.Date(x = data$Date) #estraggo dati dal 2011 ma poi dopo il test di Nyblom li prendo dal 2015
data  <- data[ind, , drop = FALSE]


###ESTRAGGO DA ANLISI PRELIMINARE 2 TIPI DI PREZZI: YC CLOSE, Y ADJUSTED
time  <- as.Date(x = data$Date)
yc    <- data$Close
yclog <- log(yc)
y     <- data$Adjusted
ylog  <- log(y)




################################################################################
###########          ANALISI DEI PREZZI                             ############
################################################################################
#PLOT CLOSE E ADJUSTED
nobs <- NROW(y)

par(mfrow = c(2,2))
plot(x = time, y = yc,    main = "Close",        xlab = "Time", ylab = "Price", type = "l")   
plot(x = time, y = yclog, main = "Ln(close)",    xlab = "Time", ylab = "Price", type = "l")  
plot(x = time, y = y,     main = "AdjClose",     xlab = "Time", ylab = "Price", type = "l")   
plot(x = time, y = ylog,  main = "Ln(AdjClose)", xlab = "Time", ylab = "Price", type = "l")


#ANALISI CORRELAZIONI
par(mfrow = c(2,1))
Acf(x = ylog, lag.max = 100, type = "correlation", main = "Price")
Acf(x = ylog, lag.max = 100, type = "partial", main = "Price")


# ACF E PACF CARATTERISTICHE DELLE SERIE GIORNALIERE DEI DATI FINANZIARI
#Acf e Pacf, non stazionarietà e decadimento lineare di Acf e picco in lag1 in Pacf, CHIARO INDIZIO DI NON STAZIONARITA'
#Non c'è stagionalità




#################################################
#########UNIT ROOT###############################
#################################################

cat("\n-----------------------------------------------------------------
  Unit root analysis following the Di Fonzo-Lisi procedure\n")
adf.1 <- ur.df(y = ylog, type = "trend", lags = 20, selectlags = "AIC")

cat("\n-----\nTest1: ADF with trend\n") ## ACCETTO PASSO 1 E 2 E VADO AL PASSO 3
print( summary(adf.1) )

adf.2 <- ur.df(y = ylog, type = "drift", lags = 20, selectlags = "AIC")
cat("\n-----\nTest1: ADF with drift\n")
print( summary(adf.2) )                              ##ACCETTO AL PASSO 3, HO UNIT ROOT 


################################################################################
#PREZZI NON STAZIONARI USO I LOG RENDIMENTI
#### LOG RETURNS
yret <- xts(x = 100 * data$cc.ret, order.by = time)  #MI ESTRAGGO I RENDIMENTI SU SCALA LOG
#yret <- 100 * data$cc.ret



###########ANALISI PRELIMINARE 
cat("\n-----------------------------------------------------------------
  Preliminary analysis of log-returns\n")
par(mfrow = c(1,1))
plot(x = time, y = yret, main = "Returns", 
     xlab = "Time", ylab = "Log-Returns", type = "l")

#I PREZZI HANNO UN ANDAMENTO STILE RW, ANCHE I PLOT PRIMA CI DICONO CHE NON SONO TROPPO DIVERSI DA UN RW
#QUESTO VALE ANCHE PER I LOG, E QUINDI SE FACCIAMO LE DIFFERENZE PRIME SARA STILE WN
#QUESTO E' QUELLO CHE SEMBRA DEI RENDIMENTI, SI MUOVONO INTORNO ALLO 0 COME WN A MEDIA 0,
#MA CE' UNA DIFFERENZA, IL WN E' OMOSCHEDASTICO QUINDI LA VARIABILITA INTRONO ALLA MEDIA E' SEMPRE
#LA STESSA MENTRE QUI NO, FA UN PO DA FISARMONICA, CI SONO DEI PICCHI IN 2012 E 2020 MOLTO MARCATI.




par(mfrow = c(2,1))
Acf(x = yret, lag.max = 100, type = "correlation", main = "Log-Returns Acf")
Acf(x = yret, lag.max = 100, type = "partial", main = "Log-Returns Pacf")
#ESSENDO MOLTO SIMILI SIAMO VICINI A UN WN 

cat("\nLjung-Box statistics on log-returns\n")
npar <- 0
lag <- c(2, 5, 10, 15, 20, 30, 50) + npar
lb <- mapply(FUN = Box.test, lag = lag, 
             MoreArgs = list(x = yret, type = "Ljung-Box", fitdf = npar))[1:3,]
print(rbind(lag = lag, lb))

#RIFIUTO HO QUINDI IL FATTO CHE LE PRIME 50 AUTOCORRELAZIONE SONO TUTTE DIVERSE 



######TEST UNIT ROOT SU LOG RETURNS
cat("\n-----------------------------------------------------------------
  Unit root analysis following the Di Fonzo-Lisi procedure\n")
adf.1 <- ur.df(y = yret, type = "trend", lags = 20, selectlags = "AIC")   #yret[-1] quando ho tutte le osservazioni

cat("\n-----\nTest1: ADF with trend\n") 
print( summary(adf.1) )

adf.2 <- ur.df(y = yret, type = "drift", lags = 20, selectlags = "AIC")    #yret[-1] quando ho tutte le osservazioni
cat("\n-----\nTest1: ADF with drift\n")
print( summary(adf.2) )  

#RIFIUTO PASSO 1, NON HO UNIT ROOT




#TEST INDIPENDENZA E IDENTICA DITRIBUZIONE
x1 <- yret  #yret[-1] quando ho tutte le osservazioni
bds <- bds.test(x = x1, m = 4, 
                eps = seq(from = 0.5 * sd(x1), to = 2 * sd(x1), length = 4),
                trace = FALSE)
cat("BDS test on returns\n")
print(bds)

#I LOG RETURNS NON SONO IID 




###########ARCH TEST

cat("\n-----------------------------------------------------------------
  ARCH based preliminary analyses\n")
cat("ARCH test on demeaned log-returns\n")
lag <- c(4, 8, 12, 16)
at <- mapply(FUN = ArchTest, lags = lag, 
             MoreArgs = list(x = yret, demean = TRUE))
print(at[1:3,])

#TUTTI SIGNIFICATIVI E QUINDI SI RIFIUTA HO E I RENDIMENTI SONO ETEROSCHEDASTICI






#FACCIAMO ACF RENDIMENTI, RENDIMENTI VALORE ASSOLUTO E RENDIMENTI AL QUADRATO

par(mfrow = c(3,1))
Acf(x = yret, lag.max = 100, type = "correlation", main = "Returns")
Acf(x = abs(yret), lag.max = 100, type = "correlation", main = "|Returns|")
Acf(x = yret^2, lag.max = 100, type = "correlation", main = expression(Returns^2))

#QUESTO CERTIFICA DA UN ALTRO PUNTO DI VISTA IL VOLATILITY CLUSTERING MISURATO IN TERMINI DI CORRELAZIONI
#SI VEDE MOLTO MEGLIO CHEI RENDIMENTI IN VALORE ASSOLUTO CHE IN QUELLI AL QUADRATO



#Jarque-Bera LOG RETURNS ISTOGRAMMA    #yret[-1] quando ho tutte le osservazioni
par(mfrow = c(1,2))  
.hist(x = yret, xlim = c(-10, 10), n = 200, breaks = 200, main = "Returns")
qqnorm(y = scale(yret))
abline(a = 0, b = 1, col = "red")

cat("\nJarque-Bera statistics on log-returns")
print( jarque.bera.test(x = yret) )

#RIFIUTO IL TEST
#ALTRA FACCIA DEL VOLATILITY CLUSTERING
#IL FATTO CHE ESCONO E OSCILLANO E' QUELLO CHE CREA CODE PIU GRNADI DELLA NORMALE, IPERNORMALITA' E QUESTA E' UNA 
#CERTIFICAZIONE DEL VOLATILITY CLUSTERING DA UN ALTRO PUNTO DI VISTA, PIU SIMILE LA T DI STUDENT CHE 
#LA NORMALE


################################################################################
## ARMA
#################################################################################
#arma(0,0)   Akaike       3.365064 Bayes        3.374389


spec0 <- arfimaspec(
  mean.model = list(armaOrder = c(0,0), 
                    include.mean = TRUE, external.regressors = NULL),  #VISTO CHE ABBIAMO TANTA LEPTOCURTOSI METTIAMO LA T DI STUDENT
  distribution.model = "std") 
fit0 <- arfimafit(spec = spec0, data = yret, 
                  solver = "solnp")

#yret[-1] quando ho tutte le osservazioni

np0 <- NROW(fit0@fit$coef)

cat( "\nInformation Criteria" ) #VALORI MEDI AIC BIC
print( infocriteria(fit0) )


cat("\nMatcoef\n")             #MATRICE COEFFICENTI
print( fit0@fit$matcoef )


cat("\nRobust matcoef\n")      #MATRICE COEFFICENTI ROBUSTI
print( fit0@fit$robust.matcoef )

#PRENDO MODELLO ARMA(0,0) CON AIC E BIC MEDI PIU BASSI


########DIAGNOSTICA
res <- as.numeric( residuals(fit0) )
par(mfrow = c(3,1))
Acf(x = res, lag.max = 100, type = "correlation", main = "Returns")
Acf(x = abs(res), lag.max = 100, type = "correlation", main = "|res|")
Acf(x = res^2, lag.max = 100, type = "correlation", main = expression(res^2))

#NON SI FA LA PACF, PRIMA BISOGNAVA DECIDERE AR, MA IN QUESTO CASO NON CE' MOLTO DA DECIDERE
# E LA PACF SIMILE ALLA ACF
#3 GRAFICI RESIDUI, AL QUADRATO E VALORE ASSOLUTO,A LIVELLO DI ETEROSCHEDASTICITA' NON ABBBIAMO FATTO NULLA 
#ARMA NON RIESCE A CORREGGERLA


################################################################################
## ARCH/GARCH modeling
################################################################################


####SGARCH
cat("\n-----------------------------------------------------------------
  GARCH on log-returns\n")

spec1 <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1), 
                        submodel = NULL, external.regressors = NULL, variance.targeting = FALSE), 
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,  
                    external.regressors = NULL), 
  distribution.model = "std")
fit1 <- ugarchfit(spec = spec1, data = yret, solver = "solnp")
#yret[-1] quando ho tutte le osservazioni

np1 <- NROW(fit1@fit$coef)


cat( "\nInformation Criteria" )
print( infocriteria(fit1) )


cat("\nMatcoef\n")
print( fit1@fit$matcoef )


cat("\nRobust matcoef\n")
print( fit1@fit$robust.matcoef )


#garch(1,0)  Akaike      3.417801 Bayes        3.428508
#garch(2,0)  Akaike      3.399040 Bayes        3.411888
#garch(2,1)  Akaike      3.374612 Bayes        3.389602
#garch(1,1)  Akaike      3.373783 Bayes        3.386632



#DIAGNOSTICHE SUL MODELLO
fit <- fit1
par(mfrow = c(3,1))
Acf(x = fit@fit$z,      lag.max = 100, type = "correlation", main = "z")
Acf(x = abs(fit@fit$z), lag.max = 100, type = "correlation", main = "|z|")
Acf(x = fit@fit$z^2,    lag.max = 100, type = "correlation", main = expression(z^2))
lag1 <- np1 + c(1, 2, 5, 10, 15, 20)

cat("\nLjung-Box on standardized residuals:\n")
lb1 <- mapply(FUN = Box.test, lag = lag1, 
              MoreArgs = list(x = fit@fit$z, type = "Ljung-Box", fitdf = np1) )
print(rbind(lag = lag1, lb1[1:3,]))

cat("\nLjung-Box on |standardized residuals|\n")
lb1 <- mapply(FUN = Box.test, lag = lag1, 
              MoreArgs = list(x = abs(fit@fit$z), type = "Ljung-Box", fitdf = np1) )
print(rbind(lag = lag1, lb1[1:3,]))

cat("\nLjung-Box on standardized residuals^2:\n")
lb1 <- mapply(FUN = Box.test, lag = lag1, 
              MoreArgs = list(x = fit@fit$z^2, type = "Ljung-Box", fitdf = np1) )
print(rbind(lag = lag1, lb1[1:3,]))








#### ARCH test
cat("\nARCH test on standardized residuals\n")
lag <- c(4, 8, 12, 16)
at <- mapply(FUN = ArchTest, lags = lag, 
             MoreArgs = list(x = fit1@fit$z, demean = TRUE))
print(at[1:3,])

#ACCETTO HO

#DISTRIBUZIONE DEI RESIDUI
par(mfrow = c(1,2))
xlim <- c(-5, 5)
.hist.fit(fit = fit1, xlim = xlim, ylim = c(0,0.75), n = 200, breaks = 100, 
          plot.norm = TRUE, main = "")
.qqplot.fit(fit = fit1)

#LA DISTRIBUZIONE T DI STUDENT SI ADATTA MOLTO MEGLIO ALL'ISTOGRAMMA


#### Leverage check
cat("\nSign bias test\n")
print( signbias(fit1) )
#ACCETTO HO, QUINDI NON HO ASIMMETRIE DI QUEI TIPI


###GJR GARCH
spec2 <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1), 
                        submodel = NULL, external.regressors = NULL, variance.targeting = FALSE), 
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE, 
                    external.regressors = NULL), 
  distribution.model = "std")
fit2 <- ugarchfit(spec = spec2, data = yret, solver = "solnp")
np2 <- NROW(fit2@fit$coef)
#yret[-1] quando ho tutte le osservazioni

cat( "\nInformation Criteria" )
print( infocriteria(fit2) )


cat("\nMatcoef\n")
print( fit2@fit$matcoef )


cat("\nRobust matcoef\n")
print( fit2@fit$robust.matcoef )







#LA DIFFERENZA TRA I DUE MODELLI IN UN GRAFICO NEWS IMPACT CURVE


ni1 <- newsimpact(z = NULL, fit1)
ni2 <- newsimpact(z = NULL, fit2)
legend <- c("Simple-GARCH", "GJR-GARCH")
col  <- c("black", "red")
ylim <- range( ni1$zy, ni2$zy )
par(mfrow = c(1,1), mar = c(4, 4.5, 3, 1) + 0.1, lwd = 2)
plot(x = ni1$zx, y = ni1$zy, ylab = ni1$yexpr, xlab = ni1$xexpr, type = "l", 
     ylim = ylim, main = "News Impact Curve", col = col[1])
lines(x = ni2$zx, y = ni2$zy, col = col[2], lwd = 2)
legend(x = "topright", y = NULL, legend = legend, border = FALSE, col = col, 
       lty = 1, text.col = col)


#CURVA NIC OVVERO QUELLO CHE HA FATTO IERI IL MERCATO
#E' ANDATO BENE ALLORA >0 SENNO <0
#SIMPLE GARCH NON E' ASIMMETRICA MA A FORMA DI PARABOLA
#GJR GARCH C'E' ASIMMETRIA
#NEL MIO CASO NEI NEGATVI REASCISCE IN MODO PIU RAPIDO RISPETTO AI POSITIVI
#CON UT-1=-3 LE CURVE HANNO UN VALORE VICINO A 5


#### TEST DI NYBLOM SULLA STABILITA'
cat("\nStability check (Nyblom test)\n")
print( nyblom(fit2) )

#L'UNICO COEFFICENTE SIGNFICATIVO E' MU E GLI ALTRI SI ACCETTANO TUTTI 
#QUINDI VUOL DIRE CHE IL VETTORE E' STABILE E SI PROCEDE A DIMINUIRE LA FONTE DI DATI 
#DAL 2011 AL 2015
#ORA CHE HO DIMINUITO LE OSSERVAZIONI I PARAMETRI SONO QUASI TUTTI SIGNIFICATIVI E QUINDI HO STABILITA'



################################################################################
## Stationarity and Variance Targeting
################################################################################

#### Include variance targeting in GJR
spec3 <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1), 
                                          submodel = NULL, external.regressors = NULL, variance.targeting = TRUE),  
                    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE, arfima = FALSE, 
                                      external.regressors = NULL, archex = FALSE), distribution.model = "std")
fit3 <- ugarchfit(spec = spec3, data = yret, solver = "solnp")
## Store the number of parameters
np3 <- NROW(fit3@fit$coef)

#### Compare the two fits
cat("\nGJR without VT\n")
print( infocriteria(fit2) )
print(fit2@fit$robust.matcoef)
cat("\nGJR with VT\n")
print( infocriteria(fit3) )
print(fit3@fit$robust.matcoef)






################################################################################
## Alternative GARCH formulations via fGARCH: GJRGARCH and T-GARCH
################################################################################

#### GJR using fGARCH
spec4 <- ugarchspec(
  variance.model = list(model = "fGARCH", garchOrder = c(1, 1), 
                        submodel = "GJRGARCH", external.regressors = NULL, variance.targeting = FALSE),  
  mean.model = list(armaOrder = c(1, 0), include.mean = TRUE,  
                    external.regressors = NULL), 
  distribution.model = "std")
fit4 <- ugarchfit(spec = spec4, data = yret, solver = "solnp")
## Store the number of parameters
np4 <- NROW(fit4@fit$coef)
#### Conversion to the "traditional" GJR form
fit4c <- .fgarch.2.gjr(fit = fit4)
#### Comparison
cat("\n\ngjrGARCH vs fGARCH(GJRGARCH)\n")
cat("Direct GJR\n")
print( infocriteria(fit2) )
print(fit2@fit$robust.matcoef)
cat("GJR via fGARCH\n")
print( infocriteria(fit4) )
print(fit4c$robust.matcoef)




#### TGARCH (using fGARCH) 
spec5 <- ugarchspec(
  variance.model = list(model = "fGARCH", garchOrder = c(1, 1), 
                        submodel = "TGARCH", external.regressors = NULL, variance.targeting = FALSE),  
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE, 
                    external.regressors = NULL), 
  distribution.model = "std")
fit5 <- ugarchfit(spec = spec5, data = yret, solver = "solnp")
## Store the number of parameters
np5 <- NROW(fit5@fit$coef)
#### Conversion to the "traditional" GJR form
fit5c <- .fgarch.2.gjr(fit = fit5)
#### Coefficient comparison
cat("\n\nGJR-GARCH vs T-GARCH\n")
cat("GJR-GARCH\n")
print( infocriteria(fit4) )
print(fit4c$robust.matcoef)
cat("T-GARCH\n")
print( infocriteria(fit5) )
print(fit5c$robust.matcoef)





#### Compare the News Impact Curves (sGARCH vs gjrGARCH vs TGARCH) 
ni1 <- newsimpact(z = NULL, fit1)
ni2 <- newsimpact(z = NULL, fit2)
ni5 <- newsimpact(z = NULL, fit5)
legend <- c("Simple-GARCH", "GJR-GARCH", "T-GARCH")
col  <- c("black", "red", "blue")
ylim <- range( ni1$zy, ni2$zy, ni5$zy)
par(mfrow = c(1,1), mar = c(4, 4.5, 3, 1) + 0.1)
plot(x = ni1$zx, y = ni1$zy, ylab = ni1$yexpr, xlab = ni1$xexpr, type = "l", 
     ylim = ylim, main = "News Impact Curve", col = col[1])
lines(x = ni2$zx, y = ni2$zy, col = col[2])
lines(x = ni5$zx, y = ni5$zy, col = col[3])
legend(x = "topright", y = NULL, legend = legend, border = FALSE, col = col, 
       lty = 1, text.col = col)

#### Use standardized residuals!
fit <- fit4
par(mfrow = c(3,1))
Acf(x = fit@fit$z, lag.max = 100, type = "correlation", main = "z")
Acf(x = abs(fit@fit$z), lag.max = 100, type = "correlation", main = "|z|")
Acf(x = fit@fit$z^2, lag.max = 100, type = "correlation", main = expression(z^2))
cat("\nLjung-Box statistics on z residuals\n")
lag1 <- np5 + c(1, 2, 5, 10, 15, 20)
cat("\nLjung-Box on standardized residuals:\n")
lb1 <- mapply(FUN = Box.test, lag = lag1, 
              MoreArgs = list(x = fit@fit$z, type = "Ljung-Box", fitdf = np5) )
print(rbind(lag = lag1, lb1[1:3,]))
cat("\nLjung-Box statistics on |z residuals|\n")
lb1 <- mapply(FUN = Box.test, lag = lag1, 
              MoreArgs = list(x = abs(fit@fit$z), type = "Ljung-Box", fitdf = np5) )
print(rbind(lag = lag1, lb1[1:3,]))
cat("\nLjung-Box statistics on (z residuals)^2\n")
lb1 <- mapply(FUN = Box.test, lag = lag1, 
              MoreArgs = list(x = fit@fit$z^2, type = "Ljung-Box", fitdf = np5) )
print(rbind(lag = lag1, lb1[1:3,]))

#### Independence test on ln(|standardized residuals|)
x1 <- log( abs(fit@fit$z) )
bds <- bds.test(x = x1, m = 4, 
                eps = seq(from = 0.5 * sd(x1), to = 2 * sd(x1), length = 4),
                trace = FALSE)
cat("BDS test on log(abs(z residuals))\n")
print(bds)
#SI RIFIUTA HO E LE VARIABILI SONO IID









################################################################################
## Alternative GARCH specifications: iGARCH
################################################################################

#### IGARCH 
spec6 <- ugarchspec(variance.model = list(model = "iGARCH", garchOrder = c(1, 1), 
                                          submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),  
                    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE, arfima = FALSE, 
                                      external.regressors = NULL, archex = FALSE), distribution.model = "std")
fit6 <- ugarchfit(spec = spec6, data = yret, solver = "solnp")
## Store the number of parameters
np6 <- NROW(fit6@fit$coef)
## Some statistics
cat( "\nInformation Criteria" )
print( infocriteria(fit6) )
cat("\nMatcoef\n")
print( fit6@fit$matcoef )
cat("\nRobust matcoef\n")
print( fit6@fit$robust.matcoef )






################################################################################
## Forecasting ability
##
## Details: 
##  1. Volatility estimates from GARCH are compared with an external benchmark: 
##    the Garman-Klass volatility, a measure of volatility less noisy than 
##    squared returns (note that a GARCH, at the end, is a model of squared or 
##    absolute returns).
##  2. The section below shows IS (in-sample) forecasts, in the sense that the
##    forecasting period is included in the estimation period. More genuine OOS
##    (out-of-sample) forecasts should be considered in the comparison
##   
################################################################################





#### External benchmark
y  <- data$gkVol * 100
#GKVOL LO AVEVO GIA CALCOLATO ALL'INIZIO PRENDENDO COME RIFERIMENTO LA DEVIAZIONE STANDARD

####  PLOT PER CONFRONTARE GARMANKLASS RISPETTO AI RENDIMENTI IN VALORE ASSOLUTO 
#SI USA IL VALORE ASSOLUTO PERCHE CON AL QUADRATO VENGONO PEGGIO
par(mfrow = c(2,1), lwd = 1)
plot(x = time, y = y, type = "l", main = "Garman-Klass volatility measure", ylab="")
plot(x = time, y = sqrt(2/pi)*abs(yret), type = "l", ylab = "", main="Absolute returns")


par(mfrow = c(1,1), lwd = 1)
plot(x = time, y = y, type = "l", ylab = "Garman-Klass volatility measure")
lines(x = time, y = fit5@fit$sigma, col = "red")
lines(x = time, y = fit6@fit$sigma, col = "blue")
legend("topright", inset=0, title="Volatility measure:", 
       legend=c("Garman-Klass","T-GARCH", "I-GARCH"), 
       col=c("black","red","blue"), lty=1, cex=0.70)
#GK HA DEI COMPORTAMENTI SIMILI AL VALORE ASSOLUTO PERO HA UN ANDAMENTO PIU NITIDO, MENTRE
#I RENDIMENTI IN VALORE ASSOLUTO  HANNO ANDAMENDTO PIU SPORCO




####  FACCIO NAIVE SIA DELLA VOLATILITA CHE DELLA VARIANZA
naive.vol <- sd(yret) #PRENDO DEVIAZIONE STANDARD DEI RENDIMENTI
naive.var <- naive.vol^2 #PRENDO LA VARIANZA CAMPIONARIA PERCHE SONO STAZIONARI

#### Error measures
cat("---------------------------------------------------------------------", 
    "\nError measures\n")

#MISURE DI ERRORE NON TANTO AFFIDABILI PERCHE USANO TUTTI I DATI
ErrorMeas <- data.frame(
  measure = c("Volatility", "Volatility", "Volatility", "Volatility", 
              "Variance", "Variance", "Variance", "Variance"), 
  model = c("GARCH", "GJR-GARCH", "T-GARCH", "IGARCH", 
            "GARCH", "GJR-GARCH", "T-GARCH", "IGARCH"), 
  rbind( 
    .ErrorMeasures(y = y,   fit = fit1@fit$sigma,   naive = naive.vol), 
    .ErrorMeasures(y = y,   fit = fit2@fit$sigma,   naive = naive.vol), 
    .ErrorMeasures(y = y,   fit = fit5@fit$sigma,   naive = naive.vol), 
    .ErrorMeasures(y = y,   fit = fit6@fit$sigma,   naive = naive.vol), 
    .ErrorMeasures(y = y^2, fit = fit1@fit$sigma^2, naive = naive.var), 
    .ErrorMeasures(y = y^2, fit = fit2@fit$sigma^2, naive = naive.var), 
    .ErrorMeasures(y = y^2, fit = fit5@fit$sigma^2, naive = naive.var), 
    .ErrorMeasures(y = y^2, fit = fit6@fit$sigma^2, naive = naive.var) ) ) 
print( ErrorMeas )

#SI GUARDANO LE MISURE SCALATE, SULLA VOLATILITA VINCE TGARCH DI SOLITO, PERCHE MODELLO COSTRUITO APPOSTA
#SULLA VOLATILITA. SULLA VARIANZA E' DIVERSO. SUL SCMAE VINCE TGARCH, SU SCRMSE SONO TUTTI LI LI
#SE ME SI AVVICINA AL MAE DIVENTA IL BIAS E ABBIAMO TOPPATO

#### Mincer-Zarnowitz forecasting diagnostics
cat("---------------------------------------------------------------------", 
    "\nMincer-Zarnowitz\n" )
x1 <- .MincerZarnowitz(y = y, fit = fit1@fit$sigma, msg = "GARCH\n") #rifiuto sia i coefficenti(HAC pvalue significativi) che statistica test f 
x1 <- .MincerZarnowitz(y = y, fit = fit3@fit$sigma, msg = "GJR-GARCH\n") #rifiuto  i coefficenti(HAC pvalue significativi) ma accetto che statistica test f 
x1 <- .MincerZarnowitz(y = y, fit = fit5@fit$sigma, msg = "T-GARCH\n") ##rifiuto sia i coefficenti(HAC pvalue significativi) che statistica test f 
x1 <- .MincerZarnowitz(y = y, fit = fit6@fit$sigma, msg = "IGARCH\n") #rifiuto sia i coefficenti(HAC pvalue significativi) che statistica test f 



#TEST DIEBOLD MARIANO
#### Diebold-Mariano forecasting comparison
cat("---------------------------------------------------------------------", 
    "\nDiebold-Mariano comparison\n\n")
## Volatility
cat("Volatility\n")#PRIMO CONTROLLO E' SULLE VOLATILITA
#Y INDICA GARMAN KLASS IN PERCENTUALI SCALATO SULLE VOLATILITA
h <- 1 #ORIZZONTE SERVE PER IL CALCOLO SE
e1 <- y - fit1@fit$sigma #I FIT SONO LE VOLATILITA STIMATE DAL MODELLO 1 SIMPLE GARCH
e2 <- y - fit2@fit$sigma #MODELLO 2 GJR
e5 <- y - fit5@fit$sigma #MODELLO 5 TGARCH
e6 <- y - fit6@fit$sigma #MODELLO INTEGRATO
.DieboldMariano(e1 = e1, e2 = e2, h = h, power = 1, msg = "GARCH vs GJR-GARCH ->") 
.DieboldMariano(e1 = e1, e2 = e2, h = h, power = 2, msg = "GARCH vs GJR-GARCH ->") 
.DieboldMariano(e1 = e2, e2 = e5, h = h, power = 1, msg = "GJR-GARCH vs T-GARCH   ->") 
.DieboldMariano(e1 = e2, e2 = e5, h = h, power = 2, msg = "GJR-GARCH vs T-GARCH   ->")
.DieboldMariano(e1 = e1, e2 = e6, h = h, power = 1, msg = "GARCH vs IGARCH    ->")
.DieboldMariano(e1 = e1, e2 = e6, h = h, power = 2, msg = "GARCH vs IGARCH    ->")

#POWER STA PER POTENZA A CUI ELEVO L'ERRORE IN VALORE ASSOLUTO
#IL TEST HA ALTERNATIVA BIDIREZIONALE 

## Conditional variance
cat("Conditional variance\n")
h <- 1
e1 <- y^2 - fit1@fit$sigma^2
e2 <- y^2 - fit2@fit$sigma^2
e5 <- y^2 - fit5@fit$sigma^2
e6 <- y^2 - fit6@fit$sigma^2
.DieboldMariano(e1 = e1, e2 = e2, h = h, power = 1, msg = "GARCH vs GJR-GARCH ->") 
.DieboldMariano(e1 = e1, e2 = e2, h = h, power = 2, msg = "GARCH vs GJR-GARCH ->") 
.DieboldMariano(e1 = e2, e2 = e5, h = h, power = 1, msg = "GJR-GARCH vs T-GARCH   ->") 
.DieboldMariano(e1 = e2, e2 = e5, h = h, power = 2, msg = "GJR-GARCH vs T-GARCH  -> ") 
.DieboldMariano(e1 = e1, e2 = e6, h = h, power = 1, msg = "GARCH vs IGARCH    ->") 
.DieboldMariano(e1 = e1, e2 = e6, h = h, power = 2, msg = "GARCH vs IGARCH    ->")


################################################################################
## Forecasts using rugarch
################################################################################

#### Settings
H <- 10  #2 settimane


## data
## out.sample
## n.roll
#### Rule:
## last t in info = NROW(data) - out.sample

#### 1) ex-ante, h = 1:H
forc1 <- ugarchforecast(fitORspec = fit1, n.ahead = H, 
                        data = NULL, out.sample = 0, n.roll = 0)
#forc1@forecast
#plot(forc1)

####Â 2) ex-post, h = 1:H at t = 1661 
spec1x <- getspec(fit1)
setfixed(spec1x) <- as.list(coef(fit1))
# forc2a <- ugarchforecast(fitORspec = spec1x, n.ahead = H, 
#  data = yret[1:1761, , drop = FALSE])
forc2 <- ugarchforecast(fitORspec = spec1x, n.ahead = H, 
                        data = yret, out.sample = NROW(yret)-1661, n.roll = 0)
#forc2@forecast
#plot(forc2)

####Â 3) ex-post, h = 1:H at t = 1661 
forc3 <- ugarchforecast(fitORspec = spec1x, n.ahead = H,  
                        data = yret, out.sample = NROW(yret) - 1661, n.roll = 10)
#forc3@forecast


