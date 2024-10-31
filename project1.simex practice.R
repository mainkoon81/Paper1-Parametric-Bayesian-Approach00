############################## simex practice ##################################
library(simex)

#################### Use testset......
############################# family default ###################################
# - gaussian GLM:	     gaussian(link = "identity")
# - Gamma GLM:	       Gamma(link = "inverse")
# - binomial GLM:	     binomial(link = "logit")
# - quasi-binomialGLM: quasibinomial(link = "logit")
# - poisson GLM:       poisson(link = "log")
# - quasi-poissonGLM:	 quasipoisson(link = "log")
# - quasi-GLM:	       quasiglm(link = "identity", variance = "constant")
################################################################################

##################################[note]########################################
# ## error for SD @: Homoscedastic? --------------------------------------------
# model_real <- lm(y ~ x_real)
# 
# model_naiv <- lm(y ~ x_measured, x = T)
# 
# model_simex <- simex(model_naiv, 
#                      SIMEXvariable = "x_measured", 
#                      measurement.error = [single value],
#                      asymptotic = T)
# plot(model_simex)
# 
# ## error for SD @: Heteroscedastic? ------------------------------------------
# for example,
#
# model_real <- lm(y ~ x_real)
#
# now what is the variance of your error?
# alphaa1 <- sort( runif( [n], min = 0, max = 0.6 ) ) # for example
# alphaa2 <- abs( sort( runif( [n], min = 0, max = 0.6 ) ) + rnorm(100, sd = 0.1) ) #for example
#
#       ************************************************************
#       **** x_measured.H = x_real + [alphaaaa]*rnorm([n], 0,1) ****
#       ************************************************************
# mod_naiv.H <- lm(y ~ x_measured.H, x = T)
# mod_simex.H <- simex(mod_naiv.H, 
#                      SIMEXvariable = "x_measured.H",
#                      measurement.error = [alphaaaa], 
#                      asymptotic = F)
# plot(model_simex.H) # 
################################################################################







#################
#### Brazil  ###################################################################
#################
test.bra <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/Brazil.paper3-1.test.csv",
                     header=T, 
                     na.strings=c("."), 
                     stringsAsFactors=F)
head(test.bra)
summary(test.bra$AggClaim) # no zero
summary(test.bra$PremTotal ) # no zero
cor(test.bra$AggClaim, test.bra$PremTotal)       # cor: 0.6
cor(test.bra$AggClaim, test.bra$ln_Expos)      # cor: 0.3
cor(log(test.bra$AggClaim), test.bra$PremTotal)  # cor: 0.3
cor(log(test.bra$AggClaim), test.bra$ln_Expos) # cor: 0.35

### > GLM set-up +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
model_real <- glm(AggClaim ~ ln_Expos + I(as.factor(Affiliation)), family=Gamma(link="inverse"), data=test.bra, x=T, y=T)
summary(model_real) #AIC: 46124
################################################################################
model_real <- glm(log(AggClaim) ~ ln_Expos + I(as.factor(Affiliation)), family=gaussian(link="identity"), data=test.bra, x=T, y=T)
summary(model_real) # **********************************************************
#AIC: 10928
################################################################################
# so... who is the winner?  "Log(AggClaim) with gaussian(link="identity")"
density(predict(model_real))
plot( x=density(predict(model_real)) )
# log(AggClaim)   predictive density            
# Min.   : 2.597   Min.   :0.0000101  
# 1st Qu.: 4.667   1st Qu.:0.0026518  
# Median : 6.737   Median :0.0385264  
# Mean   : 6.737   Mean   :0.1206688  
# 3rd Qu.: 8.807   3rd Qu.:0.2321782  
# Max.   :10.876   Max.   :0.3879934  

model_real$coefficients
  # (Intercept)      ln_Expos   I(as.factor(Affiliation))1
  #   8.5917477     1.0492582                   -2.2509698
    
summary(model_real$fitted.values)
#   Min. 1st Qu.  Median  Mean   3rd Qu.    Max. 
# 3.151   6.470   7.165   7.162   7.836  10.323 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++ Now, what's the reality? ++++++++++++++++++++++++++ #

# -> simex is so stupid, to build your "naiv model, u cannot use log(.), I(.) on your variable...

## > Build a dirty model
Y <- log(test.bra$AggClaim)
Xstar <-  test.bra$ln_Expos_err
Z <- I( as.factor(test.bra$Affiliation) )

model_naiv <- glm(Y ~ Xstar + Z, family=gaussian(link="identity"), x=T, y=T)
model_naiv$coefficients
 # (Intercept)      Xstar           Z1
 #   8.3623894  0.6450517   -2.2152368
   
model_naiv$fitted.values
plot( x=density(model_naiv$fitted.values) )



## > what is your best guess on "ln_insured_err"?
alphaa3 <- rep(1, nrow(test.bra))
for (i in 1:nrow(test.bra)) {
  if(test.bra$Affiliation[i]==0){
    alphaaa <- ifelse( test.bra$PremTotal[i]<77600 & test.bra$PremTotal[i]>225, 0.1, 0.4 )
    delta_sq <- ifelse( test.bra$ExposTotal[i]<=0.5, 1, 
                        ifelse(test.bra$ExposTotal[i]>35, 100, 3) )
    alphaa3[i] <- alphaaa*delta_sq
  }
  else{
    alphaa3[i] <- 1
  }
} # alphaa3 = alphaaa*delta_sq

### RUN ------------------------------------------------------------------------
mod_simex.H <- simex(model_naiv,
                     SIMEXvariable = "Xstar",
                     measurement.error = alphaa3,
                     fitting.method = "quadratic",
                     B = 100,
                     asymptotic = F)
par(mfrow=c(2,3))
plot(mod_simex.H)
par(mfrow=c(1,2))

summary(mod_simex.H)
print(mod_simex.H)
# (Intercept)        Xstar          Z1
#     8.5228       0.9284      -2.8060

    
predict(mod_simex.H)
density(predict(mod_simex.H))
plot( x=density(predict(mod_simex.H)) )

p.glm <- predict(mod_simex.H)
plot(x=density( p.glm ), xlab="Log(Y)", main = "predictive density by GLM (gamma)" )
polygon(x=density( p.glm ), col = "grey" )
# abline(v = mean(exp(Y)), col = "red") # so ugly.....
# lines(x=density(x=exp(Y)), col = "green", lty=5, lwd=0.5) # empirical
# lines(x=density(x=exp(Y), bw=4), col="blue", lwd=0.5) # Gaussian kernel density a bandwidth equal to 4: 

SSPE.glm <- sum( (p.glm - log(test.bra$AggClaim))^2 ); SSPE.glm            # 12400.17
SAPE.glm <- sum( abs(p.glm - log(test.bra$AggClaim)) ); SAPE.glm           # 4273.977




qqplot( x=predict(mod_simex.H), y=model_real$fitted.values, 
        xlab="SIMEX qunatile", ylab="Naiv quantile" )
abline(a=0, b=1, col="red")
















#################
#### French  ###################################################################
#################
test.fre <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BFren.paper3-1.test.csv",
                     header=T, 
                     na.strings=c("."), 
                     stringsAsFactors=F)
head(test.fre)
summary(test.fre$ClaimAmount) # no zero
summary(test.fre$Density) # no zero
cor(test.fre$ClaimAmount, test.fre$Density)       # cor: 0.97
cor(test.fre$ClaimAmount, test.fre$ln_expo)      # cor: 0.53
cor(log(test.fre$ClaimAmount), test.fre$Density)  # cor: 0.45
cor(log(test.fre$ClaimAmount), test.fre$ln_expo) # cor: 0.86

### > GLM set-up +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
model_real <- glm(ClaimAmount ~ ln_expo + I(as.factor(CarAge)), family=Gamma(link="identity"), data=test.fre, x=T, y=T)
summary(model_real)
#AIC: 36848
################################################################################
model_real <- glm(log(ClaimAmount) ~ ln_expo + I(as.factor(CarAge)), family=gaussian(link="log"), data=test.fre, x=T, y=T)
summary(model_real)
#AIC: 6566.7
model_real <- glm(log(ClaimAmount) ~ ln_expo + I(as.factor(CarAge)), family=gaussian(link="identity"), data=test.fre, x=T, y=T)
summary(model_real) # **********************************************************
#AIC: 6566.6
model_real <- glm(log(ClaimAmount) ~ ln_expo + I(as.factor(CarAge)), family=gaussian(link="inverse"), data=test.fre, x=T, y=T)
summary(model_real)
#AIC: 6566.8
################################################################################

# so... who is the winner?  "Log(ClaimAmount) with gaussian(link="identity")"
density(predict(model_real))
plot( x=density(predict(model_real)) )
# Log(ClaimAmount) predictive density            
# Min.   :6.684   Min.   : 0.000211  
# 1st Qu.:6.861   1st Qu.: 0.027266  
# Median :7.037   Median : 0.401314  
# Mean   :7.037   Mean   : 1.412522  
# 3rd Qu.:7.214   3rd Qu.: 1.675482  
# Max.   :7.391   Max.   :12.439175 

model_real$coefficients
# (Intercept)       ln_expo   I(as.factor(CarAge))1
# 6.79761531      -0.09530088           -0.08285480

# model_real$coefficients
# (Intercept)       ln_expo   I(as.factor(CarAge))1
# 6.83739987    -0.01876840             -0.07614282


summary(model_real$fitted.values)
#   Min. 1st Qu.  Median  Mean   3rd Qu.    Max. 
#  6.715   6.789   6.798   6.821   6.860   7.360 

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++ Now, what's the reality? ++++++++++++++++++++++++++ #

# -> simex is so stupid, to build your "naiv model, u cannot use log(.), I(.) on your variable...

## > Build a dirty model
Y <- log(test.fre$ClaimAmount)
Xstar <-  test.fre$ln_expo_err
Z <- I( as.factor(test.fre$CarAge) )

model_naiv <- glm(Y ~ Xstar + Z, family=gaussian(link="identity"), x=T, y=T)
model_naiv$coefficients
# (Intercept)       Xstar          Z1
# 6.83739987  -0.01876840 -0.07614282

# model_naiv$coefficients
# (Intercept)       Xstar           Z1
# 6.79761531  -0.09530088  -0.08285480


## > what is your best guess on "ln_expo_err"?
alphaa3 <- rep(1, nrow(test.fre))
for (i in 1:nrow(test.fre)) {
  if(test.fre$CarAge[i]==1){
    alphaaa <- ifelse( test.fre$Density[i]>=17000, 20, 
                       ifelse(test.fre$Density[i]<=20, 2, 10) )
    delta_sq <- ifelse( test.fre$Exposure[i]<1 & test.fre$Exposure[i]>0.07, 0.1, 0.2 )
    alphaa3[i] <- alphaaa*delta_sq
  }
  else{
    alphaa3[i] <- 1
  }
} # alphaa3 = alphaaa*delta_sq


### RUN ------------------------------------------------------------------------
mod_simex.H <- simex(model_naiv,
                     SIMEXvariable = "Xstar",
                     measurement.error = alphaa3,
                     fitting.method = "quadratic",
                     B = 100,
                     asymptotic = F)
par(mfrow=c(2,3))
plot(mod_simex.H)
par(mfrow=c(1,2))

summary(mod_simex.H)

print(mod_simex.H)
# (Intercept)       Xstar           Z1
#   6.82880      -0.03537     -0.07825

predict(mod_simex.H)
density(predict(mod_simex.H))
plot( x=density(predict(mod_simex.H)) )

p.glm <- predict(mod_simex.H)
plot(x=density( p.glm ), xlab="Log(Y)", main = "predictive density by GLM (gamma)" )
polygon(x=density( p.glm ), col = "grey" )
# abline(v = mean(exp(Y)), col = "red") # so ugly.....
# lines(x=density(x=exp(Y)), col = "green", lty=5, lwd=0.5) # empirical
# lines(x=density(x=exp(Y), bw=4), col="blue", lwd=0.5) # Gaussian kernel density a bandwidth equal to 4: 

SSPE.glm <- sum( (p.glm - log(test.fre$ClaimAmount))^2 ); SSPE.glm            # 2769.916
SAPE.glm <- sum( abs(p.glm - log(test.fre$ClaimAmount)) ); SAPE.glm           # 1599.461










#################
#### Swedish ###################################################################
#################
test.swe <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/BSwed.paper3-1.test.csv",
                     header=T, 
                     na.strings=c("."), 
                     stringsAsFactors=F)

head(test.swe)
summary(test.swe$TotalLoss) # no zero
summary(test.swe$SumInsAvg) # no zero
cor(test.swe$TotalLoss, test.swe$SumInsAvg)       # cor: 0.97
cor(test.swe$TotalLoss, test.swe$ln_insured)      # cor: 0.53
cor(log(test.swe$TotalLoss), test.swe$SumInsAvg)  # cor: 0.45
cor(log(test.swe$TotalLoss), test.swe$ln_insured) # cor: 0.86

### > GLM set-up +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
model_real <- glm(TotalLoss ~ ln_insured + I(as.factor(Experience)), family=Gamma(link="inverse"), data=test.swe, x=T, y=T)
summary(model_real)
#AIC: ?
model_real <- glm(TotalLoss ~ ln_insured + I(as.factor(Experience)), family=gaussian(link="log"), data=test.swe, x=T, y=T)
summary(model_real)
#AIC: 6644.1
model_real <- glm(TotalLoss ~ ln_insured + I(as.factor(Experience)), family=gaussian(link="inverse"), data=test.swe, x=T, y=T)
summary(model_real)
#AIC: 7202.9
################################################################################
model_real <- glm(log(TotalLoss) ~ ln_insured + I(as.factor(Experience)), family=gaussian(link="log"), data=test.swe, x=T, y=T)
summary(model_real)
#AIC: 662.72
model_real <- glm(log(TotalLoss) ~ ln_insured + I(as.factor(Experience)), family=gaussian(link="identity"), data=test.swe, x=T, y=T)
summary(model_real) # **********************************************************
#AIC: 670.11
model_real <- glm(log(TotalLoss) ~ ln_insured + I(as.factor(Experience)), family=gaussian(link="inverse"), data=test.swe, x=T, y=T)
summary(model_real)
#AIC: 687.06
################################################################################

# so... who is the winner?  "Log(TotalLoss) with gaussian(link="identity")"
density(predict(model_real))
plot( x=density(predict(model_real)) )
# log(TotalLoss) predictive density            
#Min.   : 4.566   Min.   :4.851e-05  
#1st Qu.: 7.822   1st Qu.:9.121e-03  
#Median :11.077   Median :2.892e-02  
#Mean   :11.077   Mean   :7.672e-02  
#3rd Qu.:14.333   3rd Qu.:1.572e-01  
#Max.   :17.588   Max.   :2.211e-01  

model_real$coefficients
#   (Intercept)       ln_insured   I(as.factor(Experience))1 
#    5.9971451        0.9564930              -0.5362616    

model_real$coefficients
  (Intercept)   ln_insured   I(as.factor(Experience))1
   1.7552229     1.7791422                  -1.3157922





summary(model_real$fitted.values)
#   Min. 1st Qu.  Median  Mean   3rd Qu.    Max. 
# 6.085   9.295  10.362  10.460  11.560  16.070 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++ Now, what's the reality? ++++++++++++++++++++++++++ #

# -> simex is so stupid, to build your "naiv model, u cannot use log(.), I(.) on your variable...

## > Build a dirty model
Y <- log(test.swe$TotalLoss)
Xstar <-  test.swe$ln_insured_err
Z <- I( as.factor(test.swe$Experience) )

model_naiv <- glm(Y ~ Xstar + Z, family=gaussian(link="identity"), x=T, y=T)
model_naiv$coefficients
#  (Intercept)        Xstar           Z1 
#   1.7552229     1.7791422   -1.3157922

# model_naiv$coefficients
#  (Intercept)        Xstar           Z1 
#    5.9971451    0.9564930   -0.5362616   

## > what is your best guess on "ln_insured_err"?
alphaa3 <- rep(1, nrow(test.swe))
for (i in 1:nrow(test.swe)) {
  if(test.swe$Experience[i]==0){
    alphaaa <- ifelse( test.swe$SumInsAvg[i]<9800 & test.swe$SumInsAvg[i]>7.5, 0.1, 0.2 )
    delta_sq <- ifelse( test.swe$Zone[i]<=2, 1, 
                        ifelse(test.swe$Zone[i]>4, 100, 3) )
    alphaa3[i] <- alphaaa*delta_sq
  }
  else{
    alphaa3[i] <- 1
  }
} # alphaa3 = alphaaa*delta_sq




### RUN ------------------------------------------------------------------------
mod_simex.H <- simex(model_naiv,
                     SIMEXvariable = "Xstar",
                     measurement.error = alphaa3,
                     fitting.method = "quadratic",
                     B = 100,
                     asymptotic = F)
par(mfrow=c(2,3))
plot(mod_simex.H)
par(mfrow=c(1,2))

summary(mod_simex.H)

print(mod_simex.H)
# (Intercept)       Xstar           Z1
#      1.922        1.943       -1.872

print(mod_simex.H)


predict(mod_simex.H)
density(predict(mod_simex.H))
plot( x=density(predict(mod_simex.H)) )

p.glm <- predict(mod_simex.H)
plot(x=density( p.glm ), xlab="Log(Y)", main = "predictive density by GLM (gamma)" )
polygon(x=density( p.glm ), col = "grey" )
# abline(v = mean(exp(Y)), col = "red") # so ugly.....
# lines(x=density(x=exp(Y)), col = "green", lty=5, lwd=0.5) # empirical
# lines(x=density(x=exp(Y), bw=4), col="blue", lwd=0.5) # Gaussian kernel density a bandwidth equal to 4: 

SSPE.glm <- sum( (p.glm - log(test.swe$TotalLoss))^2 ); SSPE.glm            # 402.3893
SAPE.glm <- sum( abs(p.glm - log(test.swe$TotalLoss)) ); SAPE.glm           # 240




























# [residual VS covariate]: plots for predictors (partial residuals, coeff)
termplot(model_real, partial.resid = T, se = T, col.res="blue")


# [Deviance residual VS theoretical quantile] 
# Under certain conditions, the deviance residuals from a GLM are approximately N(0, 1). 
# If those conditions are satisfied, then the median of deviance residuals should be close to 0 and 
# the minimum and maximum values of the deviance residuals should be close to -3 and +3, respectively.
AIC(model_real)
library(statmod)
qresid(model_real, dispersion=NULL)
qqnorm(residuals(model_real, type="deviance"), col="blue")
abline(0,1, col="red")


# [surface plot for two continuous predictor]
head(test.swe)

fit.cont <- glm(log(TotalLoss) ~ ln_insured_err + SumInsAvg, 
                family=gaussian(link="identity"), 
                data=test.swe, x=T, y=T)
summary(fit.cont)

# First, let's define functions
pglm <- function(x1,x2){
  return( predict(fit.cont, newdata=data.frame(ln_insured_err=x1, SumInsAvg=x2), type="response") )
  }

M <- dim(test.swe)[1]
cx1 <- seq(from=min(test.swe$SumInsAvg), to=max(test.swe$SumInsAvg), length=M); cx1 #continuous
cx2 <- seq(from=min(test.swe$ln_insured_err), to=max(test.swe$ln_insured_err), length=M); cx2 #continuous
z=outer(X=cx1, Y=cx2, FUN=pglm)

persp( x=cx1, y=cx2, z=outer(X=cx1, Y=cx2, FUN=pglm), col="lightblue", border=NA, theta=-30 ) # perpective plot
surf <- persp( x=cx1, y=cx2, z=outer(X=cx1, Y=cx2, FUN=pglm), col="lightblue", border=NA, theta=-30 )   
# draw lines parallel to x axis. seq(...) depends on your data's length(y)
for(i in seq(10, M, length=10)) lines(trans3d(cx1, cx2[i], z[,i], pmat = surf), col = "blue")
# draw lines parallel to y axis. seq(...) depends on your data's length(x)
for(i in seq(10, M, length=10)) lines(trans3d(cx1[i], cx2, z[i,], pmat = surf), col = "blue")
contour( x=cx1, y=cx2, z=outer(X=cx1, Y=cx2, FUN=pglm) ) # contour plot















































