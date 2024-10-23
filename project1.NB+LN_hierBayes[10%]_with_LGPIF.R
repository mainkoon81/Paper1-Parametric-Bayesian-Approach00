########### Bayes Parametric: Hierarchical LN-GLM [BASIC] with NDB ############# see 1155 for zeta





##### LGPIF Data with overall error: [10%]-[10%]-[10%]-[10%]-[10%]-[10%]-[10%]-[10%]-[10%]
LGPIF <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/LGPIF.paper3-10.csv")
#LGPIF <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/LGPIF.paper3-10+.csv")
head(LGPIF)
LGPIF$Classes <- as.integer(as.factor(LGPIF$Classes))
train.df <- subset(x=LGPIF, subset=Year<2010); head(train.df)
#summary(train.df)
test.df <- subset(x=LGPIF, subset=Year==2010); head(test.df)
#summary(test.df)

n.train <- nrow(train.df); n.train
n.test <- nrow(test.df); n.test

#N:   Freq
#Y:   yAvg
#X_F: { z:AC15, x:Premium } 
#X_S: { z:Fire5 , x:lnDeduct } 
#Extra:  Year, y, TypeCounty, LnCoverage

#%%%% train %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### Variable Definitions (S, z) for F
N = train.df$Freq 
ZF = train.df$AC15
##### Variable Definitions (S, z) for S
Y = train.df$yAvg
ZS = train.df$Fire5
##### Dirty data
XSstar = train.df$lnDeduct_err
##### Clean data
XF = train.df$Premium*0.001
XS = train.df$lnDeduct

##### COVARIATE Definition (z, x)
matXF = as.matrix(cbind(1, XF, ZF))

matXSstar = as.matrix(cbind(1, XSstar, ZS))
matXS = as.matrix(cbind(1, XS, ZS))

n.breaks.train = sqrt( nrow(train.df) ) #******Rule of thumb

#%%%% test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##### Variable Definitions (S, z) for F
N.test = test.df$Freq 
ZF.test = test.df$AC15
##### Variable Definitions (S, z) for S
Y.test = test.df$yAvg
ZS.test = test.df$Fire5
##### Dirty data
XSstar.test = test.df$lnDeduct_err

##### Clean data
XF.test = test.df$Premium*0.001
XS.test = test.df$lnDeduct

##### COVARIATE Definition (z, x)
matXF.test = as.matrix(cbind(1, XF.test, ZF.test))
matXS.test = as.matrix(cbind(1, XS.test, ZS.test))

matXSstar.test = as.matrix(cbind(1, XSstar.test, ZS.test))

n.breaks.test = sqrt( nrow(test.df) ) #******Rule of thumb

#-----JUST STUDY---------------------------------------------------------------#
# .... Negative Binomial Density
# x <- rnbinom(n=?, size = ?, mu = ?)
# hist(x, freq = FALSE, col = "gray", main = "Negative Binomial Density")
# curve( dnbinom(x, size = size, mu = mu), add = TRUE, lwd = 2, col = "red" )

rnbinom(n=5, size=2, mu=5) # for N
dnbinom(x=rnbinom(n=5, size=2, mu=5), size=2, mu=5) 

plot( dnbinom(x=rnbinom(n=5, size=2, mu=5), size=2, mu=5) )

# .... Lognormal Density
x = seq(0, 100, length.out = 1000) # Y ~ LogN
mu = 2
sigma = 5

rlnorm(n=5, meanlog=mu, sdlog=sigma) # for Y
plot(x=x, y=dlnorm(x=x, meanlog=mu, sdlog=sigma), type = "l", xlab = "x", ylab = "Density", main = "Log-Normal Distribution")
curve(dnorm(x=x, mean=mu, sd=sigma), add = TRUE, col = "red" ) # for logY ~ N
#-------------------------------------------------------------------------------

# .... Negative Binomial Density for "N"
mk.dnbinom = function(N, x, beta, psi) {
  mu = exp(sum(x*beta)) # x^T*B
  p = psi/(mu+psi) # in case numerator is too small...1/(1+exp(log(psi)-log(mu)))
  return( choose(N+psi-1, N) * (1-p)^N*(p)^psi )
}

mk.dnbinom_log = function(N, x, beta, psi) {
  mu = exp(sum(x*beta)) # x^T*B
  p = psi/(mu+psi) # in case numerator is too small...1/(1+exp(log(psi)-log(mu)))
  return( log( choose(N+psi-1, N) ) + N*log(1-p) + psi*log(p) ) 
}

# .... normal Density for "log(Y)"
dlognorm = function(y, x, beta, sig2) {
  mu = sum(x*beta) - 0.5*sig2 
  return( (1/(y*sqrt(2*pi*sig2)))*exp(-((log(y)-mu)^2 / (2*sig2))) )
}

# ... log of normal Density for "log(Y)"
dlognorm_log = function(y, x, beta, sig2) {
  mu = sum(x*beta) - 0.5*sig2
  return( log(1) - log(y*sqrt(2*pi*sig2)) - ((log(y)-mu)^2 / (2*sig2)) )
}
#-------------------------------------------------------------------------------

#] Initial Clustering Based on Classes
cl_membership <- train.df$Classes
table(cl_membership)

JF = length(table(cl_membership));JF
JS = length(table(cl_membership));JS
cl_membershipF = train.df$Classes
cl_membershipS = train.df$Classes









#-------------------------- Experiments start here! ---------------------------#
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
########################### with Dirty data [10%] ##############################
################################################################################
### Step01> Define prior, hyperprior model and hyperparameter values
################################################################################

#] ------ Outcome ------ ::: for logY ~ N(X*betaparamS, sig2paramS)
#------------------------------------------------------------------------------------------------------------------------


# PRIOR: "betaparam" ~ MVN(beta0S, SIG_b0S)
# Hyperprior               beta0S ~ MVN(m0S, 1/deltaS*SIG_b0S)

#> beta0
#> Gaussian? Gamma? regression for initialize OUTCOME parameter, using 'm0', 'SIG_b0' 
#fit1S <- glm( log(Y) ~ XS + factor(ZS), family=gaussian(link = "identity")) # AIC: 4308.6
#summary(fit1S)
#vcov(fit1S)

#fit1S <- glm( Y ~ XS + factor(ZS), family=Gamma(link = "log")) # AIC: 27056
#summary(fit1S)
#vcov(fit1S)
#fit1S.xx <- glm( Y ~ XSstar + factor(ZS), family=Gamma(link = "log"))
fit1S.xx <- glm( Y ~ XSstar + factor(ZS), family=poisson(link = "log")) # AIC: ?
summary(fit1S.xx)
vcov(fit1S.xx) # too small...


m0S.xx = coef(fit1S.xx)       # 1x3 initial reg_coeff vector (sort of Regression result as "mean")
SIG_b0S.xx = vcov(fit1S.xx)*varinf*10000   # 3x3 initial cov matrix of reg_coeff
p = length(m0S.xx) # intercept,b1,b2

#****** for VAR Inflation factor
# - hierarchical models can introduce additional sources of variability that need to be accounted for.
# - In hierarchical model, multiple levels of parameters can cause the variance of the estimated parameters to be larger 
# - compared to a simple model without these hierarchical structures. We use this [variation inflation factor] 
deltaS.xx = 0.01             # equal importance on m0 and beta ???? #(1)          #(0.1)          #(0.01)
varinf = n.train/100    # rule of thumb: n.train divide by 5?? #(n.train/10) #(n.train/1000) #(n.train/100)

options(scipen = 999)
SIG_b0S.xx

#SIG_b0... = SIG_b0*varinf ???
SIG_b0invS.xx = solve(a=SIG_b0S.xx) # inverse of cov matrix of reg_coeff for later use (for posterior on reg_coeff)!!
varinf; SIG_b0S.xx; SIG_b0invS.xx



# PRIOR: "betaparamS" ~ MVN(beta0S, SIG_b0S)
# Hyperprior                        SIG_b0S ~ IW( qu0S, LAMBS )


#> SIG_b0S
#- First, SIG_b0S = VCOV(fit1S)       ???????????? Really? 
#- Next,  SIG_b0S ~ IW( qu0S, LAMBS )  ???????????? Really?

# it requires -------
qu0S.xx = p+2
# LAMB ?????????????????????????????????????????????
LAMBS.xx = SIG_b0S.xx#*varinf # really??????????????????????????????                  #~xxxxxxxxxxxxxxxxxxxxx Q.DAVID
LAMBS.xx


###> PRIOR: "sig2param" ~ IG( u0S, v0S )
# Hyperprior               u0S ~ Fink( rho_u1S, rho_u2S )
# Hyperprior                   v0S ~ Ga( rho_v1S, rho_v2S )

# - Hyperparameter a0 for Fink
rho_u1S.xx = 1/8 #1/8  #1/8 #8 #2                                                  
rho_u2S.xx = 3/2 #2    #1/2 #1 #1

# - Hyperparameter b0 - 2024-02-09
rho_v1S.xx = 8 # 4
rho_v2S.xx = 1 # 1

# - NEED to adjust for "Proposal design" InvGa (shape=u0S, scale=v0S ) ---------#~xxxxxxxxxxxxxxxxxxxxx Q.DAVID
# Based on [Method of Moment]?
# u0_new = 2.561215 # from 2/22 
# v0_new = 2.761394 # from 2/22 

# Based on the result? mean(unlist(list_u0paramS)), mean(unlist(list_v0paramS))
# E[invGa] = scale/(shape-1) = 50
u0_newS.xx = 1.7591 #3 #2 #1.1 # from 3/7  # from # summary(unlist(list_u0paramS)) # for u0_new                            
v0_newS.xx = 7.113  #3 #5 #8.5 # from 3/7  # from # summary(unlist(list_v0paramS)) # for v0_new
#-------------------------------------------------------------------------------


#] ------ Covariate ------ ::: for ZS ~ Bin( 1, [piparamS] )  
# PRIOR                                         "piparamS" ~ Beta(g0S, h0S)
#------------------------------------------------------------------------------------------------------------------------
g0S.xx = 0.5 
h0S.xx = 0.5 

#] ------ Covariate ------ ::: for XS ~ N( XS_bar, [lambda2paramS] ) where 
# PRIOR                                            "lambda2paramS" ~ IG(c0S, d0S)
#------------------------------------------------------------------------------------------------------------------------

#> KAPPA ????
# dirty covariate x is NOT tied to z under the NDB assumption, by this logic, regression of [x|z] and [xstar|z] should return 
# almost the same regression parameter-results...other than parameter-variance...
# It is for the sake of investigation or computation or mathematical convenience ?!?! ...since we need model (x,z)=(x|z)(z)
# It is for the connection between [clean] - [dirty] exposure model ?!?!
# We d not know VAR(Xstar). but by conditioning VAR(Xstar|Z), 
# which is "\hat{lambda}^2", and by comparing VAR(Xstar|Z) vs VAR(X|Z), we can perceive.... tau^2     

#> XS_bar 
#XS_bar = mean(XS) #????

#> lambda2
c0S.xx = 0.5          
d0S.xx = 0.5


# Now, we need a posterior sample of hyperPRIOR...
##################################################################################################
### Step02> initialize major parameters (single value), using analytically driven "POSTERIOR" pdf  
#                                                 Based on natural clustering
##################################################################################################
JS = length(table(cl_membershipS))
set.seed(1)

#] ------ Covariate ZS ------ ::: for ZS ~ Bin( 1, [piparamS_j] ) 
# prior                                            "piparamS_j" ~ Beta(g0S, h0S)
#------------------------------------------------------------------------------------------------------------------------
# Covariate parameter from Posterior
piparamS.xx = numeric(JS)
for (j in 1:JS) {
  ZSj=ZS[cl_membershipS==j]
  nS_j=length(ZSj)
  piparamS.xx[j] = rbeta( n=1, shape1=g0S.xx+sum(ZSj), shape2=h0S.xx+nS_j-sum(ZSj) ) 
}

piparamS.xx #% pi parameter (posterior sample)


#] ------ Covariate XS ------ ::: for XS ~ N( XS_bar, [lambda2paramS] ) where 
# PRIOR                                            "lambda2paramS" ~ IG(c0S, d0S)
#------------------------------------------------------------------------------------------------------------------------
# Covariate parameter from Posterior
lambda2paramS.xx = numeric(JS)
for (j in 1:JS) {
  XSj=XSstar[cl_membershipS==j]
  nS_j=length(XSj)
  lambda2paramS.xx[j]=rinvgamma( n=1, shape=(c0S.xx+nS_j)/2, scale=(d0S.xx + sum( (XSj - mean(XSj))^2 ) )/2 )  
}                                                                                                     

lambda2paramS.xx #% lambda2 parameter (posterior sample)

#% empirical X_barF ? using "aggregate(.)": investigation by splitting the data into subset
X_barS.xx = aggregate(x=XSstar, by=list(cl_membershipS), FUN=mean, na.rm= T)$x; X_barS.xx
# ---------------------------------------------------------------- Perfect!!!!!!





############################## But Outcome parameter sampling is trickyyyyyyyyyy

#> This is where [new outcome parameters] go. 
betaparamS.xx = matrix(data=NA, nrow=JS, ncol=3) 
sig2param.xx = numeric(JS) # for cluster-wise
sig2paramfull.xx = 0       # for population

#] ------ Outcome Y ------ ::: [betaparamS], [sig2param]....Analytical Posterior is inaccessible, so try MH-algorithm
#------------------------------------------------------------------------------------------------------------------------

#####> First, sample hyperparameter from hyperPRIOR.
# Set in place the old hyperparameters..

### (A) beta0S, SIG_b0S For [betaparamS]
# 1) beta0S ~ MVN(m0S, 1/deltaS*SIG_b0S)  
# 2) SIG_b0S ~ ~ IW( qu0S, LAMBS )

###> initial hyperparameter values from hyperprior +++++++++++++++++++++++(no need to care cluster)
SIG_b0_prevS.xx <- riwish(v=qu0S.xx, S=LAMBS.xx)                                         #~~~~~~~~~~~~~~~(CHECK)
beta0_prevS.xx <- rmvn(n=1, mu=m0S.xx, sigma=1/deltaS.xx*SIG_b0_prevS.xx)                   #~~~~~~~~~~~~~~~(CHECK)


### (B) u0S, v0S For [sig2param]
# 1) u0S ~ Fink(rho_u1S, rho_u2S)
# 2) v0S ~ Ga(rho_v1S, rho_v2S)

# since.... Fink function for ....... "u0" sampling 
fink = function(u0, A, B) {
  A^(u0/2-1)/gamma(u0/2)^B
}
u0_proposed_shapeS = 3.75 # 3
u0_proposed_rateS = 2 # 2  based on Ga(3.75,2)...intuition..? for mean 3/2 ?

###> initial hyperparameter values from hyperprior +++++++++++++++++++++++(no need to care cluster)
u0_prevS.xx = u0_proposed_shapeS/u0_proposed_rateS     # instead of sampling from Fink function???
v0_prevS.xx = rgamma(n=1, shape=rho_v1S.xx, rate=rho_v2S.xx)

#jpeg(file="plot.test.jpeg", width=1000, height=800)
#-------------------------------------------------------------------------------
# See? Fink function is close to Ga(3.75,2) ?
curve(fink(x,rho_u1S.xx,rho_u2S.xx)/integrate(function(x) fink(x,rho_u1S.xx,rho_u2S.xx),0,Inf)$value, from = 0.01, to = 10, ylim = c(0,0.5))
curve(dgamma(x, shape = u0_proposed_shapeS, rate = u0_proposed_rateS), add=TRUE, col="blue")
#-------------------------------------------------------------------------------
#dev.off()
#ggsave("my_plot.jpeg", plot = p, width=10, height=8, units = "cm")



##### Now....ready for Single SAMPLE FROM the PRIOR
###> initial [main parameter] values from PRIOR ++++++++++++++++++++++++++++++(cluster-wise, but...)
beta_prev_fullS.xx = rmvn(n=1, mu=beta0_prevS.xx, sigma=SIG_b0_prevS.xx) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
beta_prevS.xx = matrix( rep(beta_prev_fullS.xx, JS), nrow=JS, byrow= TRUE ) # cluster-wise!!!!

###> initial [main parameter] values from PRIOR ++++++++++++++++++++++++++++++(cluster-wise, but...)
sig2_prev_full.xx = rinvgamma(n=1, shape=u0_prevS.xx, scale=v0_prevS.xx)   # population 
sig2_prev.xx = rep(sig2_prev_full.xx, JS)                               # cluster-wise!!!!


#####> Second, Prepare for Metropolis Hastings for [betaparamS], [sig2param] ~ LN(.)*MVN(_)*InvGa(_)
###########################################################################
#> It's time to use......... 
# - "beta0_prevS", "SIG_b0_prevS", "u0_prevS", "v0_prevS"         for MVN(_), InvGa(_) :based on population
# - "beta_prevS", "sig2_prev", "beta_prev_fullS", "sig2_prev_full"  for LN(.)         :based on cluster-wise

#> update hyperparameter first, then ..... update a single sampling: [main parameter]
# 1) beta0paramS, SIG_b0paramS for: [betaparamS]
# 2) u0paramS, v0paramS for: [sig2param]
#---------------------------------------------------------------------------------------------------------------
############################################ from posterior of hyperprior ###########################################

##] serious Hyperparameter Toward [betaparam]
SIG_b0paramS.xx = riwish(v = qu0S.xx+2,
                         S = t(beta0_prevS.xx - beta_prev_fullS.xx) %*% (beta0_prevS.xx - beta_prev_fullS.xx) +
                           deltaS.xx*t(beta0_prevS.xx - m0S.xx) %*% (beta0_prevS.xx - m0S.xx) + LAMBS.xx) #~~~~~~~~~~~~~~~(CHECK)

beta0paramS.xx = rmvn(n = 1, 
                      mu = deltaS.xx/(deltaS.xx+1)*m0S.xx + 1/(deltaS.xx+1)*beta_prev_fullS.xx,
                      sigma = 1/(deltaS.xx+1)*SIG_b0paramS.xx) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)



##] serious Hyperparameter Toward [sig2param]
u0_proposalS.xx = rgamma(n=1, shape=u0_proposed_shapeS, rate=u0_proposed_rateS)

ratio_u0S.xx = min(fink(u0=u0_proposalS.xx, A=sig2_prev_full.xx*rho_u1S.xx/(v0_prevS.xx/2), B=rho_u2S.xx+1)/
                     fink(u0=u0_prevS.xx, A=sig2_prev_full.xx*rho_u1S.xx/(v0_prevS.xx/2), B=rho_u2S.xx+1)*
                     dgamma(x = u0_prevS.xx, shape = u0_proposed_shapeS, rate = u0_proposed_rateS)/
                     dgamma(x = u0_proposalS.xx, shape = u0_proposed_shapeS, rate = u0_proposed_rateS), 1 )

if(runif(n = 1, min = 0, max = 1) < ratio_u0S.xx) {
  u0paramS.xx = u0_proposalS.xx
} else {
  u0paramS.xx = u0_prevS.xx
}


v0paramS.xx = rgamma(n=1, shape=rho_v1S.xx+u0paramS.xx/2, rate=rho_v2S.xx+0.5/sig2_prev_full.xx)

SIG_b0paramS.xx; beta0paramS.xx; u0paramS.xx; v0paramS.xx  
# for sharing...#### Now we have single samples to start MH engine.

#####> From here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
############################################################################### "Cluster-wise"!!!!!!!
##] Run your MH to obtain single posterior samples of serious main parameters: [betaparamS], [sig2param] 

for(j in 1:JS) {
  
  ##> Sample proposals from PRIORS
  beta_propS.xx = rmvn(n=1, mu=beta0paramS.xx, sigma=SIG_b0paramS.xx) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
  sig2_prop.xx = rinvgamma(n=1, shape=u0paramS.xx, scale=v0paramS.xx) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
  
  ##> subsetting by cluster
  YSj = Y[cl_membershipS==j]
  matXSj = matXSstar[cl_membershipS==j, ]
  
  
  
  ##> A. MH algorithm for: [betaparamS]
  # Note "beta_propS" is from PRIOR proposal (new) while "beta_prevS[j,]" is from the single POSTERIOR sample (old). 
  numerator=0
  denominator=0
  
  for (i in 1:length(YSj)){
    numerator = numerator + dlognorm_log( y=YSj[i], x=matXSj[i,], beta=t(as.matrix(beta_propS.xx)), sig2=sig2_prev.xx[j] ) + 0
    #log( dmvn(beta_propS, mu=beta0paramS, sigma=SIG_b0paramS) ) +
    #log( dmvn(beta_prevS[j,], mu=beta0paramS, sigma=SIG_b0paramS) ) # cancelled out?? #~xxxxxxxxxxxxxxxxx Q.DAVID
    
    denominator = denominator + dlognorm_log( y=YSj[i], x=matXSj[i,],
                                              beta=t(as.matrix(beta_prevS.xx[j,])), sig2=sig2_prev.xx[j] ) + 0
    #log( dmvn(beta_prevS[j,], mu=beta0paramS, sigma=SIG_b0paramS) ) +
    #log( dmvn(beta_propS, mu=beta0paramS, sigma=SIG_b0paramS) )    # cancelled out?? #~xxxxxxxxxxxxxxxxxx Q.DAVID
  }
  ratio_betaS.xx = min(exp(numerator-denominator), 1)
  
  # Example of adaptive variance adjustment
  #if (ratio_betaS < 0.2) {
  #    SIG_b0paramS <- SIG_b0paramS * 0.9
  #} else if (ratio_betaS > 0.4) {
  #    SIG_b0paramS <- SIG_b0paramS * 1.1
  
  # what if ... the ratio is likely to be so small???? no-update??????????
  # --- Adjust the Proposal Distribution: 1) Decrease the variance of the proposal distribution
  #                         2) Use an adaptive proposal distribution to adjust its variance based on the acceptance rate.
  #                         3) it can be varinf problem.....   
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio_betaS.xx) {
    betaparamS.xx[j, ] = beta_propS.xx
  } 
  else {
    betaparamS.xx[j, ] = beta_prevS.xx[j, ]
  }
  
  
  ##> B. MH algorithm for: [sig2param]
  # Note "sig2_prop" is from PRIOR proposal (new) while "sig2_prev[j,]" is from the single POSTERIOR sample (old).
  numerator=0
  denominator=0
  for (i in 1:length(YSj)){
    numerator=numerator + dlognorm_log( y=YSj[i], x=matXSj[i,], beta=t(as.matrix(betaparamS.xx[j,])), sig2=sig2_prop.xx )+0
    #log(dinvgamma(sig2_prop, shape=u0paramS, scale=v0paramS)) + 
    #log(dinvgamma(sig2_prev[j], shape=u0paramS, scale=v0paramS)) # cancelled out? #~xxxxxxxxxxxxxxxxx Q.DAVID
    
    denominator=denominator + dlognorm_log( y=YSj[i], x=matXSj[i,], beta=t(as.matrix(betaparamS.xx[j,])), sig2=sig2_prev.xx[j] )+0 
    #log(dinvgamma(sig2_prev[j], shape=u0paramS, scale=v0paramS)) + 
    #log(dinvgamma(sig2_prop, shape=u0paramS, scale=v0paramS))    # cancelled out? #~xxxxxxxxxxxxxxxxx Q.DAVID
  }
  ratio_sig2.xx = min(exp(numerator-denominator), 1)
  
  # what if ... the ratio is likely to be so small???? no-update????????????
  # --- Adjust the Proposal Distribution: 1) Decrease the variance of the proposal distribution
  #                         2) Use an adaptive proposal distribution to adjust its variance based on the acceptance rate.
  #                         3) it can be varinf problem.....   
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio_sig2.xx) {
    sig2param.xx[j] = sig2_prop.xx
  } 
  else {
    sig2param.xx[j] = sig2_prev.xx[j]
  }
}
betaparamS.xx; sig2param.xx



# NOW..................

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# We NEED a Cluster-wise set of....... Many Smaples ............................
# - outcome parameters: [betaparamS], [sig2param]
# - covariate parameters: [piparamS], [X_barS], [lambda2paramS]
# 
#..............................................................................#
#..............................................................................#
#..............................................................................#
#................................ now be ready ................................# 
#..............................................................................#
#..............................................................................#
#..............................................................................#
################################################################################
### Step03> Gibbs Sampler --------- cl_membership and param update ---with J= ?
################################################################################
set.seed(1)
total_iter = 1600
# r_convergence = 1000   # first run Gibbs sampler to determine r value for convergence
r_convergence = 100

# loglikelihood = numeric(total_iter)                  # for monitor convergence
loglikelihood.xx = matrix(0, nrow = n.train, ncol = total_iter)
# [note] why matrix? each col gives each iteration result

###] Prepare a "pool" to store the ultimate hyperparameters, parameters ++++++++++++++++++++++++++++++++++++++++
#> w/o cluster (hyperparameter)
list_beta0paramS.xx = list() # For beta0 hyperparameter
list_SIG_b0paramS.xx = list() # For SIG_b0 hyperparameter
list_u0paramS.xx = list() # for u0 hyperparameter
list_v0paramS.xx = list() # for v0 hyperparameter

#> w/o cluster (outcome)
list_betaparamfullS.xx = list() # for N - full data
list_sig2paramfull.xx = list()  # for N - full data

#> w/o cluster (covariate)
# N/A                                                                       #~xxxxxxxxxxxxxxxxxxxxx Q.DAVID
# N/A                                                                       #~xxxxxxxxxxxxxxxxxxxxx Q.DAVID
#-----------------------------------------------------------



#> cluster-wise (outcome)
list_betaparamS.xx = list()    #for N
list_sig2param.xx = list()    #for N

#> cluster-wise (covariate)
list_piparamS.xx = list()   #for ZF
list_lambda2S.xx = list()    #for XF
#list_X_barS = list()     #for XF...always same?



###] Store the current versions of single "cluster-wise" parameters as the "previous" ++++++++++++++++++++++++++++
#> hyperparameter
beta0_prevS.xx = beta0paramS.xx
SIG_b0_prevS.xx = SIG_b0paramS.xx
u0_prevS.xx = u0paramS.xx
v0_prevS.xx = v0paramS.xx

#> outcome
beta_prevS.xx = betaparamS.xx
sig2_prev.xx = sig2param.xx

#> covariate
#piparamS #always same?
lambda2param_prevS.xx = lambda2paramS.xx # do we need???? 
#X_barS # always same?





###] Miscellaneous ???

max_iter = 10   # let's say.. investigation( sampling rejection ) limit...?
counts_matS.xx <- matrix(0, nrow = JS*total_iter, ncol = 3) # count..How many rejected???

nS_j = table(cl_membershipS)

beta_p_acceptedS.xx <- 0
sig2_p_accepted.xx <- 0

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

nS_j # just temporary...for cluster-wise sampling...

########################################################## START ###############################################################

for (r in 1:total_iter) {
  
  ###[1] Updating Posterior -----------------------------------------------------------------------------------------------
  
  #>>> First, sample ***HYPERPARAMETERS*** from hyperPosterior based on full data {w/o clusters}
  
  ## hyperparam for [betaparamS]
  SIG_b0paramS.xx = riwish(v=qu0S.xx+2, 
                           S=t(beta0_prevS.xx-beta_prev_fullS.xx) %*% (beta0_prevS.xx-beta_prev_fullS.xx) +
                             deltaS.xx*t(beta0_prevS.xx - m0S.xx) %*% (beta0_prevS.xx - m0S.xx) + LAMBS.xx) #~~~~~~~~~~~~~~~~~~(CHECK)
  
  beta0paramS.xx = rmvn(n=1, mu=deltaS.xx/(deltaS.xx+1)*m0S.xx + 1/(deltaS.xx+1)*beta_prev_fullS.xx, 
                        sigma=1/(deltaS.xx+1)*SIG_b0paramS.xx) 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
  
  ## hyperparam for [sig2param]  
  #> - Sample u0 from MH using Gamma(3.75, 2) as proposal 
  u0_propS.xx = rgamma(n=1, shape=u0_proposed_shapeS, rate=u0_proposed_rateS) # Fink approx
  ratio_u0S.xx = min(fink(u0_propS.xx, A = sig2_prev_full.xx*rho_u1S.xx/(v0_prevS.xx/2), B = rho_u2S.xx+1)/
                       fink(u0_prevS.xx, A = sig2_prev_full.xx*rho_u1S.xx/(v0_prevS.xx/2), B = rho_u2S.xx+1)*
                       dgamma(x = u0_prevS.xx, shape = u0_proposed_shapeS, rate = u0_proposed_rateS)/
                       dgamma(x = u0_propS.xx, shape = u0_proposed_shapeS, rate = u0_proposed_rateS), 1)
  if(runif(n = 1, min = 0, max = 1) < ratio_u0S.xx) {
    u0paramS.xx = u0_propS.xx
  } else {
    u0paramS.xx = u0_prevS.xx
  }
  
  #> - Sample v0 from Gamma distribution 
  v0paramS.xx = rgamma(n=1, shape=rho_v1S.xx + u0paramS.xx/2, rate = rho_v2S.xx+0.5/sig2_prev_full.xx)
  
  
  
  
  
  #------------------------------------------------------------------------------------------------------------------
  #------------------------------------------------------------------------------------------------------------------
  
  #>>> Next I, Sample ****[Main parameters]**** 
  
  #> - Proposal Sampling using proposals OF PRIORS {w/o clusters} 
  beta_prop_fullS.xx = rmvn(n=1, mu=beta0paramS.xx, sigma=SIG_b0paramS.xx*varinf) #~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
  sig2_prop_full.xx = rinvgamma(n=1, shape=u0paramS.xx, scale=v0paramS.xx) 
  
  
  #> - Sifting a single [betaparam] based on FULL data 
  numerator=0
  denominator=0
  # plugging prior sample first...
  for (i in 1:n.train){
    numerator = numerator + dlognorm_log( y=Y[i], x=matXSstar[i,], 
                                          beta=t(as.matrix(beta_prop_fullS.xx)), sig2=sig2_prev_full.xx )
    denominator = denominator + dlognorm_log( y=Y[i], x=matXSstar[i,],
                                              beta=t(as.matrix(beta_prev_fullS.xx)), sig2=sig2_prev_full.xx )
  }
  ratio_betaS.xx = min(exp(numerator-denominator), 1)
  
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio_betaS.xx) {
    betaparam_fullS.xx = beta_prop_fullS.xx
  } 
  else {
    betaparam_fullS.xx = beta_prev_fullS.xx
  }
  
  
  #> - Sifting a single [sig2param] based on FULL data 
  numerator=0
  denominator=0
  # plugging prior sample first...
  for (i in 1:n.train){
    numerator = numerator + dlognorm_log( y=Y[i], x=matXSstar[i,], 
                                          beta=t(as.matrix(betaparam_fullS.xx)), sig2=sig2_prop_full.xx )
    denominator = denominator + dlognorm_log( y=Y[i], x=matXSstar[i,],
                                              beta=t(as.matrix(betaparam_fullS.xx)), sig2=sig2_prev_full.xx )    
  }
  ratio_sig2.xx = min(exp(numerator-denominator), 1)
  
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio_sig2.xx) {
    sig2param_full.xx = sig2_prop_full.xx
  } 
  else {
    sig2param_full.xx = sig2_prev_full.xx
  }
  
  #--------------------------------------------------------------------------------------------------------#
  #--------------------------------------------------------------------------------------------------------#
  #-----------------------------------Now it's time for cluster-wise---------------------------------------#
  #--------------------------------------------------------------------------------------------------------#
  #--------------------------------------------------------------------------------------------------------#
  for(j in 1:JS) {
    
    #>>> First.., Sample **[covariate parameters_j]**
    XS_j = matXSstar[cl_membershipS==j] #covariate
    ZS_j = ZS[cl_membershipS==j] #covariate
    nS_j = length(XS_j)
    
    #> - proposal Sampling [lambda2paramF_j] using proposals OF POSTERIOR {with clusters} 
    lambda2paramS.xx[j] = rinvgamma( n=1, 
                                     shape=(c0S.xx+nS_j)/2, 
                                     scale=(d0S.xx+sum( (XS_j - mean(XS_j))^2 ))/2 )  
    
    #> - proposal Sampling [piparamF_j] using proposals OF POSTERIOR {with clusters} 
    piparamS.xx[j] = rbeta( n=1, shape1=g0S.xx+sum(ZS_j), shape2=h0S.xx+nS_j-sum(ZS_j) )
    
    
    #------------------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------------------
    #>>> Next,... Sample ****[Main parameters_j]****
    
    #>>> subsetting by cluster
    Y_j = Y[cl_membershipS==j]           #outcome
    matXS_j = matXSstar[cl_membershipS==j, ] #covariate
    
    
    #>>> - Sifting [betaparam]_j
    ratio_betaS.xx = 0
    count_betaS.xx = 0
    U = 1
    print("betaS sampling")
    while((U > ratio_betaS.xx) & (count_betaS.xx < max_iter)) {
      # if new samples keep being rejected, ?????
      beta_propS.xx = rmvn(n=1, mu=beta0paramS.xx, sigma=SIG_b0paramS.xx*varinf) #~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
      numerator=0
      denominator=0
      for (i in 1:length(Y_j)){
        numerator = numerator + dlognorm_log( y=Y_j[i], x=matXS_j[i,], 
                                              beta=t(as.matrix(beta_propS.xx)), sig2=sig2_prev.xx[j] )
        
        denominator = denominator + dlognorm_log( y=Y_j[i], x=matXS_j[i,],
                                                  beta=t(as.matrix(beta_prevS.xx[j,])), sig2=sig2_prev.xx[j] )
      }
      ratio_betaS.xx = min(exp(numerator-denominator), 1)
      U = runif(n = 1, min = 0, max = 1)
      count_betaS.xx = count_betaS.xx+1
    }
    print(count_betaS.xx)
    
    counts_matS.xx[JS*(r-1)+j, 1] <- count_betaS.xx
    
    if(count_betaS.xx < max_iter) {
      betaparamS.xx[j, ] = beta_propS.xx
      beta_p_acceptedS.xx <- beta_p_acceptedS.xx+1
    } else {
      betaparamS.xx[j,] = beta_prevS.xx[j,]
    }
    
    
    #>>> - Sifting [sig2param]_j
    ratio_sig2.xx = 0
    count_sig2.xx = 0
    U = 1
    print("sig2 sampling")
    while((U > ratio_sig2.xx) & (count_sig2.xx < max_iter)) {
      sig2_prop.xx = rinvgamma(n=1, shape=u0_newS.xx, scale=v0_newS.xx) 
      numerator=0
      denominator=0
      for (i in 1:length(Y_j)){
        numerator = numerator + 
          dlognorm_log( y=Y_j[i], x=matXS_j[i,], beta=t(as.matrix(betaparamS.xx[j,])), sig2=sig2_prop.xx ) + 
          log(dinvgamma(sig2_prop.xx, shape=u0paramS.xx, scale=v0paramS.xx)) + 
          log(dinvgamma(sig2_prev.xx[j], shape=u0_newS.xx, scale=v0_newS.xx))
        
        denominator = denominator + 
          dlognorm_log( y=Y_j[i], x=matXS_j[i,], beta=t(as.matrix(betaparamS.xx[j,])), sig2=sig2_prev.xx[j] ) + 
          log(dinvgamma(sig2_prev.xx[j], shape=u0paramS.xx, scale=v0paramS.xx)) + 
          log(dinvgamma(sig2_prop.xx, shape=u0_newS.xx, scale=v0_newS.xx))
      }
      ratio_sig2.xx = min(exp(numerator-denominator), 1)
      U = runif(n = 1, min = 0, max = 1)
      count_sig2.xx = count_sig2.xx+1
    }
    print(count_sig2.xx)
    counts_matS.xx[JS*(r-1)+j,2] <- count_sig2.xx
    if(U < ratio_sig2.xx) {
      sig2param.xx[j] = sig2_prop.xx
      sig2_p_accepted.xx <- sig2_p_accepted.xx+1
    }
    else {
      sig2param.xx[j] = sig2_prev.xx[j]
    }
  }
  #################################################################################################################
  
  
  #### For NEXT!!!
  ###> Set previous set of parameter to the current for the next iteration
  ## - Manner! main param for next iteration
  lambda2param_prevS.xx = lambda2paramS.xx
  piparam_prevS.xx = piparamS.xx
  
  beta_prev_fullS.xx = betaparam_fullS.xx
  beta_prevS.xx = betaparamS.xx
  
  sig2_prev_full.xx = sig2param_full.xx 
  sig2_prev.xx = sig2param.xx
  
  ## - Manner! hyperparam for next iteration
  beta0_prevS.xx = beta0paramS.xx
  
  SIG_b0_prevS.xx = SIG_b0paramS.xx
  
  u0_prevS.xx = u0paramS.xx
  
  v0_prevS.xx = v0paramS.xx
  
  
  
  
  
  
  
  ###[2] Calculating Loglikelihood to monitor the chain convergence ---------------------------------------------
  #>>> First
  loglike_r = numeric(n.train)
  
  for(i in 1:n.train) {
    xS_i = matXSstar[i, ]      # first row..{1,x,z}
    jS = cl_membershipS[i] # membership ID
    
    loglike_r[i] = dlognorm_log( y=Y[i], x=xS_i, beta=t(as.matrix(betaparamS.xx[jS,])), sig2=sig2param.xx[jS] ) + 
      dnorm(x=xS_i[2], mean=X_barS.xx[jS], sd=sqrt(lambda2paramS.xx[jS]), log=TRUE) +
      dbinom(x=xS_i[3], size=1, prob=piparamS.xx[jS], log=TRUE) 
  }  
  loglikelihood.xx[,r] = loglike_r
  
  if(r > r_convergence) {
    
    # W/O CLUSTER
    list_beta0paramS.xx[[r-r_convergence]] = beta0paramS.xx
    list_SIG_b0paramS.xx[[r-r_convergence]] = SIG_b0paramS.xx
    list_u0paramS.xx[[r-r_convergence]] = u0paramS.xx # summary(unlist(list_u0paramS)) # for u0_new
    list_v0paramS.xx[[r-r_convergence]] = v0paramS.xx # summary(unlist(list_v0paramS)) # for v0_new
    
    list_betaparamfullS.xx[[r-r_convergence]] = betaparam_fullS.xx
    list_sig2paramfull.xx[[r-r_convergence]] = sig2param_full.xx
    
    
    # Cluster-wise
    list_piparamS.xx[[r-r_convergence]] = piparamS.xx     
    list_lambda2S.xx[[r-r_convergence]] = lambda2paramS.xx  
    
    list_betaparamS.xx[[r-r_convergence]] = betaparamS.xx  
    list_sig2param.xx[[r-r_convergence]] = sig2param.xx 
  }
  print(paste("r=",r))
  trace(what = "print", where = getNamespace("base"), exit = flush.console, print = FALSE)
}

################################################################################
####### End of Gibbs Sampler ###################################################
################################################################################

##> Check your sampling performance
beta_p_acceptedS.xx/(total_iter*JS)
sig2_p_accepted.xx/(total_iter*JS)

par(mfrow = c(1,1))
plot(colSums(loglikelihood.xx), type="l")

lppd_EX_S.Y.xx = mean( colSums(loglikelihood.xx) ); lppd_EX_S.Y.xx
#[1] -17599.1, -17731.03

# Proportion of acceptance for betaparamS, sig2param
prop.table(table(counts_matS.xx[,1]))
prop.table(table(counts_matS.xx[,2]))

#-------------------------------------------------------------------------------
###> Investigation for "sig2"
list_sig2param.xx[1]

sig2.vec_cl1.xx = sapply(list_sig2param.xx, function(x) x[1])
sig2.vec_cl2.xx = sapply(list_sig2param.xx, function(x) x[2])
sig2.vec_cl3.xx = sapply(list_sig2param.xx, function(x) x[3])
sig2.vec_cl4.xx = sapply(list_sig2param.xx, function(x) x[4])
sig2.vec_cl5.xx = sapply(list_sig2param.xx, function(x) x[5])
sig2.vec_cl6.xx = sapply(list_sig2param.xx, function(x) x[6])

sig2.xx.df <- data.frame(
  value = c(sig2.vec_cl1.xx, sig2.vec_cl2.xx, sig2.vec_cl3.xx, sig2.vec_cl5.xx, sig2.vec_cl5.xx, sig2.vec_cl6.xx),
  group = factor(rep(1:6, each = total_iter-r_convergence))
)
#option A (together)                    
#ggplot(sig2.df, aes(x = value, color = group)) +
#                     geom_density() +
#                     labs(title = "Density for Each sig2.vec", x = "sig2", y = "Density") + 
#                     theme_minimal()

#option B (separate) #----------------------------------------- R they InvGamma?
tapply(sig2.xx.df$value, sig2.xx.df$group, mean)
# 1        2        3        4        5        6 
# 3.526035 3.347704 3.610494 3.767755 3.767755 3.520474 
aggregate(value ~ group, data=sig2.xx.df, FUN=mean)

# Calculate the average sig2
sn <- length(list_sig2param.xx); sn
avg_sig2.vec1.xx <- sum(sig2.vec_cl1.xx)/sn; avg_sig2.vec1.xx #
avg_sig2.vec2.xx <- sum(sig2.vec_cl2.xx)/sn; avg_sig2.vec2.xx #
avg_sig2.vec3.xx <- sum(sig2.vec_cl3.xx)/sn; avg_sig2.vec3.xx #
avg_sig2.vec4.xx <- sum(sig2.vec_cl4.xx)/sn; avg_sig2.vec4.xx #
avg_sig2.vec5.xx <- sum(sig2.vec_cl5.xx)/sn; avg_sig2.vec5.xx #
avg_sig2.vec6.xx <- sum(sig2.vec_cl6.xx)/sn; avg_sig2.vec6.xx #



#-------------------------------------------------------------------------------
###> Investigation for "beta0", "beta1", "beta2"
list_betaparamS.xx[1]
# Sum all the matrices
sum_matrixS.xx <- Reduce("+", list_betaparamS.xx)
# Number of matrices
nb <- length(list_betaparamS.xx)
# Calculate the average matrix of betaS
avg_beta.matS.xx <- sum_matrixS.xx / nb
avg_beta.matS.xx

#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++ PREDICTION +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# Expected Value of LN: exp(sum(x*beta) # x^T*B
#                       = exp( mu + 0.5*sig2 )     
#                              mu = sum(x*beta) - 0.5*sig2

###> with Training data
expvalS.train.xx = matrix(0, nrow = n.train, ncol = total_iter-r_convergence)
mu.train.xx = matrix(0, nrow = n.train, ncol = total_iter-r_convergence)

for(r in 1:(total_iter-r_convergence)) {
  piparamS.xx = list_piparamS.xx[[r]]
  lambda2paramS.xx = list_lambda2S.xx[[r]]
  betaparamS.xx = list_betaparamS.xx[[r]]
  sig2param.xx = list_sig2param.xx[[r]] 
  
  for(i in 1:n.train) {
    j = cl_membershipS[i]
    #cleanxS = rnorm(n=1, mean=X_barS[j], sd=sqrt(lambda2paramS[j]))
    dirtyXS = XSstar[i]
    expvalS.train.xx[i,r] <- exp( betaparamS.xx[j,1] + betaparamS.xx[j,2]*dirtyXS + 
                                    betaparamS.xx[j,3]*ZS[i] - 0.5*sig2param.xx[j] ) #********** this is the key???
    mu.train.xx[i,r] <- betaparamS.xx[j,1] + betaparamS.xx[j,2]*dirtyXS + betaparamS.xx[j,3]*ZS[i]
  }
}

EX_S.train.dirty <- apply(X=expvalS.train.xx, MARGIN=1, FUN=mean); EX_S.train.dirty
Y
mu.vec.train.dirty <- apply(X=mu.train.xx, MARGIN=1, FUN=mean)

#options(scipen = 999)
format( aggregate(log(Y), by=list(cl_membershipS), FUN=summary), scientific=FALSE)
format( aggregate(log(EX_S.train.dirty), by=list(cl_membershipS), FUN=summary), scientific=FALSE)



#
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW with COREECTION NOW NOW NOW NOW NOW NOW NOW NOW NOW 
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW
# NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW NOW

###############################################################################################################################
######################################### Y ~ hierarchical LN with Gustafson ##################################################
###############################################################################################################################

################################################################################
### Step01> Define prior, hyperprior model and hyperparameter values
################################################################################

#] ------ Outcome ------ ::: for logY ~ N(X*betaparamS, sig2paramS)
#------------------------------------------------------------------------------------------------------------------------

# PRIOR: "betaparam" ~ MVN(beta0S, SIG_b0S)
# Hyperprior               beta0S ~ MVN(m0S, 1/deltaS*SIG_b0S)

#> beta0
#> Gaussian regression for initialize OUTCOME parameter, using 'm0', 'SIG_b0' with dirty data
#fit1S.st <- glm( log(Y) ~ XSstar + factor(ZS), family=gaussian(link = "identity")) # AIC: 2596.5
#summary(fit1S.st)

fit1S.st <- glm( Y ~ XSstar + factor(ZS), family=poisson(link="log")) # AIC: ?
summary(fit1S.st)
vcov(fit1S.st) # too small...

m0S.st = coef(fit1S.st)       # 1x3 initial reg_coeff vector (sort of Regression result as "mean")
SIG_b0S.st = vcov(fit1S.st)*varinf*10000   # 3x3 initial cov matrix of reg_coeff

p = length(m0S.st) # intercept,b1,b2
options(scipen = 999)
SIG_b0S.st


deltaS.st = 0.01     #1             # equal importance on m0 and beta ???? (0.1)
varinf = n.train/100 #n.train/10    # rule of thumb divide by 5    (n.train/1000)


#SIG_b0... = SIG_b0*varinf ???
SIG_b0invS.st = solve(a=SIG_b0S.st) # inverse of cov matrix of reg_coeff for later use (for posterior on reg_coeff)!!
SIG_b0invS.st

options(scipen = 999)
SIG_b0S.st*varinf

# PRIOR: "betaparamS" ~ MVN(beta0S, SIG_b0S)
# Hyperprior                        SIG_b0S ~ IW( qu0S, LAMBS )


#> SIG_b0S
#- First, SIG_b0S = VCOV(fit1S)       ???????????? Really? 
#- Next,  SIG_b0S ~ IW( qu0S, LAMBS )  ???????????? Really?


#> it requires -------
qu0S.st = p+2
# LAMB ?????????????????????????????????????????????
LAMBS.st = SIG_b0S.st#*varinf # really??????????????????????????????
LAMBS.st



# PRIOR: "sig2param" ~ IG( u0S, v0S )
# Hyperprior               u0S ~ Fink( rho_u1S, rho_u2S )
# Hyperprior                    v0S ~ Ga( rho_v1S, rho_v2S )

# Hyperparameter a0 for Fink
rho_u1S.st = 1/8 
rho_u2S.st = 3/2 

# Hyperparameter b0 - 2024-02-09
rho_v1S.st = 8 # 4
rho_v2S.st = 1 # 1
# this is an exponential distribution for the prior on v ?

# Based on [Method of Moment], the Proposal Distribution was obtained, and thus Parameters for U0 and V0
# u0_new = 2.561215 # from 2/22 
# v0_new = 2.761394 # from 2/22 

##### NEED to change --------------------------------------------------------------------------------------------------------
# E[invGa] = scale/(shape-1) = 50
u0_newS.st = 1.7591 #3 #2 #1.1 # from 3/7  # from # summary(unlist(list_u0paramS)) # for u0_new                            
v0_newS.st = 7.113  #3 #5 #8.5 # from 3/7  # from # summary(unlist(list_v0paramS)) # for v0_new

#-------------------------------------------------------------------------------
#u0_newF = ?
#v0_newF = ?
# How can we expand the sampling range??? for psi~Gamma(.) ?

plot( density(unlist(list_sig2param.st)) ) # sth we obtain...
curve( dinvgamma(x, shape=mean(unlist(list_u0paramS.st)), scale=mean(unlist(list_v0paramS.st))), add=TRUE, col="red") 
curve( dinvgamma(x, shape=u0_newS.st, scale=v0_newS.st), add=TRUE, col="blue")
#-------------------------------------------------------------------------------


#] ------ Covariate ------ ::: for ZS ~ Bin( 1, [piS] )  
# PRIOR                                         "piS" ~ Beta(g0S, h0S)
#------------------------------------------------------------------------------------------------------------------------
g0S.st = 0.5 
h0S.st = 0.5 


#] ------ Covariate ------ ::: for XSstar|ZS ~ N( K0+K1*ZS, [lambda2S] ) where 
# PRIOR                                      "Kappa" ~ MVN(kappa_v0, SIG_k0)
# PRIOR                                               "lambda2S" ~ IG(c0S, d0S)
#------------------------------------------------------------------------------------------------------------------------

#> KAPPA ????
# dirty covariate x is NOT tied to z under the NDB assumption, by this logic, regression of [x|z] and [xstar|z] should return 
# almost the same regression parameter-results...other than parameter-variance...
# It is for the sake of investigation or computation or mathematical convenience ?!?! ...since we need model (x,z)=(x|z)(z)
# It is for the connection between [clean] - [dirty] exposure model ?!?!
# We d not know VAR(Xstar). but by conditioning VAR(Xstar|Z), 
# which is "\hat{lambda}^2", and by comparing VAR(Xstar|Z) vs VAR(X|Z), we can perceive.... tau^2     

#> Kappa 
fit2S.st <- lm(XSstar ~ factor(ZS))
summary(fit2S.st) 

kappa_v0 <- coef(fit2S.st)
SIG_k0 <- vcov(fit2S.st)    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[CHECK]
SIG_inv_k0 <- solve(SIG_k0)

#> lambda2
c0S.st = 0.5          
d0S.st = 0.5

#> [note] 
# - we need KK1 matrix for all obvs w/o concerning cluster
# - KK1_j : cluster-wise ?
KK1 = cbind(1, ZS); KK1 
# for the analytical solution for Parameter-Free covariate model development later ??? (in DPM context)

#] ------ Measurement Model ------ ::: for XSstar|XStrue ~ N( XStrue, [tau2] ) 
# PRIOR:                                                              "tau2" ~ IG( c0S, [1-zeta_scale]*d0S )
#------------------------------------------------------------------------------------------------------------------------

# since..... tau2 = [1-zeta_scale]*lambda2 
# HEY HEY HEY ...We sample tau2 from the posterior of the "lambda model" in the future!!! so...we don't need prior for tau2. 
# zeta_scale = c(0.01,0.05,0.1,0.2,0.5)    # for 5 clusters
# zeta_scale = c(0.60,0.70,0.80,0.90,0.95) # for 5 clusters




# Large ########################
#zeta_scale = 0.4              #
# Medium                       #
#zeta_scale = 0.6              #
# Small                        #
#zeta_scale = 0.5               #
################################
zeta_scale = 0.7 


# Now, we need a posterior sample of main parameters
##################################################################################################
### Step02> initialize major parameters (single value), using analytically driven "POSTERIOR" pdf  
#                                                 Based on natural clustering
##################################################################################################
JS = length(table(cl_membershipS))

set.seed(1)

#] ------ Covariate ZS ------ ::: for ZS ~ Bin( 1, [piparamS_j] ) 
# prior                                            "piparamS_j" ~ Beta(g0S, h0S)
#------------------------------------------------------------------------------------------------------------------------
# Covariate parameter from Posterior
piparamS.st = numeric(JS)
for (j in 1:JS) {
  ZSj=ZS[cl_membershipS==j]
  nS_j=length(ZSj)
  piparamS.st[j] = rbeta( n=1, shape1=g0S.st+sum(ZSj), shape2=h0S.st+nS_j-sum(ZSj) ) 
}

piparamS.st #% pi parameter (posterior sample)



#] ------ Covariate X*|ZS ::: for XSstar|ZS ~ N( K0+K1*ZS, [lambda2paramS] ) where prior: "lambda2paramS" ~ IG(c0S, d0S)
#                                                                                            "Kappa" ~ MVN(kappa_v0, SIG_k0)
#------------------------------------------------------------------------------------------------------------------------
# Covariate parameter from Posterior
kappaparam.st = matrix(nrow = JS, ncol = length(kappa_v0))
lambda2paramS.st = numeric(JS)

for (j in 1:JS) {
  XSstarj=XSstar[cl_membershipS==j]
  ZSj = ZS[cl_membershipS==j]
  nS_j=length(XSstarj)
  lambda2paramS.st[j]=rinvgamma( n=1, shape=(nS_j+c0S.st)/2, scale=(sum((XSstarj-kappa_v0[1]-kappa_v0[2]*ZSj)^2) + d0S.st)/2 )
  KK1j = cbind(1, ZSj)
  KK2j = matrix(c(sum(XSstarj), sum(XSstarj*ZSj)), nrow =2, ncol = 1) # 2x1 vector
  
  SIG_inv_j = solve(SIG_inv_k0+(t(KK1j) %*% KK1j)) # solve( ..[matrix].. ) = [matrix]^(-1)
  kappaparam.st[j,] <- rmvn(n=1, mu=SIG_inv_j %*% (SIG_inv_k0 %*% kappa_v0 + KK2j), sigma = lambda2paramS.st[j]*SIG_inv_j)      
}
kappaparam.st       #% kappa parameter (posterior sample)
lambda2paramS.st #% lambda2 parameter (posterior sample)

#% [Check!] empirical kappa, lambda^2 ? using "aggregate(.)": investigation by splitting the data into subset
kappa_empirical=rbind(coef(lm(XSstar[cl_membershipS==1] ~ ZS[cl_membershipS==1])),
                      coef(lm(XSstar[cl_membershipS==2] ~ ZS[cl_membershipS==2])),
                      coef(lm(XSstar[cl_membershipS==3] ~ ZS[cl_membershipS==3])),
                      coef(lm(XSstar[cl_membershipS==4] ~ ZS[cl_membershipS==4])),
                      coef(lm(XSstar[cl_membershipS==5] ~ ZS[cl_membershipS==5])),
                      coef(lm(XSstar[cl_membershipS==6] ~ ZS[cl_membershipS==6]))); kappa_empirical

lambda2_empiricalS.st = aggregate(x=as.numeric(XSstar), by=list(cl_membershipS), FUN=var, na.rm=T)$x
lambda2_empiricalS.st # I dont' know


#] ------ Measurement XSstar|XS ------ ::: for XSstar|XS ~ N( XS, "[tau2]") where prior: "tau2" = (1-zeta_scale)*lambda2paramS
#------------------------------------------------------------------------------------------------------------------------
tau2param = (1-zeta_scale)*lambda2paramS.st


# But Outcome parameter sampling would be trickyyyyyyyyyy

#> This is where [new outcome parameters] go. 
betaparamS.st = matrix(data=NA, nrow=JS, ncol=3) 
sig2param.st = numeric(JS) # for cluster-wise
sig2paramfull.st = 0       # for population


#] ------ Outcome Y ------ ::: [betaparamS], [sig2param]....Analytical Posterior is inaccessible, so try MH-algorithm
#------------------------------------------------------------------------------------------------------------------------

#####> First, sample hyperparameter from hyperPRIOR.
# Set in place the old hyperparameters..

### (A) beta0S, SIG_b0S For [betaparamS]
# 1) beta0S ~ MVN(m0S, 1/deltaS*SIG_b0S)  
# 2) SIG_b0S ~ ~ IW( qu0S, LAMBS )

###> initial hyperparameter values from hyperprior +++++++++++++++++++++++(no need to care cluster)
SIG_b0_prevS.st <- riwish(v=qu0S.st, S=LAMBS.st)           
beta0_prevS.st <- rmvn(n=1, mu=m0S.st, sigma=1/deltaS.st*SIG_b0_prevS.st)


### (B) u0S, v0S For [sig2param]
# 1) u0S ~ Fink(rho_u1S, rho_u2S)
# 2) v0S ~ Ga(rho_v1S, rho_v2S)

# since.... Fink function for ....... "u0" sampling 
fink = function(u0, A, B) {
  A^(u0/2-1)/gamma(u0/2)^B
}
u0_proposed_shapeS.st = 3.75 # 3
u0_proposed_rateS.st = 2 # 1  based on Ga(3.75,2)...intuition..? for mean 3/2 ?

###> initial hyperparameter values from hyperprior +++++++++++++++++++++++(no need to care cluster)
u0_prevS.st = u0_proposed_shapeS.st/u0_proposed_rateS.st     # instead of sampling from Fink function???
v0_prevS.st = rgamma(n=1, shape=rho_v1S.st, rate=rho_v2S.st)


##### Now....ready for Single SAMPLE FROM the PRIOR

###> initial [main parameter] values from PRIOR ++++++++++++++++++++++++++++++(cluster-wise, but...)
beta_prev_fullS.st = rmvn(n=1, mu=beta0_prevS.st, sigma=SIG_b0_prevS.st*varinf) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
beta_prevS.st = matrix( rep(beta_prev_fullS.st, JS), nrow=JS, byrow= TRUE ) # cluster-wise!!!!

###> initial [main parameter] values from PRIOR ++++++++++++++++++++++++++++++(cluster-wise, but...)
sig2_prev_full.st = rinvgamma(n=1, shape=u0_prevS.st, scale=v0_prevS.st)   # population #?????? 
sig2_prev.st = rep(sig2_prev_full.st, JS)



#####> Second, Prepare for Metropolis Hastings for [betaparamS], [sig2param] ~ LN(.)*MVN(_)*InvGa(_)
###########################################################################
#> It's time to use......... 
# - "beta0_prevS", "SIG_b0_prevS", "u0_prevS", "v0_prevS"         for MVN(_), InvGa(_) :based on population
# - "beta_prevS", "sig2_prev", "beta_prev_fullS", "sig2_prev_full"  for LN(.)         :based on cluster-wise

#> update hyperparameter first, then ..... update a single sampling: [main parameter]
# 1) beta0paramS, SIG_b0paramS for: [betaparamS]
# 2) u0paramS, v0paramS for: [sig2param]
#---------------------------------------------------------------------------------------------------------------
############################################ from posterior of hyperprior ###########################################

##] serious Hyperparameter Toward [betaparam]
SIG_b0paramS.st = riwish(v = qu0S.st+2,
                         S = t(beta0_prevS.st - beta_prev_fullS.st) %*% (beta0_prevS.st - beta_prev_fullS.st) +
                           deltaS.st*t(beta0_prevS.st - m0S.st) %*% (beta0_prevS.st - m0S.st) + LAMBS.st) #~~~~~~~~~~(CHECK)

beta0paramS.st = rmvn(n = 1, 
                      mu = deltaS.st/(deltaS.st+1)*m0S.st + 1/(deltaS.st+1)*beta_prev_fullS.st,
                      sigma = 1/(deltaS.st+1)*SIG_b0paramS.st) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)



##] serious Hyperparameter Toward [psiparam]
u0_proposalS.st = rgamma(n=1, shape=u0_proposed_shapeS.st, rate=u0_proposed_rateS.st)

ratio_u0S.st = min(fink(u0=u0_proposalS.st, A=sig2_prev_full.st*rho_u1S.st/(v0_prevS.st/2), B=rho_u2S.st+1)/
                     fink(u0=u0_prevS.st, A=sig2_prev_full.st*rho_u1S.st/(v0_prevS.st/2), B=rho_u2S.st+1)*
                     dgamma(x = u0_prevS.st, shape = u0_proposed_shapeS.st, rate = u0_proposed_rateS.st)/
                     dgamma(x = u0_proposalS.st, shape = u0_proposed_shapeS.st, rate = u0_proposed_rateS.st), 1 )

if(runif(n = 1, min = 0, max = 1) < ratio_u0S.st) {
  u0paramS.st = u0_proposalS.st
} else {
  u0paramS.st = u0_prevS.st
}


v0paramS.st = rgamma(n = 1, shape = rho_v1S.st + u0paramS.st/2, rate = rho_v2S.st + 0.5/sig2_prev_full.st)

SIG_b0paramS.st; beta0paramS.st; 
u0paramS.st; v0paramS.st; 


###############################################################################
##] Run your MH to obtain single posterior samples of serious main parameters: [betaparamF], [psiparam] "Cluster-wise"!!!!!!!

for(j in 1:JS) {
  
  ##> Sample proposals from PRIORS
  beta_propS.st = rmvn(n=1, mu=beta0paramS.st, sigma=SIG_b0paramS.st*varinf) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
  sig2_prop.st = rinvgamma(n=1, shape=u0paramS.st, scale=v0paramS.st) # coz...??
  
  ##> subsetting by cluster
  YSj = Y[cl_membershipS==j]
  matXSjstar = matXSstar[cl_membershipS==j, ]
  
  
  
  ##> A. MH algorithm for: [betaparamS]
  # Note "beta_propS" is from PRIOR proposal (new) while "beta_prevS[j,]" is from the single POSTERIOR sample (old). 
  numerator=0
  denominator=0
  
  for (i in 1:length(YSj)){
    numerator = numerator + dlognorm_log( y=YSj[i], x=matXSjstar[i,], 
                                          beta=t(as.matrix(beta_propS.st)), sig2=sig2_prev.st[j] ) + 0
    #log( dmvn(beta_propS, mu=beta0paramS, sigma=SIG_b0paramS) ) +
    #log( dmvn(beta_prevS[j,], mu=beta0paramS, sigma=SIG_b0paramS) ) # cancelled out??
    
    denominator = denominator + dlognorm_log( y=YSj[i], x=matXSjstar[i,],
                                              beta=t(as.matrix(beta_prevS.st[j,])), sig2=sig2_prev.st[j] ) + 0
    #log( dmvn(beta_prevS[j,], mu=beta0paramS, sigma=SIG_b0paramS) ) +
    #log( dmvn(beta_propS, mu=beta0paramS, sigma=SIG_b0paramS) )    # cancelled out??
  }
  ratio_betaS.st = min(exp(numerator-denominator), 1)
  
  # Example of adaptive variance adjustment
  #if (ratio_betaS < 0.2) {
  #    SIG_b0paramS <- SIG_b0paramS * 0.9
  #} else if (ratio_betaS > 0.4) {
  #    SIG_b0paramS <- SIG_b0paramS * 1.1
  
  # what if ... the ratio is likely to be so small???? no-update??????????
  # --- Adjust the Proposal Distribution: 1) Decrease the variance of the proposal distribution
  #                         2) Use an adaptive proposal distribution to adjust its variance based on the acceptance rate.
  #                         3) it can be varinf problem.....   
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio_betaS.st) {
    betaparamS.st[j, ] = beta_propS.st
  } 
  else {
    betaparamS.st[j, ] = beta_prevS.st[j, ]
  }
  
  
  ##> B. MH algorithm for: [sig2param]
  # Note "sig2_prop" is from PRIOR proposal (new) while "sig2_prev[j,]" is from the single POSTERIOR sample (old).
  numerator=0
  denominator=0
  for (i in 1:length(YSj)){
    numerator=numerator + dlognorm_log( y=YSj[i], x=matXSjstar[i,], 
                                        beta=t(as.matrix(betaparamS.st[j,])), sig2=sig2_prop.st )+0
    #log(dinvgamma(sig2_prop, shape=u0paramS, scale=v0paramS)) + 
    #log(dinvgamma(sig2_prev[j], shape=u0paramS, scale=v0paramS)) # cancelled out?
    
    denominator=denominator + dlognorm_log(y=YSj[i], x=matXSjstar[i,], 
                                           beta=t(as.matrix(betaparamS.st[j,])), sig2=sig2_prev.st[j])+0 
    #log(dinvgamma(sig2_prev[j], shape=u0paramS, scale=v0paramS)) + 
    #log(dinvgamma(sig2_prop, shape=u0paramS, scale=v0paramS))    # cancelled out?
  }
  ratio_sig2.st = min(exp(numerator-denominator), 1)
  
  # what if ... the ratio is likely to be so small???? no-update????????????
  # --- Adjust the Proposal Distribution: 1) Decrease the variance of the proposal distribution
  #                         2) Use an adaptive proposal distribution to adjust its variance based on the acceptance rate.
  #                         3) it can be varinf problem.....   
  U = runif(n = 1, min = 0, max = 1)
  if(U < ratio_sig2.st) {
    sig2param.st[j] = sig2_prop.st
  } 
  else {
    sig2param.st[j] = sig2_prev.st[j]
  }
}

# Check if properly updated....

betaparamS.st   #(J,p) : membership(J), b0,b1,b2(p)
beta_prevS.st #(J,p)

sig2param.st   #(J)
sig2_prev.st #(J)
















# NOW.................................................................. Be ready
# NOW.................................................................. Be ready
# NOW.................................................................. Be ready
# NOW.................................................................. Be ready
# NOW.................................................................. Be ready
# NOW.................................................................. Be ready
# NOW.................................................................. Be ready
# NOW.................................................................. Be ready
# NOW.................................................................. Be ready
# NOW.................................................................. Be ready
# NOW.................................................................. Be ready
# NOW.................................................................. Be ready
# NOW.................................................................. Be ready
# NOW.................................................................. Be ready

# We NEED a Cluster-wise set of....... Many Many Many Many Many Many Many Many Smaples 
# - outcome parameters: [betaparamS], [sig2param]
# - covariate parameters: [piparamS], [X_barS], [lambda2paramS]
# 
#..............................................................................#
#..............................................................................#
#..............................................................................#
#................................... now ......................................# 
#..............................................................................#
#..............................................................................#
#..............................................................................#
################################################################################
### Step03> Gibbs Sampler --------- cl_membership and param update ---with J= ?
################################################################################
set.seed(1)
total_iter = 3000
# r_convergence = 1000   # first run Gibbs sampler to determine r value for convergence
r_convergence = 1000

# loglikelihood = numeric(total_iter)                  # for monitor convergence
loglikelihood.st = matrix(0, nrow = n.train, ncol = total_iter)
# [note] why matrix? each col gives each iteration result




#########################] Prepare a "pool"(list) to store the ultimate hyperparameters, parameters [##########################

#> w/o cluster ++++++++++++++++++++++++++++++++ (hyperparameter)
list_beta0paramS.st = list() # For beta0 hyperparameter
list_SIG_b0paramS.st = list() # For SIG_b0 hyperparameter
list_u0paramS.st = list() # for u0 hyperparameter
list_v0paramS.st = list() # for v0 hyperparameter

#> w/o cluster +++++++++++++++++++++++++++++++++++++++ (outcome)
list_betaparamfullS.st = list() # for Y - full data
list_sig2paramfull.st = list()  # for Y - full data

#> w/o cluster ++++++++++++++++++++++++++++++++++++++ (covariate)
# N/A
# N/A
# N/A

#-----------------------------------------------------------
#-----------------------------------------------------------
#-----------------------------------------------------------
#-----------------------------------------------------------
#-----------------------------------------------------------

####> with cluster-wise +++++++++++++++++++++++++++++++ (outcome) 
list_betaparamS.st = list()    #for Y
list_sig2param.st = list()    #for Y

# -------------------- after correction ---------------------- #
list_betaparamS.clean = list()    
list_sig2param.clean = list()    

####> with cluster-wise +++++++++++++++++++++++++++++ (covariate) 
list_piparamS.st = list()       #for ZS
list_lambda2paramS.st = list()  #for XS
list_kappaparam.st = list()             #for XS

list_tau2param = list()         #for measurement model

# -------------------- after correction ---------------------- #
list_piparamS.clean = list()        #for Z
list_lambda2paramS.clean = list()
list_kappaparam.clean = list()     #for X|Z

###########################################################################################################################


#--------------------------------------------------------------------- ALWAYS KEEP and NO-update !!!! --------------------------
### What we have so far................................................................................ for {initialization}
###] Prepare for the current versions of single "cluster-wise" parameters as the "previous" +++++++++++++++++++++++++++++++++

#> hyperparameter
beta0_prevS.st = beta0paramS.st
SIG_b0_prevS.st = SIG_b0paramS.st
u0_prevS.st = u0paramS.st
v0_prevS.st = v0paramS.st

#> outcome
beta_prevS.st = betaparamS.st
sig2_prev.st = sig2param.st
#beta_prev_fullS.st = rmvn(n=1, mu=beta0_prevS.st, sigma=SIG_b0_prevS.st) # prebviously given!!!!
#sig2_prev_full.st = rgamma(n=1, shape=u0_prevS.st, rate=v0_prevS.st)     # prebviously given!!!!


#> covariate ---------------------------------------------------------
piparam_prevS.st = piparamS.st           # we do not need. unaffected.
piparam_prevS.clean = piparamS.st        # we do not need. unaffected.

lambda2param_prevS.st = lambda2paramS.st # do we need???? 
lambda2param_prevS.clean = lambda2paramS.st #????

kappaparam_prev.st = kappaparam.st
kappaparam_prev.clean = kappaparam.st

tau2_prev = tau2param  # tau2param:::(1-zeta_scale)*lambda2paramS.st
#---------------------------------------------------------------------
#--------------------------------------------------------------------- ALWAYS KEEP and NO-update !!!! --------------------------


###] Miscellaneous ???

max_iter = 10   # let's say.. investigation( sampling rejection ) limit...?
counts_matS <- matrix(0, nrow = JS*total_iter, ncol = 3) # count..How many rejected???

nS_j.st = table(cl_membershipS)

beta_p_acceptedS.st <- 0
sig2_p_accepted.st <- 0

beta_p_acceptedS.clean <- 0
sig2_p_accepted.clean <- 0

diff_lambda2S <- matrix(0, nrow = n.train, ncol = total_iter)
diff_tau2 <- matrix(0, nrow = n.train, ncol = total_iter)

nS_j.st # just temporary...for cluster-wise sampling...


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
########################################################## START ###############################################################

for (r in 1:total_iter) {
  
  kappaparamclean = matrix(0, nrow = JS, ncol = ncol(kappaparam.st))     #for X|Z vs Xstar|Z
  lambda2paramcleanS = numeric(JS)                                        #for X|Z vs Xstar|Z
  
  betaparamcleanS =  matrix(0, nrow = JS, ncol = ncol(betaparamS.st))     #for Y
  sig2paramclean =numeric(JS)                                            #for Y
  
  xSclean = numeric(n.train)                                             #for X
  
  
  for(j in 1:JS) {
    
    XSstar_j=XSstar[cl_membershipS==j]
    ZSj = ZS[cl_membershipS==j]
    
    difference1 = -1
    difference2 = -1
    difference3 = -1
    
    count_iterations = 0
    
    ###> Don't stop until all the three positive values are collected.  
    while( difference1 < 0 | difference2 < 0 ) {
      
      ###[1] Updating Posterior -------------------------------------------------------------------------------
      
      #>>> First, sample ***HYPERPARAMETERS*** from hyperPosterior based on full data {w/o clusters}
      
      ## hyperparam for [betaparamS]
      SIG_b0paramS.st = riwish(v=qu0S.st+2, 
                               S=t(beta0_prevS.st-beta_prev_fullS.st) %*% (beta0_prevS.st-beta_prev_fullS.st) +
                                 deltaS.st*t(beta0_prevS.st - m0S.st) %*% (beta0_prevS.st - m0S.st) + LAMBS.st) 
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
      
      beta0paramS.st = rmvn(n=1, mu=deltaS.st/(deltaS.st+1)*m0S.st + 1/(deltaS.st+1)*beta_prev_fullS.st, 
                            sigma=1/(deltaS.st+1)*SIG_b0paramS.st) 
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
      
      ## hyperparam for [sig2param]  
      #> - Sample u0 from MH using Gamma(3,1) as proposal 
      u0_propS.st = rgamma(n=1, shape=u0_proposed_shapeS.st, rate=u0_proposed_rateS.st)
      ratio_u0S.st = min(fink(u0_propS.st, A = sig2_prev_full.st*rho_u1S.st/(v0_prevS.st/2), B = rho_u2S.st+1)/
                           fink(u0_prevS.st, A = sig2_prev_full.st*rho_u1S.st/(v0_prevS.st/2), B = rho_u2S.st+1)*
                           dgamma(x = u0_prevS.st, shape = u0_proposed_shapeS.st, rate = u0_proposed_rateS.st)/
                           dgamma(x = u0_propS.st, shape = u0_proposed_shapeS.st, rate = u0_proposed_rateS.st), 1)
      if(runif(n = 1, min = 0, max = 1) < ratio_u0S.st) {
        u0paramS.st = u0_propS.st
      } else {
        u0paramS.st = u0_prevS.st
      }
      
      #> - Sample v0 from Gamma distribution 
      v0paramS.st = rgamma(n=1, shape = rho_v1S.st + u0paramS.st/2, rate = rho_v2S.st+1/sig2_prev_full.st)
      
      
      #------------------------------------------------------------------------------------------------------------------
      #------------------------------------------------------------------------------------------------------------------
      
      #>>> Next I, Sample ****[Main parameters]**** 
      
      #> - Proposal Sampling using proposals OF PRIORS {w/o clusters} 
      beta_prop_fullS.st = rmvn(n=1, mu=beta0paramS.st, sigma=SIG_b0paramS.st*varinf) #~~~~~~~~~~~~~~~~~~~(CHECK)
      sig2_prop_full.st = rinvgamma(n=1, shape=u0paramS.st, scale=v0paramS.st)
      
      
      #> - Sifting a single [betaparam] based on FULL data 
      numerator=0
      denominator=0
      # plugging prior sample first...
      for (i in 1:n.train){
        numerator = numerator + dlognorm_log( y=Y[i], x=matXSstar[i,], 
                                              beta=t(as.matrix(beta_prop_fullS.st)), sig2=sig2_prev_full.st )
        denominator = denominator + dlognorm_log( y=Y[i], x=matXSstar[i,],
                                                  beta=t(as.matrix(beta_prev_fullS.st)), sig2=sig2_prev_full.st )
      }
      ratio_betaS.st = min(exp(numerator-denominator), 1)
      
      U = runif(n = 1, min = 0, max = 1)
      if(U < ratio_betaS.st) {
        betaparam_fullS.st = beta_prop_fullS.st
      } 
      else {
        betaparam_fullS.st = beta_prev_fullS.st
      }
      
      
      #> - Sifting a single [sig2param] based on FULL data 
      numerator=0
      denominator=0
      # plugging prior sample first...
      for (i in 1:n.train){
        numerator = numerator + dlognorm_log( y=Y[i], x=matXSstar[i,], 
                                              beta=t(as.matrix(betaparam_fullS.st)), sig2=sig2_prop_full.st )
        denominator = denominator + dlognorm_log( y=Y[i], x=matXSstar[i,],
                                                  beta=t(as.matrix(betaparam_fullS.st)), sig2=sig2_prev_full.st )    
      }
      ratio_sig2.st = min(exp(numerator-denominator), 1)
      
      U = runif(n = 1, min = 0, max = 1)
      if(U < ratio_sig2.st) {
        sig2param_full.st = sig2_prop_full.st
      } 
      else {
        sig2param_full.st = sig2_prev_full.st
      }
      
      #--------------------------------------------------------------------------------------------------------#
      #--------------------------------------------------------------------------------------------------------#
      #-----------------------------------Now it's time for cluster-wise---------------------------------------#
      #--------------------------------------------------------------------------------------------------------#
      #--------------------------------------------------------------------------------------------------------#
      
      
      #>>> First.., Sample **[covariate parameters_j]**
      
      #XSstar_j=XSstar[cl_membershipS==j] #: already defined
      #ZSj = ZS[cl_membershipS==j]        #: already defined
      #nS_j.st[j]                         #: already defined
      
      #> - proposal Sampling [piparamF_j] using proposals OF POSTERIOR {with clusters} 
      piparamS.st[j] = rbeta( n=1, shape1=g0S.st+sum(ZSj), shape2=h0S.st+nS_j.st[j]-sum(ZSj) )
      # no need...previous value...?
      
      #> - proposal Sampling [lambda2paramF_j] using proposals OF POSTERIOR {with clusters} 
      lambda2paramS.st[j] = rinvgamma( n=1, 
                                       shape=(c0S.st+nS_j.st[j])/2, 
                                       scale=(d0S.st + 
                                                sum( (XSstar_j-kappaparam_prev.st[j,1]-kappaparam_prev.st[j,2]*ZSj)^2 ))/2
      )
      
      #> - proposal Sampling [tau2param_j] using proposals OF POSTERIOR {with clusters}                               
      tau2param[j] = rinvgamma( n=1, 
                                shape=(c0S.st+nS_j.st[j])/2, 
                                scale=(1-zeta_scale)*(d0S.st + 
                                                        sum( (XSstar_j-kappaparam_prev.st[j,1]-kappaparam_prev.st[j,2]*ZSj)^2 ))/2
      )
      #> - proposal Sampling [kappaparam]_j
      KK1j = cbind(1, ZSj)
      KK2j = matrix(c(sum(XSstar_j), sum(XSstar_j*ZSj)), nrow=2, ncol=1) # 2x1 vector
      SIG_inv_j = solve(SIG_inv_k0+(t(KK1j) %*% KK1j)) # solve( ..[matrix].. ) = [matrix]^(-1)
      kappaparam.st[j,] <- rmvn(n=1, mu=SIG_inv_j %*% (SIG_inv_k0 %*% kappa_v0 + KK2j), 
                                sigma = lambda2paramS.st[j]*SIG_inv_j) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
      
      
      #---------------------------------------------------------------------------------------------------------------
      #---------------------------------------------------------------------------------------------------------------
      #>>> Next,... Sample ****[Main parameters_j]****
      
      #>>> subsetting by cluster
      Y_j = Y[cl_membershipS==j]                  #outcome
      matXS_j.st = matXSstar[cl_membershipS==j, ] #covariate
      
      #>>> Sample "new" outcome paramters from PROPOSALS (priors)
      beta_propS.st = rmvn(n=1, mu=beta0paramS.st, sigma=SIG_b0paramS.st*varinf) #~~~~~~~~~~~~~~~~~~~~~~~~~(CHECK)
      sig2_prop.st = rinvgamma(n=1, shape=u0_newS.st, scale=v0_newS.st)
      
      
      
      #>>> - Sifting [betaparam]_j
      numerator=0
      denominator=0
      for (i in 1:length(Y_j)){
        numerator = numerator + dlognorm_log( y=Y_j[i], x=matXS_j.st[i,], 
                                              beta=t(as.matrix(beta_propS.st)), sig2=sig2_prev.st[j] )
        
        denominator = denominator + dlognorm_log( y=Y_j[i], x=matXS_j.st[i,],
                                                  beta=t(as.matrix(beta_prevS.st[j,])), sig2=sig2_prev.st[j] )
      }
      ratio_betaS.st = min(exp(numerator-denominator), 1)
      U = runif(n = 1, min = 0, max = 1)
      if(ratio_betaS.st > U) {
        betaparamS.st[j, ] = beta_propS.st
        beta_p_acceptedS.st <- beta_p_acceptedS.st+1
      } else {
        betaparamS.st[j,] = beta_prevS.st[j,]
      }
      
      
      #>>> - Sifting [sig2param]_j
      numerator=0
      denominator=0
      for (i in 1:length(Y_j)){
        numerator = numerator + 
          dlognorm_log( y=Y_j[i], x=matXS_j.st[i,], beta=t(as.matrix(betaparamS.st[j,])), sig2=sig2_prop.st ) + 
          log(dinvgamma(sig2_prop.st, shape=u0paramS.st, scale=v0paramS.st)) + 
          log(dinvgamma(sig2_prev.st[j], shape=u0_newS.st, scale=v0_newS.st))
        
        denominator = denominator + 
          dlognorm_log( y=Y_j[i], x=matXS_j.st[i,], beta=t(as.matrix(betaparamS.st[j,])), sig2=sig2_prev.st[j] ) + 
          log(dinvgamma(sig2_prev.st[j], shape=u0paramS.st, scale=v0paramS.st)) + 
          log(dinvgamma(sig2_prop.st, shape=u0_newS.st, scale=v0_newS.st))
      }
      ratio_sig2.st = min(exp(numerator-denominator), 1)
      U = runif(n = 1, min = 0, max = 1)
      if(U < ratio_sig2.st) {
        sig2param.st[j] = sig2_prop.st
        sig2_p_accepted.st <- sig2_p_accepted.st+1
      }
      else {
        sig2param.st[j] = sig2_prev.st[j]
      }
      
      #################################################################################################################
      #################################################################################################################
      #################################################################################################################
      #################################################################################################################
      ###> [2] Calculate differences to ensure VALIDITY of [System Equation]
      #################################################################################################################
      #>>> [main] I. For clean beta_j_1 (the problematic***slope on clean X)
      betaparamclean_j = betaparamS.st[j,2]*lambda2paramS.st[j]/(lambda2paramS.st[j]-tau2param[j])
      
      #>>> [main] II. For clean sig2_j (the problematic***variance param for LSN on clean X)
      sig2paramclean_j = sig2param.st[j] - 
        ( tau2param[j] * betaparamclean_j^2 * (lambda2paramS.st[j]-tau2param[j]) )/lambda2paramS.st[j]
      
      #>>> [supplementary] How to ensure They are positive ? 
      difference1 = lambda2paramS.st[j]-tau2param[j] # ensure clean lambda2 is positive
      difference2 = sig2paramclean_j             # ensure clean sig2 are positive  
      
      #> -while looping count- iteration for each j
      print(c(difference1, difference2)) #, sig2param.st[j], u0paramS.st, v0paramS.st))
      count_iterations = count_iterations+1
    } # End of while (Don't stop until all the three differences are positive!!!!)
    print(paste("count=", count_iterations))
    
    # Values that meet the condition are obtained. Now....
    ###> Collect the params and compute the clean parameters using [System Equation].
    #################################################################################################################
    
    #kappaparamclean = matrix(0, nrow = JS, ncol = ncol(kappaparam.st))     # already defined
    #lambda2paramcleanS = numeric(JS)                                       # already defined
    #betaparamcleanS =  matrix(0, nrow = JS, ncol = ncol(betaparamS.st))    # already defined
    #sig2paramclean =numeric(JS)                                           # already defined
    
    #xSclean = numeric(n.train)                                             # already defined
    
    
    lambda2paramcleanS[j] = lambda2paramS.st[j] - tau2param[j]    
    
    kappaparamclean[j,] = kappaparam.st[j,]
    
    betaparamcleanS[j,2] = betaparamS.st[j,2]*lambda2paramS.st[j]/(lambda2paramS.st[j]-tau2param[j]) #for Beta1
    
    betaparamcleanS[j,1] = betaparamS.st[j,1] - 
      betaparamS.st[j,2]*kappaparam.st[j,1]*tau2param[j]/(lambda2paramS.st[j]-tau2param[j]) #for Beta0
    
    betaparamcleanS[j,3] = betaparamS.st[j,3] - 
      betaparamS.st[j,2]*kappaparam.st[j,2]*tau2param[j]/(lambda2paramS.st[j]-tau2param[j]) #for Beta2
    
    sig2paramclean[j] = sig2param.st[j] - 
      ( tau2param[j] * betaparamcleanS[j,2]^2 * (lambda2paramS.st[j]-tau2param[j]) )/lambda2paramS.st[j]   #for sig2
    
    #################################################################################################################
  } # End of [j]
  
  
  #### For NEXT!!!
  ###> 1. Store the clean data parameters
  ###> 2. Set previous set of parameter to the [[[current]]] for the next iteration
  ## - Manner! main param for next iteration
  
  ##]]]]]] - Let's SEE..
  #-beta0_prevS.st = beta0paramS.st #-------------------------------------------------already defined
  #-SIG_b0_prevS.st = SIG_b0paramS.st #-----------------------------------------------already defined
  #-u0_prevS.st = u0paramS.st #-------------------------------------------------------already defined
  #-v0_prevS.st = v0paramS.st #-------------------------------------------------------already defined
  
  #-beta_prevS.st = betaparamS.st #---------------------------------------------------already defined
  #-sig2_prev.st = sig2param.st #-----------------------------------------------------already defined
  #-beta_prev_fullS.st = rmvn(n=1, mu=beta0_prevS.st, sigma=SIG_b0_prevS.st) # prebviously given!!!!
  #-sig2_prev_full.st = rgamma(n=1, shape=u0_prevS.st, rate=v0_prevS.st)     # prebviously given!!!!
  
  #piparam_prevS.st = piparamS.st           # we do not need. unaffected.
  #piparam_prevS.clean = piparamS.st        # we do not need. unaffected.
  
  #-lambda2param_prevS.st = lambda2paramS.st # do we need???? #-----------------------already defined
  #-lambda2param_prevS.clean = lambda2paramS.st #???? #-------------------------------already defined
  #-kappaparam_prev.st = kappaparam.st #----------------------------------------------already defined
  #-kappaparam_prev.clean = kappaparam.st #-------------------------------------------already defined
  
  #-tau2_prev = tau2param  # tau2param:::(1-zeta_scale)*lambda2paramS.st #------------already defined
  
  
  ##]]]]]] - The truth is.......                                     
  #>>>> WHY? ++++++++++++++++++++++++++++++++++++++++++ to sample [betaparamS] and [sig2param] next time....[O]: ONLY USEFUL
  beta0_prevS.st = beta0paramS.st #:::::::::::::::: useful for proposal sample of "SIG_b0paramS.st" 
  SIG_b0_prevS.st = SIG_b0paramS.st
  u0_prevS.st = u0paramS.st #::::::::::::::::::::::::::::::::::::::::: useful to sift "u0paramS.st"
  v0_prevS.st = v0paramS.st #::::::::::::::::::::::::::::::::::::::::: useful to sift "u0paramS.st"
  
  #>>>> WHY? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ to sample [Y] next time....but [X]
  beta_prevS.st = betaparamS.st # ::::::::::::::::::::::::::::::::::: useful to sift "betaparamS.st"
  sig2_prev.st = sig2param.st # ::::::::::::::::::::::::::::::::::::: useful to sift "betaparamS.st" and "sig2param.st"
  #beta_prev_fullS.st  # this is not updated!!!
  #sig2_prev_full.st   # this is not updated!!!
  
  #>>>> WHY? +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ to sample [ZS] next time....but [X]
  #piparam_prevS.st = piparamS.st           
  #piparam_prevS.clean = piparamS.st                            
  
  #>>>> WHY? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ to sample [XS|ZS] next time....but [X]
  lambda2param_prevS.st = lambda2paramS.st
  lambda2param_prevS.clean = lambda2paramcleanS #***********
  kappaparam_prev.st = kappaparam.st # ::::::::::::: useful for proposal sample of "lambda2paramS", and thus "tau2param"     
  kappaparam_prev.clean = kappaparamclean #*****************
  
  #>>>> WHY? ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ to sample [XSstar|XS] next time....but [X]                                       
  tau2_prev = tau2param                                  
  
  
  
  
  ###############################################################################################################
  ###[3] Calculating Loglikelihood to monitor the chain convergence ---------------------------------------------
  ###############################################################################################################                                        
  #>>> First
  loglike_r = numeric(n.train)
  
  for(i in 1:n.train) {
    xS_i = matXSstar[i, ]      # first row..{1,x,z}
    jS = cl_membershipS[i] # membership ID
    
    # loglike_r[i] = dlognorm_log( y=Y[i], x=xS_i, beta=t(as.matrix(betaparamS.st[jS,])), sig2=sig2param.st[jS] ) + 
    #   dnorm(x=xS_i[2], mean=kappaparam.st[j,1] + kappaparam.st[j,2]*ZS[i], sd=sqrt(lambda2paramS.st[jS]), log=TRUE) +
    #   dbinom(x=xS_i[3], size=1, prob=piparamS.st[jS], log=TRUE) 
    loglike_r[i] = dlognorm_log( y=Y[i], x=xS_i, beta=t(as.matrix(betaparamcleanS[jS,])), sig2=sig2paramclean[jS] ) + 
      dnorm(x=xS_i[2], mean=kappaparamclean[j,1] + kappaparamclean[j,2]*ZS[i], sd=sqrt(lambda2paramcleanS[jS]), log=TRUE) +
      dbinom(x=xS_i[3], size=1, prob=piparamS.st[jS], log=TRUE)
  }  
  loglikelihood.st[,r] = loglike_r
  
  if(r > r_convergence) {
    
    # W/O CLUSTER
    list_beta0paramS.st[[r-r_convergence]] = beta0paramS.st
    list_SIG_b0paramS.st[[r-r_convergence]] = SIG_b0paramS.st
    list_u0paramS.st[[r-r_convergence]] = u0paramS.st
    list_v0paramS.st[[r-r_convergence]] = v0paramS.st
    
    list_betaparamfullS.st[[r-r_convergence]] = betaparam_fullS.st
    list_sig2paramfull.st[[r-r_convergence]] = sig2param_full.st
    
    
    # With Cluster-wise
    list_piparamS.st[[r-r_convergence]] = piparamS.st     
    list_lambda2paramS.st[[r-r_convergence]] = lambda2paramS.st
    list_kappaparam.st[[r-r_convergence]] = kappaparam.st
    list_tau2param[[r-r_convergence]] = tau2param
    
    list_betaparamS.st[[r-r_convergence]] = betaparamS.st  
    list_sig2param.st[[r-r_convergence]] = sig2param.st 
    
    ################# STORE the CLEAN parameters ##################
    list_kappaparam.clean[[r-r_convergence]] = kappaparamclean     #for Xstar
    list_lambda2paramS.clean[[r-r_convergence]] = lambda2paramcleanS    #for Xstar
    list_betaparamS.clean[[r-r_convergence]] = betaparamcleanS    #for S
    list_sig2param.clean[[r-r_convergence]] = sig2paramclean    #for S
    
    # list_xclean[[r-r_convergence]] = xclean #~xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Q.DAVID (can we impute?)
    
  }
  print(paste("r=",r))
  trace(what = "print", where = getNamespace("base"), exit = flush.console, print = FALSE)
}

########################################
####### End of Gibbs Sampler ###########
########################################

beta_p_acceptedS.st/(total_iter*JS)
sig2_p_accepted.st/(total_iter*JS)
par(mfrow=c(1,1))
plot(colSums(loglikelihood.st),type="l")

lppd_EX_S.Y.st = mean( colSums(loglikelihood.st) ); lppd_EX_S.Y.st
# -16503.5 (0.7), -16501.98 (0.6),  -16507.48 (0.5), -16508.09 (0.4)

###> Investigation for "sig2"
list_sig2param.clean[1]
# investigation for [sig2param.clean]
sig2.vec_cl1.st = sapply(list_sig2param.clean, function(x) x[1])
sig2.vec_cl2.st = sapply(list_sig2param.clean, function(x) x[2])
sig2.vec_cl3.st = sapply(list_sig2param.clean, function(x) x[3])
sig2.vec_cl4.st = sapply(list_sig2param.clean, function(x) x[4])
sig2.vec_cl5.st = sapply(list_sig2param.clean, function(x) x[5])
sig2.vec_cl6.st = sapply(list_sig2param.clean, function(x) x[6])

sig2.st.df <- data.frame(
  value = c(sig2.vec_cl1.st, 
            sig2.vec_cl2.st, sig2.vec_cl3.st, sig2.vec_cl5.st, sig2.vec_cl5.st, sig2.vec_cl6.st),
  group = factor(rep(1:6, each = 5))
)
#option A (together)                    
#ggplot(sig2.df, aes(x = value, color = group)) +
#                     geom_density() +
#                     labs(title = "Density for Each sig2.vec", x = "sig2", y = "Density") + 
#                     theme_minimal()

#option B (separate)
ggplot(sig2.st.df, aes(x = value)) +
  geom_density(size=1) +
  facet_wrap(~ group, scales = "free") +
  labs(title = "Density for Each sig2.st.vec", x="sig2", y="Density") +
  theme_minimal()

tapply(sig2.st.df$value, sig2.st.df$group, mean)


aggregate(value ~ group, data=sig2.st.df, FUN=mean)

# # Calculate the average sig2
sn <- length(list_sig2param.clean); sn
avg_sig2.vec1.st <- sum(sig2.vec_cl1.st)/sn; avg_sig2.vec1.st #1.119387
avg_sig2.vec2.st <- sum(sig2.vec_cl2.st)/sn; avg_sig2.vec2.st #1.010708
avg_sig2.vec3.st <- sum(sig2.vec_cl3.st)/sn; avg_sig2.vec3.st #1.168137
avg_sig2.vec4.st <- sum(sig2.vec_cl4.st)/sn; avg_sig2.vec4.st #1.168856
avg_sig2.vec5.st <- sum(sig2.vec_cl5.st)/sn; avg_sig2.vec5.st #1.099336
avg_sig2.vec6.st <- sum(sig2.vec_cl6.st)/sn; avg_sig2.vec6.st #1.085018




###> investigation for [betaparamS]
# Sum all the matrices
sum_matrixS.st <- Reduce("+", list_betaparamS.clean)

# Number of matrices
r <- length(list_betaparamS.clean)

# Calculate the average matrix
avg_beta.matS.st <- sum_matrixS.st / r
avg_beta.matS.st


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++ PREDICTION +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# Expected Value of LN: exp(sum(x*beta) # x^T*B
#                       = exp( mu + 0.5*sig2 ) 

###> with Training data
expvalS.train.st = matrix(0, nrow = n.train, ncol = total_iter-r_convergence)
mu.train.st = matrix(0, nrow = n.train, ncol = total_iter-r_convergence)

### - prediction 
for(r in 1:(total_iter-r_convergence)) {
  piparamcleanS = list_piparamS.st[[r]]
  lambda2paramcleanS = list_lambda2paramS.clean[[r]]
  kappaparamclean = list_kappaparam.clean[[r]]
  betaparamcleanS = list_betaparamS.clean[[r]]
  sig2paramclean = list_sig2param.clean[[r]]
  
  for(i in 1:n.train) {
    j = cl_membershipS[i]
    #dirtyXS = rnorm(n=1, mean=kappaparamclean[j,1] + kappaparamclean[j,2]*ZS[i], sd=sqrt(lambda2paramcleanS[j]))
    dirtyXS = XSstar[i]
    expvalS.train.st[i,r] <- exp( betaparamcleanS[j,1] + betaparamcleanS[j,2]*dirtyXS + 
                                    betaparamcleanS[j,3]*ZS[i] - 0.5*sig2paramclean[j] -0.6)
    mu.train.st[i,r] <- betaparamcleanS[j,1] + betaparamcleanS[j,2]*dirtyXS + betaparamcleanS[j,3]*ZS[i] 
  }
}

EX_S.train.st <- apply(X=expvalS.train.st, MARGIN=1, FUN=mean); EX_S.train.st
Y
mu.vec.train.st <- apply(X=mu.train.st, MARGIN=1, FUN=mean); mu.vec.train.st

#options(scipen = 999)
format( aggregate(Y, by=list(cl_membershipS), FUN=summary), scientific=FALSE)
format( aggregate(EX_S.train.st, by=list(cl_membershipS), FUN=summary), scientific=FALSE)


##############################################################################################################
##############################################################################################################
#-----------------------------------Zeta-LAB-----------------------------------#

### - zeta:  Large (0.9)
EX_S.train.st.Large <- apply(X=expvalS.train.st, MARGIN=1, FUN=mean)
#options(scipen = 999)
format( aggregate(Y, by=list(cl_membershipS), FUN=summary), scientific=FALSE)
format( aggregate(EX_S.train.st.Large, by=list(cl_membershipS), FUN=summary), scientific=FALSE)


# ### - zeta:  Medium (0.8)
# EX_S.train.st.Medium <- apply(X=expvalS.train.st, MARGIN=1, FUN=mean)
# #options(scipen = 999)
# format( aggregate(Y, by=list(cl_membershipS), FUN=summary), scientific=FALSE)
# format( aggregate(EX_S.train.st.Medium, by=list(cl_membershipS), FUN=summary), scientific=FALSE)
# 
# 
# ### - zeta:  Small (0.7)
# EX_S.train.st.Small <- apply(X=expvalS.train.st, MARGIN=1, FUN=mean)
# #options(scipen = 999)
# format( aggregate(Y, by=list(cl_membershipS), FUN=summary), scientific=FALSE)
# format( aggregate(EX_S.train.st.Small, by=list(cl_membershipS), FUN=summary), scientific=FALSE)


################### Clusterwise predictive density plots #######################
par(mfrow = c(3,2))

for(j in 1:JS) {
  hist(log(Y)[cl_membershipS==j], freq=FALSE, breaks=50, col="white", xlim=c(5,15), ylim=c(0, 2.8),
       xlab = paste("Cluster", j))
  lines(density(log(EX_S.train)[cl_membershipS==j]), col="red", lwd=2)     # clean 
  lines(density(log(EX_S.train.dirty)[cl_membershipS==j]), col="black", lty = "dotted") # dirty
  lines(density(log(EX_S.train.st.Large)[cl_membershipS==j]), col="blue", lwd=1) # dirty with correction
  #lines(density(log(EX_S.train.st.Medium)[cl_membershipS==j]), col="darkmagenta", lwd=1) # dirty with correction
  #lines(density(log(EX_S.train.st.Small)[cl_membershipS==j]), col="darkgreen", lwd=1) # dirty with correction
  
  if (j == 1) {
    legend("topright", legend = c("Clean", "Dirty", "Correction"), #, "Correction Medium", "Correction Small"),
           col = c("red", "black", "blue"), #, "darkmagenta", "darkgreen"), 
           lty = c(1, 2, 1), # 1, 1), 
           lwd = c(3, 2, 3), #, 1, 1), 
           cex = 0.4
    )
  }
}    # plot: 1500 x 1500  

##############################################################################################################
##############################################################################################################
#####>>>>>>>>>>>>>>>>>> Out-of-sample test for SSAE <<<<<<<<<<<<<<<<<<<<<#######
# ?
#####> with [Clean]
expvalS.test = matrix(0, nrow = n.test, ncol = total_iter-r_convergence)

for(r in 1:(total_iter-r_convergence)) {
  piparamS = list_piparamS[[r]]
  betaparamS = list_betaparamS[[r]]
  sig2param = list_sig2param[[r]]
  
  for(i in 1:n.test) {
    j = test.df$Classes[i]
    #cleanXS = rnorm(n=1, mean=X_barS[j], sd=sqrt(lambda2paramS[j]))
    cleanXS = XS.test[i]
    expvalS.test[i,r] <- exp( betaparamS[j,1] + betaparamS[j,2]*cleanXS +
                                betaparamS[j,3]*ZS.test[i] - 0.5*sig2param[j] )
  }
}
EX_S.Y.test <- apply(X=expvalS.test, MARGIN=1, FUN=mean)

#options(scipen = 999)
cl_membership.test = test.df$Classes
table(cl_membership.test)

format( aggregate(log(test.df$yAvg), by=list(cl_membership.test), FUN=summary), scientific=FALSE)
format( aggregate(log(EX_S.Y.test), by=list(cl_membership.test), FUN=summary), scientific=FALSE)

EX_S.Y.test
SSPE.LN <- sum( (log(EX_S.Y.test) - log(test.df$yAvg))^2 ); SSPE.LN   #784.5236
SAPE.LN <- sum( abs(log(EX_S.Y.test) - log(test.df$yAvg)) ); SAPE.LN  #415.2161 


# ?
#####> with [Dirty]
expvalS.test.xx = matrix(0, nrow = n.test, ncol = total_iter-r_convergence)

for(r in 1:(total_iter-r_convergence)) {
  piparamS.xx = list_piparamS.xx[[r]]
  betaparamS.xx = list_betaparamS.xx[[r]]
  sig2param.xx = list_sig2param.xx[[r]]
  
  for(i in 1:n.test) {
    j = test.df$Classes[i]
    #cleanXS = rnorm(n=1, mean=X_barS[j], sd=sqrt(lambda2paramS[j]))
    dirtyXS = XSstar.test[i]
    expvalS.test.xx[i,r] <- exp( betaparamS.xx[j,1] + betaparamS.xx[j,2]*dirtyXS +
                                   betaparamS.xx[j,3]*ZS.test[i] - 0.5*sig2param.xx[j] )
  }
}
EX_S.test.Y.xx <- apply(X=expvalS.test.xx, MARGIN=1, FUN=mean)

#options(scipen = 999)
cl_membership.test = test.df$Classes
table(cl_membership.test)

format( aggregate(log(test.df$yAvg), by=list(cl_membership.test), FUN=summary), scientific=FALSE)
format( aggregate(log(EX_S.test.Y.xx), by=list(cl_membership.test), FUN=summary), scientific=FALSE)

EX_S.test.Y.xx
SSPE.LN.xx <- sum( (log(EX_S.test.Y.xx) - log(test.df$yAvg))^2 ); SSPE.LN.xx   #840.0246
SAPE.LN.xx <- sum( abs(log(EX_S.test.Y.xx) - log(test.df$yAvg)) ); SAPE.LN.xx  #430.2199  


# ?
#####> with [Gustafson] with zeta:[0.9] 
expvalS.test.st = matrix(0, nrow = n.test, ncol = total_iter-r_convergence)

for(r in 1:(total_iter-r_convergence)) {
  piparamS = list_piparamS.st[[r]]
  betaparamS = list_betaparamS.clean[[r]]
  sig2param = list_sig2param.clean[[r]]
  
  for(i in 1:n.test) {
    j = test.df$Classes[i]
    #cleanXS = rnorm(n=1, mean=X_barS[j], sd=sqrt(lambda2paramS[j]))
    cleanXS = XSstar.test[i]
    expvalS.test.st[i,r] <- exp( betaparamS[j,1] + betaparamS[j,2]*cleanXS +
                                   betaparamS[j,3]*ZS.test[i] - 0.5*sig2param[j] )
  }
}
EX_S.test.Y.st <- apply(X=expvalS.test.st, MARGIN=1, FUN=mean)

#options(scipen = 999)
cl_membership.test = test.df$Classes
table(cl_membership.test)

format( aggregate(log(test.df$yAvg), by=list(cl_membership.test), FUN=summary), scientific=FALSE)
format( aggregate(log(EX_S.test.Y.st), by=list(cl_membership.test), FUN=summary), scientific=FALSE)

EX_S.test.Y.st
SSPE.LN.st <- sum( (log(EX_S.test.Y.st) - log(test.df$yAvg))^2 ); SSPE.LN.st   #839.8009
SAPE.LN.st <- sum( abs(log(EX_S.test.Y.st) - log(test.df$yAvg)) ); SAPE.LN.st  #427.6068  




#
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#
###############################################################################
###############################################################################
###############################################################################
############################ S(t) with Clean ##################################
###############################################################################
###############################################################################
###############################################################################
set.seed(1)  # For reproducibility
n_simulations <- n.train
#n_simulations <- n.test


#####> Step 1: Define the parameters for each group
group_params <- list(
  list(xi_F=mean(EX_F.train[cl_membershipF==1]), psi_F=avg_psi.vec1, 
       meanlog_S=mean(mu.vec.train[cl_membershipF==1]), sdlog_S=sqrt(avg_sig2.vec1), 
       n=table(cl_membershipF)[1]),                                             #for Grp 01
  list(xi_F=mean(EX_F.train[cl_membershipF==2]), psi_F=avg_psi.vec2, 
       meanlog_S=mean(mu.vec.train[cl_membershipF==2]), sdlog_S=sqrt(avg_sig2.vec2), 
       n=table(cl_membershipF)[2]),                                             #for Grp 02
  list(xi_F=mean(EX_F.train[cl_membershipF==3]), psi_F=avg_psi.vec3, 
       meanlog_S=mean(mu.vec.train[cl_membershipF==3]), sdlog_S=sqrt(avg_sig2.vec3), 
       n=table(cl_membershipF)[3]),                                             #for Grp 03
  list(xi_F=mean(EX_F.train[cl_membershipF==4]), psi_F=avg_psi.vec4, 
       meanlog_S=mean(mu.vec.train[cl_membershipF==4]), sdlog_S=sqrt(avg_sig2.vec4), 
       n=table(cl_membershipF)[4]),                                             #for Grp 04
  list(xi_F=mean(EX_F.train[cl_membershipF==5]), psi_F=avg_psi.vec5, 
       meanlog_S=mean(mu.vec.train[cl_membershipF==5]), sdlog_S=sqrt(avg_sig2.vec5), 
       n=table(cl_membershipF)[5]),                                             #for Grp 05
  list(xi_F=mean(EX_F.train[cl_membershipF==6]), psi_F=avg_psi.vec6, 
       meanlog_S=mean(mu.vec.train[cl_membershipF==6]), sdlog_S=sqrt(avg_sig2.vec6), 
       n=table(cl_membershipF)[6])                                              #for Grp 06
); group_params





#####> Step 2: Simulate the total claim amount for each group
E.total_claims <- numeric(n_simulations)
E.grp_claims <- numeric(n_simulations)

# - without cluster:S(t)
for (grp in 1:6) {
  params <- group_params[[grp]]
  claim_counts <- rnbinom(n=n_simulations, size=params$psi_F, mu=params$xi_F)
  
  for (i in 1:n_simulations) {
    if (claim_counts[i] > 0) {
      severities <- rlnorm( n=claim_counts[i], meanlog=params$meanlog_S, sdlog=params$sdlog_S )
      E.total_claims[i] <- E.total_claims[i] + sum(severities)
    }
  }
}
summary(E.total_claims)
# Min.  1st Qu.   Median     Mean   3rd Qu.     Max. 
#   0   506881   1793785  4394145  5201678 64686438 

summary(log(E.total_claims))


###> Approach 1.
par(mfrow=c(1,1))
hist(log(E.total_claims), breaks=80, xlab="log(S(t))", ylab="Density", prob=T, col="white", 
     xlim=c(0,20), ylim=c(0,0.25))
lines( density(x=log(E.total_claims)), col="blue", lwd=1 )
#plot( density(x=E.total_claims), col="blue", lwd=1 )

# - with clusters: cluster-wise S(t)
par(mfrow=c(2,3))
for (grp in 1:6) {
  params <- group_params[[grp]]
  claim_counts <- rnbinom(n=n_simulations, size=params$psi_F, mu=params$xi_F)
  
  for (i in 1:n_simulations) {
    if (claim_counts[i] > 0) {
      severities <- rlnorm(claim_counts[i], meanlog=params$meanlog_S, sdlog=params$sdlog_S)
      E.grp_claims[i] <- sum(severities)
    }
  }
  #> summary
  cat("\nSummary for Group", grp, "\n")
  print(summary(E.grp_claims))
  summary(E.grp_claims)
  #> plot
  #hist(log(E.grp_claims), breaks=80, prob=T, 
  #     xlab=paste("log(S(t)) by Cluster", grp), ylab="Density", ylim=c(0,0.2))
  #lines( density(x=log(E.grp_claims)), col="blue", lwd=1 )
  plot( x=density(x=log(E.grp_claims)), xlab=paste("log(S(t)) for Cluster", grp), 
        ylab="Density", col="blue", lwd=2, xlim=c(0,20), ylim=c(0,0.25)  )
}

# Summary for Group 1 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0        0        0   370199   181971 11103825 
# 
# Summary for Group 2 
# Min.  1st Qu.   Median     Mean   3rd Qu.     Max. 
# 0      20373    413766  3120182  2969400 95137620 
# 
# Summary for Group 3  
# Min.  1st Qu.   Median     Mean   3rd Qu.     Max. 
# 0      20439    202712  1979686  1060599 95137620 
# 
# Summary for Group 4 
# Min.  1st Qu.   Median     Mean   3rd Qu.     Max. 
# 0      42289    245342  1538978  1020689 95137620 
# 
# Summary for Group 5 
# Min.  1st Qu.   Median     Mean   3rd Qu.     Max. 
# 0      24041    143685  1013176   568862 95137620 
# 
# Summary for Group 6 
# Min.  1st Qu.   Median     Mean   3rd Qu.     Max. 
# 0      22410    102796   781703   428263 95137620 



###> Approach 2.
# Initialize a list to store the density data
par(mfrow=c(1,1))
# Initialize a list to store the density data
densities <- list()

# Perform the simulation for each group
for (grp in 1:6) {
  params <- group_params[[grp]]
  claim_counts <- rnbinom(n=n_simulations, size=params$psi_F, mu=params$xi_F)
  
  E.grp_claims <- numeric(n_simulations)
  for (i in 1:n_simulations) {
    if (claim_counts[i] > 0) {
      severities <- rlnorm(claim_counts[i], meanlog=params$meanlog_S, sdlog=params$sdlog_S)
      E.grp_claims[i] <- sum(severities)
    }
  }
  
  # Compute the density of the log of E.grp_claims
  densities[[grp]] <- density(log(E.grp_claims))
}

# Set up the plot
#plot(NULL, xlim=c(min(sapply(densities, function(d) min(d$x))), max(sapply(densities, function(d) max(d$x)))), 
#     ylim=c(0, max(sapply(densities, function(d) max(d$y)))), 
#     xlab="log(S(t))", ylab="Density", main="Density plots for all clusters")

hist(log(E.total_claims), breaks=80, xlab="log(S(t))", ylab="Density", prob=T, col="white", 
     xlim=c(0,20), ylim=c(0,0.25))

# Define line types and colors for the plots
line_types <- 1:6
colors <- c("red", "blue", "green", "purple", "orange", "brown")

# Plot each density in the same frame
for (grp in 1:6) {
  lines(densities[[grp]], lty=line_types[grp], col=colors[grp], lwd=2, xlim=c(0,20), ylim=c(0,0.25))
}
# Add a legend
legend("topleft", legend=paste("Cluster", 1:6), col=colors, lty=line_types, lwd=2)




###> Approach 3.
# 
# # Calculate total sample size
# #total_n <- sum(sapply(group_params, function(p) p$n))
# 
# # Simulate data points for each group
# data <- numeric(0)
# 
# for (p in group_params) {
#   group_data <- rlnorm(n=p$n, meanlog=p$meanlog_S, sdlog=p$sdlog_S)
#   data <- c(data, group_data)
# }
# 
# # Plot the density of the mixture distribution
# par(mfrow=c(1,1))
# plot( density(log(data)), xlab="S(t)", ylab="Density", ylim=c(0,0.4) )
# 
# 
# # Add individual densities for each group (scaled by their mixing coefficient)
# colors <- rainbow(6)
# for (i in 1:6) {
#   group_data <- rlnorm(group_params[[i]]$n, 
#                        meanlog = group_params[[i]]$meanlog_S, 
#                        sdlog = group_params[[i]]$sdlog_S)
#   lines( density(log(group_data)), col=colors[i], lwd=1, lty=1 )
# }
# legend("topright", legend=paste("Group", 1:6), lwd=2, lty=1, col=colors )






###############################################################################
###############################################################################
###############################################################################
############################ S(t) with dirty ##################################
###############################################################################
###############################################################################
###############################################################################

#####> Step 1: Define the parameters for each group
group_params.xx <- list(
  list(xi_F=mean(EX_F.train[cl_membershipF==1]), psi_F=avg_psi.vec1, 
       meanlog_S=mean(mu.vec.train.dirty[cl_membershipF==1]), sdlog_S=sqrt(avg_sig2.vec1.xx), 
       n=table(cl_membershipF)[1]),                                             #for Grp 01
  list(xi_F=mean(EX_F.train[cl_membershipF==2]), psi_F=avg_psi.vec2, 
       meanlog_S=mean(mu.vec.train.dirty[cl_membershipF==2]), sdlog_S=sqrt(avg_sig2.vec2.xx), 
       n=table(cl_membershipF)[2]),                                             #for Grp 02
  list(xi_F=mean(EX_F.train[cl_membershipF==3]), psi_F=avg_psi.vec3, 
       meanlog_S=mean(mu.vec.train.dirty[cl_membershipF==3]), sdlog_S=sqrt(avg_sig2.vec3.xx), 
       n=table(cl_membershipF)[3]),                                             #for Grp 03
  list(xi_F=mean(EX_F.train[cl_membershipF==4]), psi_F=avg_psi.vec4, 
       meanlog_S=mean(mu.vec.train.dirty[cl_membershipF==4]), sdlog_S=sqrt(avg_sig2.vec4.xx), 
       n=table(cl_membershipF)[4]),                                             #for Grp 04
  list(xi_F=mean(EX_F.train[cl_membershipF==5]), psi_F=avg_psi.vec5, 
       meanlog_S=mean(mu.vec.train.dirty[cl_membershipF==5]), sdlog_S=sqrt(avg_sig2.vec5.xx), 
       n=table(cl_membershipF)[5]),                                             #for Grp 05
  list(xi_F=mean(EX_F.train[cl_membershipF==6]), psi_F=avg_psi.vec6, 
       meanlog_S=mean(mu.vec.train.dirty[cl_membershipF==6]), sdlog_S=sqrt(avg_sig2.vec6.xx), 
       n=table(cl_membershipF)[6])                                              #for Grp 06
); group_params.xx




#####> Step 2: Simulate the total claim amount for each group
E.total_claims.xx <- numeric(n_simulations)
E.grp_claims.xx <- numeric(n_simulations)

# - without cluster:S(t)
for (grp in 1:6) {
  params <- group_params.xx[[grp]]
  claim_counts <- rnbinom(n=n_simulations, size=params$psi_F, mu=params$xi_F)
  
  for (i in 1:n_simulations) {
    if (claim_counts[i] > 0) {
      severities <- rlnorm( n=claim_counts[i], meanlog=params$meanlog_S, sdlog=params$sdlog_S )
      E.total_claims.xx[i] <- E.total_claims.xx[i] + sum(severities)
    }
  }
}
summary(E.total_claims.xx)
summary(log(E.total_claims.xx))

par(mfrow=c(1,1))
hist(log(E.total_claims.xx), breaks=80, xlab="log(S(t))", ylab="Density", prob=T, col="white")#, xlim=c(1:50000))
lines( density(x=log(E.total_claims.xx)), col="blue", lwd=1 )
#plot( density(x=E.total_claims), col="blue", lwd=1 )

# - with clusters: cluster-wise S(t)
par(mfrow=c(2,3))
for (grp in 1:6) {
  params <- group_params.xx[[grp]]
  claim_counts <- rnbinom(n=n_simulations, size=params$psi_F, mu=params$xi_F)
  
  for (i in 1:n_simulations) {
    if (claim_counts[i] > 0) {
      severities <- rlnorm(claim_counts[i], meanlog=params$meanlog_S, sdlog=params$sdlog_S)
      E.grp_claims.xx[i] <- sum(severities)
    }
  }
  #> summary
  cat("\nSummary for Group", grp, "\n")
  print(summary(E.grp_claims.xx))
  summary(E.grp_claims.xx)
  #> plot
  #hist(log(E.grp_claims), breaks=80, prob=T, 
  #     xlab=paste("log(S(t)) by Cluster", grp), ylab="Density")
  #lines( density(x=log(E.grp_claims)), col="blue", lwd = 1 )
  plot( x=density(x=log(E.grp_claims.xx)), xlab=paste("log(S(t)) for Cluster", grp), 
        ylab="Density", col="blue", lwd=2  )
}



#####> Step 3: Show the mixture results ???? ------------------------------------------------ # 

# Calculate total sample size
#total_n <- sum(sapply(group_params, function(p) p$n))

# Simulate data points for each group
data.xx <- numeric(0)

for (p in group_params.xx) {
  group_data <- rlnorm(n=p$n, meanlog=p$meanlog_S, sdlog=p$sdlog_S)
  data.xx <- c(data.xx, group_data)
}

# Plot the density of the mixture distribution
par(mfrow=c(1,1))
plot( density(log(data.xx)), xlab="S(t)", ylab="Density", ylim=c(0,0.4) )


# Add individual densities for each group (scaled by their mixing coefficient)
colors <- rainbow(6)
for (i in 1:6) {
  group_data <- rlnorm(group_params.xx[[i]]$n, 
                       meanlog = group_params.xx[[i]]$meanlog_S, 
                       sdlog = group_params.xx[[i]]$sdlog_S)
  lines( density(log(group_data)), col=colors[i], lwd=1, lty=1 )
}
legend("topright", legend=paste("Group", 1:6), lwd=2, lty=1, col=colors )


################################################################################
################################################################################
############################## simex practice ##################################
################################################################################
################################################################################
library(simex)
# binomial(link = "logit")
# gaussian(link = "identity")
# Gamma(link = "inverse")
# poisson(link = "log")
# quasibinomial(link = "logit")
# quasipoisson(link = "log")

###> GLM (Gold standard) -------------------------------------------------------
fit_m1.S <- glm( N*Y ~ XS + factor(ZS), family=Gamma(link = "log") )
fit_m1.Y <- glm( Y ~ XS + factor(ZS), family=Gamma(link = "log") )
#fit_m1 <- glm( log(N.test*Y.test) ~ XS.test + factor(ZS.test), family=gaussian(link="identity") )
#fit_m1 <- glm( log(N*Y) ~ XS + factor(ZS), family=gaussian(link="identity") )
summary(fit_m1.S)
summary(fit_m1.Y)
fit_m1.S$fitted.values
fit_m1.Y$fitted.values
summary(fit_m1.S$fitted.values)
summary(fit_m1.Y$fitted.values)
#str(fit_m1.S$fitted.values)
#str(fit_m1.Y$fitted.values)
E.claims.S.glm <- fit_m1.S$fitted.values
E.claims.Y.glm <- fit_m1.Y$fitted.values

new_data <- data.frame(XS=XS.test, ZS=ZS.test); new_data
pred.val.S.glm = predict(fit_m1.S, newdata = new_data, type = "response")
plot( x=density(pred.val.S.glm) )
pred.val.Y.glm = predict(fit_m1.Y, newdata = new_data, type = "response")
plot( x=density(pred.val.Y.glm) )


###> SIMEX with NDB ------------------------------------------------------------
#matXS.test
#matXSstar.test
# SIMEX does not have Predict!!!SIMEX does not have Predict!!!so...no training!!!!!! go test.....so stupid
#####  SIMEX...need to use test set? SIMEX does not have ....."predict(.)" so..

# ---> start with "test data"
fit_m1.S.sim <- glm( N.test*Y.test ~ XSstar.test + factor(ZS.test), family=poisson(link="log") )
#fit_m1.Y.sim <- glm( Y.test ~ XSstar.test + factor(ZS.test), family=Gamma(link="log") )
fit_m1.Y.sim <- glm( Y.test ~ XSstar.test + factor(ZS.test), family=poisson(link="log") )
## > what is your best guess on "XSstar"?
err.S <- rep(1, nrow(test.df))

for (i in 1:nrow(test.df)) {
  if(test.df$Classes[i]==2){
    uuu <- ifelse( train.df$LnCoverage[i]<0 & test.df$LnCoverage[i]>5, 0.2, 0.8 )
    delta_sq <- ifelse( test.df$lnDeduct[i]<=6, 1.8, 
                        ifelse(test.df$lnDeduct[i]>10, 0.2, 0.8) )
    err.S[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err.S[i] <- 0.01
  }
};err.S

###> RUN SIMEX.test for S
mod_simex.S.test <- simex(fit_m1.S.sim,
                          SIMEXvariable = "XSstar.test",
                          measurement.error = abs(err.S),
                          fitting.method = "quadratic",
                          B = 10,
                          asymptotic = F)
par(mfrow=c(2,2))
plot(mod_simex.S.test)
summary(mod_simex.S.test)
pred.val.S.simex <- mod_simex.S.test$fitted.values

#E.claims.S.simex <- pred.val.S.simex

###> RUN SIMEX.test for Y
mod_simex.Y.test <- simex(fit_m1.Y.sim,
                          SIMEXvariable = "XSstar.test",
                          measurement.error = abs(err.S),
                          fitting.method = "quadratic",
                          B = 10,
                          asymptotic = F)
par(mfrow=c(2,2))
plot(mod_simex.Y.test)
summary(mod_simex.Y.test)
pred.val.Y.simex <- mod_simex.Y.test$fitted.values




# - start again
fit_m1.S.xx <- glm( N*Y ~ XSstar + factor(ZS), family=poisson(link = "log") )
## > what is your best guess on "XSstar"?
err.SS <- rep(1, nrow(train.df))

for (i in 1:nrow(train.df)) {
  if(train.df$Classes[i]==2){
    uuu <- ifelse( train.df$LnCoverage[i]<0 & train.df$LnCoverage[i]>5, 0.2, 0.8 )
    delta_sq <- ifelse( train.df$lnDeduct[i]<=6, 1.8,
                        ifelse(train.df$lnDeduct[i]>10, 0.2, 0.8) )
    err.SS[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
  }
  else{
    err.SS[i] <- 0.01
  }
};err.SS
### RUN SIMEX.train ------------------------------------------------------------
mod_simex.S.train <- simex(fit_m1.S.xx,
                           SIMEXvariable = "XSstar",
                           measurement.error = abs(err.SS),
                           fitting.method = "quadratic",
                           B = 10,
                           asymptotic = F)
par(mfrow=c(2,2))
plot(mod_simex.S.train)
summary(mod_simex.S.train)
summary(fit_m1.S.xx)
E.claims.S.simex <- mod_simex.S.train$fitted.values




###############################################################################
###############################################################################
###############################################################################
########################## S(t) with Gustafson ################################
###############################################################################
###############################################################################
###############################################################################


#####> Step 1: Define the parameters for each group
group_params.st <- list(
  list(xi_F=mean(EX_F.train[cl_membershipF==1]), psi_F=avg_psi.vec1, 
       meanlog_S=mean(mu.vec.train.st[cl_membershipF==1]), 
       sdlog_S=sqrt(avg_sig2.vec1.st), 
       n=table(cl_membershipF)[1]),
  list(xi_F=mean(EX_F.train[cl_membershipF==2]), psi_F=avg_psi.vec2, 
       meanlog_S=mean(mu.vec.train.st[cl_membershipF==2]), 
       sdlog_S=sqrt(avg_sig2.vec2.st), 
       n=table(cl_membershipF)[2]),
  list(xi_F=mean(EX_F.train[cl_membershipF==3]), psi_F=avg_psi.vec3, 
       meanlog_S=mean(mu.vec.train.st[cl_membershipF==3]), 
       sdlog_S=sqrt(avg_sig2.vec3.st), 
       n=table(cl_membershipF)[3]),
  list(xi_F=mean(EX_F.train[cl_membershipF==4]), psi_F=avg_psi.vec4, 
       meanlog_S=mean(mu.vec.train.st[cl_membershipF==4]), 
       sdlog_S=sqrt(avg_sig2.vec4.st), 
       n=table(cl_membershipF)[4]),
  list(xi_F=mean(EX_F.train[cl_membershipF==5]), psi_F=avg_psi.vec5, 
       meanlog_S=mean(mu.vec.train.st[cl_membershipF==5]), 
       sdlog_S=sqrt(avg_sig2.vec5.st), 
       n=table(cl_membershipF)[5]),
  list(xi_F=mean(EX_F.train[cl_membershipF==6]), psi_F=avg_psi.vec6, 
       meanlog_S=mean(mu.vec.train.st[cl_membershipF==6]), 
       sdlog_S=sqrt(avg_sig2.vec6.st), 
       n=table(cl_membershipF)[6])
); group_params



#####> Step 2: Simulate the total claim amount for each group
E.total_claims.st <- numeric(n_simulations)
E.grp_claims.st <- numeric(n_simulations)

# - without cluster
for (grp in 1:6) {
  params <- group_params.st[[grp]]
  claim_counts <- rnbinom(n=n_simulations, size=params$psi_F, mu=params$xi_F)
  
  for (i in 1:n_simulations) {
    if (claim_counts[i] > 0) {
      severities <- rlnorm(n=claim_counts[i], meanlog=params$meanlog_S, sdlog=params$sdlog_S)
      E.total_claims.st[i] <- E.total_claims[i] + sum(severities)
    }
  }
}
summary(E.total_claims.st)
summary(log(E.total_claims.st))


# - with clusters
par(mfrow=c(2,3))
for (grp in 1:6) {
  params <- group_params.st[[grp]]
  claim_counts <- rnbinom(n=n_simulations, size=params$psi_F, mu=params$xi_F)
  
  for (i in 1:n_simulations) {
    if (claim_counts[i] > 0) {
      severities <- rlnorm(claim_counts[i], meanlog=params$meanlog_S, sdlog=params$sdlog_S)
      E.grp_claims.st[i] <- sum(severities)
    }
  }
  #> summary
  cat("\nSummary for Group", grp, "\n")
  print(summary(E.grp_claims.st))
  
  #> plot
  #hist(log(E.grp_claims.st), breaks=1000, prob=T, xlab=paste("S(t) for Group", grp), ylab="Frequency", xlim=c(0,300000))
  #lines( density(x=log(E.grp_claims.st)), col="blue", lwd = 1 )
  plot( x=density(x=log(E.grp_claims)), xlab=paste("log(S(t)) for Group", grp), 
        ylab="Density", col="blue", lwd=2  )
}


# Step 3: Show the mixture results ???? ------------------------------------------------ # 

# Calculate total sample size
#total_n <- sum(sapply(group_params, function(p) p$n))

# Simulate data points for each group
data.st <- numeric(0)

for (p.st in group_params.st) {
  group_data.st <- rlnorm(n=p.st$n, meanlog=p.st$meanlog_S, sdlog=p.st$sdlog_S)
  data.st <- c(data.st, group_data.st)
}

# Plot the density of the mixture distribution
par(mfrow=c(1,1))
plot( density(log(data.st)), xlab="S(t)", ylab="Density", ylim=c(0,0.25) )

# Add individual densities for each group (scaled by their mixing coefficient)
colors <- rainbow(6)
for (i in 1:6) {
  group_data.st <- rlnorm(group_params.st[[i]]$n, 
                          meanlog = group_params.st[[i]]$meanlog_S, 
                          sdlog = group_params.st[[i]]$sdlog_S)
  lines( density(log(group_data.st)), col=colors[i], lwd=1, lty=1) #, xlim=c(1:50000))
}
legend("topright", legend=paste("Group", 1:6), lwd=2, lty=1, col=colors )





# ------------------------------------------------------------------------------
# ---- COMPARE all together !!!! -----------------------------------------------
# ------------------------------------------------------------------------------

###> LPPD for Bayes model (for Y)
lppd_EX_S.Y   #-16155.9
lppd_EX_S.Y.xx#-17731.03
lppd_EX_S.Y.st#-16470.44



###> SSE (for Y)
SSPE.LN <- sum( (log(EX_S.Y.test) - log(test.df$yAvg))^2 ); SSPE.LN   #784.5236
SAPE.LN <- sum( abs(log(EX_S.Y.test) - log(test.df$yAvg)) ); SAPE.LN  #415.2161 
SSPE.glm <- sum( (log(pred.val.Y.glm) - log(test.df$yAvg))^2 ); SSPE.glm   #1269.91
SAPE.glm <- sum( abs(log(pred.val.Y.glm) - log(test.df$yAvg)) ); SAPE.glm  #594.6237 

SSPE.LN.xx <- sum( (log(EX_S.test.Y.xx) - log(test.df$yAvg))^2 ); SSPE.LN.xx   #840.0246
SAPE.LN.xx <- sum( abs(log(EX_S.test.Y.xx) - log(test.df$yAvg)) ); SAPE.LN.xx  #430.2199 
SSPE.LN.st <- sum( (log(EX_S.test.Y.st) - log(test.df$yAvg))^2 ); SSPE.LN.st   #839.8009
SAPE.LN.st <- sum( abs(log(EX_S.test.Y.st) - log(test.df$yAvg)) ); SAPE.LN.st  #427.6068
SSPE.simex <- sum( (log(pred.val.Y.simex) - log(test.df$y))^2 ); SSPE.simex   #1546.584
SAPE.simex <- sum( abs(log(pred.val.Y.simex) - log(test.df$y)) ); SAPE.simex  #635.0063



###> CTE (for S)
summary( E.total_claims )
summary( E.claims.S.glm )

summary( E.total_claims.xx )
summary( E.total_claims.st )
summary( E.claims.S.simex )

# - Calculate the quantile threshold
alalpha = 0.95
alalpha = 0.9
alalpha = 0.5
alalpha = 0.1
quantile_threshold <- quantile(E.total_claims, alalpha); quantile_threshold
quantile_threshold.glm <- quantile(E.claims.S.glm, alalpha); quantile_threshold.glm
quantile_threshold.xx <- quantile(E.total_claims.xx, alalpha); quantile_threshold.xx
quantile_threshold.st <- quantile(E.total_claims.st, alalpha); quantile_threshold.st
quantile_threshold.simex <- quantile(E.claims.S.simex, alalpha); quantile_threshold.simex
# - Filter the data to include only the values in the tail
tail_data <- E.total_claims[E.total_claims > quantile_threshold]; tail_data # NA?????
tail_data.glm <- E.claims.S.glm[E.claims.S.glm > quantile_threshold.glm]; tail_data.glm # NA?????
tail_data.xx <- E.total_claims.xx[E.total_claims.xx > quantile_threshold.xx]; tail_data.xx # NA?????
tail_data.st <- E.total_claims.st[E.total_claims.st > quantile_threshold.st]; tail_data.st # NA?????
tail_data.simex <- E.claims.S.simex[E.claims.S.simex > quantile_threshold.simex]; tail_data.simex # NA?????

CTE <- mean(tail_data, na.rm = TRUE); CTE                    # 27499631 (0.95), 20976138 (0.90), 8159358 (0.5)  4878240 (0.1), 
CTE.glm <- mean(tail_data.glm, na.rm = TRUE); CTE.glm        # meaningless....  

CTE.xx <- mean(tail_data.xx, na.rm = TRUE); CTE.xx           # 16985933 (0.95), 12659919 (0.90), 4910311 (0.5)  2959975 (0.1), 
CTE.st <- mean(tail_data.st, na.rm = TRUE); CTE.st           # 28066650 (0.95), 21627316 (0.90), 8721810 (0.5)  5299564 (0.1),
CTE.simex <- mean(tail_data.simex, na.rm = TRUE); CTE.simex  # 320487.2 (0.95), 231342.2 (0.90), 83736.26 (0.5)  59237.36 (0.1),
#-------------------------------------------------------------------------------

# [1] 91930.49
# > CTE.glm <- mean(tail_data.glm, na.rm = TRUE); CTE.glm        # 
# [1] 96329.76
# > CTE.xx <- mean(tail_data.xx, na.rm = TRUE); CTE.xx           # 53319.54, 52937.02, 51968.08, 49727.14 
# [1] 51968.08
# > CTE.st <- mean(tail_data.st, na.rm = TRUE); CTE.st           # 85482.6, 80047.48, 70590.93, 62906.14
# [1] 70590.93
# > CTE.simex <- mean(tail_data.simex, na.rm = TRUE); CTE.simex  # 53038.49, 53015.64, 54451.99, 53038.49
# [1] 53038.49








#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Versus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#>>> Similarity?
# (1) Kolmogorov-Smirnov Test (KS Test): 
# This test compares the cumulative distributions of two datasets and is useful for determining 
# if they differ significantly.
ks.test(E.total_claims, E.total_claims.st) # D = 0.068966
ks.test(E.total_claims, E.total_claims.xx) # D = 0.11991
ks.test(E.total_claims, E.claims.S.simex) # D = 0.83542
# D: The maximum distance between the empirical cumulative distribution functions (ECDFs) 
# of the two samples. A smaller value indicates that the distributions are more similar.
# p-value: Indicates whether the difference between the distributions is statistically significant.
# A low p-value (typically < 0.05) suggests that the two distributions are significantly different.





# (2) Kullback-Leibler (KL) Divergence:
# This measures how one probability distribution diverges from a second, expected probability 
# distribution. For sample data, this requires estimating the probability density functions.
library(entropy)

##> Estimate the probability density functions
##> 
# #- with aggregate claims ?
# pdf <- density(E.total_claims)$y
# pdf.xx <- density(E.total_claims.xx)$y
# pdf.st <- density(E.total_claims.st)$y
# pdf.simex <- density(E.claims.S.simex)$y

# # - with test claims ?
pdf <- density(EX_S.Y.test)$y
pdf.xx <- density(EX_S.test.Y.xx)$y
pdf.st <- density(EX_S.test.Y.st)$y
pdf.simex <- density(pred.val.Y.simex)$y

# # - with train claims ?
# pdf <- density(EX_S.train)$y
# pdf.xx <- density(EX_S.train.dirty)$y
# pdf.st <- density(EX_S.train.st)$y
# pdf.simex <- density(E.claims.S.simex)$y


# Calculate the KL divergence
kl.divergence.xx <- entropy::KL.empirical(pdf, pdf.xx)
print(kl.divergence.xx) # 0.6979195
kl.divergence.st <- entropy::KL.empirical(pdf, pdf.st)
print(kl.divergence.st) #0.1938633
kl.divergence.simex <- entropy::KL.empirical(pdf, pdf.simex)
print(kl.divergence.simex) #1.619171



























