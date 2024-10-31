############ Bayes Parametric: Hierarchical LN-GLM [BASIC] with NDB ############

library(mvnfast)
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
#library(extraDistr)
library(MCMCpack)
library('dplyr')


full.LGPIF <- read.csv("C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/PropertyFundInsample_original.csv")

str(full.LGPIF)


full.LGPIF <- full.LGPIF[full.LGPIF$y>0, ]

head(full.LGPIF)
str(full.LGPIF) # 1679 obs. of  23 variables:

full.LGPIF <- full.LGPIF %>% mutate(Classes = case_when(
  TypeCity==1 & TypeCounty==0 & TypeMisc==0 & TypeSchool==0 & TypeTown==0 & TypeVillage==0 ~ "city",
  TypeCity==0 & TypeCounty==1 & TypeMisc==0 & TypeSchool==0 & TypeTown==0 & TypeVillage==0 ~ "county",
  TypeCity==0 & TypeCounty==0 & TypeMisc==1 & TypeSchool==0 & TypeTown==0 & TypeVillage==0 ~ "misc",
  TypeCity==0 & TypeCounty==0 & TypeMisc==0 & TypeSchool==1 & TypeTown==0 & TypeVillage==0 ~ "school",
  TypeCity==0 & TypeCounty==0 & TypeMisc==0 & TypeSchool==0 & TypeTown==1 & TypeVillage==0 ~ "town",
  TypeCity==0 & TypeCounty==0 & TypeMisc==0 & TypeSchool==0 & TypeTown==0 & TypeVillage==1 ~ "village",
  TRUE ~ "none"
)) %>%
  mutate(Classes = as.factor(Classes))

head(full.LGPIF)


#N:   Freq
#Y:   yAvg
#X_F: { z:AC15, x:Premium } 
#X_S: { z:Fire5 , x:lnDeduct } 

#Extra:  Year, y, TypeCounty, LnCoverage
df = dplyr::select(.data=full.LGPIF, Year, y, LnCoverage, AC15, Premium, Fire5, lnDeduct, 
                   Freq, yAvg, Classes); head(df)

# Add NDB error : 1%, 10%, 40%

quantile(df$lnDeduct, c(0.025, 0.975))
#  2.5%  97.5% 
#  2.5%  97.5% 



#-------------------------------------------------------------------------- [1%]
set.seed(1)
df.n <- length(df$yAvg)
err.S <- numeric(df.n)
n.breaks <- sqrt( nrow(df) ) #******Rule of thumb

# for (i in 1:df.n) {
#   if(df$Classes[i]=="county"){
#     # From unknown factor
#     uuu <- ifelse( df$LnCoverage[i]<0 & df$LnCoverage[i]>5, 0.2, 0.8 )
#     
#     # From clean covariate
#     delta_sq <- ifelse( df$lnDeduct[i]<=6, 1.8, 
#                         ifelse(df$lnDeduct[i]>10, 0.2, 0.8) )
#     
#     err.S[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
#   }
#   else{
#     err.S[i] <- 0
#   }
# }

for (i in 1:df.n) {
  if(df$Fire5[i]==1){
    rho2 <- ifelse( df$Classes[i]=="city", 0.91,
                    ifelse( df$Classes[i]=="county", 0.99,
                            ifelse( df$Classes[i]=="misc", 0.89,
                                    ifelse( df$Classes[i]=="school", 0.98,
                                            ifelse( df$Classes[i]=="town", 0.92,
                                                    0.95 ) ) ) ) ) 
    vvv <- var(df$lnDeduct)

    err.S[i] <- rnorm( n=1, mean=0, sd=sqrt( vvv/rho2 - vvv ) )
  }
  else{
    err.S[i] <- 0
  }
}
hist(err.S,  breaks=n.breaks)


df$Error_S <- err.S
df$lnDeduct_err <- df$lnDeduct + df$Error_S 
head(df)

overall_err_S = sum(abs(df$Error_S)) /  sum(abs(df$lnDeduct)); overall_err_S ######### error: 1%

write.csv(df, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/LGPIF.paper3-1+.csv", row.names=FALSE)


#-------------------------------------------------------------------------- [10%]
set.seed(1)
df.n <- length(df$yAvg)
err.S <- numeric(df.n)
n.breaks <- sqrt( nrow(df) ) #******Rule of thumb

# for (i in 1:df.n) {
#   if(df$Classes[i]=="county"){
#     # From unknown factor
#     uuu <- ifelse( df$LnCoverage[i]<0 & df$LnCoverage[i]>5, 2.2, 8.5 )
#     
#     # From clean covariate
#     delta_sq <- ifelse( df$lnDeduct[i]<=6, 2.2, 
#                         ifelse(df$lnDeduct[i]>10, 8.5, 6) )
#     
#     err.S[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
#   }
#   else{
#     err.S[i] <- 0
#   }
# }
for (i in 1:df.n) {
  if(df$Fire5[i]==1){
    rho2 <- ifelse( df$Classes[i]=="city", 0.7,
                    ifelse( df$Classes[i]=="county", 0.2,
                            ifelse( df$Classes[i]=="misc", 0.6,
                                    ifelse( df$Classes[i]=="school", 0.4,
                                            ifelse( df$Classes[i]=="town", 0.8,
                                                    0.3 ) ) ) ) ) 
    vvv <- var(df$lnDeduct)
    
    err.S[i] <- rnorm( n=1, mean=0, sd=sqrt( vvv/rho2 - vvv ) )
  }
  else{
    err.S[i] <- 0
  }
}
hist(err.S,  breaks=n.breaks)


df$Error_S <- err.S
df$lnDeduct_err <- df$lnDeduct + df$Error_S 
tail(df)


overall_err_S = sum(abs(df$Error_S)) /  sum(abs(df$lnDeduct)); overall_err_S ######### error: 10%

write.csv(df, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/LGPIF.paper3-10+.csv", row.names=FALSE)


#-------------------------------------------------------------------------- [40%]
set.seed(1)
df.n <- length(df$yAvg)
err.S <- numeric(df.n)
n.breaks <- sqrt( nrow(df) ) #******Rule of thumb

# for (i in 1:df.n) {
#   if(df$Classes[i]=="county"){
#     # From unknown factor
#     uuu <- ifelse( df$LnCoverage[i]<0 & df$LnCoverage[i]>5, 8.5, 55 )
#     
#     # From clean covariate
#     delta_sq <- ifelse( df$lnDeduct[i]<=6, 8.5, 
#                         ifelse(df$lnDeduct[i]>10, 55, 10) )
#     
#     err.S[i] <- rnorm( n=1, mean=0, sd=sqrt( uuu*delta_sq ) )
#   }
#   else{
#     err.S[i] <- 0
#   }
# }
for (i in 1:df.n) {
  if(df$Fire5[i]==1){
    rho2 <- ifelse( df$Classes[i]=="city", 0.07,
                    ifelse( df$Classes[i]=="county", 0.05,
                            ifelse( df$Classes[i]=="misc", 0.06,
                                    ifelse( df$Classes[i]=="school", 0.04,
                                            ifelse( df$Classes[i]=="town", 0.08,
                                                    0.02 ) ) ) ) ) 
    vvv <- var(df$lnDeduct)
    
    err.S[i] <- rnorm( n=1, mean=0, sd=sqrt( vvv/rho2 - vvv ) )
  }
  else{
    err.S[i] <- 0
  }
}
hist(err.S,  breaks=n.breaks)

df$Error_S <- err.S
df$lnDeduct_err <- df$lnDeduct + df$Error_S 
tail(df)

overall_err_S = sum(abs(df$Error_S)) /  sum(abs(df$lnDeduct)); overall_err_S ######### error: 40%

write.csv(df, "C:/Users/kimm4/Desktop/WORKSPACE/2nd 3rd paper/DATASET/LGPIF.paper3-40+.csv", row.names=FALSE)



