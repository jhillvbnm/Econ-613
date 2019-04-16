# Jeff Hill, ECON 613 Assignment 4
library(plyr)
library(dplyr)
library(reshape2)
library(fastDummies)
library(nlme)
library(plm)
library(nnet)
# set seed
set.seed(613)

# load data
datapath <- "/Users/admin/Documents/Econ_613/Data/Assignment 4/"
koop <- read.csv(file=paste(datapath,"Koop-Tobias.csv",sep=""), header=TRUE, sep=",")

#************************************************
# Exercise 1 Data 
#************************************************
# randomly select 5 of the 2178 individuals
sample <- subset(koop[,c(1,3,5)], subset = PERSONID %in% c(sample(1:2178,5,replace=F)))
# with seed = 613, the individuals should be: 226, 707, 1790, 2081, and 2132 

# since our sample df is in long form, convert it to wide form for clarity, and rename it wagepanel
wagepanel <- dcast(sample, PERSONID ~ TIMETRND, value.var="LOGWAGE")

# wagepanel is missing a column for time period 1, as no individuals have data for that period.
# add in that column and then rename rows and columns
wagepanel <- wagepanel[,c(2,1,3:15)]
colnames(wagepanel) <- c('t0','t1','t2','t3','t4','t5','t6','t7','t8','t9','t10','t11','t12','t13','t14')
rownames(wagepanel) <- c(226,707,1790,2081,2132)
wagepanel$t1 <- c(NA,NA,NA,NA,NA)
# wagepanel now contains the panel dimension of wages for 5 randomly selected individuals.
wagepanel



#************************************************
# Exercise 2 Random Effects 
#************************************************

# Random effects model using gls
REfit <- gls(LOGWAGE ~ EDUC + POTEXPER, data=koop)

# notice that the results from this regression yield the same coefficients as a linear OLS
# however the standard errors may differ.


#************************************************
# Exercise 3 Fixed Effects Model 
#************************************************


# Between Estimator------------------------------

# turn our koop-tobias data from long to wide for the sake of the between estimator but also because
# wide data is easier to manipulate in certain regards
koopwide <- recast(koop, PERSONID ~ TIMETRND + variable, measure.var = c("LOGWAGE", "POTEXPER","EDUC"))
# reorder columns
koopwide <- koopwide[,c(1,2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,
                        3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,
                        4,7,10,13,16,19,22,25,28,31,34,37,40,43,46)]

# define which columns are LOGWAGE, POTEXPER, and EDUC for convenience:
lwagecol <- c(2:16)
pexpcol <- c(17:31)
educcol <- c(32:46)
# generate averages for logwage and potential experience across time for each individual.
koopwide$lwageavg <- rowMeans(koopwide[,lwagecol], na.rm = TRUE)
koopwide$pexpavg <- rowMeans(koopwide[,pexpcol], na.rm = TRUE)
koopwide$educavg <- rowMeans(koopwide[,educcol], na.rm = TRUE)

# Now with the data in this wide format, we can just run a linear model on these average values to
# calculate the between estimator
betweenfit <- lm(lwageavg ~ educavg + pexpavg, koopwide)
betweenfit$coefficients


# Within Estimator------------------------------

# Now we will need to merge our averages for lwage, potexper, and educ we calculated for the between 
# estimator back into our long dataframe (namely koop), since we need data for each individual in 
# each time period for the within estimator
koop <- merge(koop,koopwide[,c("PERSONID","lwageavg","pexpavg","educavg")], by.x = "PERSONID", by.y = "PERSONID")

# now that we have the average lwage, education, and potential experience for each individual, we need to generate 
# the differences from the means for each individual in each time period
koop$lwagewithin <- koop$LOGWAGE - koop$lwageavg
koop$pexpwithin <- koop$POTEXPER - koop$pexpavg
koop$educwithin <- koop$EDUC - koop$educavg

# run the model
withinfit <- lm(lwagewithin ~ educwithin + pexpwithin -1, koop)
withinfit$coefficients


# First Difference Estimator --------------------

# I subset the columns of interest in koops into a dataframe called firstd which i will manipulate for the first
# difference estimator. Also of note is that the time periods missing are not systematic. Looking at individual 1
# for example, they have data for time periods 0,7,9, and 10. For our first difference model, we only consider 
# adjacent time periods, and we accomplish that using trend_diff below

firstd <- koop[,c("PERSONID","LOGWAGE","EDUC","POTEXPER","TIMETRND")]
# create lagged columns for lwage, educ, and potexper that yield NAs for the first time period of each individual
# we then take the difference between periods and place the difference values in columns lwagelag, educlag, etc.
# this is all done in the commands below
firstd$lwagelag <- unlist(by(firstd$LOGWAGE , list(firstd$PERSONID) , function(i) c(NA,diff(i))))
firstd$educlag <- unlist(by(firstd$EDUC , list(firstd$PERSONID) , function(i) c(NA,diff(i))))
firstd$pexplag <- unlist(by(firstd$POTEXPER , list(firstd$PERSONID) , function(i) c(NA,diff(i))))
firstd$trend_diff <- unlist(by(firstd$TIMETRND, list(firstd$PERSONID), function(i) c(NA,diff(i))))
# this last column now calculates the difference between time periods. we only want to regress on data where
# the time difference is 1.


# run the model
firstdfit <- lm(lwagelag ~ educlag + pexplag -1, subset(firstd, trend_diff == 1))
firstdfit$coefficients


# Create a table for comparison of coefficients ------------

comp_df <- data.frame(matrix(NA, ncol = 3, nrow = 3))
colnames(comp_df) <- c("Between","Within","First Difference")
rownames(comp_df) <- c("Intercept","Education","Experience")

comp_df[,1] <- betweenfit$coefficients
comp_df[,2] <- c(NA,withinfit$coefficients)
comp_df[,3] <- c(NA,firstdfit$coefficients)

# comp_df now is a dataframe comparing the coefficents of the 3 estimators
comp_df

#************************************************
# Exercise 4 Understanding Fixed Effects
#************************************************
# subset our koop dataframe to randomly select 100 individuals
FEkoop <- subset(koop, subset = PERSONID %in% c(sample(1:2178,100,replace=F)))

# betas contain our two betas for education and experience, and then 100 alphas.

# subset our x data to just pull the columns we want
FEx <- FEkoop[,c(2,4,1)]
FEy <- FEkoop[,3]

test <- as.matrix(FEx[,c(1,2)])%*%c(3,3)

FE_neg_ll <- function (param,x,y) {
  # x will contain education and experience, as well as idividual id since we need to match the alphas to 
  # the individuals
  # parameters will be vector of length 103. the first 100 will be the 100 alphas, parameters 101 and 102 will 
  # be the two betas for educ and exper, and parameter 103 will be the sigma
  sigma <- param[103]
  alpha <- param[1:100]
  
  # i need to map the alphas for each individual to the observations for each individual. I do that below
  mapdf <- data.frame(old=c(x[with(x, c(PERSONID[-1]!= PERSONID[-nrow(x)], TRUE)),3]),new=param[1:100]) # this creates a map between unique values of person ID
  # (which we have 100 of) and the 100 alphas (just 1 to 100)
  alphalong <- mapdf$new[match(FEx[,3],mapdf$old)] # this maps the alphas to the person IDs, and creates the 
  # alpha vector which is of correct length, namely 798, which matches the dim of our X matrix.
  
  cut_x <- x[,c(1,2)] # dropping the person ID column of x since it was only used for matching the alphas
  xb <- as.matrix(cut_x)%*%param[101:102] + alphalong # create XB
  u <- y - xb # create u, the difference between y and yhat
  -sum(log(dnorm(u,sd = sigma))) # sum the log of each individual likelihood, then make negative
}

FEfit <- nlm(FE_neg_ll,rep(1,103),x=FEx,y=FEy)

# this optimization is very sensitive to starting values. with starting values that are 0 or negative, the optimization
# yields incorrect results (significantly higher likelihood). with starting values above 5, the optimization yields 
# incorrect results (significiantly higher likelihood). With starting values between .1 and 5, the optimization yields the 
# lowest likelihood (lik minimum = 145.3604). Since the estimates are consistent across this range, and higher or lower 
# starting values we get less optimal likelihood, I will use the results from this range as my estimates. also note
# that I varied the 100 alphas, 2 betas, and the standard deviation separately.



alphas <- FEfit$estimate[1:100]
betas <- FEfit$estimate[101:102]
sigma <- FEfit$estimate[103]

# parameters reported below:
alphas
betas
sigma

# subset our FE data to get the time-invariant data for our 100 individuals, and add alphas to the data.
FE_time_inv <- unique(FEkoop[,c(6:10)])
FE_time_inv$alpha <- alphas

# here we regress the alphas on the time invariant variables.
time_inv_fit <- lm(alpha ~ ABILITY + MOTHERED + FATHERED + BRKNHOME + SIBLINGS, FE_time_inv)


# The standard errors are potentially not correctly estimated, because the alphas are generated from our optimization,
# and so they carry some amount of error in their estimation. But we then treat them as our endogenous variable in the regression,
# not accounting for the error in their values. so to correct for these standard errors, we bootstrap.



# FEkoop going in

fe_bootstrap <- function(data, reps) {
  beta_df <- data.frame(matrix(NA, nrow = reps, ncol = 6)) # here we create the empty matrix to temporarily store betas in
  colnames(beta_df) <- c("Intercept", "Ability","MotherEd","FatherEd","Broken Home","Siblings") # name the columns
  for (i in 1:reps) {          
    ind_list <- sample <-sample(unique(data$PERSONID), 100, replace = T) # sample a list of the 100 individuals selected with replacement
    boot_sample <- data.frame() # boot_sample will be FEkoop resampled with replacement of individuals each iteration
    for(j in ind_list) {
      boot_sample <- rbind(boot_sample,data[which(data$PERSONID==j),])
    }
    boot_ss <- boot_sample[,1:4] # take only the columns of interest from boot_sample
    boot_ss <- with(boot_ss, data.frame(class.ind(PERSONID), EDUC, POTEXPER, LOGWAGE)) # generates all the dummy variables for the model
    boot_fit <- lm(LOGWAGE~.,boot_ss) # run the FE model with dummies using OLS
    # using ols here vs the log-lik function for speed and becuase of the potential inconsistency of the Log-lik with different 
    # starting parameters.
    
    alphas <- as.vector(boot_fit$coefficients[2:(length(boot_fit$coefficients)-2)]) # pull out all the alphas
    boot_time_inv <- unique(boot_sample[,c(6:10)]) # gather our time-invariant variables still in boot_sample
    boot_time_inv$alpha <- alphas 
    boot_inv_fit <- lm(alpha ~ ABILITY + MOTHERED + FATHERED + BRKNHOME + SIBLINGS, boot_time_inv) # run alpha ~ time invariants
    beta_df[i,] <- t(as.vector(boot_inv_fit$coefficients)) # store the betas
  }  
  apply(beta_df, 2, sd, na.rm = T)
}
fe_bootstrap(FEkoop,199)
# the end result is 5 standard errors (6 if you include constant) for our 5 time invariant variables on alpha generated through bootstrap.
