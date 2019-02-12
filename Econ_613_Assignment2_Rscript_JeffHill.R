# Jeff Hill, ECON 613 Assignment 2

library(stringr)
library(dplyr) 
library(StatMeasures)
library(ggplot2)
library(Hmisc)
library(boot)
library(numDeriv)
#*******************************************************
# Exercise 1 Data Creation
set.seed(613)

# generate variables pulling from different distributions
# uniform [1:3]
x1 <- runif(10000, min = 1, max = 3)
# gamma shape:3 scale: 2
x2 <- rgamma(10000, 3, scale = 2)
# binomial
x3 <- rbinom(10000, 1, 0.3)
# normal mean 2
eps <- rnorm(10000, mean = 2, sd = 1)

# creating y and ydum
y <- 0.5 + 1.2*x1 -0.9*(x2) + 0.1*x3 +eps
ydum <- as.numeric(y > mean(y))


#******************************************************
# Exercise 2 OLS
# correlation between y and X1
rcorr(y, x1)
# the correlation between y and x1 is 0.22, which is different from 1.2 by 0.98. since correlation
# is bound between -1 and 1, it would be problematic if we got correlation close to 1.2 (above 1)

# Regression of Y on x = [1,x1,x2,x3]
# creating X matrix
intercept <- rep(1, 10000)
x <- matrix(c(intercept,x1,x2,x3), nrow = 10000, ncol = 4, byrow = FALSE)

#solving for OLS betas using B = (X'X)^(-1)X'Y
xpxinv <- solve(t(x)%*%x)
betas <- xpxinv%*%t(x)%*%y
# the betas are:
betas[1]
# intercept: 2.466498
betas[2]
# beta for x1: 1.234452
betas[3]
# beta for x2: -.9058537
betas[4]
# beta for x3: 0.1240774
# these values make sense as they match quite closely to the expected values for the betas from the formula for y.
# the only notable difference is the intercept term which is 2 larger than expected (2.5 vs. 0.5). This is due to the error term,
# eps, having a mean of 2, which is captured by the constant term.

# calculating standard errors
# using standard OLS formula: variance(Betas) = (X'X)^(-1)*sigma^2 where sigma^2 is the variance of eps
sigmasq <- var(eps)
varbetas <- sigmasq*xpxinv
seOLS <- sqrt(diag(varbetas))
# the standard error of the betas are:
seOLS[1]
# for intercept: 0.04094909
seOLS[2]
# for beta_x1: 0.01739891
seOLS[3]
# for beta_x2: 0.002882897
seOLS[4]
# for beta_x3: 0.02200445

# Now using bootstrap with 49 and 499 replications respectively.
# create data frame suitable for bootstrapping
dat <- data.frame(y,intercept,x1,x2,x3,eps,ydum)

#replications = 49
# create empty dataframe to put beta values from for loop into
betadf49 <- data.frame(matrix(NA, nrow = 49, ncol = 4))
colnames(betadf49) <- c("intercept","beta_x1","beta_x2","beta_x3")

# bootstrap in for loop, sampling y 10000 times with replacement, then calculating betas. store those betas in betadf
for (i in 1:49) {
  bootd <- dat[sample(nrow(dat), 10000, replace = T), ]
  booty <- as.matrix(bootd[,1])
  bootx <- as.matrix(bootd[,2:5])
  bootinv <- as.matrix(solve(t(bootx)%*%bootx))
  bootbeta <- as.matrix(bootinv%*%t(bootx)%*%booty)
  betadf49[i,] <- t(bootbeta)
}

# calculate standard errors of these betas
se49 <- apply(betadf49, 2, sd)
# the standard error of the betas from a trial using bootstrap 49 replications are:
se49[1]
# se for intercept: 0.03701252    ***     Compared to OLS formula SE's      ***     for intercept: 0.04094909
se49[2]
# se for beta_x1:   0.01581407                                                      for beta_x1: 0.01739891
se49[3]
# se for beta_x2:   0.002948432                                                     for beta_x2: 0.002882897
se49[4]
# se for beta_x3:   0.01596766                                                      for beta_x3: 0.02200445


# Replications = 499
# create another empty dataframe to put beta values from for loop into
betadf499 <- data.frame(matrix(NA, nrow = 499, ncol = 4))
colnames(betadf499) <- c("intercept","beta_x1","beta_x2","beta_x3")
# bootstrap in for loop, sampling y 10000 times with replacement, then calculating betas. store those betas in betadf
for (i in 1:499) {
  bootd <- dat[sample(nrow(dat), 10000, replace = T), ]
  booty <- as.matrix(bootd[,1])
  bootx <- as.matrix(bootd[,2:5])
  bootinv <- as.matrix(solve(t(bootx)%*%bootx))
  bootbeta <- as.matrix(bootinv%*%t(bootx)%*%booty)
  betadf499[i,] <- t(bootbeta)
}
# calculate the standard errors of the betas
se499 <- apply(betadf499, 2, sd)
# the standard error of the betas using bootstrap 499 replications are:
se499[1]
# se for intercept: 0.040901293    ***     Compared to OLS formula SE's    ***     for intercept: 0.04094909
se499[2]
# se for beta_x1:   0.017351702                                                    for beta_x1: 0.01739891
se499[3]   
# se for beta_x2:   0.002780576                                                    for beta_x2: 0.002882897
se499[4]
# se for beta_x3:   0.022487883                                                    for beta_x3: 0.02200445



#************************************************* 
# Exercise 3 Numerical Optimization
# consider the probit estimation of ydum on X
# we begin by writing a function that returns the likelihood of the probit
# we want to be sure to minimize the negative log-likelihood, instead of maximizing the log-likelihood
# this function will take inputs betas (named beta), since those are what we want to optimize ultimately.
# input beta will be a vector of length 4
pro_neg_ll <- function (pro_beta) {
  pro_x <- as.matrix(dat[,2:5]) #brings in data matrix X
  pro_ydum <- as.matrix(dat[,7]) # bring in ydum
  xb <- pro_x%*%pro_beta # create XB
  p <- pnorm(xb) # apply F() to XB, where F() is cumulative standard normal distribution function
  -sum( pro_ydum*log(p) + (1 - pro_ydum)*log(1 - p) ) # sum the log of each individual likelihood, then make negative
}
# test pro_neg_ll on some betas
pro_neg_ll(betas)
pro_neg_ll(c(2.4,1.23,-0.9,.1))

# now impliment steepest ascent optimization to maximize this function
# begin by setting intial values for beta vector, alpha, old and new likelihood.
beta <- c(0,0,0,0)
e <- 10^-10 # the amount we are incrementing our 4 betas by
e_mat <- matrix(c(e,0,0,0,0,e,0,0,0,0,e,0,0,0,0,e),nrow = 4, ncol = 4) # the matrix we add to mat1 in the while loop
# that increments the betas
alpha = .00001
likold = 1
liknew = 0

while(abs(liknew - likold) > 10^-6) { # this while loop continues to iterate as long as liknew is different enough from
  # lik old. 
  likold <- liknew # store the old likelihood from the last iteration for use in calculation the partial derivatives
  mat1 <- matrix(c(beta,beta,beta,beta), nrow = 4, ncol = 4, byrow = F) # create a 4x4 matrix of betas 1 through 4
  mat2 <- mat1 + e_mat # create a second matrix that is equal to mat1 with "e" added to positions 11, 22, 33, and 44
  d1 <- (pro_neg_ll(mat2[,1]) - pro_neg_ll(mat1[,1]))/e # here we define each direction as the difference in likelihood divided by e
  d2 <- (pro_neg_ll(mat2[,2]) - pro_neg_ll(mat1[,2]))/e
  d3 <- (pro_neg_ll(mat2[,3]) - pro_neg_ll(mat1[,3]))/e
  d4 <- (pro_neg_ll(mat2[,4]) - pro_neg_ll(mat1[,4]))/e
  d <- c(d1,d2,d3,d4) # collect all the directions into a vector
  beta <- beta - alpha*d # construct our new beta = the old beta minus alpha times our direction vector
  liknew <- pro_neg_ll(beta) # compute the likelihood of the new beta, which will be compared to the old likelihood at the
    # beginning of the next iteration.
}
print(beta) # optimal beta values resulting from this optimization are: 
# intercept: 3.02055601 beta_x1: 1.16047017 beta_x2: -0.88794791 beta_x3: 0.04261497
# this while loop was tested at intial beta = c(0,0,0,0), c(-1,-1,-1,-1), c(1,1,1,1), c(2,2,2,2)
# and it yielded the same optimal betas to within 1*10^-5

# testing again at alpha = .000001
beta <- c(0,0,0,0)
e <- 10^-10 # the amount we are incrementing our 4 betas by
e_mat <- matrix(c(e,0,0,0,0,e,0,0,0,0,e,0,0,0,0,e),nrow = 4, ncol = 4)
alpha = .000001
likold = 1
liknew = 0
while(abs(liknew - likold) > 10^-6) { 
  likold <- liknew # store the old likelihood from the last iteration for use in calculation the partial derivatives
  mat1 <- matrix(c(beta,beta,beta,beta), nrow = 4, ncol = 4, byrow = F) # create a 4x4 matrix of betas 1 through 4
  mat2 <- mat1 + e_mat # create a second matrix that is equal to mat1 with "e" added to positions 11, 22, 33, and 44
  d1 <- (pro_neg_ll(mat2[,1]) - pro_neg_ll(mat1[,1]))/e
  d2 <- (pro_neg_ll(mat2[,2]) - pro_neg_ll(mat1[,2]))/e
  d3 <- (pro_neg_ll(mat2[,3]) - pro_neg_ll(mat1[,3]))/e
  d4 <- (pro_neg_ll(mat2[,4]) - pro_neg_ll(mat1[,4]))/e
  d <- c(d1,d2,d3,d4)
  beta <- beta - alpha*d
  liknew <- pro_neg_ll(beta)
}
beta # optimal beta values resulting from this optimization at alpha = .000001 are: 
# intercept: 3.02055383 beta_x1: 1.16046993 beta_x2: -0.88794748 beta_x3: 0.04261515
# which as you can see, are very close to the betas gathered when alpha = .00001. However
# as expected, this while loop took SIGNIFICANTLY longer since we were moving in much shorter
# increments.

# the coefficients for exercise 3 are quite different from the coefficients in exercise 2. Exercise 2 ran an OLS of Y on X
# where in exercise 3 we run a probit of ydum on X. We would not expect the coefficients to be identical.


#****************************************
# exercise 4 discrete choice

# using optim for probit. need to create gradient function, and optim will take in intital values for betas,
# the probit function we created before, and then the gradient function for probit, and output the optimal
# results.

# write the probit gradient function, based on the F.O.C.s
pro_grad <- function (pro_grad_betas) {
  grad_x <- as.matrix(dat[,2:5]) #brings in data matrix X
  grad_ydum <- as.matrix(dat[,7]) # bring in ydum
  xb <- grad_x%*%pro_grad_betas # create XB
  p <- pnorm(xb) # apply F() to XB, where F() is cumulative standard normal distribution function
  d <- dnorm(xb) # apply F'() to XB, where F'() is the probability standard normal distribution function
  E <- ((grad_ydum - p) * d) / (p*(1-p)) # sum of the first part of the F.O.C.
  crossprod(E,grad_x)
}

# compose the optim command calling on intial betas, probit NLL function, and gradient.
pro_fit <- optim(betas, pro_neg_ll,pro_grad)
# pro_fit$par contains the estimated coefficents, which match the betas resulting from our manual MLE optimization
# in exercise 3. This is a great victory.
# Intreptation of coefficients will be at the end once we have the coefficients from all 3 models.


# logit optimization
# write logit negative log-likelihood function.
logit_neg_ll <- function (logit_beta) {
  logit_x <- as.matrix(dat[,2:5]) #brings in data matrix X
  logit_ydum <- as.matrix(dat[,7]) # bring in ydum
  xb <- logit_x%*%logit_beta # create XB
  p <- exp(xb)/(1+exp(xb)) # apply F() to XB, where F() is the appropriate logit function
  -sum( logit_ydum*log(p) + (1 - logit_ydum)*log(1 - p) ) # sum the log of each individual likelihood, then make negative
}

# write logit gradient function
logit_grad <- function (logit_grad_betas) {
  grad_x <- as.matrix(dat[,2:5]) #brings in data matrix X
  grad_ydum <- as.matrix(dat[,7]) # bring in ydum
  xb <- grad_x%*%logit_grad_betas # create XB
  p <- exp(xb)/(1+exp(xb)) # apply F() to XB, where F() is the appropriate logit function
  d <- exp(xb)/((1+exp(xb))^2) # apply F'() to XB, where F'() is derivate of the logit function
  E <- ((grad_ydum - p) * d) / (p*(1-p)) # sum of the first part of the F.O.C.
  crossprod(E,grad_x)
}

# compose the optim command calling on intial betas, logit NLL function, and gradient.
logit_fit <- optim(betas, logit_neg_ll,logit_grad)
# logit_fit$par contains the estimated coefficents.
# Intreptation of coefficients will be at the end once we have the coefficients from all 3 models.


# linear probability optimization.
# here no need for optim, just "manually" calcuate betas for OLS regression on ydum.
lin_x <- as.matrix(dat[,2:5])
lin_ydum <- as.matrix(dat[,7])
lin_xpxinv <- solve(t(lin_x)%*%lin_x)
lin_betas <- lin_xpxinv%*%t(lin_x)%*%lin_ydum

# create a small dataframe containing coefficients from these 3 optimizations for convenient comparison.
comp_df <- as.data.frame(cbind(pro_fit$par,logit_fit$par,lin_betas))
colnames(comp_df) <- c('probit','logit','linear prob')

# now we have our coefficents for each model. to interpret we will need significance, so lets use bootstrap to get
# the standard errors we need.

#>>>>>>>>>>>>>>>>>>>
# probit bootstrap standard errors.
probit_betadf499 <- data.frame(matrix(NA, nrow = 499, ncol = 4))
colnames(probit_betadf499) <- c("intercept","beta_x1","beta_x2","beta_x3")
# bootstrap in for-loop, sampling y 10000 times with replacement, then calculating betas. store those betas in probit_betadf

for (i in 1:499) {
  boot_d <- dat[sample(nrow(dat), 10000, replace = T), ]
  boot_ydum <- as.matrix(boot_d[,7])
  boot_x <- as.matrix(boot_d[,2:5])
  # probit NLL function
  pro_neg_ll <- function (pro_beta) {
    pro_x <- boot_x #brings in data matrix X
    pro_ydum <- boot_ydum # bring in ydum
    xb <- pro_x%*%pro_beta # create XB
    p <- pnorm(xb) # apply F() to XB, where F() is cumulative standard normal distribution function
    -sum( pro_ydum*log(p) + (1 - pro_ydum)*log(1 - p) ) # sum the log of each individual likelihood, then make negative
  }
  # probit gradient
  pro_grad <- function (pro_grad_betas) {
    grad_x <- boot_x #brings in data matrix X
    grad_ydum <- boot_y # bring in ydum
    xb <- grad_x%*%pro_grad_betas # create XB
    p <- pnorm(xb) # apply F() to XB, where F() is cumulative standard normal distribution function
    d <- dnorm(xb) # apply F'() to XB, where F'() is the probability standard normal distribution function
    E <- ((grad_ydum - p) * d) / (p*(1-p)) # sum of the first part of the F.O.C.
    crossprod(E,grad_x)
  }
  boot_pro_fit <- optim(betas, pro_neg_ll,pro_grad)
  probit_betadf499[i,] <- t(boot_pro_fit$par) # entire loop takes about a minute and a half.
}
probit_se499 <- apply(probit_betadf499, 2, sd)


#>>>>>>>>>>>>>
# logit bootstrap standard errors.
logit_betadf499 <- data.frame(matrix(NA, nrow = 499, ncol = 4))
colnames(logit_betadf499) <- c("intercept","beta_x1","beta_x2","beta_x3")
# bootstrap in for-loop, sampling y 10000 times with replacement, then calculating betas. store those betas in logit_betadf

for (i in 1:499) {
  boot_d <- dat[sample(nrow(dat), 10000, replace = T), ]
  boot_ydum <- as.matrix(boot_d[,7])
  boot_x <- as.matrix(boot_d[,2:5])
  # logit NLL function
  logit_neg_ll <- function (logit_beta) {
    logit_x <- boot_x #brings in data matrix X
    logit_ydum <- boot_ydum # bring in ydum
    xb <- logit_x%*%logit_beta # create XB
    p <- pnorm(xb) # apply F() to XB, where F() is cumulative standard normal distribution function
    -sum( logit_ydum*log(p) + (1 - logit_ydum)*log(1 - p) ) # sum the log of each individual likelihood, then make negative
  }
  # logit gradient
  logit_grad <- function (logit_grad_betas) {
    grad_x <- boot_x #brings in data matrix X
    grad_ydum <- boot_ydum # bring in ydum
    xb <- grad_x%*%logit_grad_betas # create XB
    p <- exp(xb)/(1+exp(xb)) # apply F() to XB, where F() is the appropriate logit function
    d <- exp(xb)/((1+exp(xb))^2) # apply F'() to XB, where F'() is derivate of the logit function
    E <- ((grad_ydum - p) * d) / (p*(1-p)) # sum of the first part of the F.O.C.
    crossprod(E,grad_x)
  }
  boot_logit_fit <- optim(betas, logit_neg_ll,logit_grad)
  logit_betadf499[i,] <- t(boot_logit_fit$par) # entire loop takes about a minute and a half.
}
logit_se499 <- apply(logit_betadf499, 2, sd)



#>>>>>>>>>>>>>
# linear probability bootstrap standard errors.
lin_betadf499 <- data.frame(matrix(NA, nrow = 499, ncol = 4))
colnames(lin_betadf499) <- c("intercept","beta_x1","beta_x2","beta_x3")
# bootstrap in for-loop, sampling y 10000 times with replacement, then calculating betas. store those betas in lin_betadf
for (i in 1:499) {
  boot_d <- dat[sample(nrow(dat), 10000, replace = T), ]
  boot_ydum <- as.matrix(boot_d[,7])
  boot_x <- as.matrix(boot_d[,2:5])
  boot_inv <- as.matrix(solve(t(boot_x)%*%boot_x))
  boot_beta <- as.matrix(boot_inv%*%t(boot_x)%*%boot_ydum)
  lin_betadf499[i,] <- t(boot_beta)
}

lin_se499 <- apply(lin_betadf499, 2, sd)

# now create similar comparison dataframe for the standard errors of the betas of probit,logit, and lin prop models
se_comp_df <- as.data.frame(cbind(probit_se499,logit_se499,lin_se499))
colnames(se_comp_df) <- c('probit','logit','linear prob')

t_stat_df <- comp_df /se_comp_df #using the matrix of prob/logit/lin betas (namely comp_df) and the matrix of their standard errors,
# create matrix of t_statistics.

p_val_df <- data.frame(matrix(NA, nrow = 4, ncol = 3)) #create matrix of p_values of t_statistics using pt('beta', df = 10000)
p_val_df[,1] <- pt(abs(t_stat_df[,1]),10000)
p_val_df[,2] <- pt(abs(t_stat_df[,2]),10000) # fill the matrix
p_val_df[,3] <- pt(abs(t_stat_df[,3]),10000)

# p_val_df tells us the significance of the betas found in comp_df for probit, logit, and lin prob. models. for logit and probit,
# beta_x1 and beta_x2 are both highly significant. Thus they can be interpreted only as far as their sign. The x1 and x2betas for 
# logit and probit are all positive, thus we can say they have a positive effect on the probability of ydum = 1 that is statistically
# significant. the betas for x1 for logit and probit are not statistically significant at a 5% level, so we cannot say their effect
# on ydum is statistically significantly different from zero.

# interpretation for lin prob. beta_x3 is not statistically significant at even a 10% level, thus its effect on ydum is not significantly
# diffferent from zero. beta_x1 is significant at a 1% level, so we say for a 1 unit increase in x1, the estimated increase of ydum
# is .143. Beta_x2 is significant at a 1% level, so we say for a 1 unit increase in x2, the estimated increase of ydum is -0.103.
# note for the linear probability model, fitted values may exceed 1 or be lower than 0. since fitted values can be interpreted as
# probability of y = 1, probability > 1 or < 0 do not make sense.


#******************************************
# exercise 5 Marginal Effects
# fetch betas and variance/covariance matrices from the glm probit and logit regressions
logit_glm <- glm(ydum~x1+x2+x3, family = binomial(logit))
probit_glm <-glm(ydum~x1+x2+x3, family = binomial(probit))
# store betas
logit_betas <- logit_glm$coefficients
probit_betas <- probit_glm$coefficients
# store vcov matrices 
logit_vcov <- vcov(logit_glm)
probit_vcov <- vcov(probit_glm)

#calculate marginal effects: NOTE: these marginal effects are for the representative individual, xbar
#probit marginal effects
probit_me <- dnorm(c(colMeans(x))%*%probit_betas)*probit_betas 
#logit marginal effects
logit_me <- exp(c(colMeans(x))%*%logit_betas)/((1+exp(c(colMeans(x))%*%logit_betas))^2)*logit_betas


# now compute standard deviations of marginal effects via delta method.
# here we will use those vcov matrices.
# probit
# using the jacobian function from numDeriv package, we can feed it a probit marginal effect function and the probit coefficients
# and it will return the jacobian matrix. then multiply t(jacobian)*vcov(probit)*jacobian to yield var covariance matrix of probit
# Marginal effects. same process is repeated for logit.

probit_me_fun <- function(probit_b) { #function that takes probit betas and returns marginal effects
  dnorm(c(colMeans(x))%*%probit_b)*probit_b
}

jac_prob <-jacobian(probit_me_fun,probit_betas) # probit jacobian
delt_p <-t(jac_prob)%*%probit_vcov%*%jac_prob
jac_prob_sd <-sqrt(diag(delt_p)) #taking the square root of diagonal terms
# probit standard deviations produced by delta method:  0.040020902   0.017199120   0.006647109   0.018844260

#logit delta method
logit_me_fun <- function(logit_b) { #function that takes logit betas and returns marginal effects of average individual
  exp(c(colMeans(x))%*%logit_b)/((1+exp(c(colMeans(x))%*%logit_b))^2)*logit_b
}

jac_logit <-jacobian(logit_me_fun,logit_betas) # logit jacobian
delt_l <-t(jac_logit)%*%logit_vcov%*%jac_logit
jac_logit_sd <-sqrt(diag(delt_l)) #taking the square root of diagonal terms
# logit standard deviations produced by delta method:  0.025312935   0.010731260   0.003287703   0.011800993
 
# standard deviations by bootstrap
# probit bootstrat
me_probit_df499 <- data.frame(matrix(NA, nrow = 499, ncol = 4))
colnames(me_probit_df499) <- c("intercept","beta_x1","beta_x2","beta_x3")
for (i in 1:499) {
  boot_d <- dat[sample(nrow(dat), 10000, replace = T), ]
  boot_ydum <- as.matrix(boot_d[,7])
  boot_x <- as.matrix(boot_d[,2:5])
  boot_probit_glm <-glm(boot_ydum~boot_x[,2]+boot_x[,3]+boot_x[,4], family = binomial(probit))
  boot_probit_betas <- boot_probit_glm$coefficients
  boot_probit_me <- dnorm(c(colMeans(x))%*%boot_probit_betas)*boot_probit_betas
  me_probit_df499[i,] <- boot_probit_me
}
me_probit_se499 <- apply(me_probit_df499, 2, sd)
# probit standard deviations produced by bootstrap:  0.038586497   0.016186590   0.007057369   0.018398691 


# logit bootstrap
me_logit_df499 <- data.frame(matrix(NA, nrow = 499, ncol = 4))
colnames(me_logit_df499) <- c("intercept","beta_x1","beta_x2","beta_x3")
for (i in 1:499) {
  boot_d <- dat[sample(nrow(dat), 10000, replace = T), ]
  boot_ydum <- as.matrix(boot_d[,7])
  boot_x <- as.matrix(boot_d[,2:5])
  boot_logit_glm <-glm(boot_ydum~boot_x[,2]+boot_x[,3]+boot_x[,4], family = binomial(logit))
  boot_logit_betas <- boot_logit_glm$coefficients
  boot_logit_me <- dlogis(c(colMeans(x))%*%boot_logit_betas)*boot_logit_betas
  me_logit_df499[i,] <- boot_logit_me
}
me_logit_se499 <- apply(me_logit_df499, 2, sd)
# logit standard deviations produced by bootstrap: 0.07499633 0.02984987 0.01308334 0.03177344 

# comparison                                            intercept     ME_beta_x1    ME_beta_x2    ME_beta_x3     
# probit standard deviations produced by delta method:  0.040020902   0.017199120   0.006647109   0.018844260
# probit standard deviations produced by bootstrap:     0.038586497   0.016186590   0.007057369   0.018398691 
# logit standard deviations produced by delta method:   0.047285764   0.019974584   0.006438529   0.021140521
# logit standard deviations produced by bootstrap:      0.045828773   0.020022542   0.008503187   0.021023087 



