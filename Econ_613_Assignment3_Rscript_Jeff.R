# Jeff Hill, ECON 613 Assignment 3
library(bayesm)
library(plyr)
library(dplyr)

# load data
datapath <- "/Users/admin/Documents/Econ_613/Data/Assignment 3/"
demos <- read.csv(file=paste(datapath,"demos.csv",sep=""), header=TRUE, sep=",")
product <- read.csv(file=paste(datapath,"product.csv",sep=""), header=TRUE, sep=",")

#************************************************
# Exercise 1 Data Description
#************************************************
# average and variance of each product:
# mean
prodmean <- apply(product[,4:13], 2, mean)
prodmean
# standard deviation
prodsd <- apply(product[,4:13], 2, sd)
prodsd

# market share by product
mshare <-table(product$choice)/(sum(table(product$choice)))
mshare 

# market share by brand
# first combine stick and tub for brands: PPk, PFl, and PHse as they are the only brands offering both stick and tub.
brand_choice <- product$choice
brand_choice[brand_choice == 8] <- 1
brand_choice[brand_choice == 9] <- 3
brand_choice[brand_choice == 10] <- 4
# now choices are grouped only by brand in brand_choice, so again use table
brand_share <-table(brand_choice)/(sum(table(brand_choice)))
brand_share

# market share by observed characteristic, i.e. income level, family size, colllege, etc.
# need to link demographic data to product data, so will create new df, and merge via hhid.
df <- data.frame(merge(product, demos, by="hhid"))
df$intercept <- 1 # add a column to df, containing all 1's for an intercept. this will be used in exercise 2 and 3.

# create a vector of income levels
inc_levels <- c(seq(2.5, 47.5, by=5), 55, 67.5, 87.5, 130) 

# create an empty dataframe to insert income share data into.
inc_share_df <- data.frame(prod1 = NA, prod2 = NA, prod3 = NA, prod4 = NA, prod5 = NA, prod6 = NA, prod7 = NA, prod8 = NA, prod9 = NA, prod10 = NA)
colnames(inc_share_df) <- c(1:10) # rename column for sake of convenience, they will be renamed at the end.

# here we use a nested for loop, to cycle through first income levels and then cycle through the 10 products inside that.
# we insert the market share of each product at a given income level into a 14 x 10 dataframe.
for(i in inc_levels) {
  for(j in c(1:10)) {
    inc_share_df[i,j] <- sum(df[df["Income"] == i,'choice']==j)/length(df[df["Income"] == i,'choice'])
  }
}

# our nested for loop generated a bunch of empty rows because we indexed on income level which goes up to 130,
# so we just drop the empty rows here.
inc_share_df <- inc_share_df[rowSums(is.na(inc_share_df)) != ncol(inc_share_df),]

# rename rows and columns
rownames(inc_share_df) <- c('inc2.5','inc7.5','inc12.5','inc17.5','inc22.5','inc27.5','inc32.5',
                            'inc37.5','inc42.5','inc47.5','inc55','inc67.5','inc87.5','inc130')
colnames(inc_share_df) <- c('product1','product2','product3','product4','product5',
                            'product6','product7','product8','product9','product10')
inc_share_df # this is the final dataframe, containing product share broken down at each income level. for certain income levels,
             # the number of individuals was low, so there are some 0% market shares, for example at income level = 2.5 no one 
             # bought product 3, 6, or 10.


#************************************************
# Exercise 2 First Model
#************************************************
# here we create the conditional logit negative log likelihood function. The proposed model specification includes all 10 prices for
# the 10 products.
clogit_nll <- function (clogit_parameters,x,y) { # here clogit_parameters will be a 10x1 vector containing 1 beta(price) 
                                                 # and 9 alphas, since we have 10 product choices and we force the first alpha
                                                 # to be zero.
  clogit_beta <- clogit_parameters[1] # separating out beta and alphas from clogit_parameters
  alpha <- clogit_parameters[2:10]
  a_mat <- cbind(0,matrix(alpha,nrow=4470,ncol=9,byrow=T)) # create alpha matrix to account for alternative specific constants.
                                                           # we cbind in a column of zeros to set the first alpha to our reference
  xb <- x*clogit_beta + a_mat # create XB, size:n x 10. Beta is just the single beta for price, and there are 9 alphas

  rs <- rowSums(exp(xb)) # this is the denominator of P_ij
  xb_vec <- rep(0, length(y)) # create an empty vector from which we will select the correct xb value, 
                              # based on which choice individual i made.
  for (i in 1:length(y)) {  # this for loop goes through the matrix xb by row, and pulls out the correct xb value based on what choice
    xb_vec[i] <- xb[i,y[i]] # that individual made, and stores the value in xb_vec. e^(xb_vec) is now the numerator of P_ij
  }
  -sum(xb_vec - log(rs))  # here we calculate the negative log likelihood. 
}

cond_x <- as.matrix(df[,c(4:13)]) # here we subset the x for the conditional logit out of df. we only pull out the columns of
                                  # interest, namely the price columns for the 10 products.

#clogit_estimates contains the optimized values for Beta, and the 9 alphas. note these are alphas for products 2 through 10,
# as alpha for product 1 was bound = 0 to set it as our reference point.
clogit_estimates <- nlm(clogit_nll,rep(0,10),x=cond_x,y=df$choice)$estimate
clogit_price_beta <-clogit_estimates[1]
clogit_alphas <- clogit_estimates[c(2:10)]
clogit_price_beta # this is the beta for price
clogit_alphas # these are the alphas for products 2 thorugh 10, relative to product 1 (alpha 1 = 0)

#interpret the coefficient on price:
# an increase in price of a product results in decreasing the probability of choosing that product

#************************************************
# Exercise 3 Second Model
#************************************************

# here we create the multinomial logit log likelihood function. Much of the function is identical to the clogit_nll function above.
# our model specification includes income, the two family size dummy variables (with family size = 1 or 2 being the reference point),
# college, whtcollar and retired.
mlogit_nll <- function (mlogit_beta,x,y) {
  mat <- matrix(mlogit_beta,nrow=7,ncol=9, byrow=T) # turn out mlogit_beta vector into a matrix of the correct size.
  mat <- cbind(0,mat) # add on a column of zeros on front, to account for product 1 being the reference point.
  xb <- x%*%mat # create XB size:n x 10. We do this through matrix multiplication of x_i * Beta_j, which differs from cond. logit.
  rs <- rowSums(exp(xb)) # this is the denominator of multinomial logit probability P_ij
  xb_vec <- rep(0, length(y))
  for (i in 1:length(y)) { # again this for loop pulls out the correct xb term, based on the choice of that individual.
      xb_vec[i] <- xb[i,y[i]]
  }
  -sum(xb_vec - log(rs)) # sum the log of each individual likelihood, then make negative
}

multi_x <- as.matrix(df[,c(22,15,16,17,19,20,21)]) # here subset x out of df, pulling out only the columns of interest: intercept,
                                                   # income, Fs3_4, Fs5.,whtcollar,retired, and college

mlogit_opt <- nlm(mlogit_nll,rep(0,63),x=multi_x,y=df$choice)

#mlogit_estimates contains the optimized values for our 63 Betas. there are 9 for each independent variable, and 9 for the intercept.
mlogit_estimates <- mlogit_opt$estimate
m_ind_betas <- as.data.frame(matrix(mlogit_estimates,nrow=7,ncol=9,byrow=T))
rownames(m_ind_betas) <- c('intercept','income','Fs3_4','Fs5.','college','whtcollar','retired')
colnames(m_ind_betas) <- c('product2','product3','product4','product5',
                           'product6','product7','product8','product9','product10')
m_ind_betas # m_ind_betas contains the betas for the 6 individual variables, for products 2 thorugh 10 relative to product 1.

#interpret the coefficient on family: I was not sure if I was supposed to interpret the coefficient on
# family (Fs3_4,Fs5.) or family income, so I do both.

# interpretation of income:
m_ind_betas[2,]
# the coefficient from income for product 2 says an increase in income results in having a lower probability
# of choosing product 2 relative to product 1. This is the same interpretation for all income coefficients
# based on the SIGN of the coefficient (negative or positive). (This is also assuming significance)

# I included the dummy variables for Family size 3-4 and family size >5
m_ind_betas[3:4,] # these are the coefficients for the family dummys.
# the coefficient for Fs3_4 for Product 2 says Having a family size of 3-4 compared to having a family size of 1-2
# has a negative impact on your likelihood of choosing product 2 relative to product 1. (i.e you are more likely
# to switch product 1 from product 2 if you "switch" from family size 1-2 to size 3-4, holding all else constant.
# this holds the same for "switching" from family size 1-2 to size >5 based on the SIGN of the coefficient
# (This is assuming significance)

#*****************************************
# Exercise 4 Marginal Effects
#*****************************************

#NOTE: I did not complete this section before the midnight deadline, I only made it through part of the 
# multinomial logit marginal effects.

mlogit_p <- function (mlogit_beta,x,y) {
  mat <- matrix(mlogit_beta,nrow=7,ncol=9, byrow=T) 
  mat <- cbind(0,mat) 
  xb <- x%*%mat 
  rs <- rowSums(exp(xb)) 
  exp(xb)/rs
}
p <- mlogit_p(mlogit_estimates,multi_x,df$choice) #following the formula for multinomial
# marginal effects, here we calculate the matrix p containing all p_ij
# the p matrix must then be multiplied by (B_j - B_i_bar), where B_i_bar = the sum of p_il*B_l
me_mlogit_betas <- t(cbind(0,m_ind_betas))
B_bar <- p%*%me_mlogit_betas


#*****************************************
# Exercise 5 IIA
#*****************************************
# here we construct the mixed logit negative log likelihood function, which combines product and individual characteristics.
mixedlogit_nll <- function (mixedlogit_parameters,x,y) { # mixedlogit_parameters is now a vector of length 64, 1 price beta and
                                                         # 63 individual characteristic betas.
  c_x <- x[,c(1:10)]  # here since x will now contain product data AND demographic data, we separate them into
  m_x <- x[,c(11:17)] # c_x containing product data and m_x containing demographic data.
  beta <- mixedlogit_parameters[1] # separating out beta, and gammas from mixedlogit_parameters
  gamma <- mixedlogit_parameters[2:64]
  # there is no need for an alpha matrix in this specification, since the intercept is already included in m_x
  # we cbind in a column of zeros to set the first alpha to our reference
  g_mat <- matrix(gamma,nrow=7,ncol=9, byrow=T) #create gamma matrix for demographic data
  g_mat <- cbind(0,g_mat)
  xb <- c_x*beta + m_x%*%g_mat # create XB now with c_x and m_x
  rs <- rowSums(exp(xb)) # this is the denominator of P_ij
  xb_vec <- rep(0, length(y)) # create an empty vector from which we will select the correct xb value, 
                              # based on which choice individual i made.
  for (i in 1:length(y)) {  # this for loop goes through the matrix xb by row, and pulls out the correct xb value based on what choice
    xb_vec[i] <- xb[i,y[i]] # that individual made, and stores the value in xb_vec. e^(xb_vec) is now the numerator of P_ij
  }
  -sum(xb_vec - log(rs))  # here we calculate the negative log likelihood. 
}

mixed_x <- as.matrix(df[,c(4:13,22,15,16,17,19,20,21)]) #this is the subset of df, containing all relevant x data.

mixed_logit_opt <- nlm(mixedlogit_nll,rep(0,64),x=mixed_x,y=df$choice) #the optimization
mixed_logit_estimates <- mixed_logit_opt$estimate
mixed_price_beta <- mixed_logit_estimates[1]
mixed_price_beta # this is the beta for price in our mixed model.
mixed_ind_betas <- as.data.frame(matrix(mixed_logit_estimates[2:64],nrow=7,ncol=9,byrow=T))
rownames(mixed_ind_betas) <- c('intercept','income','Fs3_4','Fs5.','college','whtcollar','retired')
colnames(mixed_ind_betas) <- c('product2','product3','product4','product5',
                           'product6','product7','product8','product9','product10')
mixed_ind_betas # mixed_ind_betas contains the betas for the 6 individual variables, for products 2 thorugh 10 relative to product 1
                # from the mixed logit regression. This is Beta^f in the assignment.

# Now an alternative specification: I will remove data for choice 4. This involves removing associated column, but also all
# observations that selected choice 4.
df_cut <- df[df$choice!=4,] # removed all observations where someone chose product 4
df_cut <- subset(df_cut, select=-c(PHse_Stk))
# now since 4 has been removed, I will shift choice 5 to 4, 6 to 5, etc. 
for (i in 1:3877) { # this quick forloop finds all values of choice above 4, and reduces them by 1.
  if (df_cut[['choice']][i] > 4) {
    df_cut[['choice']][i] <- df_cut[['choice']][i]-1
  }
}

# now modify the mixed logit log likelihood function to run on the reduced dataset:
cut_mixedlogit_nll <- function (mixedlogit_parameters,x,y) { # It is exactly the same, with a couple vector lengths shifted.
  # 56 individual characteristic betas. (7x8 now)
  c_x <- x[,c(1:9)] 
  m_x <- x[,c(10:16)]
  beta <- mixedlogit_parameters[1] 
  gamma <- mixedlogit_parameters[2:57]
  g_mat <- matrix(gamma,nrow=7,ncol=8, byrow=T) 
  g_mat <- cbind(0,g_mat)
  xb <- c_x*beta + m_x%*%g_mat
  rs <- rowSums(exp(xb)) 
  xb_vec <- rep(0, length(y)) 
  for (i in 1:length(y)) {  
    xb_vec[i] <- xb[i,y[i]] 
  }
  -sum(xb_vec - log(rs)) 
}

cut_x <- as.matrix(df_cut[,c(4:12,21,14,15,16,18,19,20)])

cut_mixed_logit_opt <- nlm(cut_mixedlogit_nll,rep(0,57),x=cut_x,y=df_cut$choice) #the optimization
cut_mixed_logit_estimates <- cut_mixed_logit_opt$estimate
cut_mixed_price_beta <- cut_mixed_logit_estimates[1]
cut_mixed_price_beta # this is the beta for price in our mixed model.
cut_mixed_ind_betas <- as.data.frame(matrix(cut_mixed_logit_estimates[2:57],nrow=7,ncol=8,byrow=T))
rownames(cut_mixed_ind_betas) <- c('intercept','income','Fs3_4','Fs5.','college','whtcollar','retired')
colnames(cut_mixed_ind_betas) <- c('product2','product3','product5',
                               'product6','product7','product8','product9','product10')
cut_mixed_ind_betas # mixed_ind_betas contains the betas for the 6 individual variables, for products 2 through 10 (minus 4 now) 
# relative to product 1 from the mixed logit regression.

# here we calculate the negative log likelihood of the restricted model (less product 4) of the restricted betas.
restricted_lik <- cut_mixedlogit_nll(cut_mixed_logit_estimates,cut_x,df_cut$choice)

# to calculate the same for the unrestricted betas (beta^f), we need to select the proper betas out of our beta vector
# selecting all betas except those associated with product 4:
unrestricted_betas <- mixed_logit_estimates[-c(4,13,22,31,40,49,58)]
unrestricted_lik <- cut_mixedlogit_nll(unrestricted_betas,cut_x,df_cut$choice)

# construct test statistic:
MTT <- 2*(unrestricted_lik-restricted_lik) #since our log likelihood function is a negative log likelihood function,
# multiplying by another negative negates the minus in front.
MTT 
p
qchisq(.001,df= 57)
# since our test statistic MTT is less than the critical value for the chi- squared distribution at alpha = .95
# we conclude that IIA does hold. IIA may not hold depending on the product choice. Here product 4 was chosen roughly
# 500 times out of the 4470 choices, but selecting a more popular choice to remove from the data might have a
# more powerful effect on altering IAA.