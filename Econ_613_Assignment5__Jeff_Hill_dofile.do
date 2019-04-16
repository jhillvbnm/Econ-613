/* Jeff Hill
Econ 613
Assignment 5 */

cd "/Users/admin/Documents/Econ_613/Assignment Output"
set seed 613
log using "/Users/admin/Documents/Econ_613/Assignment Output/Assignment5_JeffHill.smcl", replace

**********************************************
*************** Assignment 2 *****************
**********************************************

**************** Exercise 1  *****************
set obs 10000
* generate our data
gen x1 = runiform(1,3)
gen x2 = rgamma(3,2)
gen x3 = rbinomial(1,.3)
gen eps = rnormal(2,1)
gen y = 0.5 + 1.2*x1 - 0.9*x2 + 0.1*x3 + eps
egen y_mean = mean(y)
gen y_dum = 0
replace y_dum = 1 if y>y_mean

**************** Exercise 2  *****************

corr y x1
* the correlation is .1948, which is not close to 1.2. It should not be, as this
* is correlation, bound between -1 and 1

reg y x1 x2 x3
eststo ols
* the standard OLS standard errors are contained within the regression output.

* bootstrap se's
* 49 replications
bootstrap, reps(49) : reg y x1 x2 x3
eststo boot49
* 499 replications
bootstrap, reps(499) : reg y x1 x2 x3
eststo boot499

* create a table
esttab ols boot49 boot499, se(3) title(Assignment 2: OLS, Bootstrap 49 and 499) ///
 nonumbers mtitles("OLS" "Bootstrap 49" "Bootstrap499")
 

**************** Exercise 3  *****************
* probit
probit y_dum x1 x2 x3
eststo probit2

/* the estimates for the probit model are all significant at the 99% 
significance level. Interpreting the coefficent on x2 (-.8879) says that an
increase in x2 has a statisticaly significant negative impact on y_dum. */

**************** Exercise 4  *****************
* logit
logit y_dum x1 x2 x3
eststo logit2

/* the estimates for the logit model are all significant at the 99% 
significance level. Interpreting the coefficent on x2 (-1.604) says that an
increase in x2 has a statisticaly significant negative impact on y_dum. */


* linear model
reg y_dum x1 x2 x3
eststo lin2

/* the estimates for x1 and x2 in the linear model are significant at the 99% 
significance level. x3 is significant at the 90% level. Interpreting the 
coefficent on x2 (-.10431) says that an increase in x2 by 1 unit results  y_dum
decreasing by -.1043. */

* table comparing estimates from the 3 regressions
esttab probit2 logit2 lin2, se(3) title(Assignment 2: Probit, Logit, & Linear) ///
 nonumbers mtitles("Probit" "Logit" "Linear Prob")

**************** Exercise 5  *****************
* marginal effects for probit
probit y_dum x1 x2 x3
margins, dydx(*) post
eststo probit_me_delta2



* marginal effects for logit
logit y_dum x1 x2 x3
margins, dydx(*) post
eststo logit_me_delta2

/* the delta method standard errors are calculated and reported in the marginal
effect tables above. now all that is left is to calculate them via bootstrap. */ 

bootstrap, reps(499) : probit y_dum x1 x2 x3
margins, dydx(*) post
eststo probit_me_boot2

bootstrap, reps(499) : logit y_dum x1 x2 x3
margins, dydx(*) post
eststo logit_me_boot2

esttab probit_me_delta2 logit_me_delta2 probit_me_boot2 logit_me_boot2, se(6) ///
 title(Assignment 2: Probit, Logit, & Linear) nonumbers ///
 mtitles("Probit Delta" "Logit Delta" "Probit Boot" "Logit Boot")
clear
**********************************************
*************** Assignment 3 *****************
**********************************************
cd "/Users/admin/Documents/Econ_613/Data/Assignment 3"
insheet using "product.csv",clear
cd "/Users/admin/Documents/Econ_613/Assignment Output"

**************** Exercise 1  *****************
sum
* the product means and standard deviations are listed in the sum table which
* describes the data.

tab choice
* tab choice breaks down the market share of each product in percentages

gen brand = choice
replace brand = 1 if brand == 8
replace brand = 3 if brand == 9
replace brand = 4 if brand == 10
tab brand
* brand now groups sticks and tubs made by the same company together, and tab
* displays the market share.

* import data that was merged outside stata, since we need to use data for demos 
* and product
cd "/Users/admin/Documents/Econ_613/Data/Assignment 3"
insheet using "margarine data.csv",clear
cd "/Users/admin/Documents/Econ_613/Assignment Output"

* two way tabluate breaks down choice by income group
tab choice income, row nofreq

**************** Exercise 2  *****************
* need to reshape the data

reshape long sel p, i(id) j(selection 1 2 3 4 5 6 7 8 9 10) 
* asmixlogit is easier to use than clogit, and yields the same results.
asmixlogit sel p, case(id) alternatives(selection)
eststo clogit_betas_3

* need tables for each regression here since the betas aren't comparable.
esttab clogit_betas_3, noomitted nostar unstack compress se(3) ///
nonumbers title(Assignment 3: Conditional Logit Estimates)


**************** Exercise 3  *****************

mlogit choice income fs3_4 fs5 college whtcollar retired
eststo mlogit_betas_3

* again need a separate table
esttab mlogit_betas_3, noomitted nostar unstack compress se(3) ///
nonumbers title(Assignment 3: Multinomial Logit Estimates)

**************** Exercise 4  *****************
* quietly because no need to repeat the results twice, we just want the marginal
* effects.
quietly clogit sel p, group(id)
margins, dydx(*) post
eststo clogit_me_3

esttab clogit_me_3, noomitted nostar unstack compress se(3) ///
nonumbers title(Assignment 3: Conditional Logit Marginal Effects)

quietly mlogit choice income fs3_4 fs5 college whtcollar retired
mfx

**************** Exercise 5  *****************
asclogit sel p , case(id) alternatives(selection) ///
  casevars(income fs3_4 fs5 college whtcollar retired)
eststo mixed_logit_betas_3
esttab mixed_logit_betas_3, noomitted nostar unstack compress se(3) ///
nonumbers title(Assignment 3: Mixed Logit Regression)

* drop choice 4 and rerun the regression
drop if choice == 4
drop if selection == 4

asclogit sel p , case(id) alternatives(selection) ///
  casevars(income fs3_4 fs5 college whtcollar retired)
eststo mixed_logit_no4_3
esttab mixed_logit_no4_3, noomitted nostar unstack compress se(3) ///
nonumbers title(Assignment 3: Mixed Logit Regression Without Choice = 4)

* hausman test
hausman mixed_logit_betas_3 mixed_logit_no4_3, alleqs constant

**********************************************
*************** Assignment 4 *****************
**********************************************
cd "/Users/admin/Documents/Econ_613/Data/Assignment 4"
insheet using "Koop-Tobias.csv",clear
cd "/Users/admin/Documents/Econ_613/Assignment Output"

**************** Exercise 1  *****************
keep personid logwage timetrnd
reshape wide logwage, i(personid) j(timetrnd)
sample 5, count
count
list
clear

**************** Exercise 2  *****************
cd "/Users/admin/Documents/Econ_613/Data/Assignment 4"
insheet using "Koop-Tobias.csv",clear
cd "/Users/admin/Documents/Econ_613/Assignment Output"
xtset personid timetrnd
xtreg logwage educ potexper
eststo re_betas_4
esttab re_betas_4, se(3) nonumbers title(Assignment 4: Random Effects Model)

**************** Exercise 3  *****************
* between estimator
xtreg logwage educ potexper, be
eststo between_betas_4

* within estimator
xtreg logwage educ potexper, fe
eststo within_betas_4

* first time difference estimator
reg d.(logwage educ potexper), noconstant
eststo firstd_betas_4
esttab between_betas_4 within_betas_4 firstd_betas_4, se(3) nonumbers ///
title(Assignment 4: Fixed Effects Estimators) mtitles("Between" "Within" "First Difference")

**************** Exercise 4  *****************
keep personid educ potexper logwage timetrnd ability mothered fathered ///
brknhome siblings
reshape wide educ logwage potexper, i(personid) j(timetrnd)
sample 100, count
reshape long

xtreg logwage educ potexper, fe
predict alphas, u
drop if missing(alphas)
duplicates drop alphas, force

* regress the alphas on the time invariant variables 
reg alphas ability mothered fathered brknhome siblings
eststo alpha_betas_4

bootstrap, reps(199): reg alphas ability mothered fathered brknhome siblings
eststo alpha_betas_4_boot
esttab alpha_betas_4 alpha_betas_4_boot, se(4) title(Assignment 4: Individual Fixed Effect ///
  Regressions) nonumbers ///
  mtitles("Non-Corrected" "Bootstrap 199")

*********************************************
log close
