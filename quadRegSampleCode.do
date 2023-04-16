
****************************************************
****************************************************
*												   *
*   Stata File to Implement Quadratic Regression   *
*   to test for Predictive Adaptive Response, 	   *
*   Developmental Constraints, and Adult           *
*   Environmental Quality theories                 *
*												   *
****************************************************
****************************************************

/***************************************************
* Background

y = adult outcome
e0 = early life environment
e1 = late life environment
ed = |e1 - e0| = absolute value of change in environment

In general, adult outcome may depend on early environment and change in environment:

y = f(e0, ed)

where f is continuous and assumed to be continuously differentiable and we assume that the data generating process is that nature gives e0, then ed, adn e1 is determined by ed.  

We don't know the functional form for f so we use a second order Taylor approxmation.  The second order allows an interaction term, which is what the prior literature used to test for Predictive Adaptice Response (PAR).  Since we can demean early life environment, we take an approximation around 0, which is the  average value of the demeaner early environment.  Moreover we assume the  environment is trendless, with occasional or cross-sectional shocks, so 
E(e1)=e0.  The Taylor approximation is

y = [f0] + [d/d0 f0] e0 + [d/d1 f0] ed 
    + [d^2/de0^2 f0] e0^2 + [d^2/d^ed f0] ed^2 + [d^2/(de0 ded) f0] e0*ed + o  

  = [b] + [b0] e0 + [bd] ed 
    + [b00] e0^2 + [bdd] ed^2 + [b0d] e0*ed + u  

	
where f0 = f(0,0), d are artial derivatives, the expressions in square brackets are regression cefficients, and u is the (regression) error.  When we estimate the equation above, we assume the researcher has access to exogenous variation in e0 and ed.  Equivalently, we assume away problems with causation.
	
We use the following definitions of the developmental constraints (DC) and predictive adaptive response (PAR) theories:

* Developmental constraints (DC) implies: dy / de0 > 0
* PAR: dy / ded < 0
* AEQ: dy / de1 > 0

The definitions imply the following tests using estimated coefficients from the regression above:

* DC: b0 + 2 b00 e0 + b0d ed > 0
* PAR: bd + 2 bdd ed + b0d e0 < 0

In other words, e.g., the test for whether the data show DC is whether the non-linear equation above produces a result that is statistically significantly greater than 0.   

(We do not list the adult environmental quality (AEQ) test because one cannot modify e1 independently of ed or e0. To test AEQ, we need to assume a different data generating process, e.g., nature gives e0 and e1, and ed comes from |e1 - e0|, and suppose that y = f(e0,e1).  In this specification, one can test DC and AEQ, but not PAR.)   
	
***************************************************/

* Load data

// use "filename.dta", clear

* Clean data

// We assume the data are labeled y, e0, and ed using the definintions in the background section.  We also assume there are no additional controls.  THese include group or individual fixed effects.  I.e., we are assuming the analysis is cross-sectional.

// Generate quadratic terms
gen e02 = e0^2
gen ed2 = ed^2
gen e0d = e0*ed

* Basic quadratic regression

reg y e0 ed e02 ed2 e0d

* Test for DC and PAR at mean of data

// Gather means of each input variable
su e0
local m0 = r(mean)

su ed
local md = r(mean)

// Test for DC: check if the estimate below is positive and the p-value is below your desired critical value 
testnl _b[e0] + (2 * _b[e02] * `m0') + (_b[e0d] * `md')


// Test for PAR: check if the estimate below is negative and the p-value is below your desired critical value 
testnl _b[ed] + (2 * _b[ed2] * `md') + (_b[e0d] * `m0')

