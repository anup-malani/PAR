* This is a Stata do file (compatible with Stata 17.0) that implements the 
* simulation in Malani, Rosenbaum, Archie, Alberts (2022) from scratch.

********************
********************
* Preliminaries
********************
********************

// Continuous scrolling, clear data, set memory, set seed for random numbers
	set more 1
	drop _all
	set max_memory 128g
	set seed 32679 // Happy bday, Stacy!

// Log file
	cap log close
	cap log using sim_210629, replace

// Libraries
	local datalib "/afs/crc.nd.edu/user/a/amalani/PAR"
	local homelib "/afs/crc.nd.edu/user/a/amalani/PAR"
	local parlib "/afs/crc.nd.edu/user/a/amalani/PAR"
	*local datalib "/scratch/midway2/amalani"
	*local homelib "/home/amalani/"
	*local parlib "/home/amalani/rosenbaum/par"
	*local datalib "\\midwaysmb.rcc.uchicago.edu\homes\rosenbaum\parcap 

//Change directory
	cap cd "`datalib'"

// Notes to self
	*Q: Do we want to aggregate groups of monkeys?
	
// Toggles --> allows you to run parts of code at a time
	// Toggle macros to 1 to run, 0 to not run
	// Toggle 1 sim at a time so you don't generate excess files and violate quota

	local gtoolkit 0
	local model 0 // permanently 0 so we can collapse text
	local func_strat 0 // permanently 0
	local func_ninp 0 // permanently 0

	local indvar 0 // generate unique covariate vectors
	local make_coef_10 0
	local test_coef_10x16 0
	local test_coef_10x16_merge 0
	local screen_coef 0
	
	local sim_data 0 // generate sim data sets
	local test_regs 0 // test for DC and PAR with sim data sets

// Load gtoolkit to work with large data sets more quickly
	if `gtoolkit'==1{

		* Load necessary software
			
			cap ssc install gtools
			gtools, upgrade

	}

********************
********************
* Simulation primitives
********************
********************

	// We generate realities.  Then we generate a data set for each reality to test for PAR, DC in that reality.  
	// Each data set has `regressor' covariate vectors x `u' error terms. Total `ressors'*`u' obs.
	local p = 3 // number of powers for Taylor series expansion of for truth in each reality = 3
	local regressors 100 // building number of covariate vectors for a reality
	local u 20 // multiples of each caovariate vector, i.e., generate u errors for each covariate vector  
	
********************
********************
* INDEPENDENT VARS
********************
********************

if `indvar'==1{ 

	/* 
	The goal is to create a data set of `regressors' obs on var for purpose of simulating.
	Later we can add obs to make the distribution of errors look nicer.  
	*/

	// Create data set with `regressors' empty obs 	
		clear *
		set obs `regressors' // sim size --> keep low to reduce size of data sets for Stata

	// Create e0 (initial rank): uniformly distributed from 0 to 1
		gen e0 = (_n-1)/(_N-1) // higher is better
		la var e0 "Childhood rank"

	// Create \Delta e or echange (change in rank): normally distributed with mean and SD from baboon data	
		gen echange = rnormal(-0.03,0.21) // drawn from su of momsubjectrankdiff 
		// i.e., (rank of baboon as adult) MINUS (rank that baboon's mother held at time infant was born)
		la var echange "Random change in rank"

	// Create final rank: is sum of e0 and echange
	// This is NOT an important step for the reality that we create 
	// because our reality is y_1 = f(e0,echange). e1 plays no role.
		gen e1 = e0 + echange // calculate final rank as sum of e0 and echange
		la var e1 "Adult rank"
		
		/* DO NOT FLATTEN OUT e1 AND MODIFY echange
		// Do not run the following code unless you want to flatten e1 and care about e1 more than echange
		gsort +e1 // sort the final rank
		gen e1_order = _n
		replace e1 = (e1_order - 1)/(_N-1) // create final rank based on 
		drop e1_order	
		// Generate change in rank if we use the flattened e1.  
		// Calculate the mean and SD to see if it aligns with baboon data after flattening.	
			gen echange_final = e1 - e0
			la var echange_final "Final change in rank (after making e1 uniform)"
			su echange echange_final
		*/
	
	// Calculate absolute value of change using flattened e1
		gen ed = ((e1 - e0)^2)^0.5
		la var ed "|e1 - e0|"
	
	// Generate obs number variable for each covariate vector	
		sort e0
		gen obsID = _n // sorted by e0
		la var obsID "Regressor combination ID"
	
	// Generate more regressors --> combinations of regressors for null, DC, par

		// 2nd order
		gen e02 = e0^2
		la var e02 "e0^2"

		gen e12 = e1^2
		la var e12 "e1^2"

		gen ed2 = ed^2
		la var ed2 "|e1 - e0|^2"

		gen e0_e1 = e0 * e1
		la var e0_e1 "e1 x e0"
		
		gen e0_ed = e0 * ed
		la var e0_ed "e0 x ed"
		
		gen e1_ed = e1 * ed
		la var e1_ed "e1 x ed"
		
		// 3rd order
		gen e03 = e0^3
		la var e03 "e0^3"

		gen e13 = e1^3
		la var e13 "e1^3"

		gen ed3 = ed^3
		la var ed3 "|e1 - e0|^3"

		gen e0_e12 = e0 * (e1^2)
		la var e0_e12 "e0 x (e1^2)"
		
		gen e02_e1 = (e0^2) * e1
		la var e02_e1 "(e0^2) x e1"
		
		gen e0_ed2 = e0 * (ed^2)
		la var e0_ed2 "e0 * (|e1 - e0|^2)"

		gen e02_ed = (e0^2) * ed
		la var e02_ed "(e0^2) * |e1 - e0|"

		gen e1_ed2 = e1 * (ed^2)
		la var e1_ed2 "e1 * (|e1 - e0|^2)"

		gen e12_ed = (e1^2) * ed
		la var e12_ed "(e1^2) * |e1 - e0|"

		su e*
	
	* Save data with unique covariate vectors
	
		save vars_indep, replace

} // End toggle; checks out

********************
********************
* REALITIES
********************
********************

	********************
	* Generate coef files
	********************

	if `make_coef_10'==1 {

		* Primitives

		// Clear all variables
		clear *
		
		// Set range and increments of each coefficient
			local mincoef = -1 // min coefficient in power function
			local maxcoef = 1 // max coefficient in power function
			local increm = 0.5 // Increment for coefficients in power function
			local steps `=((`maxcoef'-`mincoef')/`increm')+1'

		// Specify number of and create counter for coefficeients
			local k 10 // 10 coef for 2 input & 3 powers or for 3 inputs and 2 powers
			local s = `k'-1 // useful when starting index at 0
			local o `=`steps'^`k''
			di "steps = `steps', k=`k', s=`s', o=`o'"

		* Generate coefficients
		
		// Initialise coefs (i.e., the vars representing coefs) at missing	
			set obs `steps'
			forval m = 0/`s' {
				gen a`m' = . // set equal to missing to start
				la var a`m' "Coef on x^`p'" // coefficient values
			} 
		
		// Generate variation in coefs
		// This is done in a sequance of expand commands; be sure to sort before each one!
			gen check = .
			gen obsID0 = _n

			forval i=1/`steps' {
				replace a0 = `mincoef'+((`i'-1)*`increm') if obsID0==`i'
			}

			forval j = 1/`s' {	
				local j_1 = `j'-1
				qui: expand `steps'
				hashsort obsID`j_1'
				qui: gen obsID`j' = _n
				qui: replace check = (obsID`j'-((obsID`j_1'-1)*`steps'))
				forval i=1/`steps' {
					qui: replace a`j' = `mincoef'+((`i'-1)*`increm') if (obsID`j'-((obsID`j_1'-1)*`steps'))==`i'
				}
			}

			su // check result

		* CLean up and save
			
		// Eliminate excess vars; add ID variable	
			drop obsID* check
			gen coefID = _n
			la var coefID "ID for coef vector"
		
		// Compress to save space, then save	
			compress 
			save coef_`k', replace

	} // end toggle; checks out: 9.765m obs (i.e., coef vectors or possible realities)

	*stop

	********************
	* To screen realities (coefficient vectors) for PAR/DC, create variation in covariates e0, ed
	********************

	// Generate 4 values of each of 2 covariates (total 16)
	// These covariates will be equally spaced because this is to test for PAR/DC in each reality
	// The goal is not (yet) to generate a realistic data set that examines whether IR, QR regression work 

	if `test_coef_10x16' == 1{ 

		// Load coefs file	
			use coef_10, clear
				
		* I have to create combinations of e0 and ed, so need 4*4 (if increment by 0.3333)
		
		// Code to reduce space for debugging		
			*keep a0 // just for testing code, keep subset of vars
			*cap drop *flag* *ID
			*keep if _n < 1001
		
		// Keep just the coef id to save space; you can merge in coefs later
			keep coefID 
		
		// Create 1st variable, call it "x", with 4 stop (0,.3,.6,1)
			timer on 1 // track time since these steps take time
			expand 4
			hashsort coefID // hashsort is faster (it's from gtools, so make sure that's installed and updated)
			gen long xID = _n
			gen byte x = ((xID - ((coefID-1)*4))-1) // x goes from 0(0.334)1 for each coef vector		
			timer off 1
		
		// Create 2nd variable, call it "z", with 4 stops 
			timer on 2
			expand 4 
			hashsort xID
			gen long zID = _n
			gen byte z = ((zID - ((xID-1)*4))-1) // z goes from 0(0.334)1 for each coef vector		
			timer off 2
			
		// Check times for each step
			timer list 
		
		// We worked in integers to save space and mem; now we divide by 3 to get fractions, i.e., change 1,...,4 --> 0,.3,.6,1		
			recast float x z // change type to allow fractions		
			gen e0 = x/3
			gen ed = z/3
			drop x z 
			rename xID e0ID
			rename zID edID
		
		// Summarize, and drop unnecessary vars		
			su coefID e0ID edID e0 ed		
			drop e0ID edID  

		// Save file of e0 and ed w/ coefID's, but not coefs themselves		
			compress
			hashsort coefID
			save coef_10x16, replace

	} // end toggle; checks out: 156m obs = 5^10 x 16

	*stop

	********************
	* Merge coefs and coariates for testing
	********************
		
	if `test_coef_10x16_merge' == 1 {
		
			// Load coef data			
				use coef_10, clear
			
			// Merge in e0, ed from coef_10x16			
				timer on 1
				hashsort coefID
				*gen w = runiform() // use these 2 lines to take a subsample fo data for debugging
				*keep if w < 0.01
				merge 1:m coefID using coef_10x16, keep(3) // NB: coef_10x16 sorted before saving
				
			// Conpress and save data
				compress
				save coef_10x16_merged, replace // we save file in case program crashes 
				timer off 1
			
			// Check times for each step
				timer list
			
	} // End toggle; Checks out 156m obs w inputs

	*stop

	********************
	* Screen out infeasible realities, test realities for PAR, DC
	********************
	
	if `screen_coef==1' {	

		// Load data
			use coef_10x16_merged, clear
		
		// Gen flags for range, monotonicity
			cap gen byte range_flagy = 0
			cap gen byte monot_flage0 = 0
			cap gen byte monot_flaged = 0
		
		// Clear timers
			timer clear

		// Generate outcome, test range // works for binary or fraction/share/percent outcome vars		
			timer on 1 // start timer for feasible y check
			gen float y = a0 + (a1*e0) + (a2*ed) + (a3*e0*e0) + (a4*ed*ed) + (a5*e0*ed) + (a6*(e0^3)) + (a7*(ed^3)) + (a8*(e0*ed*ed)) + (a9*(e0*e0*ed))
				// This simulates y (fertility)
			replace range_flagy = 1 if (y<0 | y>1) // flag any y that is outside acceptable range [0,1]
			hashsort coefID // use gtools version of sort bec faster with big data
			by coefID: gegen byte range_flag = max(range_flagy) 
				// flag any coef vector if any e0 or ed causes out of range error for coef vector
			su e0 ed if range_flag==0
			drop range_flagy
			timer off 1
			drop if range_flag==1 // save time on subsequent calcs
				// drop any coef vector if any e0 or ed causes out of range error for coef vector
				// Drop 154 million obs (before drop, total # obs = # coef vectors x 16 covariate)!
		
		// Test for DC, PAR		
			timer on 2 // start timer for DC and PAR flags
			gen float yde0 = a1 + (2*a3*e0) + (a5*ed) + (3*a6*e0*e0) + (a8*ed*ed) + (2*a9*e0*ed) // generate derivative wrt e0
			gen float yded = a2 + (2*a4*ed) + (a5*e0) + (3*a7*ed*ed) + (2*a8*e0*ed) + (a9*e0*e0) // generate derivative wrt ed
			replace monot_flage0 = 1 if sign(yde0) > 0 // flag any funcs that are monotone increasing in e0 (DC)
			replace monot_flaged = 1 if sign(yded) < 0 // flag any funcs that are monotone decreasing in ed (PAR)
			hashsort coefID 
			gegen byte dc_flag_meantrue = mean(monot_flage0), by(coefID) 
			gegen byte par_flag_meantrue = mean(monot_flaged), by(coefID) 
			gegen byte dc_flag_true = min(monot_flage0), by(coefID) // always true 
			gegen byte par_flag_true = min(monot_flaged), by(coefID) 
			gegen byte dc_flag_false = max(monot_flage0), by(coefID) 
			gegen byte par_flag_false = max(monot_flaged), by(coefID) 
			replace dc_flag_false = 1 - dc_flag_false // always false 
				// Have to do this in 2 steps because gegen does not allow this in 1 step
			replace par_flag_false = 1 - par_flag_false
			la var dc_flag_meantrue "DC model is true (on ave.)"
			la var par_flag_meantrue "PAR model is true (on ave.)"
			la var dc_flag_true "DC model is true (all regressor values)"
			la var par_flag_true "PAR model is true (all regressor values)"
			la var dc_flag_false "DC model is false (all regressor values)"
			la var par_flag_false "PAR model is false (all regressor values)"
			su e0 ed yde0 yded if dc_flag_meantrue==1
			su e0 ed yde0 yded if par_flag_meantrue==1
			drop monot_flage0 monot_flaged 
			timer off 2

		// Check times for each step
			timer list 
		
		// Keep coefs and DC PAR flags only
			local keepvars "a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 dc_flag* par_flag*"
			keep `keepvars'
		
		// Drop duplicates.  This drops all the covariate variation	used to test range and PAR/DC.
			gduplicates report `keepvars'
			gduplicates drop `keepvars', force // drops 1.95 m obs
			gen long coefID=_n
		
		// Check variables you keep, check how many realities have DC, PAR
			su `keepvars'
			tab dc_flag_meantrue par_flag_meantrue
		
		// Save coefs		
			save coef_flagged, replace
			*save /home/amalani/rosenbaum/par/coef_flagged, replace
			*cap rm coef_10x16_merged

	}  // end toggle 

	// 130,201 survive out of range restriction
	// 2% pass DC on ave, 6% fail DC on ave
	// 2% pass PAR on ave, 6% fail PAR on ave
	// They are symmetric because regression is symmetric in e0 and ed and we span e0 and ed symmetrically
	// But it's not the same set of coeff that pass each one:

	/*
	. tab dc_flag_meantrue par_flag_meantrue

	  DC model | PAR model is true (on
	   is true |         ave.)
	 (on ave.) |         0          1 |     Total
	-----------+----------------------+----------
			 0 |   124,749      2,697 |   127,446 
			 1 |     2,697         58 |     2,755 
	-----------+----------------------+----------
		 Total |   127,446      2,755 |   130,201 
	*/

	*stop
		
********************
********************
* SIMULATE IR REGS
********************
********************

********************
* Flesh out data set with error terms
********************

	if `sim_data'==1 {

		// Load screened coefs (ie realities)
			use "`datalib'/coef_flagged.dta", clear // 130k obs
		
		// Create space to merge in (realistic) indep variables for data set.
			recast int coefID
			sort coefID // recall this is a bew coefID that counts at coef vectors that meet tests
			*su coefID
			expand `regressors' // create space for all possible e0's -> 130k x 20 = 2.6m
			hashsort coefID // sort in coefID because expand tacks on end
			gen int obsID = _n
			replace obsID = (_n - ((coefID-1)*`regressors')) // create variable to m:1 merge on	
			*su // obsID should be from 1-`regressors'

		// Merge in indep vars for simulated regs
			merge m:1 obsID using "`datalib'/vars_indep.dta"
			tab _merge
			drop _merge
			drop e03 e13 ed3 e0_e12 e02_e1 e0_ed2 e02_ed e1_ed2 e12_ed echange* 
			compress	
			*su coefID obsID
			
		// Save data for use later
			hashsort coefID obsID // use to check data editor
			save "`datalib'/simdatafortest_noerror.dta", replace

		// Create space for more obs for each x; we will add in errors later
			hashsort coefID obsID // use to check data editor
			expand `u' // so this makes 2000 = 100x20 obs per coef vector
			// now at 52m obs
		
		// Take stock
			hashsort coefID obsID 	
			*su coefID obsID
			*su

		// Generate hypothetical y's
			gen float y = a0 + (a1*e0) + (a2*ed) + (a3*e0*e0) + (a4*ed*ed) + (a5*e0*ed) + (a6*(e0^3)) + (a7*(ed^3)) + (a8*(e0*ed*ed)) + (a9*(e0*e0*ed))
				// deterministic

		// Generate draws on an error term for each coef vector
		// We take a average organism, who starts at rank 0 and sees no change.  What is it's outcome?
		// We can calculate the variance of the error from that mean outcome.
		// We draw normal errors with mean 0 and the variance just calculated.
		// Why normal?  Number of positve draws from n binomials with mean p converges to normal with mean np and variance np(1-p).   
			gen A = a0 + (a1*0.5) + (a3*0.25) + (a6*0.125) // calculate the mean at the average value of x 
			gen s2 = A*(1-A) // calculate the binomial varance 
			su s2
			gen u = rnormal(0,s2) // generate a random variable with mean 0 and variance equal to the bonomial variance
			/* We do not do this:
				* Use variance of residual from reg of fertility on rank 1 to add randomness
				* Mean (SD) of residual is 0 (.3524) if repstats_state=1
			*/
		
		// Calculate outcomes with this variance
			replace y = y + u
			la var y "Fertility"
		
		// Save data
			save "`datalib'/simdatafortest.dta", replace

	}

	*stop

********************
* Test interaction and quadratic regression models
********************
	
	if `test_regs'==1 { // 
		
		*******************
		* Run regressions *
		*******************

		// Load sim data	
			use "`datalib'/simdatafortest.dta", clear // n = 260,402,000

		// Run IR reg, save results
			*regressby y e0 e1 e0_e1, by(coefID) save(ir_results)
		
		// Do QR regression, save results
			*regressby y e0 ed e02 ed2 e0_ed, by(coefID) save(qr_results)
		
		
		// Run visualization test regardless of whether IR regression yield sig coef on interaction term
			// We'll screen on sig coef on interaction later
			// Check how outcomes change with adult environment
			// First, for those who start in low quality environments
			*keep if e0 < 0.5 // regressby does not allow if statements
			*regressby y e1, by(coefID) save(vislow_results) 
			//Second, for thsoe who start in high quality environments
			*suse "`datalib'/simdatafortest.dta", clear // n = 260,402,000
			keep if e0 > 0.5
			regressby y e1, by(coefID) save(vishigh_results)
		
		*stop
		
		
		***********************
		* Analyze regressions *
		***********************

		// Interaction regression
		
			// Load IR results 			
			use "`datalib'/ir_results.dta", clear
			
			// Conventional test for PAR with IR
			// STEP ONE: Is interaction sig?
			gen ir_par_inter_flag = abs(_b_e0_e1 / _se_e0_e1) > 1.96 // flag for interaction term being sig
			la var ir_par_inter_flag "Positive test - sig interaction term with IR reg"
			// STEP TWO comes later
			
			// Conventional test for DC with IR
			gen ir_dc_conv_flag = abs(_b_e0 / _se_e0) > 1.96 // flag for coef on e0 term being sig
			la var ir_dc_conv_flag "Positive test - sig coef on e0 with IR reg"
			
			// Definition aligned test for PAR with IR reg:
			// Is derivative of IR reg wrt ed (|echange|) negative at mean values of e0 and ed?
			gen ir_par_def_mean = (_b_e1 + (_b_e0_e1 * 0.5)) * (-1) 
				// multiply by -1 because mean value of echange is -0.03
			gen ir_par_def_var = (_se_e1^2) + ((_se_e0_e1 * 0.5)^2) + 2*_cov_e0_e1_e1*0.5
			gen ir_par_def_se = sqrt(ir_par_def_var)
			gen ir_par_def_flag = ir_par_def_mean / ir_par_def_se < -1.96 
				// negative and sig
			la var ir_par_def_flag "Positive test - definition aligned PAR test with IR reg"

			//Definition aligned test for DC with IR				
			gen ir_dc_def_mean = _b_e0 + (_b_e0_e1 * 0.47) // 0.47 because mean of echange is -0.03 
			gen ir_dc_def_var = (_se_e0^2) + ((_se_e0_e1 * 0.47)^2) + 2*_cov_e0_e1_e1*0.47
			gen ir_dc_def_se = sqrt(ir_dc_def_var)
			gen ir_dc_def_flag = ir_dc_def_mean / ir_dc_def_se > 1.96 
				// positive and sig
			la var ir_dc_def_flag "Positive test - definition aligned DC test with IR reg"
			
			// Save flags only in results file			
			keep coefID ir_*
			sort coefID
			save "`datalib'/test_results.dta", replace 
		
		// Visualization test
		
			// Load lower visualization datalib, rename coefs to mark them as lower			
			use "`datalib'/vislow_results.dta", clear
			rename _b_cons _b_cons_low
			rename _b_e1 _b_e1_low

			// Merge higher visualization datalib, rename coefs to mark them as higher 			
			sort coefID
			merge 1:1 coefID using "`datalib'/vishigh_results.dta"
			tab _merge
			drop _merge
			rename _b_cons _b_cons_high
			rename _b_e1 _b_e1_high
					
			// Flag whether high-birth organisms starts lower and ends higher
			gen vis_startlower_flag = _b_cons_high < _b_cons_low
			gen vis_endhigher_flag = _b_cons_high + _b_e1_high > _b_cons_low + _b_e1_low
			gen vis_flag = (vis_startlower_flag == 1 & vis_endhigher_flag == 1)
			la var vis_flag "Positive test - visualization test" 
			
			// Flag whether high (low) birth organisms have positive (negative) slope
			gen vis_strict_flag = (_b_e1_high > 0 & _b_e1_low < 0)
			la var vis_strict_flag "Positive test - strict visualization test"
			
			// Flag whether slope of 4th quartile is higher (relaxed visualization test)
			gen vis_relaxed_flag = (_b_e1_high > _b_e1_low)
			la var vis_relaxed_flag "Positive test - relaxed visualization test" 
			
			// Merge results with results file			
			keep coefID vis_flag vis_strict_flag vis_relaxed_flag
			sort coefID
			merge 1:1 coefID using "`datalib'/test_results.dta"
			table _m
			drop _merge
			save "`datalib'/test_results.dta", replace 
			
			// Generate conventional PAR test for IR reg			
			gen ir_par_conv_flag = (ir_par_inter_flag == 1 & vis_flag == 1)
			la var ir_par_conv_flag "Positive test - conventional PAR test with IR reg"
		
			// Generate strict conventional PAR test for IR reg			
			gen ir_par_conv_strict_flag = (ir_par_inter_flag == 1 & vis_strict_flag == 1)
			la var ir_par_conv_strict_flag "Positive test - strict conventional PAR test with IR reg"

			// Generate relaxed conventional PAR test for IR reg			
			gen ir_par_conv_relaxed_flag = (ir_par_inter_flag == 1 & vis_relaxed_flag == 1)
			la var ir_par_conv_relaxed_flag "Positive test - relaxed conventional PAR test with IR reg"

			// Save results from IR tests
			save "`datalib'/test_results_temp.dta", replace 			
			
			
		// Quadratic regression
		
			// Load QR results 			
			use "`datalib'/qr_results.dta", clear	
			
			// Definition aligned test for PAR with QR				
			gen qr_par_def_mean = _b_ed + (2 * _b_ed2 * 0.03) + (_b_e0_ed * 0.5) 
			gen qr_par_def_var = (_se_ed^2) + ((2 * _se_ed2 * 0.03)^2) ///
				+ ((_se_e0_ed * 0.5)^2) ///
				+ (2 * _cov_ed2_ed * 0.03) ///
				+ (2 * _cov_e0_ed_ed2 * 0.5) ///
				+ (2 * _cov_e0_ed_ed2 * 0.03 * 0.5)
			gen qr_par_def_se = sqrt(qr_par_def_var)
			gen qr_par_def_flag = (qr_par_def_mean / qr_par_def_se < -1.96) 
				// negative and sig
			la var qr_par_def_flag "Positive test - def aligned PAR test with QR reg"
			
			// Def aligned test for DC with QR				
			gen qr_dc_def_mean = _b_e0 + (2 * _b_e02 * 0.5) + (_b_e0_ed * 0.03) 
			gen qr_dc_def_var = (_se_e0^2) + ((2 * _se_e02 * 0.5)^2) ///
				+ ((_se_e0_ed * 0.03)^2) ///
				+ (2 * _cov_e02_e0 * 0.5) ///
				+ (2 * _cov_e0_ed_e0 * 0.03) ///
				+ (2 * _cov_e0_ed_e02 * 0.5 * 0.03)
			gen qr_dc_def_se = sqrt(qr_dc_def_var)
			gen qr_dc_def_flag = (qr_dc_def_mean / qr_dc_def_se > 1.96) 
				// positive and sig
			la var qr_dc_def_flag "Positive test - def aligned DC test with QR reg"
			
			// Save flags only			
			keep coefID qr_*
			
			// Merge in IR test results
			sort coefID
			merge 1:1 coefID using "`datalib'/test_results_temp.dta"
			table _m
			drop _merge
			
			// Save test results with both IR and QR tests
			rm "`datalib'/test_results_temp.dta"
			save test_results_temp, replace 
			
		// Coin flip

			// Generate coin-flip test			
			gen coinflip_flag = rbinomial(1,0.5)
			la var coinflip_flag "Positive test - from coin flip"
		
		
		***********************
		* Sensitivity and specificity
		***********************
		
			// Merge in truth, keeping only the DC and PAR true flags (at mean)			
			merge 1:1 coefID using "`datalib'/coef_flagged.dta", keepusing(dc_flag_meantrue par_flag_meantrue)
			tab _m
			drop _merge
			
			// Tabulate true and false positives			
			gen total = 1 // counter
			table par_flag_meantrue, stat(fvpercent total coinflip_flag ir_par_conv_strict_flag ) 
			table par_flag_meantrue, stat(fvpercent total ir_par_conv_flag ir_par_conv_relaxed_flag ) 
			table par_flag_meantrue, stat(fvpercent total  ir_par_def_flag qr_par_def_flag) 
			table dc_flag_meantrue, stat(fvpercent total coinflip_flag ir_dc_conv_flag)
			table dc_flag_meantrue, stat(fvpercent total ir_dc_def_flag qr_dc_def_flag)

	}

	*stop

