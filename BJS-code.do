* Set directory to location of file
*cd "~/Dropbox/HetEventStudies/"

use "output/df.dta", clear

** BJS event-study
did_imputation y i t g, allhorizons pretrends(15)

event_plot, default_look graph_opt(xtitle("Periods since the event") ytitle("Average causal effect") ///
	title("Borusyak et al. (2021) imputation estimator") xlabel(-15(1)9))

graph export "figures/bjs.png", replace

** dCDH event-study in Stata
* Note: does not allow for more placebos than dynamic
did_multiplegt y i t d, robust_dynamic dynamic(9) placebo(9)
