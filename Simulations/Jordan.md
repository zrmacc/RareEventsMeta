# TODO for week of 10/15/22

* Run the new comparison_sim.R file
* For each setting do 2000 simulations (reps) with 500 iterations (mc).
* Settings to consider:
  * Case I: Very low base rate: rate = 0.003
        * Case Ia: No effect, high heterogeneity with alpha =  3, beta = 3
        * Case Ib: No effect, moderate heterogeneity with alpha =  8.5, beta = 8.5
        * Case Ic: No effect, low heterogeneity with alpha = 20, beta = 20
        * Case Id: Effect, high heterogeneity with alpha =  6, beta = 3
        * Case Ie: Effect, moderate heterogeneity with alpha =  17, beta = 8.5
        * Case If: Effect, low heterogeneity with alpha =  40, beta = 20
        
  * Case II: Low base rate: rate = 0.006
        * Case IIa: No effect, high heterogeneity with alpha =  3, beta = 3
        * Case IIb: No effect, moderate heterogeneity with alpha =  8.5, beta = 8.5
        * Case IIc: No effect, low heterogeneity with alpha = 20, beta = 20
        * Case IId: Effect, high heterogeneity with alpha =  6, beta = 3
        * Case IIe: Effect, moderate heterogeneity with alpha =  17, beta = 8.5
        * Case IIf: Effect, low heterogeneity with alpha =  40, beta = 20
        
  * Case III: Moderate base rate: rate = 0.015  
        * Case IIIa: No effect, high heterogeneity with alpha =  3, beta = 3
        * Case IIIb: No effect, moderate heterogeneity with alpha =  8.5, beta = 8.5
        * Case IIIc: No effect, low heterogeneity with alpha = 20, beta = 20
        * Case IIId: Effect, high heterogeneity with alpha =  6, beta = 3
        * Case IIIe: Effect, moderate heterogeneity with alpha =  17, beta = 8.5
        * Case IIIf: Effect, low heterogeneity with alpha =  40, beta = 20
        
  * Case IV: "High" base rate: rate = 0.040 
        * Case IVa: No effect, high heterogeneity with alpha =  3, beta = 3
        * Case IVb: No effect, moderate heterogeneity with alpha =  8.5, beta = 8.5
        * Case IVc: No effect, low heterogeneity with alpha = 20, beta = 20
        * Case IVd: Effect, high heterogeneity with alpha =  6, beta = 3
        * Case IVe: Effect, moderate heterogeneity with alpha =  17, beta = 8.5
        * Case IVf: Effect, low heterogeneity with alpha =  40, beta = 20

* The no effect settings will give you type I error
* The effect settings will give you power
* See the Results_120922.Rmd in the "Results" folder to see how to summarize results for type I error and power
* Make an Rmd file with the plots and summaries described in the bullets below:
* For each case, make a faceted bar plot for the type I error (eg. a-c) with faceting by hetereogetiy, each bar corresponds to the type I error of a given method
* For each case, make a faceted bar plot for the power (eg. d-f) with faceting by hetereogetiy, each bar corresponds to the power of a given method.  For methods with type I error > 0.05 (ie. they are anti-conservative), make the bars gray as it isnt far to compare power
* Summarize the effective K for each setting, this is the first column of "all_res", just take the mean across the simulations and print this out in the rmd

# TODO for week of 10/18/22

* Finish running Setting 1a, balanced for all values of K (K = 16, 24, 48)
* Then move to Setting 1a, unbalanced and run for all values of K
* Check if results are similar to Table 2 in the paper
* If so, then complete all the settings in Table 2
* From the paper, you can see how to change the parameters
* If all goes well, repeat and replicate Table 3

# TODO for week of 10/4/22

This repo contains the R package for out method.

There are two main files to run analyses on this repo: cp_sim.R and ci_sim.R.  ci_sim.R is an updated version of the file you were running previously while cp_sim.R is an abbreviated version that only calculates coverage probability.  We created this is designing simulation studies.

For next week go to Simulations/Cluster Scripts:

* This is an example of how to run cp_sim.R 

* Start by getting this to run on the cluster

* Then create analogous files for ci_sim.R and repeat Setting 1a (balanced) that you performed last week using the old version of the code

* Add the files you create to this repo in Simulations/Cluster Scripts

* Create a markdown file summarizing coverage probability and confidence interval length that contains these results as well as the old version you ran last week
