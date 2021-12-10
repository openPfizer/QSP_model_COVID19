# QSP model of COVID-19 (December 2021) - A Quantitative Systems Pharmacology Model of the Pathophysiology and Treatment of COVID-19 Predicts Optimal Timing of Pharmacological Interventions

## Description
The model mechanistically accounts for the influence of key mediators relevant to COVID-19 pathophysiology including, interactions between viral dynamics, the major host immune response mediators, and alveolar tissue damage and regeneration.

## Primary results
The code generates figures in the manuscript titled *A Quantitative Systems Pharmacology Model of the Pathophysiology and Treatment of COVID-19 Predicts Optimal Timing of Pharmacological Interventions*

medRxiv doi: https://doi.org/10.1101/2021.12.07.21267277

## Prerequisites
MATLAB

This code was written in MATLAB 2019b

## Setup
Add all files/directories in this repository to the MATLAB working directory/path.

## Contents
The repository should contain the following required files:

1. **covid19_dxdt.m -> *model file with ODEs***
2. **BE_Blaze1Ph3.m -> *main script to generate figures from Blaze-1 Ph3 trial***
3. **dVL_RRR_timing.m -> *main script to generate figure to determine sensitivity of viral load reduction and RRR to timing of intervention*** 
4. **molnupiravir_sfig.m -> *script to generate supplementary figures for molnupiravir virtual population***
5. **plausible_figure.m -> *main script to generate figure for plausible population***
6. **regen_cov.m -> *main script to generate REGEN-COV figures*** 
7. **severity_plot.m -> *main script to generate figure for RRR***
8. update_parameters_ext.m -> *updates model dictionary during virtual population simulations*
9. SolveBalances.m -> *calls ODE solver*
10. Initialize.xlsx -> *parameter and initial condition file*
11. covidEventFcn.m -> *checks whether virus is below physiological levels*
12. function_run_model_noplots.m -> *returns healthy and COVID-19 ODE solutions*
13. get_data_dictionary.m -> *loads model dictionary with parameter values, initial conditions & additional simulations settings*
14. adjust_tfso.m -> adjust time from infection to time from symptom onset for virtual population
15. merge_optimdata.m -> format observational COVID-19 studies for plotting plausible population
16. generate_figures.m -> script to run driver files and save figures
17. blaze1.mat -> prerequisite mat file for Blaze-1 simulations
18. plausible_fig.mat -> prerequisite mat file for plausible figure
19. regen_cov.mat -> prerequisite mat file for REGEN-COV simulations
20. eidd.mat - > prerequisite mat file for molnupiravir simulations
21. delta_variant.mat -> prerequisite mat file for supplementary preliminary delta variant simulations 
22. Add folder Violinplot-Matlab-master to path to plot figures with violin plots in manuscript

## Usage
Running the command below generates figures and saves them as .png files.

```matlab
run generate_figures.m
```

Running all scripts to generate the figures in the paper takes approximately 2h on a 2019 Macbook Pro (2.4 GHz 8-Core Intel Core i9)
## Authors
Rohit Rao*, CJ Musante, Richard Allen*

\*Correspondence to: rohit.rao@pfizer.com or richard.allen@pfizer.com

## Acknowledgements
We sincerely thank Annaliesa Anderson, Arthur Bergman, Britton Boras, Phylinda Chan, Wei Dai, Bharat Damle, Sandeep Menon, Gianluca Nucci, Theodore Rieger, Ravi Singh, Nessy Tania and RES group for their comments and feedback on the manuscript and during the development of the model. 
